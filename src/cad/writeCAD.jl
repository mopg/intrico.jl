# ---------------------------------------------------------------------------- #
#
#   writeCAD.jl
#
#   Couple to jegads.jl to write a STEP file of the geometry.
#
#   sterno
#   Spring 2018
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

function writeCAD( mesh::MeshF,  lat::Lattice, flname::String; eps = 5.e-3, Δeps = 5.e-3 )

    # open egads
    (context, status) = jegads.EG_open( )
    if (status < 0) error("Can't open, failure code: %i", status) end

    @printf( "CAD: %d nodes, %d edges\n", mesh.n, size(mesh.e,1) )

    # 1. connect all nodes to their respective half-length edges
    println("CAD:   generating nodes")
    nodes = genNodesCAD( mesh, lat, context, eps, Δeps )

    # 2. connect those nodal hubs together to generate model
    println("CAD:   connecting nodes")
    finmodel = genModelCAD( mesh, nodes, context )

    # 3. shave off nodes interfering with boundary
    # TODO

    # 4. save model
    println("CAD:   saving model")
    if isfile( flname )
        println("   File exists, overwrite? (y/N)")
        user_input = chomp(readline())
        if isempty(user_input)
            # don't overwrite
            println("   Not overwriting, data lost")
        elseif user_input[1] == 'Y' || user_input[1] == 'y'
            rm( flname )
            status = jegads.EG_saveModel( finmodel, flname )
            if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end
        else
            println("   Not overwriting, data lost")
        end

    else
        status = jegads.EG_saveModel( finmodel, flname )
        if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end
    end

    # close everything
    status = jegads.EG_deleteObject( finmodel )
    if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end
    # delete nodes
    for kk in 1:length( nodes )
        status = jegads.EG_deleteObject( nodes[kk] )
        if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end
    end
    jegads.cleanup( context )

end

function genNodesCAD( mesh::MeshF3D, lat::Lattice, context::jegads.ego, eps::Float64, Δeps::Float64 )

    nodes = fill( jegads.ego(0), mesh.n )

    for ii in 1:mesh.n
        @printf(".")
        ar_max = maximum( lat.ar[ mesh.n2e[ii] ] )

        status_boolean = Cint(1)

        emodel = jegads.ego(0)
        ebody  = jegads.ego(0)

        inod = 0

        while (status_boolean != jegads.EGADS_SUCCESS)

            # Generate sphere
            datasph = vcat( mesh.p[ii,:], sqrt(ar_max/π)*(1.0 + eps) )

            (ebody, status) = jegads.EG_makeSolidBody(context, jegads.SPHERE, datasph)
            if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

            # Generate half length rods
            for jj in 1:length(mesh.n2e[ii])
                ie   = mesh.n2e[ii][jj]
                nods = mesh.e[ie,1:2]
                inod = find( nods .== ii )[]
                ave = mesh.p[nods[1],:] + 0.50 * (mesh.p[nods[2],:] - mesh.p[nods[1],:])
                datacyl = vcat( mesh.p[ nods[1] ,:], # always generate this half cylinder, but then move it into the correct pos
                                ave,
                                sqrt( lat.ar[ie] / π ) )

                (etemp,status) = jegads.EG_makeSolidBody(context, jegads.CYLINDER, datacyl)
                if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

                if inod == 2
                    Δx = 0.50 * (mesh.p[nods[2],:] - mesh.p[nods[1],:])

                    matrix = [ 1.0, 0.0, 0.0, Δx[1],
                               0.0, 1.0, 0.0, Δx[2],
                               0.0, 0.0, 1.0, Δx[3] ]

                    exform_ptr = Ref{jegads.ego}()
                    status = jegads.EG_makeTransform(context, matrix, exform_ptr)
                    if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

                    (etemp2,status) = jegads.EG_copyObject(etemp, exform_ptr[])
                    if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

                    status = jegads.EG_deleteObject(exform_ptr[])
                    if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

                    status = jegads.EG_deleteObject(etemp)
                    if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

                    (etemp,status) = jegads.EG_copyObject(etemp2, C_NULL)
                    if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

                    status = jegads.EG_deleteObject(etemp2)
                    if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

                end

                # join to rest of node
                (emodel, status) = jegads.EG_solidBoolean(ebody, etemp, jegads.FUSION)
                # if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end
                status_boolean = status

                # delete temporary bodies
                status = jegads.EG_deleteObject(ebody)
                if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end
                status = jegads.EG_deleteObject(etemp)
                if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

                if (status_boolean != jegads.EGADS_SUCCESS)
                    eps += Δeps
                    @printf("   NOTE: increasing radius of node %d by %3.2f%%\n", ii, Δeps*100.)
                    break
                end

                # get body from model
                (ebody, status) = jegads.getBodyFromModel( emodel )
                if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

            end

        end

        # save model
        (nodes[ii], status) = jegads.EG_copyObject(emodel, C_NULL)
        if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

        # delete body and model
        status = jegads.EG_deleteObject(ebody)
        if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end
        status = jegads.EG_deleteObject(emodel)
        if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

    end

    @printf("\n")

    return nodes # NOTE: nodes are models

end

function genModelCAD( mesh::MeshF3D, nodes::Vector{jegads.ego}, context::jegads.ego )

    # nod1 =  7
    # nod2 = 18
    #
    # # get bodies from model
    # (ebody1,status) = jegads.getBodyFromModel( nodes[ nod1 ] )
    # if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end
    # (ebody2,status) = jegads.getBodyFromModel( nodes[ nod2 ] )
    # if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end
    #
    # emodel, status = jegads.EG_join( ebody1, ebody2, Cdouble(0.0) )

    # NOTE: This doesn't work because the faces are rotated with respect to each other. How to fix this????

    efaces_ptr = Vector{ Ptr{jegads.ego} }( mesh.n )
    nface_vec  = Vector{Int64}( mesh.n )

    # get faces of all nodes first
    toler = Cdouble(0.0)
    for kk in 1:mesh.n

        # get body from model
        (ebody,status) = jegads.getBodyFromModel( nodes[ kk ] )
        if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

        # get faces from body
        nface_vec[kk], ef_ptr, status = jegads.EG_getBodyFaces( ebody, toler )
        if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

        efaces_ptr[kk] = ef_ptr[]

    end

    # get matching faces
    tolerf = Cdouble(0.0); cntmax = 14
    matches_vec = Vector{ Vector{Int64} }( mesh.n )
    for ii in 1:mesh.n
        matches_vec[ii] = fill( 0, 0 )
    end
    for jj in 1:size(mesh.e,1)

        nod1 = mesh.e[jj,1]
        nod2 = mesh.e[jj,2]

        # get bodies from model
        (ebody1,status) = jegads.getBodyFromModel( nodes[ nod1 ] )
        if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end
        (ebody2,status) = jegads.getBodyFromModel( nodes[ nod2 ] )
        if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

        # match faces
        nmatch  = 0
        matches = fill( 0, 0 )

        cnt = 0
        while nmatch != 1 && cnt < cntmax
            (nmatch,matches,status) = jegads.EG_matchBodyFaces( ebody1, ebody2, tolerf )
            if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end
            if cnt == 0
                tolerf += 1e-15
            else
                tolerf *= 10.0
                warn("WriteCAD: Was not able to find matching faces between nodes ",
                    nod1, " and ", nod2, " on edge ",jj, " increasing tolerance to ", tolerf)
                println("nmatch: ",nmatch)

            end
            cnt += 1
        end
        if cnt == cntmax
            error("WriteCAD: Was not able to find matching faces between nodes ",
                nod1, " and ", nod2, " on edge ",jj)
        end
        tolerf *= 0.0

        append!( matches_vec[nod1], matches[1] )
        append!( matches_vec[nod2], matches[2] )

        # delete temporary bodies
        # status = jegads.EG_deleteObject(ebody1)
        # if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end
        # status = jegads.EG_deleteObject(ebody2)
        # if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

    end

    # generate list of faces
    efaces = Vector{ jegads.ego }( sum(nface_vec) )
    fcnt = 1
    for kk in 1:mesh.n

        for ii in 1:nface_vec[kk]
            if !any( ii .== matches_vec[kk] )
                efaces[fcnt] = unsafe_load(efaces_ptr[kk],ii); fcnt += 1
            end
        end

    end
    fcnt -= 1

    # join model
    tolerj = Cdouble(0.0)
    (emodel, status) = jegads.EG_sewFaces( Cint(fcnt), efaces, tolerj, Cint(0) )
    if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

    return emodel

end
