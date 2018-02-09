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

function writeCAD( mesh::MeshF,  lat::Lattice, flname::String )

    # open egads
    (context, status) = jegads.EG_open( )
    if (status < 0) error("Can't open, failure code: %i", status) end

    @printf( "CAD: %d nodes, %d edges\n", mesh.n, size(mesh.e,1) )

    # 1. connect all nodes to their respective half-length edges
    println("CAD:   generating nodes")
    nodes = genNodesCAD( mesh, lat, context )

    # 2. connect those nodal hubs together to generate model
    println("CAD:   connecting nodes")
    finmodel = genModelCAD( mesh, nodes, context )

    # 3. shave off nodes interfering with boundary
    # TODO

    # 4. save model
    println("CAD:   saving model")
    status = jegads.EG_saveModel( finmodel, flname )
    # if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end # commented because otherwise context dies

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

function genNodesCAD( mesh::MeshF3D, lat::Lattice, context::jegads.ego )

    nodes = fill( jegads.ego(0), mesh.n )

    for ii in 1:mesh.n
        @printf(".")
        ar_max = maximum( lat.ar[ mesh.n2e[ii] ] )

        status_boolean = Cint(1)

        eps = 5.0e-3 # 1.0e-2
        emodel = jegads.ego(0)
        ebody  = jegads.ego(0)

        while (status_boolean != jegads.EGADS_SUCCESS)

            # Generate sphere
            datasph = vcat( mesh.p[ii,:], sqrt(ar_max/π)*(1.0 + eps) )

            (ebody, status) = jegads.EG_makeSolidBody(context, jegads.SPHERE, datasph)
            if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

            # Generate half length rods
            for jj in 1:length(mesh.n2e[ii])
                ie  = mesh.n2e[ii][jj]
                ave = (mesh.p[mesh.e[ie,1],:] + mesh.p[mesh.e[ie,2],:]) / 2.0
                datacyl = vcat( mesh.p[ii,:],
                                ave,
                                sqrt( lat.ar[ie] / π ) )

                (etemp, status) = jegads.EG_makeSolidBody(context, jegads.CYLINDER, datacyl)
                if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

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
                    eps += 5.0e-3
                    println("   NOTE: increasing radius by 0.5\% of node ", ii)
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

    nod = 1
    acnodes = [nod]

    (emodel, status) = jegads.EG_copyObject( nodes[nod], C_NULL )
    if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

    nrem = mesh.n - nod

    while nrem > 0

        nodall = fill( 0, 0 )

        for ii in 1:length(acnodes), jj in 1:length( mesh.n2e[ acnodes[ii] ] )
            ie = mesh.n2e[ acnodes[ii] ][jj]
            append!( nodall, mesh.e[ie,1:2] )
        end

        nodall = unique( nodall )
        deleteat!(nodall, findin(nodall, acnodes) )

        nrem -= length(nodall)

        acnodes = append!( acnodes, nodall )

        for kk in 1:length(nodall)

            # get body from model
            (ebody, status) = jegads.getBodyFromModel( emodel )
            if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end
            (ebody2,status) = jegads.getBodyFromModel( nodes[ nodall[kk] ] )
            if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

            # join to rest of node
            (emodel, status) = jegads.EG_solidBoolean(ebody, ebody2, jegads.FUSION)
            if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

            # delete temporary bodies
            status = jegads.EG_deleteObject(ebody)
            if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end
            status = jegads.EG_deleteObject(ebody2)
            if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

            @printf(".")
        end

    end

    @printf("\n")

    return emodel

end
