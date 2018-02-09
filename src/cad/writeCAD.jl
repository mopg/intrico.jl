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

    # 1. connect all nodes to their respective half-length edges

    nodes = genNodesCAD( mesh, lat, context )

    # status = jegads.EG_saveModel( nodes[1], flname )
    # if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end
    # status = jegads.EG_saveModel( nodes[2], "test2.egads" )
    # if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

    # 2. connect those nodal hubs together to generate model
    finmodel = genModelCAD( mesh, nodes, context )

    # 3. shave off nodes interfering with boundary

    # 4. save model
    println("   saving model")
    status = jegads.EG_saveModel( finmodel, flname )
    if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

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

    eps = 1.0e-2

    for ii in 1:mesh.n

        println("node ",ii)

        ar_max = maximum( lat.ar[ mesh.n2e[ii] ] )

        # Generate sphere
        datasph = vcat( mesh.p[ii,:], sqrt(ar_max/π)*(1.0 + eps) )
        println("   datasph ",datasph)
        (ebody, status) = jegads.EG_makeSolidBody(context, jegads.SPHERE, datasph)
        if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

        # Generate half length rods
        emodel = jegads.ego(0)
        for jj in 1:length(mesh.n2e[ii])
            ie  = mesh.n2e[ii][jj]
            ave = (mesh.p[mesh.e[ie,1],:] + mesh.p[mesh.e[ie,2],:]) / 2.0
            datacyl = vcat( mesh.p[ii,:],
                            ave,
                            sqrt( lat.ar[ie] / π ) )
            println("   datacyl ", datacyl)
            (etemp, status) = jegads.EG_makeSolidBody(context, jegads.CYLINDER, datacyl)
            if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end
            println("   solid body ", jj)
            # join to rest of node
            (emodel, status) = jegads.EG_solidBoolean(ebody, etemp, jegads.FUSION)
            if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end
            println("   solid boolean ", jj)
            # delete temporary bodies
            status = jegads.EG_deleteObject(ebody)
            if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end
            status = jegads.EG_deleteObject(etemp)
            if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

            # get body from model
            (ebody, status) = jegads.getBodyFromModel( emodel )
            if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end

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

    return nodes # NOTE: nodes are models

end

function genModelCAD( mesh::MeshF3D, nodes::Vector{jegads.ego}, context::jegads.ego )

    println("joining model ")

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

        println("   acnodes ", acnodes )
        println("   nodall  ", nodall )

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
            println("   solid boolean ", kk)

            # delete temporary bodies
            status = jegads.EG_deleteObject(ebody)
            if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end
            status = jegads.EG_deleteObject(ebody2)
            if (status != jegads.EGADS_SUCCESS) jegads.cleanup(status, context) end
            println("   done")
        end

    end

    return emodel

end
