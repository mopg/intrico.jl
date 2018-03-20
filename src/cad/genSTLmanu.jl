
function genSTLmanu( mesh::MeshF,  lat::Lattice, flname::String, manu::Vector{Bool}; n = 8, name = "object" )

    # 1. compute correct distance from node
    fdist = compDist( mesh, lat )

    # open STL file
    fid  = open( string(flname[1:end-4], "_manu.stl"), "w" )
    fid2 = open( string(flname[1:end-4], "_nmanu.stl"), "w" )
    @printf( fid,  "solid %s\n", name )
    @printf( fid2, "solid %s\n", name )

    # generate STL for cylinders
    genSTLcyls( mesh, lat, fdist, n, fid, fid2, manu )

    # generate STL for nodes
    genSTLAnods( mesh, lat, fdist, n, fid, fid2, manu )

    # close STL
    @printf( fid,  "endsolid %s\n", name )
    @printf( fid2, "endsolid %s\n", name )
    close( fid )

end

function genSTLcyls( mesh::MeshF,  lat::Lattice, fdist::Matrix{Float64}, n::Int64,
                     fid::IOStream, fid2::IOStream, manu::Vector{Bool} )

    θ = linspace( 0, 2*π, n )

    for ee in 1:size( mesh.e, 1 )

        if lat.ar[ee] < 0.0 || sqrt( lat.ar[ee] / π ) < 1e-14
            continue
        end

        nod1 = mesh.e[ee,1]
        nod2 = mesh.e[ee,2]

        vec = mesh.p[nod2,:] - mesh.p[nod1,:]

        l = norm( vec )

        # begin and start nodes of cylinder
        x1 = mesh.p[ nod1, : ] + fdist[ee,1] / l * vec
        x2 = mesh.p[ nod2, : ] - fdist[ee,2] / l * vec

        # generate points on end-faces
        r = sqrt( lat.ar[ee] / π )
        xx = r * cos.( θ )
        yy = r * sin.( θ )

        nvec = vec / l

        xx, yy, zz = rotTransCirc( xx, yy, nvec )
        #   translate from origin to end points
        xx1 = xx + x1[1]
        yy1 = yy + x1[2]
        zz1 = zz + x1[3]
        xx2 = xx + x2[1]
        yy2 = yy + x2[2]
        zz2 = zz + x2[3]

        if manu[ee]
            writeSTLcyl( Val{1}, xx1, yy1, zz1, xx2, yy2, zz2, fid )
        else
            writeSTLcyl( Val{1}, xx1, yy1, zz1, xx2, yy2, zz2, fid2 )
        end

    end

end

function genSTLAnods( mesh::MeshF, lat::Lattice, fdist::Matrix{Float64}, n::Int64,
                      fid::IOStream, fid2::IOStream, manu::Vector{Bool} )

    θ = linspace( 0, 2*π, n ).^1.0

    for nn in 1:mesh.n

        edges = copy(mesh.n2e[nn])
        inddel = Vector{Int64}(0)
        for ii in 1:length(edges)
            if lat.ar[ edges[ii] ] < 0.0 || sqrt( lat.ar[ edges[ii] ] / π ) < 1e-14
                append!(inddel,ii)
            end
        end

        manuNod = all( manu[ mesh.n2e[nn] ] )

        deleteat!(edges,inddel)

        (verts,nvecs) = compNodSingle( mesh, lat, edges, fdist, nn, n, θ )

        if length(edges) > 1
            # compute convex hull
            ch = QHull.chull( verts )
            faces = ch.simplices

            centroid = mean( verts, 1 )[:] # need to compute centroid to determine whether or not normal is correct

            ifaces = cleanupCHull( faces, verts, nvecs, centroid )

            # write STLs
            for ff in ifaces
                if manuNod
                    writeFacetSTLA( verts[ faces[ff],: ], fid::IOStream )
                else
                    writeFacetSTLA( verts[ faces[ff],: ], fid2::IOStream )
                end
            end

        end

    end

end
