# ---------------------------------------------------------------------------- #
#
#   genSTL.jl
#
#   Write STL from lattice
#
#   sterno
#   Spring 2018
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

import QHull

"""
    genSTL( mesh::MeshF,  lat::Lattice, flname::String; n = 8, name = "object" )

Generates a binary .stl file for the lattice in `mesh` with the areas defined in `lat`.
"""
function genSTL( mesh::MeshF,  lat::Lattice, flname::String; n = 8, name = "object" )

    # compute correct distance from node
    fdist = compDist( mesh, lat )
    fdist .*= 0.0

    # open STL file
    fid = open( flname, "w" )

    # generate vertices for nodes
    #   need to do that here because we need to know how many faces we have
    factot, nface = genFacesNods( mesh, lat, fdist, n, fid )

    nnzcyl = length( find(sqrt.( lat.ar / π ) .> 1e-14) ) # number of cylinders with nonzero radius
    nface = nnzcyl * (n-1) * 2

    # write header
    for jj in 1:80
        write(fid,' ')
    end
    write(fid,UInt32(nface))

    # generate STL for cylinders
    writefc = genSTLcyls( mesh, lat, fdist, n, fid, Val{2} )

    # write STL for nodes
    # writefn = genSTLBnods( factot, fid )

    # close STL
    close( fid )

end

"""
    genSTLA( mesh::MeshF,  lat::Lattice, flname::String; n = 8, name = "object" )

Generates an ASCII .stl file for the lattice in `mesh` with the areas defined in `lat`.
"""
function genSTLA( mesh::MeshF,  lat::Lattice, flname::String; n = 8, name = "object" )

    # 1. compute correct distance from node
    fdist = compDist( mesh, lat )

    # open STL file
    fid = open( flname, "w" )
    @printf( fid, "solid %s\n", name )

    # generate STL for cylinders
    genSTLcyls( mesh, lat, fdist, n, fid, Val{1} )

    # generate STL for nodes
    genSTLAnods( mesh, lat, fdist, n, fid )

    # close STL
    @printf( fid, "endsolid %s\n", name )
    close( fid )

end

"""
    compDist( mesh::MeshF,  lat::Lattice )

Computes the distance from the end face to the node
"""
function compDist( mesh::MeshF,  lat::Lattice )

    fdist = fill( 0.0, size(mesh.e,1), 2 )

    for kk in 1:mesh.n
        # loop over nodes

        # generate pairs of edges
        edg   = [i for i in 1:length(mesh.n2e[kk])]
        #   check for min area
        inddel = Vector{Int64}(0)
        for ii in 1:length(mesh.n2e[kk])
            if lat.ar[ mesh.n2e[kk][ii] ] < 0.0 || sqrt( lat.ar[ mesh.n2e[kk][ii] ] / π ) < 1.e-14
                append!(inddel,ii)
            end
        end
        edg0  = copy(edg)
        deleteat!(edg0,inddel)
        edg1  = edg0'
        eedg  =   edg0 .+ 0*edg1
        eedg1 = 0*edg0 .+   edg1
        pairs = hcat( eedg[:], eedg1[:] )

        upairs = unique(sort(pairs,2), 1) # unique pairs

        # compute required distance
        mindist = fill( 0.0, length(edg), length(edg)-1)
        indcnt = fill( 1, length(edg) )
        for jj in 1:size(upairs,1)
            edg_i1 = upairs[jj,1]
            edg_i2 = upairs[jj,2]
            if edg_i1 == edg_i2
                continue
            end
            e1 = mesh.n2e[kk][ edg[ edg_i1 ] ]
            e2 = mesh.n2e[kk][ edg[ edg_i2 ] ]

            # get different radii
            r1 = sqrt(lat.ar[ e1 ]/π)
            r2 = sqrt(lat.ar[ e2 ]/π)

            # compute angle between the edges
            nod1 = mesh.e[e1, find( mesh.e[e1,:] .!= kk )[] ]
            nod2 = mesh.e[e2, find( mesh.e[e2,:] .!= kk )[] ]
            vec1 = mesh.p[ nod1, : ] - mesh.p[ kk, : ]
            vec2 = mesh.p[ nod2, : ] - mesh.p[ kk, : ]
            arg = dot(vec1,vec2) / ( norm(vec1) * norm(vec2) )
            arg = minimum( ( 1.0, arg) ) # guard for round-off errors
            arg = maximum( (-1.0, arg) ) # guard for round-off errors
            θ = acos( arg )

            # compute min distance
            l2 = ( r1 + cos(θ)*r2 ) / sin(θ)
            l1 = sin(θ) * r2 + l2 * cos(θ)
            if θ > π - 1.e-6
                l1 = 0.0
                l2 = 0.0
            end

            mindist[ edg_i1, indcnt[edg_i1] ] = 1.001 * l1 # add 0.1% buffer
            mindist[ edg_i2, indcnt[edg_i2] ] = 1.001 * l2 # add 0.1% buffer
            indcnt[edg_i1] += 1
            indcnt[edg_i2] += 1

        end

        # save maximum distance
        for jj in 1:length(edg)
            e1 = mesh.n2e[kk][ edg[ jj ] ]
            ind1 = find( mesh.e[e1,:] .== kk )[]
            maxdist = maximum( mindist[jj,:] )
            if maxdist < 1e-14
                # if two edges are exactly coplanar, maxdist will be zero.
                # In that case, the convex hull computation will fail because both faces are the same
                maxdist = 0.01*maximum( sqrt.(abs.(lat.ar[ mesh.n2e[kk] ])/π) )
            end
            fdist[e1,ind1] = maxdist
        end

    end

    return fdist

end

"""
    genSTLcyls( mesh::MeshF,  lat::Lattice, fdist::Matrix{Float64}, n::Int64, fid::IOStream )

Writes the STL facets for all cylinders.
"""
function genSTLcyls( mesh::MeshF,  lat::Lattice, fdist::Matrix{Float64}, n::Int64, fid::IOStream, stltype::DataType )

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

        writeSTLcyl( stltype, xx1, yy1, zz1, xx2, yy2, zz2, fid )

    end

end

"""
    genFacesNods( mesh::MeshF, lat::Lattice, fdist::Matrix{Float64}, n::Int64, fid::IOStream )

Writes the STL facets for all nodes for ASCII files.
"""
function genFacesNods( mesh::MeshF, lat::Lattice, fdist::Matrix{Float64}, n::Int64, fid::IOStream )

    θ = linspace( 0, 2*π, n ).^1.0

    factot = Vector{Vector{Matrix{Float64}}}( mesh.n )
    nfac = 0

    for nn in 1:mesh.n

        edges = copy(mesh.n2e[nn])
        inddel = Vector{Int64}(0)
        for ii in 1:length(edges)
            if lat.ar[ edges[ii] ] < 0.0 || sqrt( lat.ar[ edges[ii] ] / π ) < 1e-14
                append!(inddel,ii)
            end
        end
        deleteat!(edges,inddel)

        (verts,nvecs) = compNodSingle( mesh, lat, edges, fdist, nn, n, θ )

        if length(edges) > 1
            # compute convex hull
            ch = try
                QHull.chull( verts )
            catch y
                println( "output ", y )
                println( "node ", nn )
                println( "edges ", edges )
                println( "distances ", fdist[edges] )
                error("Problems!!")
            end

            faces = ch.simplices

            centroid = mean( verts, 1 )[:] # need to compute centroid to determine whether or not normal is correct

            ifaces = cleanupCHull( faces, verts, nvecs, centroid )

            # save face vertices
            factot[nn] = Vector{Matrix{Float64}}( length(ifaces) )
            kk = 1
            for ff in ifaces
                factot[nn][kk] = verts[ faces[ff],: ]
                kk   += 1
                nfac += 1
            end

        end

    end

    return (factot, nfac)

end

"""
    genSTLBnods( factot::Vector{Vector{Matrix{Float64}}}, fid::IOStream )

Writes the STL facets for all nodes for binary files.
"""
function genSTLBnods( factot::Vector{Vector{Matrix{Float64}}}, fid::IOStream )

    for f1 in factot, verts in f1
        writeFacetSTLB( verts, fid::IOStream )
    end

end


"""
    genSTLAnods( mesh::MeshF, lat::Lattice, fdist::Matrix{Float64}, n::Int64, fid::IOStream )

Writes the STL facets for all nodes for ASCII files.
"""
function genSTLAnods( mesh::MeshF, lat::Lattice, fdist::Matrix{Float64}, n::Int64, fid::IOStream )

    θ = linspace( 0, 2*π, n ).^1.0

    for nn in 1:mesh.n

        edges = copy(mesh.n2e[nn])
        inddel = Vector{Int64}(0)
        for ii in 1:length(edges)
            if lat.ar[ edges[ii] ] < 0.0 || sqrt( lat.ar[ edges[ii] ] / π ) < 1e-14
                append!(inddel,ii)
            end
        end
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
                writeFacetSTLA( verts[ faces[ff],: ], fid::IOStream )
            end

        end

    end

end

function compNodSingle( mesh::MeshF, lat::Lattice, edges::Vector{Int64},
                        fdist::Matrix{Float64}, nn::Int64, n::Int64, θ::Vector{Float64} )

    # find vertices for convex hull
    verts = Matrix{Float64}( length(edges)*n, 3 )
    nvecs = Vector{Vector{Float64}}( length(edges) )
    cnt = 1
    cntnv = 1

    for ee in edges

        inod = find(mesh.e[ee,:] .== nn )[]

        nod1 = mesh.e[ee,1]
        nod2 = mesh.e[ee,2]

        vec = mesh.p[nod2,:] - mesh.p[nod1,:]

        l = norm( vec )

        # begin and start nodes of cylinder
        xface = mesh.p[ nn, : ]
        if inod == 1
            xface += fdist[ee,1] / l * vec
        else
            xface -= fdist[ee,2] / l * vec
        end

        # generate points on end-faces
        r = sqrt( lat.ar[ee] / π )
        xx = r * cos.( θ )
        yy = r * sin.( θ )

        nvec = vec / l

        xx, yy, zz = rotTransCirc( xx, yy, nvec )
        #   translate from origin to end points
        xx += xface[1]
        yy += xface[2]
        zz += xface[3]

        verts[cnt:(cnt+n-1),1] = xx
        verts[cnt:(cnt+n-1),2] = yy
        verts[cnt:(cnt+n-1),3] = zz

        cnt += n

        nvecs[cntnv] = nvec
        if inod == 2
            nvecs[cntnv] = -nvec
        end
        cntnv += 1

    end

    # is this a sharp node?
    sharp = true
    nvecsum = [0.0,0.0,0.0]
    for el in 1:length(edges)
        nvecsum += nvecs[el]
    end
    for el in 1:length(edges)
        dotp = dot( nvecsum, nvecs[el] )
        if dotp <= 0.0
            sharp = false
            break
        end
    end
    if all(abs.(nvecsum) .< 1e-14)
        sharp = false
    end

    if sharp
        # NOTE: NOT SURE IF THIS WILL AT SOME POINT COLLIDE WITH REST OF VERTICES
        ρ = 0.5 # factor of radius
        nnvecsum = nvecsum/norm(nvecsum)

        rmax = sqrt( maximum( lat.ar[edges] ) / π )
        r = ρ * rmax
        xx = r * cos.( θ )
        yy = r * sin.( θ )

        xx, yy, zz = rotTransCirc( xx, yy, -nnvecsum )
        #   translate from origin to end points
        offset = sqrt( rmax^2 - r^2 )
        xface = mesh.p[ nn, : ] - offset * nnvecsum
        xx += xface[1]
        yy += xface[2]
        zz += xface[3]

        verts = [verts; hcat(xx, yy, zz)]
    end

    return (verts,nvecs)

end

"""
    writeSTLcyl( xx1::Vector{Float64}, yy1::Vector{Float64}, zz1::Vector{Float64},
                 xx2::Vector{Float64}, yy2::Vector{Float64}, zz2::Vector{Float64},
                 fid::IOStream )

Writes the STL facets for one cylinder for an ASCII STL file.
"""
function writeSTLcyl( ::Type{Val{1}}, xx1::Vector{Float64}, yy1::Vector{Float64}, zz1::Vector{Float64},
                                      xx2::Vector{Float64}, yy2::Vector{Float64}, zz2::Vector{Float64},
                                      fid::IOStream )

    for kk in 1:size(xx1,1)-1

        # write facet 1
        vert = [ xx1[kk]   yy1[kk]   zz1[kk];
                 xx2[kk+1] yy2[kk+1] zz2[kk+1];
                 xx2[kk]   yy2[kk]   zz2[kk]]
        writeFacetSTLA( vert, fid )

        # write facet 2
        vert = [ xx1[kk+1]  yy1[kk+1] zz1[kk+1];
                  xx2[kk+1] yy2[kk+1] zz2[kk+1];
                  xx1[kk]   yy1[kk]   zz1[kk] ]
        writeFacetSTLA( vert, fid )

    end

end

"""
    writeSTLcyl( xx1::Vector{Float64}, yy1::Vector{Float64}, zz1::Vector{Float64},
                 xx2::Vector{Float64}, yy2::Vector{Float64}, zz2::Vector{Float64},
                 fid::IOStream )

Writes the STL facets for one cylinder for a binary STL file.
"""
function writeSTLcyl( ::Type{Val{2}}, xx1::Vector{Float64}, yy1::Vector{Float64}, zz1::Vector{Float64},
                                      xx2::Vector{Float64}, yy2::Vector{Float64}, zz2::Vector{Float64},
                                      fid::IOStream )

    for kk in 1:size(xx1,1) - 1

        # write facet 1
        vert = [ xx1[kk]   yy1[kk]   zz1[kk];
                 xx2[kk+1] yy2[kk+1] zz2[kk+1];
                 xx2[kk]   yy2[kk]   zz2[kk]]
        writeFacetSTLB( vert, fid )

        # write facet 2
        vert = [ xx1[kk+1]  yy1[kk+1] zz1[kk+1];
                  xx2[kk+1] yy2[kk+1] zz2[kk+1];
                  xx1[kk]   yy1[kk]   zz1[kk] ]
        writeFacetSTLB( vert, fid )

    end

end

"""
    cleanupCHull( faces::Vector{Vector{Int64}},   verts::Matrix{Float64},
                  nvecs::Vector{Vector{Float64}}, centroid::Vector{Float64} )

Flips the normals outward for the convex hulls and deletes faces that are
coincident with the end-faces.
"""
function cleanupCHull( faces::Vector{Vector{Int64}},   verts::Matrix{Float64},
                       nvecs::Vector{Vector{Float64}}, centroid::Vector{Float64} )

    ifaces = [i for i in 1:length(faces)]; cnt = 1

    for ff in 1:length(faces)

        # determine if normal points inward or outward
        cent   = mean( verts[ faces[ff], : ], 1 )[:] # face center
        veccnt = centroid - cent

        vec1   = verts[ faces[ff][2], : ] - verts[ faces[ff][1], : ]
        vec2   = verts[ faces[ff][3], : ] - verts[ faces[ff][1], : ]
        vec1   = cross( vec1, vec2 )
        normal = vec1 / norm(vec1)

        if dot(normal,veccnt) > 0.0
            ifac = faces[ff][3]
            faces[ff][3] = faces[ff][2]
            faces[ff][2] = ifac
            normal .*= -1.0
        end

        # TODO: check faces that align with cylinders

        cnt += 1

        for jj in 1:length(nvecs)
            if dot( normal, nvecs[jj] ) > 1.0 - 1e-14
                deleteat!(ifaces,cnt-1)
                cnt -= 1
                break
            end

        end
    end

    return ifaces

end

"""
    writeFacetSTLA( vert::Matrix{Float64}, fid::IOStream )

Writes one STL facet for ASCII file.
"""
function writeFacetSTLA( vert::Matrix{Float64}, fid::IOStream )

    vec1   = vert[2,:] - vert[1,:]
    vec2   = vert[3,:] - vert[1,:]
    vec1   = cross( vec1, vec2 )
    normal = vec1 / norm(vec1)

    @printf( fid, "facet normal %11.7E %11.7E %11.7E\n", normal[1], normal[2], normal[3] )
    @printf( fid, " outer loop\n" )
    @printf( fid, "  vertex %11.7E %11.7E %11.7E\n", vert[1,1], vert[1,2], vert[1,3] )
    @printf( fid, "  vertex %11.7E %11.7E %11.7E\n", vert[2,1], vert[2,2], vert[2,3] )
    @printf( fid, "  vertex %11.7E %11.7E %11.7E\n", vert[3,1], vert[3,2], vert[3,3] )
    @printf( fid, " endloop\n" )
    @printf( fid, "endfacet\n" )

end

"""
    writeFacetSTLB( vert::Matrix{Float64}, fid::IOStream )

Writes one STL facet for binary file.
"""
function writeFacetSTLB( vert::Matrix{Float64}, fid::IOStream )

    vec1   = vert[2,:] - vert[1,:]
    vec2   = vert[3,:] - vert[1,:]
    vec1   = cross( vec1, vec2 )
    normal = vec1 / norm(vec1)

    # normals
    write( fid, Float32( normal[1] ) )
    write( fid, Float32( normal[2] ) )
    write( fid, Float32( normal[3] ) )

    # Vertex 1
    write( fid, Float32( vert[1,1] ) )
    write( fid, Float32( vert[1,2] ) )
    write( fid, Float32( vert[1,3] ) )

    # Vertex 2
    write( fid, Float32( vert[2,1] ) )
    write( fid, Float32( vert[2,2] ) )
    write( fid, Float32( vert[2,3] ) )

    # Vertex 3
    write( fid, Float32( vert[3,1] ) )
    write( fid, Float32( vert[3,2] ) )
    write( fid, Float32( vert[3,3] ) )

    # End facet
    write( fid, UInt16(0) )

end

"""
    rotTransCirc( xx::Vector{Float64}, yy::Vector{Float64}, b::Vector{Float64} )

Rotates a circle defined in `xx` and `yy` such that it is normal to the vector `b`.
"""
function rotTransCirc( xx::Vector{Float64}, yy::Vector{Float64}, b::Vector{Float64} )
    # https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d

    a = [0.0, 0.0, 1.0]

    c = dot( a, b )

    tol = 1e-15

    if c > (-1.0-tol) && c < (-1.0+tol)
        # don't rotate, because it is just mirrored, but that should not matter
        res = hcat( xx, yy, 0.0*yy )'
    else

        v = cross( a, b )
        s = norm( v )


        vx = fill( 0.0, 3, 3 )
        vx[1,2] = -v[3]
        vx[1,3] =  v[2]
        vx[2,1] =  v[3]
        vx[2,3] = -v[1]
        vx[3,1] = -v[2]
        vx[3,2] =  v[1]

        Q = eye(3) + vx + vx * vx * 1 / (1 + c)

        res = Q * hcat( xx, yy, 0.0*yy )'

    end

    return ( res[1,:], res[2,:], res[3,:] )

end
