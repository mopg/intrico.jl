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

"""
    genSTL( mesh::MeshF,  lat::Lattice, flname::String; n = 8, name = "object" )

Generates a binary .stl file for the lattice in `mesh` with the areas defined in `lat`.
"""
function genSTLnew( mesh::MeshF,  lat::Lattice, flname::String; n = 8, name = "object" )

    # # # compute correct distance from node
    # fdist = compDist( mesh, lat )

    # open STL file
    fid = open( flname, "w" )

    # write header
    for jj in 1:80
        write(fid,' ')
    end
    write(fid,UInt32(0)) # number of faces is unknown now, so write 0

    # generate all faces (and write)
    factot = genFacesNodsCyl( mesh, lat, n, fid )

    # close STL
    close( fid )

    ## reopen to write the number of faces
    fid = open( flname, "a+" )
    seekstart(fid) # go back to top of file

    # rewrite header
    for jj in 1:80
        write(fid,' ')
    end
    # write correct number of faces
    write(fid,UInt32(factot))
    # close STL
    close( fid )
end

"""
    genSTLA( mesh::MeshF,  lat::Lattice, flname::String; n = 8, name = "object" )

Generates an ASCII .stl file for the lattice in `mesh` with the areas defined in `lat`.
"""
function genSTLAnew( mesh::MeshF,  lat::Lattice, flname::String; n = 8, name = "object" )

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
    genFacesNods( mesh::MeshF, lat::Lattice, fdist::Matrix{Float64}, n::Int64, fid::IOStream )

Writes the STL facets for all nodes for ASCII files.
"""
function genFacesNodsCyl( mesh::MeshF, lat::Lattice, n::Int64, fid::IOStream )

    nfac = 0

    ptsCylEdge = Vector{Vector{Vector{Vector{Float64}}}}( size(mesh.e,1) )
    for kk in 1:size(mesh.e,1)
        ptsCylEdge[kk] = Vector{Vector{Vector{Float64}}}(2)
    end

    # generate nodes
    for nn in 1:mesh.n

        nfac += genFacesNod( nn, mesh, lat, n, ptsCylEdge, fid )

    end

    # generate cylinders
    for kk in 1:size(mesh.e,1)
        nfac += genFacesCyl( kk, mesh, ptsCylEdge[kk], fid )
    end

    return nfac

end

"""
    compDist( mesh::MeshF,  lat::Lattice )

Computes the distance from the end face to the node
"""
function genFacesNod( kk::Int64, mesh::MeshF,  lat::Lattice, nn::Int64,
                      ptsCylEdge::Vector{Vector{Vector{Vector{Float64}}}},
                      fid::IOStream )

    nfac = 0

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
    edg1  =   edg0'
    eedg  =   edg0 .+ 0*edg1
    eedg1 = 0*edg0 .+   edg1
    pairs = hcat( eedg[:], eedg1[:] )

    upairs = unique(sort(pairs,2), 1) # unique pairs

    rmax = sqrt( maximum( lat.ar[edg] ) / π )

    # compute crossing point and normal vector
    intersec = Vector{Vector{Matrix{Float64}}}( length(edg) )
    normvec  = Vector{Vector{Vector{Float64}}}( length(edg) )
    tangvec  = Vector{Vector{Vector{Float64}}}( length(edg) )
    indedges = Vector{Vector{Int64}}( length(edg) )
    lv = Vector{Vector{Float64}}( length(edg) )
    for jj in 1:length(edg)
        intersec[jj] = Vector{Matrix{Float64}}( length(edg)-1 )
        normvec[jj]  = Vector{Vector{Float64}}( length(edg)-1 )
        tangvec[jj]  = Vector{Vector{Float64}}( length(edg)-1 )
        indedges[jj] = Vector{Int64}( length(edg)-1 )
        lv[jj] = Vector{Float64}( length(edg)-1 )
    end
    indcnt = fill(1, length(edg) )
    println("size upairs ", size(upairs))
    println("length edg ", length(edg))
    for jj in 1:size(upairs,1)
        edg_i1 = upairs[jj,1]
        edg_i2 = upairs[jj,2]
        if edg_i1 == edg_i2
            continue
        end
        e1 = mesh.n2e[kk][ edg[ edg_i1 ] ]
        e2 = mesh.n2e[kk][ edg[ edg_i2 ] ]

        # get different radii
        # r1 = sqrt(lat.ar[ e1 ]/π)
        # r2 = sqrt(lat.ar[ e2 ]/π)

        # compute angle between the edges
        nod1 = mesh.e[e1, find( mesh.e[e1,:] .!= kk )[] ]
        nod2 = mesh.e[e2, find( mesh.e[e2,:] .!= kk )[] ]
        vec1 = mesh.p[ nod1, : ] - mesh.p[ kk, : ]
        vec2 = mesh.p[ nod2, : ] - mesh.p[ kk, : ]
        vec1 ./= norm(vec1)
        vec2 ./= norm(vec2)
        arg = dot(vec1,vec2)
        arg = minimum( ( 1.0, arg) ) # guard for round-off errors
        arg = maximum( (-1.0, arg) ) # guard for round-off errors
        θ = acos( arg )

        println("e1 ", e1)
        println("e2 ", e2)
        println("vec1 ", vec1)
        println("vec2 ", vec2)
        println("θ ", θ)
        println("rmax ", rmax)

        # compute absolute distance between node center and intersection
        lv[edg_i1][indcnt[edg_i1]] = rmax / tan( θ/2 ) # offset such that rods do not collide
        lv[edg_i2][indcnt[edg_i2]] = rmax / tan( θ/2 ) # offset such that rods do not collide
        l = sqrt( lv[edg_i1][indcnt[edg_i1]]^2 + rmax^2 ) # absolute distance

        magdiff = norm(vec1 - vec2)
        magave  = norm(vec1 + vec2)

        if magave < 1e-14
            # two rods are exactly opposite

            a = [0.0, 0.0, 1.0]

            c = dot( a, vec1 )

            tol = 1e-15

            if c > (-1.0-tol) && c < (-1.0+tol)
                # don't rotate, because it is just mirrored, but that should not matter
                crossvec1 = [-l, 0.0, 0.0]
                crossvec2 = [ l, 0.0, 0.0]
            else
                v = cross( a, vec1 )
                s = norm( v )

                vx = fill( 0.0, 3, 3 )
                vx[1,2] = -v[3]
                vx[1,3] =  v[2]
                vx[2,1] =  v[3]
                vx[2,3] = -v[1]
                vx[3,1] = -v[2]
                vx[3,2] =  v[1]

                Q = eye(3) + vx + vx * vx * 1 / (1 + c)

                crossvec1 = Q * [-l, 0.0, 0.0]
                crossvec2 = Q * [ l, 0.0, 0.0]
            end
            # tangent vector
            tangvec[edg_i1][ indcnt[edg_i1] ] = crossvec1 / norm(crossvec1)
            tangvec[edg_i2][ indcnt[edg_i2] ] = crossvec2 / norm(crossvec2)
        else
            crossvec1   =  (vec1+vec2) / magave * l
            crossvec2   = -(vec1+vec2) / magave * l

            # tangent vector
            tangvec[edg_i1][ indcnt[edg_i1] ] = (vec1+vec2) / magdiff
            tangvec[edg_i2][ indcnt[edg_i2] ] = (vec1+vec2) / magdiff
        end

        crossvec1  += mesh.p[ kk, : ] # intersection location
        crossvec2  += mesh.p[ kk, : ] # intersection location

        println("crossvec1 ", crossvec1)
        println("crossvec2 ", crossvec2)
        println("normvec ", (vec1-vec2) / magdiff )
        println(" ")

        intersec[edg_i1][ indcnt[edg_i1] ] = hcat( crossvec1, crossvec2 )
        indedges[edg_i1][ indcnt[edg_i1] ] = edg_i2

        intersec[edg_i2][ indcnt[edg_i2] ] = hcat( crossvec1, crossvec2 )
        indedges[edg_i2][ indcnt[edg_i2] ] = edg_i1

        # normal vectors
        normvec[edg_i1][ indcnt[edg_i1] ] = (vec1-vec2) / magdiff
        normvec[edg_i2][ indcnt[edg_i2] ] = (vec2-vec1) / magdiff

        indcnt[edg_i1] += 1
        indcnt[edg_i2] += 1

    end

    # find maximum distance of rods
    lvmax = fill( 0.0, length(lv) )
    for ii in 1:length(lvmax)
        lvmax[ii] = maximum( lv[ii] ) +
            0.05 * norm( mesh.p[ mesh.e[mesh.n2e[kk][edg[ii]],1],: ] - mesh.p[ mesh.e[mesh.n2e[kk][edg[ii]],2],: ] )
    end
    println("lvmax ", lvmax)
    println("edg ", mesh.n2e[kk][edg] )

    # loop over each edge
    for ee in 1:length(edg)

        ## check which crossing points are active
        e1     = mesh.n2e[kk][ edg[ ee ] ]
        nod    = mesh.e[e1, find( mesh.e[e1,:] .!= kk )[] ]
        evec   = mesh.p[ nod, : ] - mesh.p[ kk, : ]
        testpt = mesh.p[ nod, : ] # distance to this point is measured to determine which is closer
        evec ./= norm(evec)

        actpts  = Vector{Vector{Float64}}( 2*length(intersec[ee]) )
        normact = Vector{Vector{Float64}}( 2*length(intersec[ee]) )
        nact = 0
        # println("intersec ",intersec[ee])
        println( "edge ", e1, " for node ", kk )
        # point on convex side
        for jj in 1:length(intersec[ee])#, kk in 1:2
            # loop over intersection points
            currpt  = intersec[ee][jj][:,1]
            println("currpt ", currpt)
            dist    =  norm( currpt - testpt )
            distorg = dist
            for ii in 1:length(edg)-1
                # loop over intersecting planes
                d = dot( mesh.p[kk,:] - currpt, normvec[ee][ii] ) / dot( evec, normvec[ee][ii] ) # scalar value for which line intersects plane
                ipt = currpt + d * evec

                tvec = tangvec[ee][ii] # vector along intersecting plane in plane between two edges

                cdist = norm( ipt - testpt )
                # println("tvec ", tvec)
                println("ipt ", ipt)
                println("cdist ", cdist)
                # if dot(tvec, ipt - mesh.p[kk,:]) < 0.0 # on concave side
                #     if cdist > dist
                #         # this point is further
                #         dist = cdist
                #     end
                # else # on convex side
                    if cdist < dist
                        # this point is closer
                        dist = cdist
                    end
                # end
            end
            println("fin dist ", dist)
            println("distorg  ", distorg)
            println(" ")
            # is point active? Only if new distance is same as original distance
            if abs(dist - distorg) < 1e-14
                nact += 1
                actpts[nact]  = currpt
                normact[nact] = copy(normvec[ee][jj])
            end
        end

        ## order points counterclockwise
        # pick first point as ϕ = 0
        println("actpts ",actpts)
        currvec   = actpts[1] - mesh.p[kk,:]
        zerovec   = currvec - dot( currvec, evec ) * evec
        zerovec ./= norm(zerovec)
        ϕ = fill(0.0,nact)
        for jj in 2:nact
            currvec  = actpts[jj] - mesh.p[kk,:]
            currnvec = currvec - dot( currvec, evec ) * evec
            # println( "currnvec ", currnvec )
            # println( "arg ", dot(zerovec,currnvec) / norm(currnvec) )
            ϕ[jj] = acos( max( min( dot(zerovec,currnvec) / norm(currnvec), 1.0), -1.0 ) )
            if dot( cross( zerovec, currnvec ), evec ) < 0.0
                ϕ[jj] = 2*π - ϕ[jj]
            end
        end

        # sort based on ϕ
        angleind = sortperm( ϕ )
        actpts[1:nact]  = actpts[ angleind ]
        normact[1:nact] = normact[ angleind ]

        # add last point to array to close the loop
        actpts[nact+1]  = copy(actpts[1])
        normact[nact+1] = copy(normact[1])
        append!( ϕ, 2*π )
        nact += 1

        # TEMP print everything
        println("active pts for ", kk, " ", e1)
        for jj in 1:nact
            println(jj, " ", actpts[jj])
        end
        println(" ")

        ## Generate points along cuts
        perpvec = cross(evec,zerovec)
        perpvec ./= norm(perpvec)
        # println("normact ", normact)
        currvec ./= norm(currvec)

        iedge = find( mesh.e[e1,:] .== kk )[]
        ptsCylEdge[e1][iedge] = Vector{Vector{Float64}}(0)
        for jj in 1:nact-1
            # compute the intersection between the planes
            lvec  = cross( normact[jj], normact[jj+1] ) # always going through the node center (because both planes go through that point)
            lvec -= dot( lvec, evec )
            if norm(lvec) < 1e-14 # vector is crossed with itself (i.e. there is only one plane)
                ϕc1 = π
            else
                ϕc1 = acos( max( min(dot(zerovec,lvec) / norm(lvec), 1.0), -1.0) )
                if dot( cross( zerovec, lvec ), evec ) < 0.0
                    ϕc1 = 2*π - ϕc1
                end
            end
            println("ϕ beg ", ϕ[jj])
            println("ϕ end ", ϕ[jj+1])
            println("ϕc1 ", ϕc1)
            if ϕc1 > ϕ[jj+1]
                ϕc1 -= π
            elseif ϕc1 < ϕ[jj]
                ϕc1 += π
            end
            println("ϕc1 (new) ", ϕc1)
            Δϕ1 = ϕc1 - ϕ[jj]
            np1 = max( convert( Int64, ceil( Δϕ1/(2*π) * nn ) + 1 ), 2 ) # ensure we have at least two points
            Δϕ2 = ϕ[jj+1] - ϕc1
            np2 = max( convert( Int64, ceil( Δϕ2/(2*π) * nn )+1 ), 2 ) # ensure we have at least two points

            ϕcurr1 = linspace( ϕ[jj], ϕc1,     np1 )
            ϕcurr2 = linspace( ϕc1,   ϕ[jj+1], np2 )

            # first part
            intpts = Vector{Vector{Float64}}( np1 )
            cylpts = Vector{Vector{Float64}}( np1 )

            for ii in 1:np1
                currpt = rmax * ( zerovec * cos(ϕcurr1[ii]) + perpvec * sin(ϕcurr1[ii]) ) # current point in plane perpendicular to evec
                # find actual point near node
                d = dot( -currpt, normact[jj] ) / dot( evec, normact[jj] ) # scalar value for which line intersects plane
                intpts[ii] = currpt + d * evec + mesh.p[kk,:]
                cylpts[ii] = sqrt(lat.ar[e1]/π) * ( zerovec * cos(ϕcurr1[ii]) + perpvec * sin(ϕcurr1[ii]) ) + mesh.p[kk,:] + evec * lvmax[ee]

            end
            # write to file
            for ii in 1:np1-1
                vert = [ cylpts[ii]';
                         intpts[ii]';
                         intpts[ii+1]' ]
                writeFacetSTLB( vert, fid )
                nfac += 1

                vert = [ intpts[ii+1]';
                         cylpts[ii+1]';
                         cylpts[ii]' ]
                writeFacetSTLB( vert, fid )
                nfac += 1
            end
            append!( ptsCylEdge[e1][iedge], cylpts[1:end-1] )
            if jj == nact-1
                append!( ptsCylEdge[e1][iedge], cylpts[end:end] )
            end

            # second part
            intpts = Vector{Vector{Float64}}( np2 )
            cylpts = Vector{Vector{Float64}}( np2 )

            for ii in 1:np2
                currpt = rmax * ( zerovec * cos(ϕcurr2[ii]) + perpvec * sin(ϕcurr2[ii]) ) # current point in plane perpendicular to evec

                # find actual point near node
                d = dot( -currpt, normact[jj+1] ) / dot( evec, normact[jj+1] ) # scalar value for which line intersects plane
                intpts[ii] = currpt + d * evec + mesh.p[kk,:]
                cylpts[ii] = sqrt(lat.ar[e1]/π) * ( zerovec * cos(ϕcurr2[ii]) + perpvec * sin(ϕcurr2[ii]) ) + mesh.p[kk,:] + evec * lvmax[ee]

            end
            # write to file
            for ii in 1:np2-1
                vert = [ cylpts[ii]';
                         intpts[ii]';
                         intpts[ii+1]' ]
                writeFacetSTLB( vert, fid )
                nfac += 1

                vert = [ intpts[ii+1]';
                         cylpts[ii+1]';
                         cylpts[ii]' ]
                writeFacetSTLB( vert, fid )
                nfac += 1
            end
            append!( ptsCylEdge[e1][iedge], cylpts[1:end-1] )
            if jj == nact-1
                append!( ptsCylEdge[e1][iedge], cylpts[end:end] )
            end
        end

    end

    return nfac

end

"""
    compDist( mesh::MeshF,  lat::Lattice )

Computes the distance from the end face to the node
"""
function genFacesCyl( kk::Int64, mesh::MeshF,
                      ptsCylEdge::Vector{Vector{Vector{Float64}}},
                      fid::IOStream )

    nfac = 0
    frac = 0.25

    e1 = mesh.e[kk,1]
    e2 = mesh.e[kk,2]

    ## first edge
    # loop over first edge
    println("kk     ", kk)
    println("node 1 ", e1)
    println("node 2 ", e2)
    println("edge points 1 ", length( ptsCylEdge[1] ))
    println("edge points 2 ", length( ptsCylEdge[2] ))

    n1 = length( ptsCylEdge[1] )
    n2 = length( ptsCylEdge[2] )

    # find which side has most points
    i1 = 1
    i2 = 2
    if n2 < n1
        i2 = 1
        i1 = 2
        n1 = (n1 + n2)
        n2 = n1 - n2
        n1 = n1 - n2
    end

    ## generate the all facets, starting with the side with the fewest points
    # find first index and generate first facet
    ind0 = findMinDistInd( ptsCylEdge, i1, i2, 1, frac )
    vert = [ ptsCylEdge[i1][1]';
             ptsCylEdge[i1][2]';
             ptsCylEdge[i2][ind0]' ]
    writeFacetSTLB( vert, fid )
    ind00 = ind0 # need to save this index to close loop
    nfac += 1
    # generate rest of facets (also for other ring)
    for jj in 2:n1-1

        ind1 = findMinDistInd( ptsCylEdge, i1, i2, jj, frac )

        vert = [ ptsCylEdge[i1][jj]';
                 ptsCylEdge[i1][jj+1]';
                 ptsCylEdge[i2][ind1]' ]

        writeFacetSTLB( vert, fid )
        nfac += 1

        # generate other faces
        # NOTE: due to ordering ind1 < ind0
        if ind0 > ind1
            for kk = 0:(ind0-ind1-1)
                vert = [ ptsCylEdge[i2][ind1+kk]';
                         ptsCylEdge[i2][ind1+kk+1]';
                         ptsCylEdge[i1][jj]' ]

                writeFacetSTLB( vert, fid )
                nfac += 1
            end
        elseif ind0 < ind1
            inds   = vcat(ind1:n2, 1:ind0)
            indsp1 = vcat(ind1+1:n2, 1:ind0+1)
            for kk in 1:length(inds)
                vert = [ ptsCylEdge[i2][ inds[kk] ]';
                         ptsCylEdge[i2][ indsp1[kk] ]';
                         ptsCylEdge[i1][ jj ]' ]

                writeFacetSTLB( vert, fid )
                nfac += 1
            end
        end

        ind0 = ind1

    end
    # last patch
    ind1 = ind00
    for kk = 0:(ind0-ind1-1)
        vert = [ ptsCylEdge[i2][ind1+kk]';
                 ptsCylEdge[i2][ind1+kk+1]';
                 ptsCylEdge[i1][1]' ]

        writeFacetSTLB( vert, fid )
        nfac += 1
    end

    return nfac

end

function findMinDistInd( pts::Vector{Vector{Vector{Float64}}},
                         ind1::Int64, ind2::Int64,
                         jj::Int64, frac::Float64)

    testpt = (1-frac) * pts[ind1][jj] + frac * pts[ind1][jj+1] # check for minimum distance to middle of edge

    # find minimum distance to node on other side
    dist = 1e4
    ind  = 0 # look for first index, but after that only allow for + or - one (because should be ordered)
    for kk in 1:length( pts[ind2] )
        cdist = norm( pts[ind2][kk] - testpt )
        if cdist < dist
            dist = cdist
            ind = kk
        end
    end

    return ind
end
