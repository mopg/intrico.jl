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

    pairs = unique(sort(pairs,2), 1) # unique pairs
    pairs = pairs[ sortperm( pairs[:,1] + pairs[:,2]*100 ), : ] # ensure highest indices are last
    # TODO: don't need to do this for each loop, can just precompute this for n = 10 and then just take part of the unique pairs
    rmax = sqrt( maximum( lat.ar[ mesh.n2e[kk] ] ) / π )

    # compute crossing point and normal vector
    normvec  = Vector{Vector{Vector{Float64}}}( length(edg) )
    tangvec  = Vector{Vector{Vector{Float64}}}( length(edg) )
    lvmax = fill( 0.0, length(edg) )
    for jj in 1:length(edg)
        normvec[jj]  = Vector{Vector{Float64}}( length(edg)-1 )
        tangvec[jj]  = Vector{Vector{Float64}}( length(edg)-1 )
    end
    indcnt = fill(1, length(edg) )
    for jj in 1:size(pairs,1)
        edg_i1 = pairs[jj,1]
        edg_i2 = pairs[jj,2]
        if edg_i1 == edg_i2
            continue
        end
        e1 = mesh.n2e[kk][ edg[ edg_i1 ] ]
        e2 = mesh.n2e[kk][ edg[ edg_i2 ] ]

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

        # compute absolute distance between node center and intersection
        lv = rmax / tan( θ/2 ) # offset such that rods do not collide
        l = sqrt( lv^2 + rmax^2 ) # absolute distance
        lvmax[edg_i1] = max( lvmax[edg_i1], lv )
        lvmax[edg_i2] = max( lvmax[edg_i2], lv )

        magdiff = norm(vec1 - vec2)
        magave  = norm(vec1 + vec2)

        if magave < 1e-14
            # two rods are exactly opposite

            # tangent vector
            tangvec[edg_i1][ indcnt[edg_i1] ] = crossvec1 / norm(crossvec1)
            tangvec[edg_i2][ indcnt[edg_i2] ] = crossvec2 / norm(crossvec2)
        else
            # tangent vector
            tangvec[edg_i1][ indcnt[edg_i1] ] = (vec1+vec2) / magdiff
            tangvec[edg_i2][ indcnt[edg_i2] ] = (vec1+vec2) / magdiff
        end

        # normal vectors
        normvec[edg_i1][ indcnt[edg_i1] ] = (vec1-vec2) / magdiff
        normvec[edg_i2][ indcnt[edg_i2] ] = (vec2-vec1) / magdiff

        indcnt[edg_i1] += 1
        indcnt[edg_i2] += 1

    end

    # add offset to rods
    for ii in 1:length(lvmax)
        lvmax[ii] += 0.05 * norm( mesh.p[ mesh.e[mesh.n2e[kk][edg[ii]],1],: ] - mesh.p[ mesh.e[mesh.n2e[kk][edg[ii]],2],: ] )
    end

    # loop over each edge
    for ee in 1:length(edg)

        e1     = mesh.n2e[kk][ edg[ ee ] ]
        nod    = mesh.e[e1, find( mesh.e[e1,:] .!= kk )[] ]
        evec   = mesh.p[ nod, : ] - mesh.p[ kk, : ]
        evec ./= norm(evec)

        # find intersections between planes and check if they are active
        testpt   = mesh.p[ nod, : ] - mesh.p[ kk, : ] # distance to this point is measured to determine which is closer, distance is relative to current node
        nunique  = uniqueNumPairs( length(edg)-1 )
        intersec = Vector{Vector{Float64}}( 2*(nunique - (length(edg)-1)) + 1)
        nact = 0
        for jj in 1:nunique # note: number of planes is length(edg)-1, so can reuse part of upairs

            if pairs[jj,1] == pairs[jj,2]
                continue
            end
            ip1 = pairs[jj,1]
            ip2 = pairs[jj,2]

            lvec  = cross( normvec[ee][ ip1 ], normvec[ee][ ip2 ] ) # always going through the node center (because both planes go through that point)
            nrmperp  = norm( lvec - dot( lvec, evec )*evec )

            # check first point
            lvec .*=  rmax / nrmperp
            # check ...

            act = checkPointAct( lvec, testpt, normvec[ee], evec )

            if act
                nact += 1
                intersec[nact] = copy(lvec)
            end
            # println("act ", act)

            # check second point
            lvec .*= -1.0
            # check ...
            act = checkPointAct( lvec, testpt, normvec[ee], evec )
            if act
                nact += 1
                intersec[nact] = copy(lvec)
            end

        end

        ## order points counterclockwise
        # pick first point as ϕ = 0
        zerovec   = intersec[1] - dot( intersec[1], evec ) * evec
        zerovec ./= norm(zerovec)
        # println( "zerovec ", zerovec )
        ϕ = fill(0.0,nact)
        for jj in 2:nact
            currnvec = intersec[jj] - dot( intersec[jj], evec ) * evec
            ϕ[jj] = acos( max( min( dot(zerovec,currnvec) / norm(currnvec), 1.0), -1.0 ) )
            if dot( cross( zerovec, currnvec ), evec ) < 0.0
                ϕ[jj] = 2*π - ϕ[jj]
            end
        end

        # sort based on ϕ
        angleind   = sortperm( ϕ )
        ϕ[1:nact]  = ϕ[angleind]
        corrind    = fill( 0, nact )
        corrind[1] = angleind[1]
        nact = 1
        for jj in 1:length(ϕ)-1
            if (ϕ[jj+1] - ϕ[jj]) > 1e-14
                nact += 1
                corrind[nact] = jj+1
            end
        end
        ϕ = ϕ[ corrind[1:nact] ]

        # add last point to array to close the loop
        append!( ϕ, 2*π )
        nact += 1

        ## Generate points along cuts
        perpvec = cross(evec,zerovec)
        perpvec ./= norm(perpvec)

        iedge = find( mesh.e[e1,:] .== kk )[]
        ptsCylEdge[e1][iedge] = Vector{Vector{Float64}}(0)
        for jj in 1:nact-1

            # figure out on which plane to project
            ϕhalf   = ( ϕ[jj] + ϕ[jj+1] ) / 2.
            currpt  = rmax * ( zerovec * cos(ϕhalf) + perpvec * sin(ϕhalf) )
            normind = findNormPlane( currpt, testpt, normvec[ee], evec )

            normcurr = normvec[ee][normind]

            Δϕ = ϕ[jj+1] - ϕ[jj]
            np = max( convert( Int64, ceil( Δϕ/(2*π) * nn ) + 1 ), 2 ) # ensure we have at least two points

            ϕcurr = linspace( ϕ[jj], ϕ[jj+1], np )

            # first part
            intpts = Vector{Vector{Float64}}( np )
            cylpts = Vector{Vector{Float64}}( np )

            for ii in 1:np
                currpt = rmax * ( zerovec * cos(ϕcurr[ii]) + perpvec * sin(ϕcurr[ii]) ) # current point in plane perpendicular to evec
                # find actual point near node
                d = dot( -currpt, normvec[ee][normind] ) / dot( evec, normvec[ee][normind] ) # scalar value for which line intersects plane
                intpts[ii] = currpt + d * evec + mesh.p[kk,:]
                cylpts[ii] = sqrt(lat.ar[e1]/π) * ( zerovec * cos(ϕcurr[ii]) + perpvec * sin(ϕcurr[ii]) ) + mesh.p[kk,:] + evec * lvmax[ee]

            end
            # write to file
            for ii in 1:np-1
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
function genFacesCyl( cyl::Int64, mesh::MeshF,
                      ptsCylEdge::Vector{Vector{Vector{Float64}}},
                      fid::IOStream )

    nfac = 0
    frac = 0.25

    e1 = mesh.e[cyl,1]
    e2 = mesh.e[cyl,2]

    ## first edge
    # loop over first edge

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
                vert = [ ptsCylEdge[i2][ ind1+kk ]';
                         ptsCylEdge[i2][ ind1+kk+1 ]';
                         ptsCylEdge[i1][ jj ]' ]

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
    if ind0 > ind1
        for kk = 0:(ind0-ind1-1)
            vert = [ ptsCylEdge[i2][ ind1+kk ]';
                     ptsCylEdge[i2][ ind1+kk+1 ]';
                     ptsCylEdge[i1][ 1 ]' ]

            writeFacetSTLB( vert, fid )
            nfac += 1
        end
    elseif ind0 < ind1
        inds   = vcat(ind1:n2, 1:ind0)
        indsp1 = vcat(ind1+1:n2, 1:ind0+1)
        for kk in 1:length(inds)
            vert = [ ptsCylEdge[i2][ inds[kk] ]';
                     ptsCylEdge[i2][ indsp1[kk] ]';
                     ptsCylEdge[i1][ 1 ]' ]

            writeFacetSTLB( vert, fid )
            nfac += 1
        end
    end

    return nfac

end

function checkPointAct( currpt::Vector{Float64}, testpt::Vector{Float64},
                        normvec::Vector{Vector{Float64}}, evec::Vector{Float64} )

    # println("currpt ", currpt)
    dist    =  norm( currpt - testpt )
    distorg = dist
    act = false
    for ii in 1:length(normvec)
        # loop over intersecting planes
        d = dot( - currpt, normvec[ii] ) / dot( evec, normvec[ii] ) # scalar value for which line intersects plane
        ipt = currpt + d * evec

        cdist = norm( ipt - testpt )

        if cdist < dist
            # this point is closer
            dist = cdist
        end
    end

    # is point active? Only if new distance is same as original distance
    if abs(dist - distorg) < 1e-14
        act = true
    end

    return act

end

function findNormPlane( currpt::Vector{Float64}, testpt::Vector{Float64},
                        normvec::Vector{Vector{Float64}}, evec::Vector{Float64} )

    dist = 3*norm( currpt - testpt )
    ind  = 0

    for ii in 1:length(normvec)
        # loop over intersecting planes
        d = dot( - currpt, normvec[ii] ) / dot( evec, normvec[ii] ) # scalar value for which line intersects plane
        ipt = currpt + d * evec

        cdist = norm( ipt - testpt )

        if cdist < dist
            # this point is closer
            dist = cdist
            ind  = ii
        end
        # end
    end

    return ind

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

function uniqueNumPairs( n::Int64 )
    #   (n + k - 1) nCr (k)
    # which is
    #   (n + k - 1)! / ( k! * (n-1)! )
    # here k = 2, so
    #   (n + 1)! / ( 2! * (n-1)! ) = (n + 1)! / ( 2 * (n-1)! )
    return convert( Int64, factorial( n+1 ) / ( 2*factorial( n-1 ) ) )
end
