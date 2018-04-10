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
function genSTLnew( mesh::MeshF,  lat::Lattice, flname::String; n::Int64 = 8 )

    edgWrite = Vector{Vector{Bool}}( size(mesh.e,1) )
    for jj in 1:size(mesh.e,1)
        edgWrite[jj] = [true]
    end

    # open STL file
    fid = open( flname, "w" )

    # write header
    for jj in 1:80
        write(fid,' ')
    end
    write(fid,UInt32(0)) # number of faces is unknown now, so write 0

    # generate all faces (and write)
    factot = [ 0 ]
    genFacesNodsCyl( mesh, lat, n, factot, fid, edgWrite )

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
    write( fid, UInt32( factot[1] ) )
    # close STL
    close( fid )
end

"""
    genSTL( mesh::MeshF,  lat::Lattice, flname::String; n = 8, name = "object" )

Generates a binary .stl file for the lattice in `mesh` with the areas defined in `lat`.
"""
function genSTLnew( mesh::MeshF,  lat::Lattice, flnames::Vector{String},
                    edgWrite::Vector{Vector{Bool}}; n::Int64 = 8 )

    # open STL file
    fid = Vector{IOStream}( length(flnames) )
    for ii in 1:length(flnames)

        fid[ii] = open( flnames[ii], "w" )

        # write header
        for jj in 1:80
            write(fid[ii],' ')
        end
        write(fid[ii],UInt32(0)) # number of faces is unknown now, so write 0

    end

    # generate all faces (and write)
    factot = fill( 0, length(flnames) )
    genFacesNodsCyl( mesh, lat, n, factot, fid, edgWrite )

    # close STL
    for ii in 1:length(flnames)
        close( fid[ii] )
    end

    ## reopen to write the number of faces
    for ii in 1:length(flnames)

        fid = open( flnames[ii], "a+" )
        seekstart(fid) # go back to top of file

        # rewrite header
        for jj in 1:80
            write(fid,' ')
        end
        # write correct number of faces
        write( fid,UInt32( factot[ii] ) )
        # close STL
        close( fid )

    end

end

"""
    genFacesNods( mesh::MeshF, lat::Lattice, fdist::Matrix{Float64}, n::Int64, fid::IOStream )

Writes the STL facets for all nodes for ASCII files.
"""
function genFacesNodsCyl( mesh::MeshF, lat::Lattice, n::Int64, nfac::Vector{Int64},
                          fid::Union{IOStream,Vector{IOStream}},
                          edgWrite::Vector{Vector{Bool}} )

    ptsCylEdge = Vector{Vector{Vector{SVector{3,Float64}}}}( size(mesh.e,1) )
    for kk in 1:size(mesh.e,1)
        ptsCylEdge[kk] = Vector{Vector{SVector{3,Float64}}}(2)
    end

    # generate nodes
    for nn in 1:mesh.n

        genFacesNod( nn, mesh, lat, n, ptsCylEdge, nfac, fid, edgWrite )

    end

    # generate cylinders
    for kk in 1:size(mesh.e,1)
        if lat.ar[ kk ] > 0.0 && sqrt( lat.ar[ kk ] / π ) > 1.e-14
            genFacesCyl( kk, mesh, ptsCylEdge[kk], nfac, fid, edgWrite[kk] )
        end
    end

end

"""
    compDist( mesh::MeshF,  lat::Lattice )

Computes the distance from the end face to the node
"""
function genFacesNod( kk::Int64, mesh::MeshF,  lat::Lattice, nn::Int64,
                      ptsCylEdge::Vector{Vector{Vector{SVector{3,Float64}}}},
                      nfac::Vector{Int64},
                      fid::Union{IOStream,Vector{IOStream}},
                      edgWrite::Vector{Vector{Bool}} )

    # generate pairs of edges
    edg   = [i for i in 1:length(mesh.n2e[kk])]
    #   check for min area
    inddel = Vector{Int64}(0)
    for ii in 1:length(mesh.n2e[kk])
        if lat.ar[ mesh.n2e[kk][ii] ] < 0.0 || sqrt( lat.ar[ mesh.n2e[kk][ii] ] / π ) < 1.e-14
            append!(inddel,ii)
        end
    end
    edg0     = [i for i in 1:length(edg)-length(inddel)]
    indedges = [jj for jj in 1:length(edg)]
    deleteat!(indedges,inddel)

    if length(edg0) == 0 # if all nodes are zero thickness, can skip the whole function
        return nothing
    end

    edg1  =   edg0'
    eedg  =   edg0 .+ 0*edg1
    eedg1 = 0*edg0 .+   edg1
    pairs = hcat( eedg[:], eedg1[:] )

    pairs = unique(sort(pairs,2), 1) # unique pairs
    pairs = pairs[ sortperm( pairs[:,1] + pairs[:,2]*100 ), : ] # ensure highest indices are last
    # TODO: don't need to do this for each loop, can just precompute this for n = 40 and then just take part of the unique pairs
    rmax = sqrt( maximum( lat.ar[ mesh.n2e[kk] ] ) / π )

    # compute crossing point and normal vector
    normvec  = Vector{Vector{SVector{3,Float64}}}( length(edg) )
    # tangvec  = Vector{Vector{Vector{Float64}}}( length(edg) )
    lvmax = fill( 0.0, length(edg) )
    for jj in 1:length(edg)
        normvec[jj]  = Vector{SVector{3,Float64}}( length(edg0)-1 )
        # tangvec[jj]  = Vector{Vector{Float64}}( length(edg)-1 )
    end
    indcnt = fill(1, length(edg) )
    for jj in 1:size(pairs,1)
        edg_i1 = indedges[ pairs[jj,1] ]
        edg_i2 = indedges[ pairs[jj,2] ]
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

        # normal vectors
        normvec[edg_i1][ indcnt[edg_i1] ] = SVector{3}( (vec1-vec2) / magdiff )
        normvec[edg_i2][ indcnt[edg_i2] ] = SVector{3}( (vec2-vec1) / magdiff )

        indcnt[edg_i1] += 1
        indcnt[edg_i2] += 1

    end

    # add offset to rods # TODO: this needs to be better such that things do not collide.
    # TODO need to do a check here for collisions
    for ii in 1:length(lvmax)
        lvmax[ii] += 0.05 * norm( mesh.p[ mesh.e[mesh.n2e[kk][edg[ii]],1],: ] - mesh.p[ mesh.e[mesh.n2e[kk][edg[ii]],2],: ] )
    end

    # loop over each edge
    for eecnt in 1:length(edg0)
        ee = indedges[ eecnt ]

        e1     = mesh.n2e[kk][ edg[ ee ] ]
        nod    = mesh.e[e1, find( mesh.e[e1,:] .!= kk )[] ]
        evec   = mesh.p[ nod, : ] - mesh.p[ kk, : ]
        evec ./= norm(evec)

        # find intersections between planes and check if they are active
        testpt   = mesh.p[ nod, : ] - mesh.p[ kk, : ] # distance to this point is measured to determine which is closer, distance is relative to current node
        nunique  = uniqueNumPairs( length(edg0)-1 )
        intersec = Vector{SVector{3,Float64}}( max( 2*(nunique - (length(edg0)-1)) + 1, 2 ) )

        nact = 0
        if nunique == 1 # if nunique is 1, there are no intersections

            indother    = find( indedges .!= ee )[]
            eother      = mesh.n2e[kk][ edg[ indedges[ indother ] ] ]
            nodother    = mesh.e[eother, find( mesh.e[eother,:] .!= kk )[] ]
            evecother   = mesh.p[ nodother, : ] - mesh.p[ kk, : ]
            evecother ./= norm( evecother )

            lvec = evec + evecother
            if norm(lvec) < 1e-14
                lvec = normvec[ee][1]
            end

            nrmperp = norm( lvec - dot( lvec, evec )*evec )
            lvec  .*=  rmax / nrmperp

            nact += 1
            intersec[nact] = SVector{3}( 1.0 * lvec)

            nact += 1
            intersec[nact] = SVector{3}(-1.0 * lvec)

        else
            for jj in 1:nunique # note: number of planes is length(edg)-1, so can reuse part of upairs

                if pairs[jj,1] == pairs[jj,2]
                    continue
                end
                ip1 = pairs[jj,1]
                ip2 = pairs[jj,2]

                lvec     = cross( normvec[ee][ ip1 ], normvec[ee][ ip2 ] ) # always going through the node center (because both planes go through that point)
                nrmperp  = norm( lvec - dot( lvec, evec )*evec )

                # check first point
                lvec1 =   lvec * rmax / nrmperp

                act = checkPointAct( lvec1, testpt, normvec[ee], evec )

                if act
                    nact += 1
                    intersec[nact] = lvec1
                end

                # check second point
                lvec2 = - lvec * rmax / nrmperp
                # check ...
                act = checkPointAct( lvec2, testpt, normvec[ee], evec )
                if act
                    nact += 1
                    intersec[nact] = lvec2
                end

            end
        end

        ## order points counterclockwise
        # pick first point as ϕ = 0
        zerovec   = (intersec[1] - dot( intersec[1], evec ) * evec) /
                    norm( intersec[1] - dot( intersec[1], evec ) * evec )
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
            if (ϕ[jj+1] - ϕ[jj]) > 1e-14 && (2*π - ϕ[jj+1]) > 1e-14
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
        ptsCylEdge[e1][iedge] = Vector{SVector{3,Float64}}(0)
        for jj in 1:nact-1

            # figure out on which plane to project
            ϕhalf   = ( ϕ[jj] + ϕ[jj+1] ) / 2.
            currpt  = rmax * ( zerovec * cos(ϕhalf) + perpvec * sin(ϕhalf) )
            normind = findNormPlane( currpt, testpt, normvec[ee], evec )

            normcurr = normvec[ee][normind]

            Δϕ = ϕ[jj+1] - ϕ[jj]
            @assert Δϕ > 1e-14
            np = max( convert( Int64, ceil( (Δϕ-1e-5)/(2*π) * nn ) + 1 ), 2 ) # ensure we have at least two points

            ϕcurr = linspace( ϕ[jj], ϕ[jj+1], np )

            # first part
            intpts = Vector{SVector{3,Float64}}( np )
            cylpts = Vector{SVector{3,Float64}}( np )

            for ii in 1:np
                currpt = rmax * ( zerovec * cos(ϕcurr[ii]) + perpvec * sin(ϕcurr[ii]) ) # current point in plane perpendicular to evec
                # find actual point near node
                d = dot( -currpt, normvec[ee][normind] ) / dot( evec, normvec[ee][normind] ) # scalar value for which line intersects plane
                intpts[ii] = SVector{3}( currpt + d * evec + mesh.p[kk,:] )
                cylpts[ii] = SVector{3}( sqrt(lat.ar[e1]/π) * ( zerovec * cos(ϕcurr[ii]) + perpvec * sin(ϕcurr[ii]) ) + mesh.p[kk,:] + evec * lvmax[ee] )

            end
            # write to file
            for ii in 1:np-1
                vert = SVector{3,SVector{3,Float64}}( cylpts[ii], intpts[ii], intpts[ii+1] )
                writeFacetSTLB( vert, nfac, fid, edgWrite[e1] )

                vert = SVector{3,SVector{3,Float64}}( intpts[ii+1], cylpts[ii+1], cylpts[ii] )
                writeFacetSTLB( vert, nfac, fid, edgWrite[e1] )

            end
            append!( ptsCylEdge[e1][iedge], cylpts[1:end-1] )
            if jj == nact-1
                append!( ptsCylEdge[e1][iedge], cylpts[end:end] )
            end

        end

    end

end

"""
    compDist( mesh::MeshF,  lat::Lattice )

Computes the distance from the end face to the node
"""
function genFacesCyl( cyl::Int64, mesh::MeshF,
                      ptsCylEdge::Vector{Vector{SVector{3,Float64}}}, nfac::Vector{Int64},
                      fid::Union{IOStream,Vector{IOStream}}, edgWrite::Vector{Bool} )

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
    vert = SVector{3,SVector{3,Float64}}( ptsCylEdge[i1][1], ptsCylEdge[i1][2], ptsCylEdge[i2][ind0] )

    writeFacetSTLB( vert, nfac, fid, edgWrite )

    ind00 = ind0 # need to save this index to close loop

    # generate rest of facets (also for other ring)
    for jj in 2:n1-1

        ind1 = findMinDistInd( ptsCylEdge, i1, i2, jj, frac )

        vert = SVector{3,SVector{3,Float64}}( ptsCylEdge[i1][jj], ptsCylEdge[i1][jj+1], ptsCylEdge[i2][ind1] )
        writeFacetSTLB( vert, nfac, fid, edgWrite )

        # generate other faces
        # NOTE: due to ordering ind1 < ind0
        if ind0 > ind1
            for kk = 0:(ind0-ind1-1)
                vert = SVector{3,SVector{3,Float64}}( ptsCylEdge[i2][ ind1+kk ],
                                                      ptsCylEdge[i2][ ind1+kk+1 ],
                                                      ptsCylEdge[i1][ jj ] )
                writeFacetSTLB( vert, nfac, fid, edgWrite )
            end
        elseif ind0 < ind1
            inds   = vcat(ind1:n2-1, 1:ind0-1)
            indsp1 = vcat(ind1+1:n2, 2:ind0)
            for kk in 1:length(inds)
                vert = SVector{3,SVector{3,Float64}}( ptsCylEdge[i2][ inds[kk] ],
                                                      ptsCylEdge[i2][ indsp1[kk] ],
                                                      ptsCylEdge[i1][ jj ] )
                writeFacetSTLB( vert, nfac, fid, edgWrite )
            end
        end

        ind0 = ind1

    end
    # last patch
    ind1 = ind00
    if ind0 > ind1
        for kk = 0:(ind0-ind1-1)
            vert = SVector{3,SVector{3,Float64}}( ptsCylEdge[i2][ ind1+kk ],
                                                  ptsCylEdge[i2][ ind1+kk+1 ],
                                                  ptsCylEdge[i1][ 1 ] )
            writeFacetSTLB( vert, nfac, fid, edgWrite )
        end
    elseif ind0 < ind1
        inds   = vcat(ind1:n2-1, 1:ind0-1)
        indsp1 = vcat(ind1+1:n2, 2:ind0)
        for kk in 1:length(inds)
            vert = SVector{3,SVector{3,Float64}}( ptsCylEdge[i2][ inds[kk] ],
                                                  ptsCylEdge[i2][ indsp1[kk] ],
                                                  ptsCylEdge[i1][ 1 ] )
            writeFacetSTLB( vert, nfac, fid, edgWrite )
        end
    end

end

function checkPointAct( currpt::SVector{3,Float64}, testpt::Vector{Float64},
                        normvec::Vector{SVector{3,Float64}}, evec::Vector{Float64} )

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
    if abs(dist - distorg)/abs(dist) < 1e-14
        act = true
    end

    return act

end

function findNormPlane( currpt::SVector{3,Float64}, testpt::Vector{Float64},
                        normvec::Vector{SVector{3,Float64}}, evec::Vector{Float64} )

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

function findMinDistInd( pts::Vector{Vector{SVector{3,Float64}}},
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
    # here k = 2, so
    #   (n + 1) nCr (k)
    return binomial( n + 1, 2 )
end

"""
    writeFacetSTLB( vert::SVector{3,SVector{3,Float64}}, fid::IOStream )

Writes one STL facet for binary file.
"""
function writeFacetSTLB( vert::SVector{3,SVector{3,Float64}}, nfac::Vector{Int64},
                         fid::IOStream, edgWrite::Vector{Bool} )

    vec1   = vert[2] - vert[1]
    vec2   = vert[3] - vert[1]
    vec1   = cross( vec1, vec2 )
    normal = vec1 / norm(vec1)

    # normals
    write( fid, Float32( normal[1] ) )
    write( fid, Float32( normal[2] ) )
    write( fid, Float32( normal[3] ) )

    # Vertex 1
    write( fid, Float32( vert[1][1] ) )
    write( fid, Float32( vert[1][2] ) )
    write( fid, Float32( vert[1][3] ) )

    # Vertex 2
    write( fid, Float32( vert[2][1] ) )
    write( fid, Float32( vert[2][2] ) )
    write( fid, Float32( vert[2][3] ) )

    # Vertex 3
    write( fid, Float32( vert[3][1] ) )
    write( fid, Float32( vert[3][2] ) )
    write( fid, Float32( vert[3][3] ) )

    # End facet
    write( fid, UInt16(0) )

    nfac[1] += 1

end

function writeFacetSTLB( vert::SVector{3,SVector{3,Float64}}, nfac::Vector{Int64},
                         fid::Vector{IOStream}, edgWrite::Vector{Bool} )

    for jj in 1:length( fid )

        if edgWrite[jj]
            writeFacetSTLB( vert, [0], fid[jj], edgWrite )
            nfac[jj] += 1
        end

    end

end

"""
    writeFacetSTLB( vert::SVector{3,SVector{3,Float64}}, fid::IOStream )

Writes one STL facet for binary file.
"""
function writeFacetSTLB( vert::SVector{3,SVector{3,Float64}}, fid::IOStream )

    vec1   = vert[2] - vert[1]
    vec2   = vert[3] - vert[1]
    vec1   = cross( vec1, vec2 )
    normal = vec1 / norm(vec1)

    # normals
    write( fid, Float32( normal[1] ) )
    write( fid, Float32( normal[2] ) )
    write( fid, Float32( normal[3] ) )

    # Vertex 1
    write( fid, Float32( vert[1][1] ) )
    write( fid, Float32( vert[1][2] ) )
    write( fid, Float32( vert[1][3] ) )

    # Vertex 2
    write( fid, Float32( vert[2][1] ) )
    write( fid, Float32( vert[2][2] ) )
    write( fid, Float32( vert[2][3] ) )

    # Vertex 3
    write( fid, Float32( vert[3][1] ) )
    write( fid, Float32( vert[3][2] ) )
    write( fid, Float32( vert[3][3] ) )

    # End facet
    write( fid, UInt16(0) )

end
