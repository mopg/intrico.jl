function genSTLslice( mesh::MeshF,  lat::Lattice, flname::String,
                      nprinter::SVector{3,Float64}, z::Float64,
                      actptsNods::Vector{Vector{SVector{3,Float64}}},
                      edgesNods::Vector{Vector{SVector{4,Int64}}};
                      n::Int64 = 8 )

    # z is defined from origin

    # TODO: - horizontal
    #       - node geometry

    # open STL file
    fid = open( flname, "w" )

    # write header
    for jj in 1:80
        write(fid,' ')
    end
    write(fid,UInt32(0)) # number of faces is unknown now, so write 0

    # z location in part coordinates
    zpart = z * nprinter

    # find intersecting edges
    edgesInt = findEdgesInt( mesh, lat, nprinter, zpart )

    # generate the slices
    factot = [ 0 ]
    writeNodSlice( mesh, lat, nprinter, zpart, z,
                   edgesInt, n, fid, factot,
                   actptsNods, edgesNods )

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

function findEdgesInt( mesh::MeshF, lat::Lattice,
                       nprinter::SVector{3,Float64}, zpart::SVector{3,Float64} )

    ne = size(mesh.e,1)

    edgesInt = fill( false, ne )

    # loop over edges
    for ee in 1:ne

        if lat.ar[ee] < 1e-10
            continue
        end

        inod1 = mesh.e[ee,1]
        inod2 = mesh.e[ee,2]
        nod1  = SVector{3,Float64}( mesh.p[inod1,:] )
        nod2  = SVector{3,Float64}( mesh.p[inod2,:] )
        ledge = norm( nod2 - nod1 )
        nvec  = SVector{3,Float64}( (nod2 - nod1) ./ ledge )

        # check if edge is intersecting or even remotely close to intersecting
        xint, d = findPlaneLineIntersec( zpart, nod1, nprinter, nvec )
        if d >= -ledge && d <= 2*ledge
            # may be intersecting
            edgesInt[ee] = true
        end

    end

    return edgesInt

end

function writeNodSlice( mesh::MeshF, lat::Lattice,
                        nprinter::SVector{3,Float64},
                        zpart::SVector{3,Float64}, z::Float64,
                        edgesInt::Vector{Bool},
                        n::Int64, fid::IOStream, factot::Vector{Int64},
                        actptsNods::Vector{Vector{SVector{3,Float64}}},
                        edgesNods::Vector{Vector{SVector{4,Int64}}} )

    zs = mesh.p * nprinter
    minz = minimum( zs )
    maxz = maximum( zs )

    φ = linspace( 0, 2*π, n+1 )
    xpts = Vector{ SVector{3,Float64} }( n+1 )

    edgWritten = fill( false, size(mesh.e,1) )

    # loop over nodes
    for nn in 1:mesh.n

        if !any(edgesInt[ mesh.n2e[nn] ])
            # no intersections
            continue
        end

        nedg = length(mesh.n2e[nn])

        # check if cutting through node
        intersecEdge = fill( false, nedg )
        for kk in 1:length(edgesNods[nn])
            nod1 = mesh.p[nn,:] + actptsNods[nn][ edgesNods[nn][kk][1] ]
            nod2 = mesh.p[nn,:] + actptsNods[nn][ edgesNods[nn][kk][2] ]
            if (dot(nod1,nprinter) - z) * ( dot(nod2,nprinter) - z ) < 0.
                # println("edgesNods[nn][kk] ", edgesNods[nn][kk])
                # println("mesh.n2e[nn] ", mesh.n2e[nn])
                # println(" bool1 ", mesh.n2e[nn] .== edgesNods[nn][kk][3] )
                # println(" bool2 ", mesh.n2e[nn] .== edgesNods[nn][kk][4] )
                intersecEdge[ mesh.n2e[nn] .== edgesNods[nn][kk][3] ] = true
                intersecEdge[ mesh.n2e[nn] .== edgesNods[nn][kk][4] ] = true
            end
        end

        # compute dist from node -- use pairings from edgesNods
        distNode = Vector{Float64}( nedg )
        radmax = sqrt( maximum(lat.ar[ mesh.n2e[nn] ]) / π )
        for kk in 1:length(edgesNods[nn])

            inde1 = edgesNods[nn][kk][3]
            inde2 = edgesNods[nn][kk][4]
            inod1 = mesh.e[inde1,1] + mesh.e[inde1,2] - nn
            inod2 = mesh.e[inde2,1] + mesh.e[inde2,2] - nn
            vec1 = mesh.p[inod1,:] - mesh.p[nn,:]
            vec2 = mesh.p[inod2,:] - mesh.p[nn,:]

            arg = dot(vec1,vec2)/( norm(vec1)*norm(vec2) )
            arg = minimum( ( 1.0, arg) ) # guard for round-off errors
            arg = maximum( (-1.0, arg) ) # guard for round-off errors
            θ = acos( arg )

            # compute absolute distance between node center and intersection
            lv = radmax / tan( θ/2 ) # offset such that rods do not collide
            l = sqrt( lv^2 + radmax^2 ) # absolute distance

            iedgloc1 = find(mesh.n2e[nn] .== edgesNods[nn][kk][3])[]
            iedgloc2 = find(mesh.n2e[nn] .== edgesNods[nn][kk][4])[]
            distNode[ iedgloc1 ] = max( distNode[ iedgloc1 ], l )
            distNode[ iedgloc2 ] = max( distNode[ iedgloc2 ], l )

        end

        # write slices for struts that do not collide with other struts
        for eloc in 1:nedg
            ee = mesh.n2e[nn][eloc]

            if intersecEdge[eloc] || edgWritten[ee]
                continue
            end

            rad = sqrt( lat.ar[ee] / π )

            inod1 = mesh.e[ee,1]
            inod2 = mesh.e[ee,2]
            nod1  = SVector{3,Float64}( mesh.p[nn,:] )
            nod2  = SVector{3,Float64}( mesh.p[inod1+inod2-nn,:] )
            ledge = norm( nod2 - nod1 )
            nvec  = SVector{3,Float64}( (nod2 - nod1) ./ ledge )

            nperp = cross( nvec, nprinter )
            nperp = nperp/norm(nperp)

            nzero = SVector{3}( nprinter[3], nprinter[1], nprinter[2] )
            if (abs( dot(nprinter,nvec) ) - 1.) < 1e-10
                nzero = cross(nzero,nprinter)
                nperp = cross(nprinter,nzero)
                nperp = nperp/norm(nperp)
            else
                nzero = SVector{3,Float64}( nullspace( Matrix([nvec nperp]') )[:] )
            end
            nzero = nzero/norm(nzero)

            d = fill(0., 3)
            intsec0, d[1] = findPlaneLineIntersec( zpart, nod1, nprinter, nvec )

            # first line
            currpt = nod1 + nzero*rad
            intsec1, d[2] = findPlaneLineIntersec( zpart, currpt, nprinter, nvec )

            # second line
            currpt = nod1 - nzero*rad
            intsec2, d[3] = findPlaneLineIntersec( zpart, currpt, nprinter, nvec )

            if all( (d .< distNode[eloc]) .| (d .> 0.5*ledge) )
                # No intersections with first half of edge
                continue
            end

            edgWritten[ee] = true

            # generate points for STL
            vec1 = rad * nperp
            vec2 = intsec1 - intsec0
            for jj in 1:n+1
                xpts[jj] = intsec0 + vec1 * cos(φ[jj]) + vec2 * sin(φ[jj])
            end

            writeSTLcontour( xpts, intsec0, nprinter, n, fid )
            factot[1] += n

        end

    end

end

function findPlaneLineIntersec( planept::SVector{3,Float64},
                                linept::SVector{3,Float64},
                                nPlane::SVector{3,Float64},
                                nLine::SVector{3,Float64} )

    d    = dot( planept-linept, nPlane ) / dot( nLine, nPlane )
    xloc = linept + d * nLine

    return xloc, d

end

function writeSTLcontour( xpts::Vector{SVector{3,Float64}},
                          intersec0::SVector{3,Float64},
                          nprinter::SVector{3,Float64},
                          n::Int64, fid::IOStream )

      for jj in 1:n

          xface = SVector{3,SVector{3,Float64}}( xpts[jj], xpts[jj+1], intersec0 )

          sterno.writeFacetSTLB( xface, nprinter, fid )

      end

end
