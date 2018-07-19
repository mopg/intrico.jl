function genSTLslice( mesh::MeshF,  lat::Lattice, flname::String,
                      nprinter::SVector{3,Float64}, z::Float64,
                      actptsNods::Vector{Vector{SVector{3,Float64}}},
                      edgesNods::Vector{Vector{SVector{5,Int64}}};
                      n::Int64 = 8 )

    # z is defined from origin

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
    writeNodStrutSlice( mesh, lat, nprinter, zpart, z,
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

function writeNodStrutSlice( mesh::MeshF, lat::Lattice,
                             nprinter::SVector{3,Float64},
                             zpart::SVector{3,Float64}, z::Float64,
                             edgesInt::Vector{Bool},
                             n::Int64, fid::IOStream, factot::Vector{Int64},
                             actptsNods::Vector{Vector{SVector{3,Float64}}},
                             edgesNods::Vector{Vector{SVector{5,Int64}}} )

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
        intersecNodEdge = fill( false, length(edgesNods[nn]) )
        for kk in 1:length(edgesNods[nn])
            nod1 = mesh.p[nn,:] + actptsNods[nn][ edgesNods[nn][kk][1] ]
            nod2 = mesh.p[nn,:] + actptsNods[nn][ edgesNods[nn][kk][2] ]
            if (dot(nod1,nprinter) - z) * ( dot(nod2,nprinter) - z ) < 0.

                intersecNodEdge[kk] = true

                intersecEdge[ mesh.n2e[nn] .== edgesNods[nn][kk][3] ] = true
                intersecEdge[ mesh.n2e[nn] .== edgesNods[nn][kk][4] ] = true
            end
        end

        # write complicated node geometry
        writeNodSlice( mesh, lat, nn, edgesNods[nn], actptsNods[nn],
                       intersecEdge, intersecNodEdge, nprinter, zpart, z, n, fid, factot )

       if any(intersecEdge)
           # intersecting the node
           continue
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
            # l = sqrt( lv^2 + radmax^2 ) # absolute distance

            iedgloc1 = find(mesh.n2e[nn] .== edgesNods[nn][kk][3])[]
            iedgloc2 = find(mesh.n2e[nn] .== edgesNods[nn][kk][4])[]
            distNode[ iedgloc1 ] = max( distNode[ iedgloc1 ], lv )
            distNode[ iedgloc2 ] = max( distNode[ iedgloc2 ], lv )

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
            if (abs( abs(dot(nprinter,nvec)) - 1.) ) < 1e-10
                nzero = cross(nzero,nprinter)
                nperp = cross(nprinter,nzero)
                nperp = nperp/norm(nperp)
            else
                nzero = cross(nvec,nperp)
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
            vec1 = - rad * nperp
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

function writeNodSlice( mesh::MeshF, lat::Lattice, nn::Int64,
                        edgesNods::Vector{SVector{5,Int64}},
                        actptsNods::Vector{SVector{3,Float64}},
                        intersecEdge::Vector{Bool}, intersecNodEdge::Vector{Bool},
                        nprinter::SVector{3,Float64},
                        zpart::SVector{3,Float64}, z::Float64,
                        n::Int64, fid::IOStream, factot::Vector{Int64} )

    nodvec = SVector{3}( mesh.p[nn,:] )

    Δn = ( z - dot(nodvec,nprinter) ) * nprinter

    Δnodvec = nodvec + Δn

    nint = count( intersecEdge )

    xpts   = Vector{ SVector{3,Float64} }( n+1 )
    xptscl = Vector{ SVector{3,Float64} }( 3 )
    xptsfl = Vector{ SVector{3,Float64} }( 4 )

    # loop over intersecting edges
    for kk in 1:length( mesh.n2e[nn] )

        e0 = mesh.n2e[nn][kk]

        if !intersecEdge[kk]
            continue
        end

        # find the other two neighboring edges
        e1 = 0
        e2 = 0
        enew = 0
        ihalf1 = 0
        ihalf2 = 0
        ct = 1
        for el in edgesNods
            if !intersecNodEdge[ct]
                ct += 1
                continue
            end
            if el[3] == e0
                enew = el[4]
            elseif el[4] == e0
                enew = el[3]
            else
                ct += 1
                continue
            end
            if e1 > 0
                e2 = enew
                ihalf2  = el[5]
            else
                e1 = enew
                ihalf1  = el[5]
            end
            ct += 1
        end

        # find which edge is active on which side by comparing ncut vectors
        nod0 = mesh.e[ e0, 1 ] + mesh.e[ e0, 2 ] - nn
        nod1 = mesh.e[ e1, 1 ] + mesh.e[ e1, 2 ] - nn
        nod2 = mesh.e[ e2, 1 ] + mesh.e[ e2, 2 ] - nn

        nvec0 = SVector{3}( mesh.p[nod0,:] - mesh.p[nn,:]); nvec0 = nvec0/norm(nvec0)
        nvec1 = SVector{3}( mesh.p[nod1,:] - mesh.p[nn,:]); nvec1 = nvec1/norm(nvec1)
        nvec2 = SVector{3}( mesh.p[nod2,:] - mesh.p[nn,:]); nvec2 = nvec2/norm(nvec2)

        ncut1 = nvec0 - nvec1; ncut1 = ncut1/norm(ncut1)
        ncut2 = nvec0 - nvec2; ncut2 = ncut2/norm(ncut2)
        nline1 = cross(ncut1,nprinter); nline1 = nline1/norm(nline1)
        nline2 = cross(ncut2,nprinter); nline2 = nline2/norm(nline2)

        # find correct orientation to make sure that nline is aligned with the halfvector of the edge (as projected from the edge)
        if dot( nline1 - dot(nline1,nvec1)*nvec1, actptsNods[ihalf1] ) < 0.
            nline1 = -nline1
        end
        if dot( nline2 - dot(nline2,nvec2)*nvec2, actptsNods[ihalf2] ) < 0.
            nline2 = -nline2
        end

        rad  = sqrt( lat.ar[ e0 ] / π )
        arg1 = dot(nvec0,nvec1)
        # arg = minimum( ( 1.0, arg) ) # guard for round-off errors
        # arg = maximum( (-1.0, arg) ) # guard for round-off errors
        θ1  = acos( arg1 )
        lv1 = rad / tan( θ1/2 ) # offset such that rods do not collide
        ℓ1  = sqrt( lv1^2 + rad^2 ) # absolute distance
        arg2 = dot(nvec0,nvec2)
        # arg = minimum( ( 1.0, arg) ) # guard for round-off errors
        # arg = maximum( (-1.0, arg) ) # guard for round-off errors
        θ2  = acos( arg2 )
        lv2 = rad / tan( θ2/2 ) # offset such that rods do not collide
        ℓ2  = sqrt( lv2^2 + rad^2 ) # absolute distance

        planeoff = dot( nprinter, Δn )
        cent1 = planeoff * cross( ncut1, cross(nprinter, ncut1 ) ) / norm( cross(ncut1, nprinter) )^2 # point of plane-plane intersection
        intersec1 = findIntersecEdge( nvec0, nline1, cent1, rad, ℓ1) + nodvec

        cent2 = planeoff * cross( ncut2, cross(nprinter, ncut2 ) ) / norm( cross(ncut2, nprinter) )^2 # point of plane-plane intersection
        intersec2 = findIntersecEdge( nvec0, nline2, cent2, rad, ℓ2) + nodvec

        # --- save this somewhere TODO for more efficient algorithm

        # --- write contour
        ledge = norm( mesh.p[nod0,:] - nodvec )
        if abs( dot(nvec0,nprinter) ) < 1e-10

            if dot( cross( intersec1 - nodvec, intersec2 - nodvec ), nprinter ) < 0.
                temp1 = intersec1
                intersec1 = intersec2
                intersec2 = temp1
            end
            # flat edge
            xptsfl[1] = intersec1
            xptsfl[2] = intersec1 + nvec0 * ( 0.5*ledge - dot(intersec1 - Δnodvec,nvec0) )
            xptsfl[3] = intersec2 + nvec0 * ( 0.5*ledge - dot(intersec2 - Δnodvec,nvec0) )
            xptsfl[4] = intersec2

            writeSTLcontour( xptsfl, Δnodvec, nprinter, 3, fid )
            factot[1] += 3

        else

            nperp = cross( nvec0, nprinter )
            nperp = nperp ./ norm(nperp)

            nzero = SVector{3}( nprinter[3], nprinter[1], nprinter[2] )
            if abs( abs(dot(nprinter,nvec0)) - 1.) < 1e-10
                nzero = cross(nzero,nprinter)
                nperp = cross(nprinter,nzero)
                nperp = nperp/norm(nperp)
            else
                nzero = cross(nvec0,nperp)
            end
            nzero = nzero/norm(nzero)

            intsec0, d0 = findPlaneLineIntersec( zpart, nodvec, nprinter, nvec0 )

            # first line
            currpt = nodvec + nzero*rad
            intsec1, = findPlaneLineIntersec( zpart, currpt, nprinter, nvec0 )

            # second line
            currpt = nodvec - nzero*rad
            intsec2, = findPlaneLineIntersec( zpart, currpt, nprinter, nvec0 )

            # find angle limits due to intersection
            vec1 = -rad * nperp
            vec2 =  intsec1 - intsec0
            if dot( cross(vec1,vec2), nprinter ) < 0.
                vec2 = -vec2
            end

            normvec1 = norm(vec1)
            normvec2 = norm(vec2)

            Δvec = intersec1 - intsec0
            φ1 = atan2( dot(vec2,Δvec)/normvec2^2, dot(vec1,Δvec)/normvec1^2 )

            Δvec = intersec2 - intsec0
            φ2 = atan2( dot(vec2,Δvec)/normvec2^2, dot(vec1,Δvec)/normvec1^2 )

            φmin = min(φ1,φ2)
            φmax = max(φ1,φ2)
            if (φmax - φmin) > π
                φ1 = φmax
                φmax = φmin
                φmin = φ1 - 2*π
            end

            φ = linspace( φmin, φmax, n+1 )
            for jj in 1:n+1
                xpts[jj] = intsec0 + vec1 * cos(φ[jj]) + vec2 * sin(φ[jj])
            end

            writeSTLcontour( xpts, intsec0, nprinter, n, fid )
            factot[1] += n

            # close the contour
            xptscl[1] = xpts[1]
            xptscl[3] = xpts[end]
            xptscl[2] = intsec0

            if d0 > 0.
                writeSTLcontour( xptscl, Δnodvec, nprinter, 2, fid )
                factot[1] += 2
            end

        end

    end

end


function findIntersecEdge( nvec1::SVector{3,Float64}, nline::SVector{3,Float64},
                           Δn::SVector{3,Float64}, r::Float64, ℓ::Float64 )

    # This uses the bisection method, so temporary.

    b1 = 0.
    b2 = 2*ℓ

    b = 0.
    rcomp = 0.
    while (b2 - b1)/b2 > 1e-5

        b = (b2 + b1)/2
        tempvec = b*nline + Δn
        rcomp = norm( tempvec - dot(tempvec,nvec1)*nvec1 )

        if rcomp > r
            b2 = b
        elseif rcomp < r
            b1 = b
        else
            break
        end

    end

    return (b*nline + Δn)

end
