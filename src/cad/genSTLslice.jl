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
        intersecNodEdge = fill( false, length(edgesNods[nn]) )
        for kk in 1:length(edgesNods[nn])
            nod1 = mesh.p[nn,:] + actptsNods[nn][ edgesNods[nn][kk][1] ]
            nod2 = mesh.p[nn,:] + actptsNods[nn][ edgesNods[nn][kk][2] ]
            if (dot(nod1,nprinter) - z) * ( dot(nod2,nprinter) - z ) < 0.
                # println("edgesNods[nn][kk] ", edgesNods[nn][kk])
                # println("mesh.n2e[nn] ", mesh.n2e[nn])
                # println(" bool1 ", mesh.n2e[nn] .== edgesNods[nn][kk][3] )
                # println(" bool2 ", mesh.n2e[nn] .== edgesNods[nn][kk][4] )
                intersecNodEdge[kk] = true
                intersecEdge[ mesh.n2e[nn] .== edgesNods[nn][kk][3] ] = true
                intersecEdge[ mesh.n2e[nn] .== edgesNods[nn][kk][4] ] = true
            end
        end

        # write complicated node geometry
        writeNodSlice( mesh, lat, nn, edgesNods[nn], actptsNods[nn],
                       intersecNodEdge, nprinter, zpart, z, n, fid, factot )

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
                # nzero = SVector{3,Float64}( nullspace( Matrix([nvec nperp]') )[:] )
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
                          n::Int64, fid::IOStream, factot::Vector{Int64} )

      for jj in 1:n

          xface = SVector{3,SVector{3,Float64}}( xpts[jj], xpts[jj+1], intersec0 )

          sterno.writeFacetSTLB( xface, nprinter, fid )

      end

end

function writeNodSlice( mesh::MeshF, lat::Lattice, nn::Int64,
                        edgesNods::Vector{SVector{4,Int64}},
                        actptsNods::Vector{SVector{3,Float64}},
                        intersecNodEdge::Vector{Bool},
                        nprinter::SVector{3,Float64},
                        zpart::SVector{3,Float64}, z::Float64,
                        n::Int64, fid::IOStream, factot )

    nodvec = SVector{3}( mesh.p[nn,:] )

    println("nn ", nn)
    println("intersecNodEdge ", intersecNodEdge)

    Δn = ( z - dot(nodvec,nprinter) ) * nprinter

    Δnodvec = nodvec + Δn

    nint = count( intersecNodEdge ); println("nint ", nint)

    intersecPts = Vector{SVector{3,Float64}}( nint )
    sintersec   = fill( SVector{2}( 0, 0 ), nint ) # which intersec points for which strut
    strutact    = fill( 0, nint ) # which struts are actually active
    isa = 1

    xpts   = Vector{ SVector{3,Float64} }( n+1 )
    xptscl = Vector{ SVector{3,Float64} }( 3 )
    xptsfl = Vector{ SVector{3,Float64} }( 4 )

    # loop over intersecting edges
    jj = 1
    for kk in 1:length(edgesNods)

        el = edgesNods[kk]

        if !intersecNodEdge[kk]
            continue
        end

        # println("el ", el[3], " ", el[4] )
        # println("strutact ", strutact)

        ind3 = find( strutact .== el[3] )
        if isempty( ind3 )
            ind3 = [isa]
            strutact[isa] = el[3]; isa += 1
        end
        ind4 = find( strutact .== el[4] )
        println("strutact ", strutact)
        if isempty( ind4 )
            ind4 = [isa]
            strutact[isa] = el[4]; isa += 1
        end

        nod1 = mesh.e[ el[3], 1 ] + mesh.e[ el[3], 2 ] - nn
        nod2 = mesh.e[ el[4], 1 ] + mesh.e[ el[4], 2 ] - nn

        nvec1 = SVector{3}( mesh.p[nod1,:] - mesh.p[nn,:]); nvec1 = nvec1/norm(nvec1)
        nvec2 = SVector{3}( mesh.p[nod2,:] - mesh.p[nn,:]); nvec2 = nvec2/norm(nvec2)

        # nmajor = nvec1 + nvec2; nmajor ./= norm(nmajor)
        # nperp  = nvec1 - nvec2; nperp  ./= norm(nperp)
        # nminor = cross(nmajor,nperp)

        ncut = nvec1 - nvec2; ncut = ncut/norm(ncut)
        nline = cross(ncut,nprinter)

        rad = sqrt( lat.ar[ el[3] ] / π )
        arg = dot(nvec1,nvec2)
        # arg = minimum( ( 1.0, arg) ) # guard for round-off errors
        # arg = maximum( (-1.0, arg) ) # guard for round-off errors
        θ = acos( arg )
        lv = rad / tan( θ/2 ) # offset such that rods do not collide
        ℓ = sqrt( lv^2 + rad^2 ) # absolute distance
        if dot(nline,nvec1) < 0.
            nline = -nline
        end

        intersecPts[jj] = findIntersecEdge( nvec1, nline, Δn, rad, ℓ) + nodvec
        println("intersecPts[jj] ", intersecPts[jj] )
        println("el[3] ", el[3], " el[4] ", el[4] )
        println("nvec1 ", nvec1)
        println("nvec2 ", nvec2)
        println(" ")

        # println("sintersec[ind3] ", sintersec[ind3], " ind3 ", ind3 )
        if sintersec[ind3[]][1] == 0
            sintersec[ind3[]] = SVector{2}( jj, 0 )
        else
            sintersec[ind3[]] = SVector{2}( sintersec[ind3[]][1], jj )
        end

        # println("sintersec[ind4] ", sintersec[ind4], " ind4 ", ind4 )
        if sintersec[ind4[]][1] == 0
            sintersec[ind4[]] = SVector{2}( jj, 0 )
        else
            sintersec[ind4[]] = SVector{2}( sintersec[ind4[]][1], jj )
        end

        jj += 1

    end

    for jj in 1:nint

        ee = strutact[jj]

        println("ee ", ee)

        ipt1 = sintersec[jj][1]
        ipt2 = sintersec[jj][2]
        println("sintersec[jj] ", sintersec[jj])

        rad = sqrt( lat.ar[ee] / π )

        inod1 = mesh.e[ee,1]
        inod2 = mesh.e[ee,2]
        nod1  = SVector{3,Float64}( mesh.p[nn,:] )
        nod2  = SVector{3,Float64}( mesh.p[inod1+inod2-nn,:] )
        ledge = norm( nod2 - nod1 )
        nvec  = SVector{3,Float64}( (nod2 - nod1) ./ ledge )

        if abs( dot(nvec,nprinter) ) < 1e-10

            # flat edge
            xptsfl[1] = intersecPts[ipt1]
            xptsfl[2] = intersecPts[ipt1] + nvec * ( ledge * dot(intersecPts[ipt1] - Δnodvec,nvec) )
            xptsfl[3] = intersecPts[ipt2] + nvec * ( ledge * dot(intersecPts[ipt2] - Δnodvec,nvec) )
            xptsfl[4] = intersecPts[ipt2]

            writeSTLcontour( xptsfl, Δnodvec, nprinter, 3, fid )
            factot[1] += 3

        else

            nperp = cross( nvec, nprinter )
            nperp = nperp ./ norm(nperp)

            nzero = SVector{3}( nprinter[3], nprinter[1], nprinter[2] )
            if abs( abs(dot(nprinter,nvec)) - 1.) < 1e-10
                nzero = cross(nzero,nprinter)
                nperp = cross(nprinter,nzero)
                nperp = nperp/norm(nperp)
            else
                nzero = cross(nvec,nperp)
                # nzero = SVector{3,Float64}( nullspace( Matrix([nvec nperp]') )[:] )
            end
            nzero = nzero/norm(nzero)
            # println("nzero ", nzero)

            intsec0, d0 = findPlaneLineIntersec( zpart, nod1, nprinter, nvec )

            # first line
            currpt = nod1 + nzero*rad
            intsec1, = findPlaneLineIntersec( zpart, currpt, nprinter, nvec )

            # second line
            currpt = nod1 - nzero*rad
            intsec2, = findPlaneLineIntersec( zpart, currpt, nprinter, nvec )

            # find angle limits due to intersection
            vec1 = -rad * nperp
            vec2 =  intsec1 - intsec0

            normvec1 = norm(vec1)
            normvec2 = norm(vec2)

            Δvec = intersecPts[ipt1] - intsec0
            println("intersecPts[ipt1] ", intersecPts[ipt1])
            println("intersecPts[ipt2] ", intersecPts[ipt2])
            println("intsec0 ", intsec0)
            println("norm(Δvec) ", norm(Δvec))
            println("norm(vec1) ", norm(vec1))
            println("norm(vec2) ", norm(vec2))
            println("rad ", rad)
            # println("Δvec ", Δvec)
            # println("arg 1 ", dot(vec2,Δvec), " arg 2 ", dot(vec1,Δvec) )
            φ1 = atan2( dot(vec2,Δvec)/normvec2^2, dot(vec1,Δvec)/normvec1^2 )
            println("φ1 ", φ1)

            Δvec = intersecPts[ipt2] - intsec0
            φ2 = atan2( dot(vec2,Δvec)/normvec2^2, dot(vec1,Δvec)/normvec1^2 )
            println("φ2 ", φ2)

            φ = linspace( min(φ1,φ2), max(φ1,φ2), n+1 )
            for jj in 1:n+1
                xpts[jj] = intsec0 + vec1 * cos(φ[jj]) + vec2 * sin(φ[jj])
            end

            println("start ", xpts[1] )
            println("end   ", xpts[end] )

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
