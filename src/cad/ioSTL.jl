# ---------------------------------------------------------------------------- #
#
#   ioSTL.jl
#
#   STL i/o routines
#
#   intrico
#   Spring 2018
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

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

"""
    writeFacetSTLB( vert::SVector{3,SVector{3,Float64}}, fid::IOStream )

Writes one STL facet for binary file.
"""
function writeFacetSTLB( vert::SVector{3,SVector{3,Float64}},
                         normal::SVector{3,Float64}, fid::IOStream )

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
