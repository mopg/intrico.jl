# ---------------------------------------------------------------------------- #
#
#   meshF2D.jl
#
#   Type for 2D frame meshes
#   Inherits from the abstract "meshF" type
#
#   sterno
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

"""
    MeshF2D

MeshF2D type:
Type for 2D meshes consisting of triangles. It holds the nodes and connectivity
information.
"""
struct MeshF2D <: MeshF

    dim::Int64              # Dimension of the problem

    n::Int64                # Number of nodes

    p::Matrix{Float64}      # Nodal locations
    t::Matrix{Int64}        # Triangle - node connectivity
    t2f::Matrix{Int64}      # Triangle - face connectivity
    f::Matrix{Int64}        # Face - node/triangle connectivity
    fb::Matrix{Int64}       # Boundary face info
    nodes::Array{Float64,3} # Nodes on which solution is evaluated -- needed to
                          #  compute Jacobian on face
    n2e::Vector{Vector{Int64}} # Node to edge connectivity

end

"""
    MeshF( mesh::luteos.Mesh2D )

Constructor that generates frame from `luteos` mesh structure in `mesh`.
"""
function MeshF( mesh::luteos.Mesh2D )

    n2e = genNEconnec2D( mesh.f, mesh.n )

    MeshF2D( 2, mesh.n, mesh.p, mesh.t, mesh.t2f, mesh.f, mesh.fb, mesh.nodes, n2e )

end

"""
    MeshF2D( name::String )

Constructor for frame from mesh written in `name`.
"""
function MeshF2D( name::String )

    if name[end-3:end] == ".su2"
        (p_, t_, bel_, tags_) = luteos.readSU2_2D( name )
    elseif name[end-3:end] == ".msh" # BAMG
        (p_, t_, bel_ ) = luteos.readBAMG( name )
    elseif name[end-4:end] == ".mesh" # FEFLOA
        (p_, t_, bel_ ) = luteos.readFEFLOA_2D( name )
    else
        error("MeshF2D: Unknown mesh type")
    end

    porder_ = luteos.P1()

    (f_, t2f_, nodes_, ploc_, tloc_, fb_) = luteos.genmesh( porder_, p_, t_, bel_ )

    n_ = size( p_, 1 )

    n2e_ = genNEconnec2D( f_, n_ )

    MeshF2D( 2, n_, p_, t_, t2f_, f_, fb_, nodes_, n2e_ )

end

"""
    genNEconnec3D( e::Matrix{Int64}, n::Int64 )

Generates node to edge connectivity given edge connectivity (`e`).
"""
function genNEconnec2D( e::Matrix{Int64}, n::Int64 )

    n2e = fill( fill(0,0), n )

    for jj in 1:size(e,1)

        ii1 = e[jj,1]
        ii2 = e[jj,2]

        if isempty( n2e[ii1] .== jj )
            append!(n2e[ii1], jj)
        end
        if isempty( n2e[ii2] .== jj )
            append!(n2e[ii2], jj)
        end

    end

    return n2e

end
