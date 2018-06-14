# ---------------------------------------------------------------------------- #
#
#   meshF3D.jl
#
#   Type for 3D frame meshes
#   Inherits from the abstract "meshF" type
#
#   sterno
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

"""
    MeshF3D

Mesh3D type:
Type for 3D meshes consisting of tetrahedrons.
"""
struct MeshF3D <: MeshF

  dim::Int64              # Dimension of the problem

  n::Int64                # Number of nodes

  p::Matrix{Float64}      # Nodal locations
  t::Matrix{Int64}        # Triangle - node connectivity
  t2f::Matrix{Int64}      # Triangle - face connectivity
  f::Matrix{Int64}        # face - node/triangle connectivity
  fb::Matrix{Int64}       # Boundary face info

  e::Matrix{Int64}        # Edge connectivity
  e2f::Vector{Vector{Int64}} # Edge to face connectivity
  f2e::Matrix{Int64} # Face to edge connectivity
  nodes::Array{Float64,3} # Nodes on which solution is evaluated -- needed to
                          #  compute Jacobian on face

  n2e::Vector{Vector{Int64}} # Node to edge connectivity

  boundtags::Vector{String} # boundary tags

end

"""
    MeshF3D( mesh::luteos.Mesh3D )

Constructor that generates frame from `mesh`.
"""
function MeshF( mesh::luteos.Mesh3D )

  e, e2f, f2e  = genEdgesF3D( mesh.f )

  n2e = genNEconnec3D( e, mesh.n )

  MeshF3D( 3, mesh.n, mesh.p, mesh.t, mesh.t2f, mesh.f, mesh.fb, e, e2f, f2e,
           mesh.nodes, n2e, mesh.boundtags )

end

"""
    MeshF3D( name::String )

Constructor for frame that reads in a mesh from `name`.
"""
function MeshF3D( name::String )

  if name[end-3:end] == ".su2"
    (p_, t_, bel_, tags_ ) = luteos.readSU2_3D( name )
  elseif name[end-4:end] == ".mesh" # FEFLOA
    (p_, t_, bel_, tags_ ) = luteos.readFEFLOA_3D( name )
  else
    error("Unknown mesh type")
  end

  porder_ = luteos.P1()

  (f_, t2f_, nodes_, ploc_, tloc_, trorder_, fb_) = luteos.genmesh3D( porder_, p_, t_, bel_ )

  n_ = size( p_, 1 )

  e_, e2f_, f2e_ = genEdgesF3D( mesh.f )

  n2e = genNEconnec3D( e_, n_ )

  MeshF3D( 3, n_, p_, t_, t2f_, f_, fb_, e_, e2f_, f2e_, nodes_, n2e, tags_ )

end

"""
    genEdgesF3D( f::Matrix{Int64} )

Generates edge connectivity given face connectivity (`f`).
"""
function genEdgesF3D( f::Matrix{Int64} )

  nf = size(f,1)

  edges = vcat( f[:,[1,2]], f[:,[2,3]], f[:,[3,1]] ) # This holds all edges (but multiple copies)
  faces  = vcat( 1:nf, 1:nf, 1:nf )

  boundsf =  f[:,5]
  boundsf[ boundsf .> 0 ] = 0
  boundsf = -boundsf
  bounds  =  vcat( boundsf, boundsf, boundsf )
  ne = size(edges,1)

  boundsUni = fill( 0, ne, 2 )

  # sort in ascending order
  edges = sort( edges, 2 )

  # This index links to the first unique index in the array edges
  ix  = groupslices( edges, 1 )

  # Find unique vector index
  # Need to find an index that links from edges to the unique edges (without gaps)
  jx = Array{Int64}( size(ix) )
  jx[1] = 0
  lind = 0 # maximum unique index
  mind = 0 # maximum index in edges that had unique index
  ne = 0  # number of unique elements

  indUni = fill( false, size(ix,1) ) # index of unique elements

  for ii = 1:length(ix)
    temp    = ix[ii] - lind - 1
    temp2   = ix[ii] - mind - 1
    jx[ii]  = temp
    mind    = max(mind,ix[ii])

    if temp2 >= 0 # Found a unique index
      lind = lind + 1
      ne  = ne + 1
      indUni[ii] = true
      # add boundary information
      boundsUni[ii,1] = bounds[ii]
    end

    if temp2 < 0 # Need to link to the unique index
      jx[ii] = jx[ix[ii]]
      # add boundary information
      if bounds[ii] > 0
        if boundsUni[ ix[ii] ,1] > 0
          boundsUni[ ix[ii] ,2] = bounds[ii]
        else
          boundsUni[ ix[ii] ,1] = bounds[ii]
        end
      end
    end

  end

  e = fill( 0::Int64, ne, 4 ) # Last two are for boundary edges
                              # We need two because some faces share two
                              # boundaries

  e[:,1:2] = edges[indUni,:]

  # add boundary information
  e[:,3:4]   = sort( boundsUni[indUni,:], 2 )

  # make e2f, f2e
  e2f = Vector{Vector{Int64}}( ne )
  for jj in 1:ne
      e2f[jj] = fill( 0, 0 )
  end
  f2e = Matrix{Int64}( nf, 3 )
  for jj = 1:3*nf

    ind  = faces[jj]
    indf = Int64(ceil(jj / nf))

    inde = ix[jj] - jx[jj]

    append!( e2f[inde], ind )

    f2e[ ind, indf ] = inde

  end

  return e, e2f, f2e

end

"""
    genNEconnec3D( e::Matrix{Int64}, n::Int64 )

Generates node to edge connectivity given edge connectivity (`e`).
"""
function genNEconnec3D( e::Matrix{Int64}, n::Int64 )

    n2e = Vector{Vector{Int64}}(n)
    for ii in 1:n
        n2e[ii] = fill( 0, 0 )
    end

    for jj in 1:size(e,1)

        ii1 = e[jj,1]
        ii2 = e[jj,2]

        if isempty( n2e[ii1] .== jj ) || !any( n2e[ii1] .== jj )
            append!(n2e[ii1], jj)
        end
        if isempty( n2e[ii2] .== jj ) || !any( n2e[ii2] .== jj )
            append!(n2e[ii2], jj)
        end

    end

    return n2e

end
