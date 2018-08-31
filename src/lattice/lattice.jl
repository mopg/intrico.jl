# ---------------------------------------------------------------------------- #
#
#   lattice.jl
#
#   Lattice type
#
#   intrico
#   Spring 2018
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

# NOTE: Ideally this information would be combined with MeshF, but typically
#       the frame structure is decided on before the actual cross-sectional areas
#       of the lattice. If Lattice would be combined with MeshF, this information
#       needs to be known at the same time, because those types are both immutable
#       (for instance the volume of the lattice could not be updated). Therefore,
#       it is more appropriate to keep this information separate in two different
#       structs.

"""
    Lattice

Lattice type:
Type for lattices. It holds the optimized result.
"""
struct Lattice

    vol::Float64            # Optimized volume

    ar::Vector{Float64}     # Cross-sectional areas of each strut
    f::Vector{Float64}      # Forces in each strut

end
