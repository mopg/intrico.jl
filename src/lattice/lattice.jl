# ---------------------------------------------------------------------------- #
#
#   lattice.jl
#
#   Lattice type
#
#   sterno
#   Spring 2018
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

"""
    Lattice

Lattice type:
Type for lattices. It holds the optimized result.
"""
struct Lattice

  vol::Float64            # Optimized volume

  ar::Vector{Float64}     # Areas of each element
  f::Vector{Float64}      # Forces in each element

end
