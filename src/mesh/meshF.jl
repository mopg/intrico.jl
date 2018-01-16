# ---------------------------------------------------------------------------- #
#
#   meshF.jl
#
#   Abstract meshF type
#   This allows for writing an n-dimensional methods for lattice optimization
#
#   sterno
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

using GroupSlices

"""
    MeshF

MeshF abstract type:
Overarching abstract type for mesh types (2D and 3D).
"""
abstract type MeshF

end

include("meshF2D.jl")
include("meshF3D.jl")
