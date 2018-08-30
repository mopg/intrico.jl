# ---------------------------------------------------------------------------- #
#
#   intrico.jl
#
#   Geometry generation for lattice structures for AM
#
#   intrico
#   Spring 2018
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

__precompile__()

"""
    intrico

Julia package to generate CAD of lattice structure for additive manufacturing.

Max Opgenoord

Spring/Fall 2018
"""

module intrico

    using Requires

    using StaticArrays

    using divido

    # Mesh datastructures
    export MeshF, MeshF2D, MeshF3D
    include("mesh/meshF.jl")

    # Lattice datastructure
    export Lattice
    include("lattice/lattice.jl")

    # Generate CAD
    export genCAD, genSTL, genSTLnew, genSTLmodel
    include("cad/genSTL.jl")
    include("cad/genSTLmanu.jl")
    include("cad/genSTLmodel.jl")
    @require egads include("cad/writeCAD.jl")

end
