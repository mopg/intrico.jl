# ---------------------------------------------------------------------------- #
#
#   genSTLmodel.jl
#
#   Write STL file for the model, based on just the mesh
#
#   intrico
#   Spring 2018
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

function genSTLmodel( mesh::Union{divido.Mesh3D,MeshF3D}, flnameBase::String )

    if length(flnameBase) > 4 && flnameBase[end-3:end] == ".stl"
        flnameBase = flnameBase[1:end-4]
    end

    nbound = length( mesh.boundtags )

    for bb in 1:nbound

        nface = count( mesh.fb[:,2] == bb )

        # open STL file
        fid = open( string(flnameBase, "_", mesh.boundtags[bb], ".stl" ), "w" )

        # write header
        for jj in 1:80
            write(fid,' ')
        end
        write(fid,UInt32( nface ))

        # generate all faces (and write)
        for jj in 1:size( mesh.fb, 1 )

            if mesh.fb[jj,2] != bb
                continue
            end

            indf = mesh.fb[jj,1]

            vert = SVector{3, SVector{3,Float64} }( SVector{3}( mesh.p[ mesh.f[indf,1], : ] ),
                                                    SVector{3}( mesh.p[ mesh.f[indf,2], : ] ),
                                                    SVector{3}( mesh.p[ mesh.f[indf,3], : ] ) )

            writeFacetSTLB( vert, fid )

        end

        # close STL
        close( fid )

    end

end
