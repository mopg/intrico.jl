using intrico
using egads

mesh   = divido.Mesh3D( "cube", 1, N = 4)
meshf  = MeshF( mesh )

lat = Lattice( 1.0, fill( 0.005, size(meshf.e,1) ), fill( 0.0, size(meshf.e,1) ) )

genCAD( meshf, lat, "test2.stp", eps = 0.0, Δeps = 0.5e-2 )
