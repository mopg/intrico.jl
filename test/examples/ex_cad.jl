using intrico
using egads

P = luteos.P1() # Polynomial order of solution

mesh   = luteos.Mesh3D( "cube", P, N = 4)
meshf  = MeshF( mesh )
master = luteos.Master3D( P )

lat = Lattice( 1.0, fill( 0.005, size(meshf.e,1) ), fill( 0.0, size(meshf.e,1) ) )

# writeCAD( meshf, lat, "test2.egads", eps = 0.0, Δeps = 0.5e-2 )
writeCAD( meshf, lat, "test2.stp", eps = 0.0, Δeps = 0.5e-2 )
