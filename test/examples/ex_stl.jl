using intrico

srand(1234) # initialize random seed

mesh   = divido.Mesh3D( "cube", 1, N = 6, scl=2.)
meshf  = MeshF( mesh )

lat = Lattice( 1.0, fill( 0.005, size(meshf.e,1) ), fill( 0.0, size(meshf.e,1) ) )

lat.ar .*= rand( length(lat.ar) )

@time genSTL( meshf,  lat, "test.stl", n=15, boundflat = [1] );
