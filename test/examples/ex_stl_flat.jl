using intrico

p = [ -1.0  0.0 0.0;
       0.0 -1.0 0.0;
       1.0 -0.1 0.0;
       1.0  1.0 1.0;
       0.0  0.0 0.0 ]

e = [ 1 2 0 0; # start horizontal z = 0
      2 3 0 0;
      3 4 0 0;
      4 1 0 0; # end horizontal z = 0
      1 5 0 0; # cross start
      2 5 0 0;
      3 5 0 0;
      4 5 0 0 ]

n2e = intrico.genNEconnec3D( e, size(p,1) )

meshf = MeshF3D( 3, size(p,1), p,
                 fill(0,0,0), fill(0,0,0),
                 fill(0,0,0), fill(0,0,0),
                 e, [ [0] ], [0 0],
                 fill(0.,0,0,0), n2e, ["bla"] )

ar = fill( 0.01, size(meshf.e,1) )
ar[1] = 0.0005
ar[4] = 0.0005
ar[5] = 0.0005
ar[6] = 0.0005
ar[8] = 0.0005

lat = Lattice( 1.0, ar , fill( 0.0, size(meshf.e,1) ) )

@time genSTL( meshf,  lat, "test_flat.stl", n=10 );

println(" ")
