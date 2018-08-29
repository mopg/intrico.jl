using intrico

srand(1234)

P = luteos.P1() # Polynomial order of solution

# mesh   = luteos.Mesh3D( "four", P, N = 2)
mesh   = luteos.Mesh3D( "cube", P, N = 2)
meshf  = MeshF( mesh )
master = luteos.Master3D( P )

# mesh.p[14,1] += -0.15
# mesh.p[14,2] +=  0.15
# mesh.p[14,3] += -0.15

lat = Lattice( 1.0, fill( 0.005, size(meshf.e,1) ), fill( 0.0, size(meshf.e,1) ) )

lat.ar .*= rand( length(lat.ar) )

# @time genSTL(    meshf,  lat, "test.stl" )
@time genSTLnew( meshf,  lat, "test_new.stl", n=15, boundflat = [1] )
# @code_warntype genSTLnew( meshf,  lat, "test_new.stl", n=10 )
# genSTLnew( meshf,  lat, "test_new.stl", n=10 )
# edgWrite = Vector{Vector{Bool}}( size(meshf.e,1) )
# for jj in 1:size(meshf.e,1)
#     edgWrite[jj] = [true,false]
# end
# edgWrite[16][2] = true
# edgWrite[16][1] = false
# @time genSTLnew( meshf,  lat, ["test_new_part1.stl","test_new_part2.stl"], edgWrite, n=10 )

# @time genSTLmodel( meshf,  "test_model" )
