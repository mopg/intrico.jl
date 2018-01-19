# ---------------------------------------------------------------------------- #
#
#   compJacob.jl
#
#   Several functions to compute jacobians of an element for both 2D and 3D
#
#   sterno
#   Spring 2018
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

"""
    compJacobFace( mesh::MeshF2D, master::luteos.Master2D, el::Int64, face::Int64 )

Returns Jacobian on on the `face` in element `el` for a 2D mesh.
"""
function compJacobFace( mesh::MeshF2D, master::luteos.Master2D, el::Int64, face::Int64 )

  indnod = 1
  rotdir = false

  if mesh.t2f[el,face] < 0
    # face DOES NOT follow counter-clockwise rotation
    indnod = 2
    rotdir = true
  end

  nod    = master.perm[:,face,indnod]

  ∂x₁∂ξ₁ = master.∇ϕ1d' * mesh.nodes[nod,1,el]
  ∂x₂∂ξ₁ = master.∇ϕ1d' * mesh.nodes[nod,2,el]

  p1d  = master.ϕ1d'  * mesh.nodes[nod,:,el]

  jac  = sqrt.( ∂x₁∂ξ₁.^2 + ∂x₂∂ξ₁.^2 )
  jcw  = master.gwts1d .* jac

  normal = [ ∂x₂∂ξ₁ -∂x₁∂ξ₁ ] ./ [jac jac]

  return (p1d, nod, normal)

end

"""
    compJacobFace( mesh::MeshF3D, master::luteos.Master3D, el::Int64, face::Int64 )

Returns Jacobian on the `face` in element `el` for a 3D mesh.
"""
function compJacobFace( mesh::MeshF3D, master::luteos.Master3D, el::Int64, face::Int64 )

  nod  = master.perm[ :, face, abs.(mesh.t2f[el,face+4]) ]

  p2d  = master.ϕ2D'  * mesh.nodes[nod,:,el]

  ∂x₁∂ξ₁ = master.∇ϕ2D[:,:,1]' * mesh.nodes[nod,1,el]
  ∂x₁∂ξ₂ = master.∇ϕ2D[:,:,2]' * mesh.nodes[nod,1,el]
  ∂x₂∂ξ₁ = master.∇ϕ2D[:,:,1]' * mesh.nodes[nod,2,el]
  ∂x₂∂ξ₂ = master.∇ϕ2D[:,:,2]' * mesh.nodes[nod,2,el]
  ∂x₃∂ξ₁ = master.∇ϕ2D[:,:,1]' * mesh.nodes[nod,3,el]
  ∂x₃∂ξ₂ = master.∇ϕ2D[:,:,2]' * mesh.nodes[nod,3,el]

  # cross product to find normal vector, normal = ∂x∂ξ₁ × ∂x∂ξ₂
  normal = hcat( ( ∂x₂∂ξ₁ .* ∂x₃∂ξ₂ - ∂x₂∂ξ₂ .* ∂x₃∂ξ₁ ),
                -( ∂x₁∂ξ₁ .* ∂x₃∂ξ₂ - ∂x₁∂ξ₂ .* ∂x₃∂ξ₁ ),
                 ( ∂x₁∂ξ₁ .* ∂x₂∂ξ₂ - ∂x₁∂ξ₂ .* ∂x₂∂ξ₁ ) )
  # normalize the normal vector
  jac      = sqrt.( normal[:,1].^2 + normal[:,2].^2 + normal[:,3].^2 )
  normal ./= jac * [1,1,1]'

  if mesh.t2f[el,face+4] < 0
    normal *= -1.0
  end

  return (p2d, nod, normal)

end
