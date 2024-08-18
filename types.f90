module types

  implicit none

  type vet1

     character(len=50) :: atom, resid, ident, code
     real ::  x, y, z
     integer :: n_atom, n_resid

  end type vet1

  type vet2

     real :: x, y, z

  end type vet2

  type vet3

     real ::  rho, phi, theta
     
  end type vet3

end module types
