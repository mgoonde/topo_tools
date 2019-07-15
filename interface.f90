module f90nautyinterf

 integer :: me

interface

subroutine c_ffnauty(n,matrix,vector1,vector2,vector3,hv1,hv2,hv3) BIND(C, name="ffnauty")
use ISO_C_BINDING
integer(kind=C_INT) :: n
integer(kind=C_INT) :: hv1,hv2,hv3
integer(kind=C_INT), dimension(n,n) :: matrix
integer(kind=C_INT), dimension(n) :: vector1, vector2, vector3
end subroutine

end interface

end module
