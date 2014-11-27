module gldata
   implicit none
!  Arrays for fields that need some extra processing and so 
!  can't be handled by readsave
   real, dimension(:,:),   allocatable, save :: psl, zs, soilt
   real, dimension(:,:,:), allocatable, save :: u, v, t, q, rh, ql, qf
   real, dimension(:),     allocatable, save :: sig, dsig, zsoil
   real, dimension(:,:,:), allocatable, save :: tgg, wbice, snowvar
   real, allocatable, dimension(:,:,:) :: zstd
   real, allocatable, dimension(:,:,:) :: hstd
   real, allocatable, dimension(:,:,:) :: tmp3d
end module gldata
