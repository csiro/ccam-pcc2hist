module gldata
   implicit none
!  Arrays for fields that need some extra processing and so 
!  can't be handled by readsave
   real, dimension(:,:),   allocatable, save :: psl, pmsl, zs, tsu, pwc, soilt
   real, dimension(:,:,:), allocatable, save :: u, v, t, q, rh, ql, qf
   real, dimension(:),     allocatable, save :: sig, dsig, zsoil
   real, dimension(:,:,:), allocatable, save :: tgg, wetfrac, wbice, snowvar, tscrn3hr
   real, dimension(:,:), allocatable, save :: taux, tauy
   real, allocatable, dimension(:,:,:) :: zstd
   integer, dimension(:,:), allocatable, save :: isflag
end module gldata
