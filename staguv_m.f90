module staguv_m
   implicit none
   public :: unstaguv
contains
   subroutine unstaguv ( u_in, v_in, uout, vout )

!     N.B. staguv & unstaguv require genuine 2D arrays as input
!     and they need to be contiguous (safest to have in common!)

!     However in the routine use 1D arrays. These are copied at the start
!     and the end.

!     This version does the reshape by hand to get around the cherax and
!     bragg performance problems.

      use newmpar_m
      use indices_m

      real, intent(in), dimension(:,:) :: u_in, v_in
      real, intent(out), dimension(:,:) :: uout, vout

      real, dimension(:), allocatable :: u, v
      real, dimension(:), allocatable, target :: ave, diff
      real, dimension(:), pointer :: ua, va, ud, vd
      integer :: i, j

      allocate ( u(1:2*ifull), v(-ifull+1:ifull) )
      allocate ( ave(-ifull+1:ifull), diff(-ifull+1:ifull) )

!      u = reshape ( u_in, (/ 2*ifull /) )
!      v = reshape ( v_in, (/ 2*ifull /) )
      do j=1,2*jl
	 do i=1,il
	    u( i + (j-1)*il ) = u_in(i,j)
	 end do
      end do
      do j=-jl+1,jl
	 do i=1,il
	    v( i + (j-1)*il ) = v_in(i,j+jl) ! Need this here because didn't
	 end do                              ! give low bound in declaration.
      end do

      ua => ave(-ifull+1:ifull)
      va => ave
      ud => diff(-ifull+1:ifull)
      vd => diff

      ua(1:ifull) = 0.5 * ( u(1:ifull) + u(i_wu2) )
      ud(1:ifull) = u(1:ifull) - u(i_wu2)
      va(1:ifull) = 0.5 * ( v(1:ifull) + v(i_sv2) )
      vd(1:ifull) = v(1:ifull) - v(i_sv2)

!     Now combine the above to give cubic interpolated values
!     Code from here on identical in staguv & unstaguv

!      uout = reshape ( ua(1:ifull) - ( ud(ieu2) - ud(iwu2) ) / 16.0,    &
!                       (/ il, jl /) )
!      vout = reshape ( va(1:ifull) - ( vd(inv2) - vd(isv2) ) / 16.0,    &
!                       (/ il, jl /) )

!     Use ua as a temp
      ua(1:ifull) = ua(1:ifull) - ( ud(i_eu2) - ud(i_wu2) ) / 16.0
      va(1:ifull) = va(1:ifull) - ( vd(i_nv2) - vd(i_sv2) ) / 16.0
      do j=1,jl
	 do i=1,il
	    uout(i,j) = ua( i + (j-1)*il )
	    vout(i,j) = va( i + (j-1)*il )
	 end do
      end do

      deallocate ( u, v, ave, diff )
   end subroutine unstaguv

   subroutine staguv ( u_in, v_in, uout, vout )

!     N.B. staguv & unstaguv require genuine 2D arrays as input
!     and they need to be contiguous (safest to have in common!)

!     However in the routine use 1D arrays. These are copied at the start
!     and the end.

!     This version does the reshape by hand to get around the cherax and
!     bragg performance problems.

      use newmpar_m
      use indices_m

      real, intent(in), dimension(:,:) :: u_in, v_in
      real, intent(out), dimension(:,:) :: uout, vout

      real, dimension(:), allocatable :: u, v
      real, dimension(:), allocatable, target :: ave, diff
      real, dimension(:), pointer :: ua, va, ud, vd
      integer :: i, j

      allocate ( u(1:2*ifull), v(-ifull+1:ifull) )
      allocate ( ave(-ifull+1:ifull), diff(-ifull+1:ifull) )

!      u = reshape ( u_in, (/ 2*ifull /) )
!      v = reshape ( v_in, (/ 2*ifull /) )
      do j=1,2*jl
	 do i=1,il
	    u( i + (j-1)*il ) = u_in(i,j)
	 end do
      end do
      do j=-jl+1,jl
	 do i=1,il
	    v( i + (j-1)*il ) = v_in(i,j+jl) ! Need this here because didn't
	 end do                              ! give low bound in declaration.
      end do

      ua => ave(-ifull+1:ifull)
      va => ave
      ud => diff(-ifull+1:ifull)
      vd => diff

      ua(1:ifull) = 0.5 * ( u(i_eu2) + u(1:ifull) )
      ud(1:ifull) = u(i_eu2) - u(1:ifull)
      va(1:ifull) = 0.5 * ( v(i_nv2) + v(1:ifull) )
      vd(1:ifull) = v(i_nv2) - v(1:ifull)

!     Now combine the above to give cubic interpolated values
!     Code from here on identical in staguv & unstaguv

!      uout = reshape ( ua(1:ifull) - ( ud(ieu2) - ud(iwu2) ) / 16.0,    &
!                       (/ il, jl /) )
!      vout = reshape ( va(1:ifull) - ( vd(inv2) - vd(isv2) ) / 16.0,    &
!                       (/ il, jl /) )

!     Use ua as a temp
      ua(1:ifull) = ua(1:ifull) - ( ud(i_eu2) - ud(i_wu2) ) / 16.0
      va(1:ifull) = va(1:ifull) - ( vd(i_nv2) - vd(i_sv2) ) / 16.0
      do j=1,jl
	 do i=1,il
	    uout(i,j) = ua( i + (j-1)*il )
	    vout(i,j) = va( i + (j-1)*il )
	 end do
      end do

      deallocate ( u, v, ave, diff )
   end subroutine staguv

end module staguv_m

