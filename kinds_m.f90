module kinds_m
   ! Define kinds for single and double precision. Note that on the NEC
   ! single precision may be 64 bit. These parameters are only for 
   ! defining overloadable functions, not for selecting particular values.
   implicit none
   ! Default real kind
   integer, parameter :: sp = selected_real_kind(6)
!   real(kind=sp), parameter :: xxx = 1.0
!  Use default real here so that it works with Intel -r8 -i8
   real, parameter :: xxx = 1.0
   ! This specification works with standard IEEE and with NEC -ew 
   integer, parameter :: dp = selected_real_kind(precision(xxx)+1)
end module kinds_m

