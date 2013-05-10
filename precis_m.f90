module precis_m
   implicit none
   integer, public, parameter :: r4 = selected_real_kind(p= 6, r=37)
   integer, public, parameter :: r8 = selected_real_kind(p=15, r=307)
   integer, public, parameter :: rx = r8
end module precis_m
