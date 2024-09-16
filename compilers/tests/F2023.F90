! Test for avalability of Fortran 2023 constructs.
!
! Will require similar approach as compilers/tests/F2018.F90 .

program F2023

   implicit none

   integer, dimension(5, 5) :: A

   A = 0

   ! According to https://wg5-fortran.org/N2201-N2250/N2212.pdf the following should be equivalent to A(3, 5)

   A(@[3,5]) = 1

   ! This way we can simplify some expressions in Piernik by using expressions with lbound() or other arrays for indexing directly.

   write(*, '(a,i0,a)')"count(A /= 0): ", count(A /= 0), (A(3,5) /= 0 ? "Yes" : "No")

end program F2023
