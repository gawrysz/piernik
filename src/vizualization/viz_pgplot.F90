!
! PIERNIK Code Copyright (C) 2006 Michal Hanasz
!
!    This file is part of PIERNIK code.
!
!    PIERNIK is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    PIERNIK is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with PIERNIK.  If not, see <http://www.gnu.org/licenses/>.
!
!    Initial implementation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.h"
!>
!! \brief Routines to perform "live" visualization using pgplot
!!
!!  Currently to use PGPLOT you need to:
!!   1. set one of i{x,y,z} to positive value and zero to the others
!!   2. use only one-element vars array in problem.par
!!
!! \deprecated BEWARE: highly experimental
!<
module viz
! pulled by PGPLOT
   implicit none
   private
   public :: draw_me

   real(kind=4), parameter :: one = 1.0
   real(kind=4), parameter :: nul = 0.0

contains

   subroutine draw_me(f, fmin, fmax)

      implicit none

      real(kind=4), dimension(:,:), intent(in) :: f
      real(kind=4), intent(in)                 :: fmin, fmax
      integer                                  :: mxi, mxj, pgopen
      real(kind=4)                             :: contra, bright
      real(kind=4), dimension(6)               :: tr
      logical, save                            :: frun = .true.

      mxi = size(f,1)
      mxj = size(f,2)

      if (frun) then
         if (pgopen('/XS') < 1) stop
         call pgpage
      endif

      tr = [ nul, one, nul, nul, nul, one ]

      ! Clear the screen. Set up window and viewport.
      call setvp
      call pgwnad(nul, one+mxi, nul, one+mxj)

      ! Set up the color map.
      bright = 0.5
      contra  = 1.0
      call palett(2, contra, bright)

      ! Annotate the plot.
      if (frun) then
         call pgsch(real(0.6,4))
         call pgbox('bcntsi',nul,0,'bcntsiv',nul,0)
         call pgmtxt('b',real(3.0,4),one,one,'pixel number')
      endif

      ! Draw the map with PGIMAG.
      call pgimag(f,mxi,mxj,1,mxi,1,mxj,fmin,fmax,tr)

      ! Draw a wedge.
      if (frun) then
         call pgwedg('bi', real(4.0,4), real(5.0,4), fmin, fmax, 'pixel value')
         call pgsch(one)
      endif
      frun = .false.

!      call pgend
!-----------------------------------------------------------------------
   end subroutine draw_me

!>
!! \brief set a "palette" of colors in the range of color indices used by pgimag.
!<
   subroutine palett(ttype, contra, bright)
      implicit none
      integer      :: ttype
      real(kind=4) :: contra, bright
!
      real(kind=4) :: gl(2), gr(2), gg(2), gb(2)
      real(kind=4) :: rl(9), rr(9), rg(9), rb(9)
      real(kind=4) :: hl(5), hr(5), hg(5), hb(5)
      real(kind=4) :: wl(10), wr(10), wg(10), wb(10)
      real(kind=4) :: al(20), ar(20), ag(20), ab(20)
!
      data gl /0.0, 1.0/
      data gr /0.0, 1.0/
      data gg /0.0, 1.0/
      data gb /0.0, 1.0/
!
      data rl /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      data rr / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      data rg / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      data rb / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
!
      data hl /0.0, 0.2, 0.4, 0.6, 1.0/
      data hr /0.0, 0.5, 1.0, 1.0, 1.0/
      data hg /0.0, 0.0, 0.5, 1.0, 1.0/
      data hb /0.0, 0.0, 0.0, 0.3, 1.0/
!
      data wl /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      data wr /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      data wg /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      data wb /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
!
      data al /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, &
               0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      data ar /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, &
               0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      data ag /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8, &
               0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      data ab /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9, &
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
!
      if (ttype == 1) then
!        -- gray scale
         call pgctab(gl, gr, gg, gb, 2, contra, bright)
      else if (ttype == 2) then
!        -- rainbow
         call pgctab(rl, rr, rg, rb, 9, contra, bright)
      else if (ttype == 3) then
!        -- heat
         call pgctab(hl, hr, hg, hb, 5, contra, bright)
      else if (ttype == 4) then
!        -- weird iraf
         call pgctab(wl, wr, wg, wb, 10, contra, bright)
      else if (ttype == 5) then
!        -- aips
         call pgctab(al, ar, ag, ab, 20, contra, bright)
      endif
   end subroutine palett

!>
!! \brief set the viewport, allowing margins around the edge for annotation.
!! \details (this is similar in effect to pgvstd, but has different margins.)
!! the routine determines the view-surface size and allocates margins
!! as fractions of the minimum of width and height.
!<
   subroutine setvp

      implicit none

      real(kind=4)            :: d, vpx1, vpx2, vpy1, vpy2
      real(kind=4), parameter :: two = 2.0, five = 5.0, eight = 8.0, one40th = 0.025
!
      call pgsvp(nul, one, nul, one)
      call pgqvp(1, vpx1, vpx2, vpy1, vpy2)
      d = min(vpx2-vpx1, vpy2-vpy1) * one40th
      vpx1 = vpx1 + five *d
      vpx2 = vpx2 - two  *d
      vpy1 = vpy1 + eight*d
      vpy2 = vpy2 - two  *d
      call pgvsiz(vpx1, vpx2, vpy1, vpy2)

   end subroutine setvp

!>
!! \brief draw the enclosing rectangle of the subarray to be contoured,
!! applying the transformation tr.
!!
!! \details for a contour map, the corners are (i1,j1) and (i2,j2);
!! for a gray-scale map, they are (i1-0.5,j1-0.5), (i2+0.5, j2+0.5).
!<
   subroutine outlin(i1,i2,j1,j2,tr)

      implicit none

      integer      :: i1, i2, j1, j2
      real(kind=4) :: tr(6)
      integer      :: k
      real(kind=4) :: xw(5), yw(5), t
!
      xw(:) = real([i1, i1, i2, i2, i1], kind=4)
      yw(:) = real([j1, j2, j2, j1, j1], kind=4)
      do k=1,5
         t = xw(k)
         xw(k) = tr(1) + tr(2)*t + tr(3)*yw(k)
         yw(k) = tr(4) + tr(5)*t + tr(6)*yw(k)
      enddo
      call pgline(5,xw,yw)
   end subroutine outlin

end module viz
