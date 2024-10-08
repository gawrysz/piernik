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

module initproblem

   use constants,    only: cbuff_len

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   character(len=cbuff_len) :: fnoise
   real                     :: rhog, eps, amp, kx, kz
   logical                  :: linear

   namelist /PROBLEM_CONTROL/  rhog, eps, amp, fnoise, kx,kz, linear

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use dataio_pub, only: nh
      use mpisetup,   only: rbuff, cbuff, lbuff, master, slave

      implicit none

      linear = .false.
      rhog    = 10.0
      eps     =  1.0
      amp    =  0.0
      kx     =  30.0
      kz     =  30.0

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PROBLEM_CONTROL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL")
         read(nh%cmdl_nml,nml=PROBLEM_CONTROL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         call nh%compare_namelist()

         cbuff(1) =  fnoise

         rbuff(1) = rhog
         rbuff(2) = eps
         rbuff(3) = amp
         rbuff(4) = kx
         rbuff(5) = kz

         lbuff(1) = linear

      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         fnoise       = cbuff(1)

         rhog         = rbuff(1)
         eps          = rbuff(2)
         amp          = rbuff(3)
         kx           = rbuff(4)
         kz           = rbuff(5)

         linear       = lbuff(1)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_leaves,    only: leaves
      use cg_list,      only: cg_list_element
      use constants,    only: dpi, xdim, ydim, zdim, LO, HI
      use dataio_pub,   only: msg, printinfo, run_id
      use domain,       only: dom
      use fluidindex,   only: flind
      use fluidtypes,   only: component_fluid
      use grid_cont,    only: grid_container
      use interactions, only: dragc_gas_dust
      use mpisetup,     only: proc
      use shear,        only: omega
#ifdef SHEAR
      use shear,        only: eta_gas, csvk
#endif /* SHEAR */
      implicit none

      real                                          :: rcx, rcy, ux, uy, wx, wy, taus, eta, vk, beta !, inv
      integer                                       :: i, j, k, n, clock
      real(kind=4), dimension(:,:,:,:), allocatable :: noise
      integer, dimension(:), allocatable            :: seed
      complex(kind=8), dimension(7)                 :: coeff
      class(component_fluid), pointer               :: dst, neu
      type(cg_list_element),  pointer               :: cgl
      type(grid_container),   pointer               :: cg

#ifdef DUST
      dst => flind%dst
#else /* !DUST */
      call warn("[initproblem]: Dust fluid not initialized. I hope you know what you are doing!"
#endif /* !DUST */
#ifdef NEUTRAL
      neu => flind%neu
#else /* !NEUTRAL */
      call warn("[initproblem]: Neutral fluid not initialized. I hope you know what you are doing!"
#endif /* !NEUTRAL */

      if (run_id == 'lnA') then
         call printinfo("Lin A")
         coeff(4) = (-0.1691398, 0.0361553 ) ! u_x
         coeff(5) = ( 0.1336704, 0.0591695 ) ! u_y
         coeff(6) = ( 0.1691389,-0.0361555 ) ! u_z
         coeff(7) = ( 0.0000224,+0.0000212 ) ! gas dens
         coeff(1) = (-0.1398623, 0.0372951 ) ! w_x
         coeff(2) = ( 0.1305628, 0.0640574 ) ! w_y
         coeff(3) = ( 0.1639549,-0.0233277 ) ! w_z
         kx = dpi/dom%L_(xdim)
         kz = dpi/dom%L_(zdim)
      else if (run_id == 'lnB') then
         call printinfo("Lin B")
         coeff(4) = (-0.0174121,-0.2770347 ) ! u_x
         coeff(5) = ( 0.2767976,-0.0187568 ) ! u_y
         coeff(6) = ( 0.0174130, 0.2770423 ) ! u_z
         coeff(7) = (-0.0000067,-0.0000691 ) ! gas dens
         coeff(1) = ( 0.0462916,-0.2743072 ) ! w_x
         coeff(2) = ( 0.2739304, 0.0039293 ) ! w_y
         coeff(3) = ( 0.0083263, 0.2768866 ) ! w_z
         kx = dpi/dom%L_(xdim)
         kz = dpi/dom%L_(zdim)
      else
         kx = dpi/dom%L_(xdim)
         kz = dpi/dom%L_(zdim)
         coeff(:) = ( 0.0, 0.0 )
      endif

      call random_seed(size=n)
      allocate(seed(n))
      call system_clock(count=clock)
      seed = clock*proc + 37 * [ (i-1, i = 1, n) ]
      call random_seed(put=seed)
      deallocate(seed)

      taus = 1./dragc_gas_dust
      vk   = neu%cs/csvk
      eta  = eta_gas
!      beta = 2.0*omega*eta*vk
!      inv  = 1./(4.0*omega**2*taus**2 + (eps+1.0)**2)
!      ux =  beta*eps*taus * inv
!      wx = -beta*    taus * inv
!      uy = -0.5*beta/omega * (4.0*omega**2*taus**2 + eps + 1.0) * inv
!      wy = -0.5*beta/omega * (                       eps + 1.0) * inv

      beta = 2.0*taus*eta*vk
      wx  = -beta / ( (1.0+eps)**2 + taus**2 )
      ux  = -eps*wx
      wy  = (1+eps)/(2.0*taus) * wx
      uy  = (1.0 + eps + taus**2)/ (2.0*taus) * wx

      write(msg,*) 'kx = ',kx,' kz = ',kz
      call printinfo(msg)
      write(msg,*) 'ux = ',ux,' uy = ',uy
      call printinfo(msg)
      write(msg,*) 'wx = ',wx,' wy = ',wy
      call printinfo(msg)
      write(msg,*) '\eta vk / \Omega = ', eta_gas * neu%cs / csvk / omega
      call printinfo(msg)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         allocate(noise(3, cg%n_(xdim), cg%n_(ydim), cg%n_(zdim)))

         do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
            rcx = cg%x(i)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               rcy = cg%y(j)
               do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
#ifdef NEUTRAL
                  cg%u(neu%idn,i,j,k) = rhog
                  cg%u(neu%imx,i,j,k) = ux * rhog
                  cg%u(neu%imy,i,j,k) = uy * rhog
                  cg%u(neu%imz,i,j,k) = 0.0
#ifndef ISO
                  cg%u(neu%ien,i,j,:) = 1.0/(neu%gam-1.0)
#endif /* !ISO */
#endif /* NEUTRAL */
#ifdef DUST
                  cg%u(dst%idn,i,j,k) = eps*rhog
                  cg%u(dst%imx,i,j,k) = wx * eps*rhog
                  cg%u(dst%imy,i,j,k) = wy * eps*rhog
                  cg%u(dst%imz,i,j,k) = 0.0

! Linear test
                  if (linear) then
!               cg%u(dst%idn,i,j,k) =  cg%u(dst%idn,i,j,k) + amp*eps*sin(kz*cg%z(k))*cos(kx*cg%x(i))
                     cg%u(dst%idn,i,j,k) =  cg%u(dst%idn,i,j,k) + amp*eps*cos(kz*cg%z(k))*cos(kx*cg%x(i))
! ...                cg%u(dst%idn,i,j,k) =  cg%u(dst%idn,i,j,k) + amp*eps*cos(kx*x(i))*cos(kz*cg%z(k))
! B              cg%u(dst%idn,i,j,k) =  cg%u(dst%idn,i,j,k) + amp * eps *&
! B                 ( real(coeff(7))*cos(kx*cg%x(i)) - &
! B                  aimag(coeff(7))*sin(kx*cg%x(i))) * cos(kz*z(k))
                     cg%u(dst%imx,i,j,k) =  cg%u(dst%imx,i,j,k) + eta*vk*amp * &
                          ( real(coeff(1))*cos(kx*cg%x(i)) - &
                          aimag(coeff(1))*sin(kx*cg%x(i))) * cos(kz*cg%z(k))
                     cg%u(dst%imy,i,j,k) =  cg%u(dst%imy,i,j,k) + eta*vk*amp * &
                          ( real(coeff(2))*cos(kx*cg%x(i)) - &
                          aimag(coeff(2))*sin(kx*cg%x(i))) * cos(kz*cg%z(k))
                     cg%u(dst%imz,i,j,k) =  cg%u(dst%imz,i,j,k) + eta*vk*(-amp) * &
                          (aimag(coeff(3))*cos(kx*cg%x(i)) + &
                          real(coeff(3))*sin(kx*cg%x(i))) * sin(kz*cg%z(k))
                     cg%u(neu%imx,i,j,k) =  cg%u(neu%imx,i,j,k) + eta*vk*amp * &
                          ( real(coeff(4))*cos(kx*cg%x(i)) - &
                          aimag(coeff(4))*sin(kx*cg%x(i))) * cos(kz*cg%z(k))
                     cg%u(neu%imy,i,j,k) =  cg%u(neu%imy,i,j,k) + eta*vk*amp * &
                          ( real(coeff(5))*cos(kx*cg%x(i)) - &
                          aimag(coeff(5))*sin(kx*cg%x(i))) * cos(kz*cg%z(k))
                     cg%u(neu%imz,i,j,k) =  cg%u(neu%imz,i,j,k) + eta*vk*(-amp) * &
                          (aimag(coeff(6))*cos(kx*cg%x(i)) + &
                          real(coeff(6))*sin(kx*cg%x(i))) * sin(kz*cg%z(k))
!               cg%u(neu%idn,i,j,k) =  cg%u(neu%idn,i,j,k) + amp * &
!                  ( real(coeff(7))*cos(kx*cg%x(i)) - &
!                   aimag(coeff(7))*sin(kx*cg%x(i))) * cos(kz*cg%z(k))
                     cg%u(neu%idn,i,j,k) =  cg%u(neu%idn,i,j,k) + (eta*vk)**2 * amp * &
                          ( real(coeff(7))*cos(kx*cg%x(i)) - &
                          aimag(coeff(7))*sin(kx*cg%x(i))) * cos(kz*cg%z(k))
                  endif
!-------

#endif /* DUST */
               enddo
            enddo
         enddo
         if (.not.linear) then
            call random_number(noise)
            cg%u(dst%imx,:,:,:) = cg%u(dst%imx,:,:,:) +amp -2.0*amp*noise(1,:,:,:) * cg%u(dst%idn,:,:,:)
            cg%u(dst%imy,:,:,:) = cg%u(dst%imy,:,:,:) +amp -2.0*amp*noise(2,:,:,:) * cg%u(dst%idn,:,:,:)
            cg%u(dst%imz,:,:,:) = cg%u(dst%imz,:,:,:) +amp -2.0*amp*noise(3,:,:,:) * cg%u(dst%idn,:,:,:)
         endif

         deallocate(noise)

         cgl => cgl%nxt
      enddo

      write(msg,*) 'linear = ',linear
      call printinfo(msg)

   end subroutine problem_initial_conditions

!-----------------------------------------------------------------------------------------------------------------------------------

#if 0
! Currently unused

   subroutine compare

      use constants,  only: dpi, xdim, zdim
      use dataio_pub, only: run_id
      use domain,     only: dom

      implicit none

      complex(kind=8), dimension(8) :: coeff

      if (run_id == 'lnA') then
         coeff(4) = (-0.1691398, 0.0361553 ) ! u_x
         coeff(5) = ( 0.1336704, 0.0591695 ) ! u_y
         coeff(6) = ( 0.1691389,-0.0361555 ) ! u_z
         coeff(7) = ( 0.0000224,+0.0000212 ) ! gas dens
         coeff(1) = (-0.1398623, 0.0372951 ) ! w_x
         coeff(2) = ( 0.1305628, 0.0640574 ) ! w_y
         coeff(3) = ( 0.1639549,-0.0233277 ) ! w_z

         coeff(8) = (-0.3480127, 0.4190204 ) ! omega
      else if (run_id == 'lnB') then
         coeff(4) = (-0.0174121,-0.2770347 ) ! u_x
         coeff(5) = ( 0.2767976,-0.0187568 ) ! u_y
         coeff(6) = ( 0.0174130, 0.2770423 ) ! u_z
         coeff(7) = (-0.0000067,-0.0000691 ) ! gas dens
         coeff(1) = ( 0.0462916,-0.2743072 ) ! w_x
         coeff(2) = ( 0.2739304, 0.0039293 ) ! w_y
         coeff(3) = ( 0.0083263, 0.2768866 ) ! w_z

         coeff(8) = ( 0.4998786, 0.0154764 ) ! omega
      endif

      kx = dpi/dom%L_(xdim)
      kz = dpi/dom%L_(zdim)

   end subroutine compare
#endif /* 0 */

end module initproblem
