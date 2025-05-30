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

   use constants, only: dsetnamelen

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   integer(kind=4)    :: norm_step
   real               :: d0, p0, bx0, by0, bz0, x0, y0, z0, r0, beta_cr, amp_cr, dtrig
#ifdef COSM_RAYS
   character(len=dsetnamelen), parameter :: aecr1_n = "aecr"
   integer, parameter :: icr = 1 !< Only first CR component is used in this test
#endif /* COSM_RAYS */

   namelist /PROBLEM_CONTROL/ d0, p0, bx0, by0, bz0, x0, y0, z0, r0, beta_cr, amp_cr, norm_step, dtrig

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use user_hooks,  only: late_initial_conditions
#ifdef COSM_RAYS
#ifdef HDF5
      use dataio_user, only: user_vars_hdf5
#endif /* HDF5 */
      use user_hooks,  only: problem_customize_solution, finalize_problem, problem_domain_update
#endif /* COSM_RAYS */

      implicit none

      late_initial_conditions    => cr_late_init
#ifdef COSM_RAYS
      problem_domain_update      => cr_dist_to_edge
      finalize_problem           => check_norm
      problem_customize_solution => check_norm_wrapper
#ifdef HDF5
      user_vars_hdf5             => crtest_analytic_ecr1
#endif /* HDF5 */
#endif /* COSM_RAYS */

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use bcast,          only: piernik_MPI_Bcast
      use constants,      only: xdim, ydim, zdim, I_ONE, I_TEN
      use dataio_pub,     only: die, nh
      use domain,         only: dom
      use func,           only: operator(.equals.)
      use mpisetup,       only: ibuff, rbuff, master, slave
#ifdef COSM_RAYS
      use constants,      only: AT_NO_B
      use cg_list_global, only: all_cg
#endif /* COSM_RAYS */

      implicit none

      d0           = 1.0e5     !< density
      p0           = 1.0       !< pressure
      bx0          =   0.      !< Magnetic field component x
      by0          =   0.      !< Magnetic field component y
      bz0          =   0.      !< Magnetic field component z
      x0           = 0.0       !< x-position of the blob
      y0           = 0.0       !< y-position of the blob
      z0           = 0.0       !< z-position of the blob
      r0           = 5.* minval(dom%L_(:)/dom%n_d) !< radius of the blob
      dtrig        = -0.01     !< fraction of top density used to trigger domain expansion

      beta_cr      = 0.0       !< ambient level
      amp_cr       = 1.0       !< amplitude of the blob

      norm_step    = I_TEN     !< how often to compute the norm (in steps)

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

         if (.not.dom%has_dir(xdim)) bx0 = 0. ! ignore B field in nonexistent direction to match the analytical solution
         if (.not.dom%has_dir(ydim)) by0 = 0.
         if (.not.dom%has_dir(zdim)) bz0 = 0.

         rbuff(1)  = d0
         rbuff(2)  = p0
         rbuff(3)  = bx0
         rbuff(4)  = by0
         rbuff(5)  = bz0
         rbuff(6)  = x0
         rbuff(7)  = y0
         rbuff(8)  = z0
         rbuff(9)  = r0
         rbuff(10) = beta_cr
         rbuff(11) = amp_cr
         rbuff(12) = dtrig

         ibuff(1)  = norm_step

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         d0        = rbuff(1)
         p0        = rbuff(2)
         bx0       = rbuff(3)
         by0       = rbuff(4)
         bz0       = rbuff(5)
         x0        = rbuff(6)
         y0        = rbuff(7)
         z0        = rbuff(8)
         r0        = rbuff(9)
         beta_cr   = rbuff(10)
         amp_cr    = rbuff(11)
         dtrig     = rbuff(12)

         norm_step = int(ibuff(1), kind=4)

      endif

      if (r0 .equals. 0.) call die("[initproblem:read_problem_par] r0 == 0")

#ifdef COSM_RAYS
      call all_cg%reg_var(aecr1_n, restart_mode = AT_NO_B)
#endif /* COSM_RAYS */

      if (norm_step <= 0) norm_step = huge(I_ONE)

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use fluidindex,     only: flind
      use fluidtypes,     only: component_fluid
      use func,           only: ekin, emag
      use grid_cont,      only: grid_container
#ifndef ISO
      use constants,      only: zdim
#endif /* !ISO */
#if defined(COSM_RAYS) || !defined(ISO)
      use constants,      only: xdim, ydim
#endif /* COSM_RAYS || !ISO */
#ifdef COSM_RAYS
      use constants,      only: HI
      use dataio_pub,     only: die
      use domain,         only: dom
      use initcosmicrays, only: gamma_cr_1, iarr_crs, ncrsp, ncrb
#endif /* COSM_RAYS */

      implicit none

      class(component_fluid), pointer :: fl
      !real                            :: cs_iso
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

#ifdef COSM_RAYS
      integer                         :: i, j, k, iecr
      real                            :: r2

      iecr = -1
      if (ncrsp+ncrb >= icr) then
         iecr = iarr_crs(icr)
      else
         call die("[initproblem:problem_initial_conditions] No CR components defined.")
      endif
#endif /* COSM_RAYS */

      fl => flind%ion

! Uniform equilibrium state

      !cs_iso = sqrt(p0/d0)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         call cg%set_constant_b_field([bx0, by0, bz0])
         cg%u(fl%idn,RNG) = d0
         cg%u(fl%imx:fl%imz,RNG) = 0.0

#ifndef ISO
         cg%u(fl%ien,RNG) = p0/fl%gam_1 + ekin(cg%u(fl%imx,RNG), cg%u(fl%imy,RNG), cg%u(fl%imz,RNG), cg%u(fl%idn,RNG)) + &
              &             emag(cg%b(xdim,RNG), cg%b(ydim,RNG), cg%b(zdim,RNG))
#endif /* !ISO */

#ifdef COSM_RAYS
         cg%u(iecr,RNG) = beta_cr * fl%cs2 * cg%u(fl%idn,RNG) / gamma_cr_1

! Explosions
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  r2 = (cg%x(i)-x0)**2+(cg%y(j)-y0)**2+(cg%z(k)-z0)**2
                  if (cg%x(i)> 2*x0-dom%edge(xdim, HI) .and. cg%y(j) > 2*y0-dom%edge(ydim, HI)) &
                       cg%u(iecr,i,j,k)= cg%u(iecr,i,j,k) + amp_cr*exp(-r2/r0**2)
               enddo
            enddo
         enddo
#endif /* COSM_RAYS */

         cgl => cgl%nxt
      enddo

#ifdef COSM_RAYS
      call check_norm
#endif /* COSM_RAYS */

   end subroutine problem_initial_conditions

!> \brief Performa late initialization of the cg added after domain expansion

   subroutine cr_late_init

      use cg_list,          only: cg_list_element
      use cg_list_dataop,   only: expanded_domain
      use constants,        only: psi_n
      use fluidindex,       only: flind
      use named_array_list, only: qna
#ifndef ISO
      use constants,        only: xdim, ydim, zdim
      use func,             only: ekin, emag
#endif /* !ISO */
#ifdef COSM_RAYS
      use initcosmicrays,   only: gamma_cr_1, iarr_crs
#endif /* COSM_RAYS */

      implicit none

      type(cg_list_element),  pointer :: cgl

      cgl => expanded_domain%first
      do while (associated(cgl))
         if (.not. cgl%cg%is_old) then
            associate (fl => flind%ion, cg => cgl%cg)
               call cgl%cg%set_constant_b_field([bx0, by0, bz0])
               cgl%cg%u(fl%idn,RNG) = d0
               cgl%cg%u(fl%imx:fl%imz,RNG) = 0.0
#ifndef ISO
               cgl%cg%u(fl%ien,RNG) = p0/fl%gam_1 + emag(cgl%cg%b(xdim,RNG), cgl%cg%b(ydim,RNG), cgl%cg%b(zdim,RNG)) + &
                    &                 ekin(cgl%cg%u(fl%imx,RNG), cgl%cg%u(fl%imy,RNG), cgl%cg%u(fl%imz,RNG), cgl%cg%u(fl%idn,RNG))
#endif /* !ISO */
#ifdef COSM_RAYS
               cgl%cg%u(iarr_crs(icr),RNG) = beta_cr * fl%cs2 * cgl%cg%u(fl%idn,RNG) / gamma_cr_1
#endif /* COSM_RAYS */
            end associate
         endif
         cgl => cgl%nxt
      enddo

      ! ToDo: move this call to some other place where it would be automagically called
      if (qna%exists(psi_n)) call expanded_domain%set_q_value(qna%ind(psi_n), 0.)

   end subroutine cr_late_init

!-----------------------------------------------------------------------------
#ifdef COSM_RAYS
   subroutine compute_analytic_ecr1

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use dataio_pub,       only: die
      use func,             only: operator(.equals.)
      use global,           only: t
      use grid_cont,        only: grid_container
      use initcosmicrays,   only: iarr_crs, ncrsp, ncrb, K_cr_paral, K_cr_perp
      use named_array_list, only: qna

      implicit none

      integer                        :: i, j, k, iecr
      real                           :: r_par2, r_perp2, delx, dely, delz, magb, ampt, r0_par2, r0_perp2, bxn, byn, bzn
      integer, parameter             :: icr = 1 !< Only first CR component
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      iecr = -1

      if (ncrsp+ncrb >= icr) then
         iecr = iarr_crs(icr)
      else
         call die("[initproblem:compute_analytic_ecr1] No CR components defined.")
      endif

      magb = sqrt(bx0**2 + by0**2 + bz0**2)
      if (magb > 0.) then
         bxn = bx0 / magb
         byn = by0 / magb
         bzn = bz0 / magb
      else
         bxn = 0.0
         byn = 0.0
         bzn = 0.0
      endif

      r0_perp2 = r0**2 + 4 * K_cr_perp(icr) * t
      r0_par2 = r0_perp2
      if (magb > 0.) r0_par2 = r0_par2 + 4 * K_cr_paral(icr) * t

      if ((r0_par2 .equals. 0.) .or. (r0_perp2 .equals. 0.)) call die("[initproblem:compute_analytic_ecr1] r0_par2 == 0. .or. r0_perp2 == 0.")

      ampt     = amp_cr * r0**2 / sqrt(r0_par2 * r0_perp2)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         if (.not. qna%exists(aecr1_n)) call die("[initproblem:compute_analytic_ecr1] cannot find aecr1")

         do k = cg%ks, cg%ke
            delz = cg%z(k) - z0
            do j = cg%js, cg%je
               dely = cg%y(j) - y0
               do i = cg%is, cg%ie
                  delx = cg%x(i) - x0

                  r_par2 = (bxn*delx + byn*dely + bzn*delz)**2 ! square of the distance form the center of the bump in direction parallel to the magnetic field
                  r_perp2 = delx**2 + dely**2 + delz**2 - r_par2
                  cg%q(qna%ind(aecr1_n))%arr(i, j, k) = ampt * exp( - r_par2/r0_par2 - r_perp2/r0_perp2)

               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine compute_analytic_ecr1

!-----------------------------------------------------------------------------

   subroutine check_norm_wrapper(forward)

      use global, only: nstep

      implicit none

      logical, intent(in) :: forward

      if (forward .and. mod(nstep, norm_step) == 0) call check_norm

   end subroutine check_norm_wrapper

!-----------------------------------------------------------------------------
   subroutine check_norm

      use allreduce,        only: piernik_MPI_Allreduce
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: pSUM, pMIN, pMAX, V_INFO
      use dataio_pub,       only: msg, die, printinfo
      use func,             only: operator(.notequals.)
      use grid_cont,        only: grid_container
      use initcosmicrays,   only: iarr_crs, ncrsp, ncrb
      use mpisetup,         only: master
      use named_array_list, only: qna

      implicit none

      integer                        :: i, j, k, iecr
      real, dimension(2)             :: norm, dev
      real                           :: crt
      integer, parameter             :: icr = 1 !< Only first CR component
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      iecr = -1

      if (ncrsp+ncrb >= icr) then
         iecr = iarr_crs(icr)
      else
         call die("[initproblem:check_norm] No CR components defined.")
      endif

      call compute_analytic_ecr1

      norm(:) = 0.
      dev(1) = huge(1.0)
      dev(2) = -dev(1)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         if (.not. qna%exists(aecr1_n)) call die("[initproblem:check_norm] cannot find aecr1")
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  crt = cg%q(qna%ind(aecr1_n))%arr(i, j, k)
                  norm(1) = norm(1) + (crt - cg%u(iecr, i, j, k))**2
                  norm(2) = norm(2) + crt**2
                  dev(1) = min(dev(1), (crt - cg%u(iecr, i, j, k)))
                  dev(2) = max(dev(2), (crt - cg%u(iecr, i, j, k)))
               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo

      call piernik_MPI_Allreduce(norm,   pSUM)
      call piernik_MPI_Allreduce(dev(1), pMIN)
      call piernik_MPI_Allreduce(dev(2), pMAX)

      if (master) then
         if (norm(2) .notequals. 0.) then
            write(msg,'(a,f12.5,a,2f12.5)')"[initproblem:check_norm] L2 error norm = ", sqrt(norm(1)/norm(2)), " min and max error = ", dev(1:2)
         else
            write(msg,'(a,2f12.5)')"[initproblem:check_norm] Cannot compute L2 error norm, min and max error = ", dev(1:2)
         endif
         call printinfo(msg, V_INFO)
      endif

   end subroutine check_norm

!-----------------------------------------------------------------------------

#ifdef HDF5
   subroutine crtest_analytic_ecr1(var, tab, ierrh, cg)

      use dataio_pub,       only: die
      use grid_cont,        only: grid_container
      use initcosmicrays,   only: iarr_crs
      use named_array_list, only: qna

      implicit none

      character(len=*),               intent(in)    :: var
      real, dimension(:,:,:),         intent(inout) :: tab
      integer,                        intent(inout) :: ierrh
      type(grid_container), pointer,  intent(in)    :: cg

      call compute_analytic_ecr1

      if (.not. qna%exists(aecr1_n)) call die("[initproblem:crtest_analytic_ecr1] cannot find aecr1")

      ierrh = 0
      select case (trim(var))
         case ("acr1")
            tab(:,:,:) = real(cg%q(qna%ind(aecr1_n))%span(cg%ijkse), kind(tab))
         case ("err1")
            tab(:,:,:) = cg%u(iarr_crs(1), RNG) - cg%q(qna%ind(aecr1_n))%span(cg%ijkse)
         case default
            ierrh = -1
      end select

   end subroutine crtest_analytic_ecr1
#endif /* HDF5 */

!> \brief Find how close is the clump to the external edges and call expansion routine if necessary

   subroutine cr_dist_to_edge

      use allreduce,      only: piernik_MPI_Allreduce
      use cg_leaves,      only: leaves
      use cg_expand_base, only: expand_base
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim, LO, HI, pMAX
      use domain,         only: dom
      use initcosmicrays, only: iarr_crs

      implicit none

      type(cg_list_element),  pointer   :: cgl
      real, dimension(xdim:zdim, LO:HI) :: ddist
      integer                           :: i
      integer, parameter                :: iprox = 2
      real                              :: cmax

      if (dtrig < 0.) return

      cmax = 0.
      cgl => leaves%first
      do while (associated(cgl))
         cmax = max(cmax, maxval(cgl%cg%u(iarr_crs(icr), cgl%cg%is:cgl%cg%ie, cgl%cg%js:cgl%cg%je, cgl%cg%ks:cgl%cg%ke), mask = cgl%cg%leafmap))
         cgl => cgl%nxt
      enddo
      call piernik_MPI_Allreduce(cmax, pMAX)

      ddist = huge(1.)
      cgl => leaves%first
      do while (associated(cgl))
         if (any(cgl%cg%ext_bnd)) then
            !> \todo roll it to a nested loop
            if (dom%has_dir(xdim)) then
               if (cgl%cg%ext_bnd(xdim, LO)) then
                  do i = cgl%cg%is, cgl%cg%ie
                     if (any(cgl%cg%u(iarr_crs(icr), i, :, :) > cmax*dtrig)) then
                        ddist(xdim, LO) = min(ddist(xdim, LO), (cgl%cg%x(i) - cgl%cg%fbnd(xdim, LO))/cgl%cg%dx)
                        exit
                     endif
                  enddo
               endif
               if (cgl%cg%ext_bnd(xdim, HI)) then
                  do i = cgl%cg%ie, cgl%cg%is, -1
                     if (any(cgl%cg%u(iarr_crs(icr), i, :, :) > cmax*dtrig)) then
                        ddist(xdim, HI) = min(ddist(xdim, HI), (cgl%cg%fbnd(xdim, HI) - cgl%cg%x(i))/cgl%cg%dx)
                        exit
                     endif
                  enddo
               endif
            endif

            if (dom%has_dir(ydim)) then
               if (cgl%cg%ext_bnd(ydim, LO)) then
                  do i = cgl%cg%js, cgl%cg%je
                     if (any(cgl%cg%u(iarr_crs(icr), :, i, :) > cmax*dtrig)) then
                        ddist(ydim, LO) = min(ddist(ydim, LO), (cgl%cg%y(i) - cgl%cg%fbnd(ydim, LO))/cgl%cg%dy)
                        exit
                     endif
                  enddo
               endif
               if (cgl%cg%ext_bnd(ydim, HI)) then
                  do i = cgl%cg%je, cgl%cg%js, -1
                     if (any(cgl%cg%u(iarr_crs(icr), :, i, :) > cmax*dtrig)) then
                        ddist(ydim, HI) = min(ddist(ydim, HI), (cgl%cg%fbnd(ydim, HI) - cgl%cg%y(i))/cgl%cg%dy)
                        exit
                     endif
                  enddo
               endif
            endif

            if (dom%has_dir(zdim)) then
               if (cgl%cg%ext_bnd(zdim, LO)) then
                  do i = cgl%cg%ks, cgl%cg%ke
                     if (any(cgl%cg%u(iarr_crs(icr), :, :, i) > cmax*dtrig)) then
                        ddist(zdim, LO) = min(ddist(zdim, LO), (cgl%cg%z(i) - cgl%cg%fbnd(zdim, LO))/cgl%cg%dz)
                        exit
                     endif
                  enddo
               endif
               if (cgl%cg%ext_bnd(zdim, HI)) then
                  do i = cgl%cg%ke, cgl%cg%ks, -1
                     if (any(cgl%cg%u(iarr_crs(icr), :, :, i) > cmax*dtrig)) then
                        ddist(zdim, HI) = min(ddist(zdim, HI), (cgl%cg%fbnd(zdim, HI) - cgl%cg%z(i))/cgl%cg%dz)
                        exit
                     endif
                  enddo
               endif
            endif
         endif
         cgl => cgl%nxt
      enddo

      call expand_base(ddist(:,:) < iprox)

   end subroutine cr_dist_to_edge
#endif /* COSM_RAYS */

end module initproblem
