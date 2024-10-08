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

   use constants, only: ndims

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real                   :: d0, p0, r0, beta_cr, amp_cr
   real, dimension(ndims) :: c_exp, c_rot, b0, sn_pos


   namelist /PROBLEM_CONTROL/ d0, p0, b0, sn_pos, r0, beta_cr, amp_cr, c_exp, c_rot

contains
!-----------------------------------------------------------------------------
   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers
!-----------------------------------------------------------------------------
   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use dataio_pub, only: die, nh
      use domain,     only: dom
      use func,       only: operator(.equals.)
      use mpisetup,   only: rbuff, master, slave

      implicit none

      d0        = 1.0e5     !< density
      p0        = 1.0       !< pressure
      b0        = [ 0., 0., 0. ] !< Magnetic field
      sn_pos    = [ 0., 0., 0. ] !< position of blob
      r0        = 5.* minval(dom%L_(:)/dom%n_d(:), mask=dom%has_dir(:))  !< radius of the blob

      beta_cr   = 0.0       !< ambient level
      amp_cr    = 1.0       !< amplitude of the blob

      c_exp = [ 0.0, 0.0, 0.0 ]
      c_rot = [ 0.0, 0.0, 0.0 ]

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

         rbuff(1)     = d0
         rbuff(2)     = p0
         rbuff(3:5)   = b0
         rbuff(6:8)   = sn_pos
         rbuff(9)     = r0
         rbuff(10)    = beta_cr
         rbuff(11)    = amp_cr
         rbuff(12:14) = c_exp
         rbuff(15:17) = c_rot

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         d0           = rbuff(1)
         p0           = rbuff(2)
         b0           = rbuff(3:5)
         sn_pos       = rbuff(6:8)
         r0           = rbuff(9)
         beta_cr      = rbuff(10)
         amp_cr       = rbuff(11)
         c_exp        = rbuff(12:14)
         c_rot        = rbuff(15:17)

      endif

      if (r0 .equals. 0.0) call die("[initproblem:read_problem_par] r0 == 0")

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine problem_initial_conditions

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim, pi, I_ONE
      use cr_data,        only: icr_H1, icr_C12, cr_index
      use crhelpers,      only: div_v
      use dataio_pub,     only: warn
      use domain,         only: dom
      use fluidindex,     only: flind
      use fluidtypes,     only: component_fluid
      use func,           only: emag, ekin, operator(.equals.), operator(.notequals.)
      use grid_cont,      only: grid_container
      use initcosmicrays, only: iarr_crn, iarr_crs, gamma_cr_1, K_cr_paral, K_cr_perp

      implicit none

      class(component_fluid), pointer  :: fl
      integer                          :: i, j, k, icr, ipm, jpm, kpm
      real                             :: cs_iso, r, r2
      type(cg_list_element),  pointer  :: cgl
      type(grid_container),   pointer  :: cg

      fl => flind%ion

! Uniform equilibrium state

      cs_iso = sqrt(p0/d0)

      where (.not. dom%has_dir)
         b0 = 0.0  ! ignore B field in nonexistent direction
      endwhere

      if ((sum(b0**2) .equals. 0.) .and. (any(K_cr_paral(:) .notequals. 0.) .or. any(K_cr_perp(:) .notequals. 0.))) then
         call warn("[initproblem:problem_initial_conditions] No magnetic field is set, K_cr_* also have to be 0.")
         K_cr_paral(:) = 0.
         K_cr_perp(:)  = 0.
      endif

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         call cg%set_constant_b_field(b0)
         cg%u(fl%idn,RNG) = d0
         cg%u(fl%imx:fl%imz,RNG) = 0.0

         do k = lbound(cg%u, zdim+I_ONE), ubound(cg%u, zdim+I_ONE)
            do j = lbound(cg%u, ydim+I_ONE), ubound(cg%u, ydim+I_ONE)
               do i = lbound(cg%u, zdim+I_ONE), ubound(cg%u, xdim+I_ONE)

                  cg%u(fl%imx,i,j,k) = c_exp(xdim) * cg%u(fl%idn,i,j,k) * cg%x(i)
                  cg%u(fl%imy,i,j,k) = c_exp(ydim) * cg%u(fl%idn,i,j,k) * cg%y(j)
                  cg%u(fl%imz,i,j,k) = c_exp(zdim) * cg%u(fl%idn,i,j,k) * cg%z(k)

                  r = sqrt( (cg%x(i))**2 + (cg%y(j))**2 + (cg%z(k))**2 )

                  cg%u(fl%imx,i,j,k) = cg%u(fl%imx,i,j,k) + c_rot(zdim) * cg%u(fl%idn,i,j,k) * (-cg%y(j)/r) * sin(pi*r)
                  cg%u(fl%imy,i,j,k) = cg%u(fl%imy,i,j,k) + c_rot(zdim) * cg%u(fl%idn,i,j,k) * ( cg%x(i)/r) * sin(pi*r)
                  cg%u(fl%imz,i,j,k) = cg%u(fl%imz,i,j,k) + 0.0
#ifndef ISO
                  cg%u(fl%ien,i,j,k) = p0 / fl%gam_1 + &
                       &               ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k)) + &
                       &               emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
#endif /* !ISO */
               enddo
            enddo
         enddo

         call div_v(flind%ion%pos, cg)

#ifdef COSM_RAYS
         do icr = lbound(iarr_crs, 1), ubound(iarr_crs, 1)
            cg%u(iarr_crs(icr),RNG) =  beta_cr * fl%cs2 * cg%u(fl%idn,RNG) / gamma_cr_1
         enddo

! Explosions
         do icr = 1, flind%crn%all
            do k = cg%ks, cg%ke
               do j = cg%js, cg%je
                  do i = cg%is, cg%ie
                     do ipm = -1, 1
                        do jpm = -1, 1
                           do kpm = -1, 1

                              r2 = (cg%x(i) - sn_pos(xdim) + real(ipm) * dom%L_(xdim))**2 + &
                                 & (cg%y(j) - sn_pos(ydim) + real(jpm) * dom%L_(ydim))**2 + &
                                 & (cg%z(k) - sn_pos(zdim) + real(kpm) * dom%L_(zdim))**2
                              if (icr == cr_index(icr_H1)) then
                                 cg%u(iarr_crn(icr),i,j,k) = cg%u(iarr_crn(icr),i,j,k) + amp_cr*exp(-r2/r0**2)
                              elseif (icr == cr_index(icr_C12)) then
                                 cg%u(iarr_crn(icr),i,j,k) = cg%u(iarr_crn(icr),i,j,k) + amp_cr*0.1*exp(-r2/r0**2) ! BEWARE: magic number
                              else
                                 cg%u(iarr_crn(icr),i,j,k) = 0.0
                              endif

                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

#endif /* COSM_RAYS */

   end subroutine problem_initial_conditions

end module initproblem
