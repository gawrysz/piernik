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
!! \brief Computation of Cosmic Ray sources and pcr gradient pcr
!!
!!
!<
module sourcecosmicrays
! pulled by COSM_RAYS
   implicit none

   private
   public :: src_gpcr_exec
#ifdef COSM_RAYS_SOURCES
   public :: src_crn_exec
#endif /* COSM_RAYS_SOURCES */

contains

!>
!! \brief Computation of Cosmic ray pressure gradient and pcr div v
!<
   subroutine src_gpcr(uu, nn, decr, grad_pcr, sweep, i1, i2, cg)

      use crhelpers,      only: set_div_v1d
      use domain,         only: dom
      use fluidindex,     only: flind
      use grid_cont,      only: grid_container
      use initcosmicrays, only: cr_active, gpcr_essential
#ifdef COSM_RAY_ELECTRONS
      use initcosmicrays, only: iarr_crn, gamma_crn
#else /* !COSM_RAY_ELECTRONS */
      use initcosmicrays, only: iarr_crs, gamma_crs
#endif /* !COSM_RAY_ELECTRONS */

      implicit none

      integer(kind=4),                    intent(in)  :: nn                 !< array size
      real, dimension(nn, flind%all),     intent(in)  :: uu                 !< vector of conservative variables
      real, dimension(nn),                intent(out) :: grad_pcr
#ifdef COSM_RAY_ELECTRONS
      real, dimension(nn, flind%crn%all), intent(out) :: decr
#else /* !COSM_RAY_ELECTRONS */
      real, dimension(nn, flind%crs%all), intent(out) :: decr
#endif /* !COSM_RAY_ELECTRONS */
      integer(kind=4),                    intent(in)  :: sweep              !< direction (x, y or z) we are doing calculations for
      integer,                            intent(in)  :: i1                 !< coordinate of sweep in the 1st remaining direction
      integer,                            intent(in)  :: i2                 !< coordinate of sweep in the 2nd remaining direction
      type(grid_container), pointer,      intent(in)  :: cg                 !< current grid piece
      real, dimension(:), pointer                     :: divv               !< vector of velocity divergence used in cosmic ray advection
      integer                                         :: icr, jcr

      call set_div_v1d(divv, sweep, i1, i2, cg)
#ifdef COSM_RAY_ELECTRONS
      do icr = 1, flind%crn%all
         decr(:, icr)      = -1. / real(dom%eff_dim) * (gamma_crn(icr)-1.0) * uu(:, iarr_crn(icr))*divv(:)  !< as above, but only for crn
      enddo
#else /* !COSM_RAY_ELECTRONS */
      do icr = 1, flind%crs%all
         ! 1/eff_dim is because we compute the p_cr*dv in every sweep (3 times in 3D, twice in 2D and once in 1D experiments)
         decr(:, icr)      = -1. / real(dom%eff_dim) * (gamma_crs(icr)-1.0) * uu(:, iarr_crs(icr))*divv(:)
      enddo
#endif /* !COSM_RAY_ELECTRONS */
      !< gpcr_essential includes electrons only if COSM_RAY_ELECTRONS not defined and cre_gpcr_ess = .true.
      grad_pcr(:) = 0.0
      do icr = 1, size(gpcr_essential)
         jcr = gpcr_essential(icr)
#ifdef COSM_RAY_ELECTRONS
         grad_pcr(2:nn-1) = grad_pcr(2:nn-1) + cr_active*(gamma_crn(jcr)-1.)*(uu(1:nn-2, iarr_crn(jcr)) - uu(3:nn, iarr_crn(jcr)))/(2.*cg%dl(sweep))
#else /* !COSM_RAY_ELECTRONS */
         grad_pcr(2:nn-1) = grad_pcr(2:nn-1) + cr_active*(gamma_crs(jcr)-1.)*(uu(1:nn-2, iarr_crs(jcr)) - uu(3:nn, iarr_crs(jcr)))/(2.*cg%dl(sweep))
#endif /* !COSM_RAY_ELECTRONS */
      enddo
      grad_pcr(1:2) = 0.0 ; grad_pcr(nn-1:nn) = 0.0

   end subroutine src_gpcr

!>
!! \brief Computation of Cosmic ray pressure gradient and pcr div v
!<
   subroutine src_gpcr_exec(uu, nn, usrc, sweep, i1, i2, cg, vx)

      use fluidindex,     only: flind, iarr_all_mx, iarr_all_en
      use grid_cont,      only: grid_container
#ifdef COSM_RAY_ELECTRONS
      use cresp_crspectrum, only: src_gpcresp
      use initcosmicrays,   only: iarr_crn, iarr_cre_e
#else /* !COSM_RAY_ELECTRONS */
      use initcosmicrays, only: iarr_crs
#endif /* !COSM_RAY_ELECTRONS */

      implicit none

      integer(kind=4),                intent(in)  :: nn                 !< array size
      real, dimension(nn, flind%all), intent(in)  :: uu                 !< vector of conservative variables
      integer(kind=4),                intent(in)  :: sweep              !< direction (x, y or z) we are doing calculations for
      integer,                        intent(in)  :: i1                 !< coordinate of sweep in the 1st remaining direction
      integer,                        intent(in)  :: i2                 !< coordinate of sweep in the 2nd remaining direction
      type(grid_container), pointer,  intent(in)  :: cg                 !< current grid piece
      real, dimension(:,:), pointer,  intent(in)  :: vx
      real, dimension(nn, flind%all), intent(out) :: usrc               !< u array update component for sources
!locals
      real, dimension(nn)                         :: grad_pcr
#ifdef COSM_RAY_ELECTRONS
      real, dimension(nn)                         :: grad_pcr_cresp
      real, dimension(nn, flind%crn%all)          :: decr
#else /* !COSM_RAY_ELECTRONS */
      real, dimension(nn, flind%crs%all)          :: decr
#endif /* !COSM_RAY_ELECTRONS */
      logical                                     :: full_dim

      full_dim = nn > 1

      usrc = 0.0
      if (.not.full_dim) return

      call src_gpcr(uu, nn, decr, grad_pcr, sweep, i1, i2, cg)
#ifdef COSM_RAY_ELECTRONS
      call src_gpcresp(uu(:,iarr_cre_e(:)), nn, cg%dl(sweep), grad_pcr_cresp)         !< cg%dl(sweep) = dx, contribution due to pressure acted upon spectrum components in CRESP via div_v
      usrc(:, iarr_crn(:)) = decr(:,:)
#else /* !COSM_RAY_ELECTRONS */
      usrc(:, iarr_crs(:)) = decr(:,:)
#endif /* !COSM_RAY_ELECTRONS */
      usrc(:, iarr_all_mx(flind%ion%pos)) = grad_pcr
#ifdef ISO
      return
#endif /* ISO */
      usrc(:, iarr_all_en(flind%ion%pos)) = vx(:, flind%ion%pos) * grad_pcr
#ifdef COSM_RAY_ELECTRONS
      usrc(:, iarr_all_en(flind%ion%pos)) = usrc(:, iarr_all_en(flind%ion%pos)) + vx(:, flind%ion%pos) * grad_pcr_cresp !< BEWARE - check it
#endif /* COSM_RAY_ELECTRONS */

   end subroutine src_gpcr_exec

!==========================================================================================

#ifdef COSM_RAYS_SOURCES
!>
!! \brief Computation of Cosmic ray particles spallation and decay
!! \deprecated BEWARE: several lines in this routine break unit consistency, move it to units.F90 and use scaling
!<
   subroutine src_crn(uu, n, decrn, rk_coeff)

      use cr_data,       only: eCRSP, cr_table, cr_tau, cr_sigma, icr_Be10, icrH, icrL
      use fluids_pub,    only: has_ion, has_neu
      use fluidindex,    only: flind

      implicit none

      integer(kind=4),                   intent(in)  :: n
      real, dimension(n, flind%all),     intent(in)  :: uu
      real,                              intent(in)  :: rk_coeff   !< coeffecient used in RK step, while computing source term
      real, dimension(n, flind%crn%all), intent(out) :: decrn

! locals
      real, dimension(n)      :: dgas
      real, dimension(n)      :: dcr
      real,         parameter :: gamma_lor = 10.0
      real(kind=8), parameter :: speed_of_light = 3e10*1e6*365.*24.*60.*60. !< cm/Myr \deprecated BEWARE: this line breaks unit consistency, move it to units.F90 and use scaling
      real,         parameter :: ndim = 2.0
      real,         parameter :: c_n = speed_of_light / ndim
      real,         parameter :: gn = 1.0 / gamma_lor / ndim

      integer               :: i, j

      dgas = 0.0
      if (has_ion) dgas = dgas + uu(:, flind%ion%idn)
      if (has_neu) dgas = dgas + uu(:, flind%neu%idn)
      dgas = c_n*dgas

      decrn(:,:) = 0.0

      if (eCRSP(icr_Be10)) then
         decrn(:, cr_table(icr_Be10)) = decrn(:, cr_table(icr_Be10)) - &
            & gn * uu(:, flind%crn%beg - 1 + cr_table(icr_Be10)) / cr_tau(cr_table(icr_Be10))
      endif

      do i = lbound(icrH, 1), ubound(icrH, 1)
         associate( Hi => icrH(i) )
            if (eCRSP(Hi)) then
               do j = lbound(icrL, 1), ubound(icrL, 1)
               associate( &
                  Lj => icrL(j), &
                  idx => flind%crn%beg - 1 + cr_table(Hi) &
               )
                  if (eCRSP(Lj)) then
                     dcr = cr_sigma(cr_table(Hi), cr_table(Lj)) * dgas * uu(:, idx)
                     dcr = min(uu(:, idx)/rk_coeff, dcr)  ! Don't decay more elements than available
                     decrn(:, cr_table(Hi)) = decrn(:, cr_table(Hi)) - dcr
                     decrn(:, cr_table(Lj)) = decrn(:, cr_table(Lj)) + dcr
                  endif
               end associate
               enddo
         endif
         end associate
      enddo

   end subroutine src_crn

!>
!! \brief Execution of src_crn procedure designed for sources module
!<
   subroutine src_crn_exec(uu, n, usrc, rk_coeff)

      use fluidindex,     only: flind
      use initcosmicrays, only: iarr_crn

      implicit none

      integer(kind=4),               intent(in)  :: n
      real, dimension(n, flind%all), intent(in)  :: uu
      real,                          intent(in)  :: rk_coeff   !< coeffecient used in RK step, while computing source term
      real, dimension(n, flind%all), intent(out) :: usrc       !< u array update component for sources
!locals
      real, dimension(n, flind%crn%all)          :: srccrn

      usrc = 0.0
      call src_crn(uu, n, srccrn, rk_coeff) ! n safe
      usrc(:, iarr_crn) = srccrn(:,:)

   end subroutine src_crn_exec
#endif /* COSM_RAYS_SOURCES */
end module sourcecosmicrays
