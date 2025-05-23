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
!! \brief Module of routines that correspond to resistivity
!!
!! In this module following namelist of parameters is specified:
!! \copydetails resistivity::init_resistivity
!<
module resistivity
! pulled by RESISTIVE
   use constants, only: dsetnamelen
   use types,     only: value

   implicit none

   private
   public  :: init_resistivity, timestep_resist, cleanup_resistivity, etamax, diffuseb, cu2max, deimin, eta1_active

   real                                  :: cfl_resist                     !< CFL factor for resistivity effect
   real                                  :: eta_0                          !< uniform resistivity
   real                                  :: eta_1                          !< anomalous resistivity
   real                                  :: j_crit                         !< critical value of current density
   real                                  :: jc2                            !< squared critical value of current density
   real                                  :: deint_max                      !< COMMENT ME
   integer(kind=4)                       :: eta_scale                      !< COMMENT ME
   real(kind=8)                          :: d_eta_factor
   type(value)                           :: etamax, cu2max, deimin
   logical, save                         :: eta1_active = .true.           !< resistivity off-switcher while eta_1 == 0.0
   character(len=dsetnamelen), parameter :: eta_n = "eta", wb_n = "wb", eh_n = "eh", dbx_n = "dbx", dby_n = "dby", dbz_n = "dbz"

contains

   subroutine cleanup_resistivity
      implicit none
   end subroutine cleanup_resistivity

!>
!! \brief Routine to set parameters values from namelist RESISTIVITY
!!
!! \n \n
!! @b RESISTIVITY
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>cfl_resist</td><td>0.4  </td><td>real value   </td><td>\copydoc resistivity::cfl_resist</td></tr>
!! <tr><td>eta_0     </td><td>0.0  </td><td>real value   </td><td>\copydoc resistivity::eta_0    </td></tr>
!! <tr><td>eta_1     </td><td>0.0  </td><td>real value   </td><td>\copydoc resistivity::eta_1    </td></tr>
!! <tr><td>eta_scale </td><td>4    </td><td>integer value</td><td>\copydoc resistivity::eta_scale</td></tr>
!! <tr><td>j_crit    </td><td>1.0e6</td><td>real value   </td><td>\copydoc resistivity::j_crit   </td></tr>
!! <tr><td>deint_max </td><td>0.01 </td><td>real value   </td><td>\copydoc resistivity::deint_max</td></tr>
!! </table>
!! The list is active while \b "RESISTIVE" is defined.
!! \n \n
!<
   subroutine init_resistivity

      use bcast,            only: piernik_MPI_Bcast
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use cg_list_global,   only: all_cg
      use constants,        only: PIERNIK_INIT_GRID, zdim, xdim, ydim, wcu_n
      use dataio_pub,       only: die, code_progress, nh
      use domain,           only: dom
      use func,             only: operator(.equals.)
      use mpisetup,         only: rbuff, ibuff, master, slave
      use named_array_list, only: qna
#ifdef ISO
      use constants,        only: zero
#endif /* ISO */

      implicit none

      real                           :: dims_twice
      type(cg_list_element), pointer :: cgl

      namelist /RESISTIVITY/ cfl_resist, eta_0, eta_1, eta_scale, j_crit, deint_max

      if (code_progress < PIERNIK_INIT_GRID) call die("[resistivity:init_resistivity] grid not initialized.")

      cfl_resist = 0.4
      eta_0      = 0.0
      eta_1      = 0.0
      eta_scale  = 4
      j_crit     = 1.0e6
      deint_max  = 0.01

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=RESISTIVITY)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=RESISTIVITY, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "RESISTIVITY")
         read(nh%cmdl_nml,nml=RESISTIVITY, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "RESISTIVITY", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=RESISTIVITY)
         close(nh%lun)
         call nh%compare_namelist()

         ibuff(1) = eta_scale

         rbuff(1) = cfl_resist
         rbuff(2) = eta_0
         rbuff(3) = eta_1
         rbuff(4) = j_crit
         rbuff(5) = deint_max

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         eta_scale  = ibuff(1)

         cfl_resist = rbuff(1)
         eta_0      = rbuff(2)
         eta_1      = rbuff(3)
         j_crit     = rbuff(4)
         deint_max  = rbuff(5)

      endif

      if (eta_scale < 0) call die("eta_scale must be greater or equal 0")

      call all_cg%reg_var(wcu_n)
      call all_cg%reg_var(eta_n)
      call all_cg%reg_var(wb_n)
      call all_cg%reg_var(eh_n)
      call all_cg%reg_var(dbx_n)
      call all_cg%reg_var(dby_n)
      call all_cg%reg_var(dbz_n)
#ifdef ISO
      if (eta_1 .equals. zero) then
         cgl => leaves%first
         do while (associated(cgl))
            cgl%cg%q(qna%ind(eta_n))%arr = eta_0
            cgl => cgl%nxt
         enddo
         etamax%val  = eta_0
         eta1_active = .false.
      endif
#endif /* ISO */

      if (eta1_active) then

         cgl => leaves%first
         do while (associated(cgl))
            if (.not. dom%has_dir(xdim)) cgl%cg%q(qna%ind(dbx_n))%arr = 0.0
            if (.not. dom%has_dir(ydim)) cgl%cg%q(qna%ind(dby_n))%arr = 0.0
            if (.not. dom%has_dir(zdim)) cgl%cg%q(qna%ind(dbz_n))%arr = 0.0
            cgl => cgl%nxt
         enddo

         jc2 = j_crit**2
         dims_twice = 2. * dom%eff_dim
         d_eta_factor = 1./(dims_twice+real(eta_scale, kind=8))
      endif

   end subroutine init_resistivity

   subroutine compute_resist

      use cg_cost_data,     only: I_MHD
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, zero, oneq, LO, HI, GEO_XYZ
      use dataio_pub,       only: die
      use domain,           only: dom
      use grid_cont,        only: grid_container
      use named_array_list, only: qna

      implicit none

      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      real, dimension(:,:,:), pointer :: eta, dbx, dby, dbz, wb, eh

      if (.not.eta1_active) return
!--- square current computing in cell corner step by step
      if (dom%geometry_type /= GEO_XYZ) call die("[resistivity:compute_resist] Unsupported geometry")

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         eta => cg%q(qna%ind(eta_n))%arr
         dbx => cg%q(qna%ind(dbx_n))%arr
         dby => cg%q(qna%ind(dby_n))%arr
         dbz => cg%q(qna%ind(dbz_n))%arr
         wb => cg%q(qna%ind(wb_n))%arr
         eh => cg%q(qna%ind(eh_n))%arr

         if (dom%has_dir(xdim)) then
            dbx(cg%lhn(xdim,LO)+1:cg%lhn(xdim,HI),:,:) = (cg%b(ydim,cg%lhn(xdim,LO)+1:cg%lhn(xdim,HI),:,:)-cg%b(ydim,cg%lhn(xdim,LO):cg%lhn(xdim,HI)-1,:,:))*cg%idl(xdim)
            dbx(cg%lhn(xdim,LO),:,:) = dbx(cg%lhn(xdim,LO)+1,:,:)
         endif
         if (dom%has_dir(ydim)) then
            dby(:,cg%lhn(ydim,LO)+1:cg%lhn(ydim,HI),:) = (cg%b(xdim,:,cg%lhn(ydim,LO)+1:cg%lhn(ydim,HI),:)-cg%b(xdim,:,cg%lhn(ydim,LO):cg%lhn(ydim,HI)-1,:))*cg%idl(ydim)
            dby(:,cg%lhn(ydim,LO),:) = dby(:,cg%lhn(ydim,LO)+1,:)
         endif
         if (dom%has_dir(zdim)) then
            dbz(:,:,cg%lhn(zdim,LO)+1:cg%lhn(zdim,HI)) = (cg%b(ydim,:,:,cg%lhn(zdim,LO)+1:cg%lhn(zdim,HI))-cg%b(ydim,:,:,cg%lhn(zdim,LO):cg%lhn(zdim,HI)-1))*cg%idl(zdim)
            dbz(:,:,cg%lhn(zdim,LO)) = dbz(:,:,cg%lhn(zdim,LO)+1)
         endif

!--- current_z **2
         eh = dbx - dby
         if (dom%has_dir(zdim)) then
            wb(:,:,cg%lhn(zdim,LO)+1:cg%lhn(zdim,HI)) =                                             oneq*(eh(:,:,cg%lhn(zdim,LO)+1:cg%lhn(zdim,HI)) + eh(:,:,cg%lhn(zdim,LO):cg%lhn(zdim,HI)-1))**2
            wb(:,:,cg%lhn(zdim,LO)) = wb(:,:,cg%lhn(zdim,LO)+1)
         else
            wb = eh**2
         endif
!--- current_x **2
         eh = dby - dbz
         if (dom%has_dir(xdim)) then
            wb(cg%lhn(xdim,LO)+1:cg%lhn(xdim,HI),:,:) = wb(cg%lhn(xdim,LO)+1:cg%lhn(xdim,HI),:,:) + oneq*(eh(cg%lhn(xdim,LO)+1:cg%lhn(xdim,HI),:,:) + eh(cg%lhn(xdim,LO):cg%lhn(xdim,HI)-1,:,:))**2
            wb(cg%lhn(xdim,LO),:,:) = wb(cg%lhn(xdim,LO)+1,:,:)
         else
            wb = wb + eh**2
         endif
!--- current_y **2
         eh = dbz - dbx
         if (dom%has_dir(ydim)) then
            wb(:,cg%lhn(ydim,LO)+1:cg%lhn(ydim,HI),:) = wb(:,cg%lhn(ydim,LO)+1:cg%lhn(ydim,HI),:) + oneq*(eh(:,cg%lhn(ydim,LO)+1:cg%lhn(ydim,HI),:) + eh(:,cg%lhn(ydim,LO):cg%lhn(ydim,HI)-1,:))**2
            wb(:,cg%lhn(ydim,LO),:) = wb(:,cg%lhn(ydim,LO)+1,:)
         else
            wb = wb + eh**2
         endif

!        eta(:,:,:) = eta_0 + eta_1 * sqrt( max(0.0,wb(:,:,:)- jc2 ))
!        the above may cause FPE because compiler may transform it to max(0.0, sqrt(wb(:,:,:)- jc2 ))
         where (wb(:,:,:) - jc2 > zero)
            eta(:,:,:) = eta_0 + eta_1 * sqrt(wb(:,:,:)- jc2)
         elsewhere
            eta(:,:,:) = eta_0
         endwhere

         eh = zero
         if (dom%has_dir(xdim)) then
            eh(cg%lhn(xdim,LO)+1:cg%lhn(xdim,HI)-1,:,:) = eh(cg%lhn(xdim,LO)+1:cg%lhn(xdim,HI)-1,:,:) + eta(cg%lhn(xdim,LO):cg%lhn(xdim,HI)-2,:,:) + eta(cg%lhn(xdim,LO)+2:cg%lhn(xdim,HI),:,:)
            eh(cg%lhn(xdim,LO),:,:) = eh(cg%lhn(xdim,LO)+1,:,:) ; eh(cg%lhn(xdim,HI),:,:) = eh(cg%lhn(xdim,HI)-1,:,:)
         endif
         if (dom%has_dir(ydim)) then
            eh(:,cg%lhn(ydim,LO)+1:cg%lhn(ydim,HI)-1,:) = eh(:,cg%lhn(ydim,LO)+1:cg%lhn(ydim,HI)-1,:) + eta(:,cg%lhn(ydim,LO):cg%lhn(ydim,HI)-2,:) + eta(:,cg%lhn(ydim,LO)+2:cg%lhn(ydim,HI),:)
            eh(:,cg%lhn(ydim,LO),:) = eh(:,cg%lhn(ydim,LO)+1,:) ; eh(:,cg%lhn(ydim,HI),:) = eh(:,cg%lhn(ydim,HI)-1,:)
         endif
         if (dom%has_dir(zdim)) then
            eh(:,:,cg%lhn(zdim,LO)+1:cg%lhn(zdim,HI)-1) = eh(:,:,cg%lhn(zdim,LO)+1:cg%lhn(zdim,HI)-1) + eta(:,:,cg%lhn(zdim,LO):cg%lhn(zdim,HI)-2) + eta(:,:,cg%lhn(zdim,LO)+2:cg%lhn(zdim,HI))
            eh(:,:,cg%lhn(zdim,LO)) = eh(:,:,cg%lhn(zdim,LO)+1) ; eh(:,:,cg%lhn(zdim,HI)) = eh(:,:,cg%lhn(zdim,HI)-1)
         endif
         eh = real((eh + eta_scale*eta)*d_eta_factor)

         where (eta > eta_0) eta = eh

         call cg%costs%stop(I_MHD)
         cgl => cgl%nxt
      enddo

   end subroutine compute_resist

!-----------------------------------------------------------------------

   subroutine timestep_resist(dt)

      use allreduce,        only: piernik_MPI_Allreduce
      use bcast,            only: piernik_MPI_Bcast
      use cg_cost_data,     only: I_OTHER
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: big, zero, pMIN, MAXL
      use grid_cont,        only: grid_container
      use func,             only: operator(.notequals.)
      use named_array_list, only: qna
      use types,            only: value
#ifndef ISO
      use constants,        only: MINL
#ifdef IONIZED
      use constants,        only: small, xdim, ydim, zdim
      use fluidindex,       only: flind
      use func,             only: ekin, emag
      use named_array_list, only: wna
#endif /* IONIZED */
#endif /* !ISO */

      implicit none

      real, intent(inout)               :: dt
      type(cg_list_element),  pointer   :: cgl
      type(grid_container),   pointer   :: cg
      real                              :: dt_eta, dt_eint
#if !defined(ISO) && defined(IONIZED)
      real, dimension(:,:,:),   pointer :: eta, wb, eh
      real, dimension(:,:,:,:), pointer :: uu, bb
#endif /* !ISO && IONIZED */

      dt_eta = big ; dt_eint = big
      call compute_resist
      call leaves%get_extremum(qna%ind(eta_n), MAXL, etamax)
      call piernik_MPI_Bcast(etamax%val)
      if (eta1_active) then
         call leaves%get_extremum(qna%ind(wb_n), MAXL, cu2max)
      else
         cu2max = value(0., 0., [0., 0., 0.], [0, 0, 0], 0_4)
      endif
      call piernik_MPI_Bcast(cu2max%val)

      if (etamax%val .notequals. zero) then
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            call cg%costs%start

            dt_eta = min(dt_eta, cfl_resist * cg%dxmn2 / (2. * etamax%val))
#ifndef ISO
#ifdef IONIZED
            eta => cg%q(qna%ind(eta_n))%span(cg%ijkse)
            wb => cg%q(qna%ind(wb_n))%span(cg%ijkse)
            eh => cg%q(qna%ind(eh_n))%span(cg%ijkse)
            uu => cg%w(wna%fi)%span(cg%ijkse)
            bb => cg%W(wna%bi)%span(cg%ijkse)
            eh = (uu(flind%ion%ien,:,:,:) - ekin(uu(flind%ion%imx,:,:,:), uu(flind%ion%imy,:,:,:), uu(flind%ion%imz,:,:,:), uu(flind%ion%idn,:,:,:)) - &
                  emag(bb(xdim,:,:,:), bb(ydim,:,:,:), bb(zdim,:,:,:)))/ (eta(:,:,:) * wb + small)
            dt_eint = min(dt_eint, deint_max * abs(minval(eh)))
#endif /* IONIZED */
#endif /* !ISO */

            call cg%costs%stop(I_OTHER)
            cgl => cgl%nxt
         enddo
      endif

      call piernik_MPI_Allreduce(dt_eta, pMIN)
#ifndef ISO
#ifdef IONIZED
      call piernik_MPI_Allreduce(dt_eint, pMIN)
#endif /* IONIZED */
      call leaves%get_extremum(qna%ind(eh_n), MINL, deimin)
      deimin%assoc = dt_eint
#endif /* !ISO */
      etamax%assoc = dt_eta ; cu2max%assoc = min(dt_eta, dt_eint)

      dt = min(dt, dt_eta, dt_eint)

   end subroutine timestep_resist

!-----------------------------------------------------------------------------
!>
!! \brief
!! \todo overload me or use class(*) if you dare
!<
   subroutine vanleer_limiter(f,a,b)

      implicit none

      real, dimension(:), intent(in)    :: a !< second order correction of left- or right- moving waves flux on the left cell boundary
      real, dimension(:), intent(in)    :: b !< second order correction of left- or right- moving waves flux on the right cell boundary
      real, dimension(:), intent(inout) :: f !< second order flux correction for left- or right- moving waves
      ! locals
      real, dimension(size(a,1))        :: c !< a*b

      c = a*b                                                                    !> \todo OPTIMIZE ME
      where (c > 0.0)
         f = f+2.0*c/(a+b)
      endwhere

   end subroutine vanleer_limiter

   subroutine tvdd_1d(b1d,eta1d,idi,dt,wcu1d)

      use constants,     only: half
      implicit none

      real, dimension(:), pointer, intent(in)    :: eta1d, b1d
      real, dimension(:), pointer, intent(inout) :: wcu1d
      real,                        intent(in)    :: idi,dt

      real, dimension(size(b1d))               :: w, wp, wm, b1
      integer                                  :: n

      n = size(b1d)
      w(2:n)    = eta1d(2:n) * ( b1d(2:n) - b1d(1:n-1) )*idi ;  w(1)  = w(2)
      b1(1:n-1) = b1d(1:n-1) + half*(w(2:n) - w(1:n-1))*dt*idi; b1(n) = b1(n-1)

      w(2:n)    = eta1d(2:n) * ( b1(2:n) - b1(1:n-1) )*idi   ; w(1)  = w(2)
      wp(1:n-1) = half*(w(2:n) - w(1:n-1))                   ; wp(n) = wp(n-1)
      wm(2:n)   = wp(1:n-1)                                  ; wm(1) = wm(2)

      call vanleer_limiter(w,wm,wp)
      wcu1d     = w*dt

   end subroutine tvdd_1d

!-------------------------------------------------------------------------------
!
! 6 routines have been substituted by one with parameters:
!   diffuseby_x  --> diffuseb(ibdir = ydim, sdir = xdim, etadir = zdim, emf = 'emfz', n1 = ydim, n2 = zdim)
!   diffusebz_x  --> diffuseb(ibdir = zdim, sdir = xdim, etadir = ydim, emf = 'emfy', n1 = ydim, n2 = zdim)
!   diffusebz_y  --> diffuseb(ibdir = zdim, sdir = ydim, etadir = xdim, emf = 'emfx', n1 = zdim, n2 = xdim)
!   diffusebx_y  --> diffuseb(ibdir = xdim, sdir = ydim, etadir = zdim, emf = 'emfz', n1 = zdim, n2 = xdim)
!   diffusebx_z  --> diffuseb(ibdir = xdim, sdir = zdim, etadir = ydim, emf = 'emfy', n1 = xdim, n2 = ydim)
!   diffuseby_z  --> diffuseb(ibdir = ydim, sdir = zdim, etadir = xdim, emf = 'emfx', n1 = xdim, n2 = ydim)

   subroutine diffuseb(ibdir, sdir)

      use cg_cost_data,     only: I_MHD
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, ndims, half, I_ONE, wcu_n, idm, INT4, LO, HI
      use domain,           only: dom
      use global,           only: dt
      use grid_cont,        only: grid_container
      use magboundaries,    only: bnd_emf
      use named_array_list, only: qna, wna

      implicit none

      integer(kind=4),  intent(in)      :: ibdir, sdir
      integer                           :: i1, i2
      integer(kind=4)                   :: n1, n2, etadir, dir, emf, wcu_i, eta_i
      integer(kind=4), dimension(ndims) :: idml, idmh
      real, dimension(:),    pointer    :: b1d, eta1d, wcu1d
      type(cg_list_element), pointer    :: cgl
      type(grid_container),  pointer    :: cg

      n1 = I_ONE + mod(sdir    ,   ndims)
      n2 = I_ONE + mod(sdir+I_ONE, ndims)
      etadir = sum([xdim,ydim,zdim]) - ibdir - sdir

      call compute_resist

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         wcu_i = qna%ind(wcu_n)
         eta_i = qna%ind(eta_n)

         idmh(:) = cg%lhn(:,HI) - idm(:,etadir)
         idml(:) = cg%lhn(:,LO) + idm(:,etadir)
         cg%q(eta_i)%arr(cg%lhn(xdim,LO):idmh(xdim),cg%lhn(ydim,LO):idmh(ydim),cg%lhn(zdim,LO):idmh(zdim)) = half*(cg%q(eta_i)%span(cg%lhn(:,LO),idmh) + cg%q(eta_i)%span(idml,cg%lhn(:,HI)))

         do i1 = cg%lhn(n1,LO), cg%lhn(n1,HI)
            do i2 = cg%lhn(n2,LO), cg%lhn(n2,HI)
               b1d   => cg%w(wna%bi)%get_sweep(sdir,ibdir,i1,i2)
               eta1d => cg%q(eta_i )%get_sweep(sdir,      i1,i2)
               wcu1d => cg%q(wcu_i )%get_sweep(sdir,      i1,i2)
               call tvdd_1d(b1d, eta1d, cg%idl(sdir), dt, wcu1d)
            enddo
         enddo

         call cg%costs%stop(I_MHD)
         cgl => cgl%nxt
      enddo

      cgl => leaves%first
      do while (associated(cgl))
         call cg%costs%start

         do dir = xdim, zdim
            emf = idm(etadir,dir) + 2_INT4
            if (dom%has_dir(dir)) call bnd_emf(wcu_i, emf, dir, cgl%cg)
         enddo

         call cg%costs%stop(I_MHD)
         cgl => cgl%nxt
      enddo

   end subroutine diffuseb

end module resistivity
