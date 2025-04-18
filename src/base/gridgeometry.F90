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
!! \brief Module that sets coordinate system
!! \date December 2010
!!
!! This module contains routines that calculates grid coefficients and source terms used during flux update in rtvd::relaxing_tvd
!!
!<
module gridgeometry

   implicit none

   private
   public   :: gc, GC1, GC2, GC3, init_geometry, set_geo_coeffs, geometry_source_terms_exec

   real, dimension(:,:,:), pointer :: gc            !< array of geometrical coefficients, \todo move it to the grid container
   enum, bind(C)
      enumerator :: GC1, GC2, GC3                   !< \todo change to a some meaningful names
   end enum

   interface
      !>
      !! \brief interface for routine setting grid coefficients
      !<
      subroutine set_gc(sweep, flind, i1, i2, cg)

         use fluidtypes, only: var_numbers
         use grid_cont,  only: grid_container

         implicit none

         integer(kind=4),  intent(in)  :: sweep         !< direction (xdim, ydim or zdim) we are doing calculations for
         type(var_numbers), intent(in) :: flind         !< \copydoc fluidindex::flind
         integer, intent(in)           :: i1            !< coordinate of sweep in the 1st remaining direction
         integer, intent(in)           :: i2            !< coordinate of sweep in the 2st remaining direction
         type(grid_container), pointer :: cg            !< current grid container

      end subroutine set_gc

      !>
      !! \brief interface for routine returning grid dependent source terms
      !!
      !! Currently, gsrc function returns accelerations
      !<
      function gsrc(u, b, sweep, i1, i2, cg) result(res)

         use fluidindex, only: flind
         use grid_cont, only: grid_container

         implicit none

         integer(kind=4),               intent(in) :: sweep !< direction (x, y or z) we are doing calculations for
         integer,                       intent(in) :: i1    !< coordinate of sweep in the 1st remaining direction
         integer,                       intent(in) :: i2    !< coordinate of sweep in the 2nd remaining direction
         real, dimension(:,:), target,  intent(in) :: u     !< sweep of fluid conservative variables
         real, dimension(:,:), target,  intent(in) :: b     !< local copy of magnetic field
         type(grid_container), pointer, intent(in) :: cg    !< current grid container

         real, dimension(size(u,1),flind%fluids)   :: res   !< output sweep of accelerations

      end function gsrc

   end interface

   procedure(set_gc), pointer :: set_geo_coeffs          !< generic pointer for routine setting geometrical coefficients
   procedure(gsrc),   pointer :: geometry_source_terms   !< generic pointer for routine calculating source terms

contains
!>
!! \brief Routine for module initialization
!!
!! \details Routine associates gridgeometry::set_geo_coeffs() and gridgeometry::geometry_source_terms() pointers.
!<
   subroutine init_geometry

      use dataio_pub, only: msg, die, code_progress
      use constants,  only: PIERNIK_INIT_DOMAIN, GEO_XYZ, GEO_RPZ
      use domain,     only: dom

      implicit none

      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[gridgeometry:init_geometry] domain not initialized")

      select case (dom%geometry_type)
         case (GEO_XYZ)
            set_geo_coeffs          => set_cart_coeffs
            geometry_source_terms   => cart_geometry_source_terms
         case (GEO_RPZ)
            set_geo_coeffs          => set_cyl_coeffs
            geometry_source_terms   => cyl_geometry_source_terms
         case default
            write(msg,'(a,i3)') "[gridgeometry:init_geometry] Unknown system of coordinates ", dom%geometry_type
            call die(msg)
      end select

   end subroutine init_geometry

!>
!! \brief Routine allocating auxiliary arrays
!!
!! We need to allocate those arrays later to have flind%all
!<
   subroutine geo_coeffs_arrays(flind, cg)

      use constants,  only: xdim, ydim, zdim
      use dataio_pub, only: die
      use fluidtypes, only: var_numbers
      use grid_cont,  only: grid_container

      implicit none

      type(var_numbers), intent(in) :: flind
      type(grid_container), pointer :: cg

      if ( any( [allocated(cg%gc_xdim), allocated(cg%gc_ydim), allocated(cg%gc_zdim)] ) ) then
         call die("[gridgeometry:geo_coeffs_arrays] double allocation")
      else
         allocate(cg%gc_xdim(GC1:GC3, cg%n_(xdim), flind%all))
         allocate(cg%gc_ydim(GC1:GC3, cg%n_(ydim), flind%all))
         allocate(cg%gc_zdim(GC1:GC3, cg%n_(zdim), flind%all))
      endif

   end subroutine geo_coeffs_arrays
!>
!! \brief A dummy routine. There is no need to set any cartesian coefficients, because all of them are equal to 1. so we don't use them in rtvd
!<
   subroutine set_cart_coeffs(sweep, flind, i1, i2, cg)

      use fluidtypes, only: var_numbers
      use grid_cont,  only: grid_container

      implicit none

      integer(kind=4),   intent(in) :: sweep
      type(var_numbers), intent(in) :: flind
      integer, intent(in)           :: i1, i2
      type(grid_container), pointer :: cg

      return
      if (.false.) write(0,*) cg%u(flind%all, sweep, i1, i2) ! suppress compiler warnings

   end subroutine set_cart_coeffs
!>
!! \brief routine setting geometrical coefficients for cylindrical grid
!<
   subroutine set_cyl_coeffs(sweep, flind, i1, i2, cg)

      use constants,  only: xdim, ydim, zdim, INV_CENTER, LEFT, RIGHT
      use dataio_pub, only: die, msg
      use fluidtypes, only: var_numbers
      use grid_cont,  only: grid_container

      implicit none

      integer(kind=4), intent(in)   :: sweep
      type(var_numbers), intent(in) :: flind
      integer, intent(in)           :: i1, i2
      type(grid_container), pointer :: cg

      integer                        :: i

      !> \todo This should be probably called from cg%init (beware of cyclic dependencies) or init_grid
      if (.not. allocated(cg%gc_xdim)) then ! assume allocated(cg%gc_xdim) .eqv. allocated(cg%gc_ydim) .eqv. allocated(cg%gc_zdim)

         call geo_coeffs_arrays(flind, cg)

         do i = 1, flind%all
            cg%gc_xdim(GC1,:,i) = cg%coord(INV_CENTER, xdim)%r(:)
            cg%gc_xdim(GC2,:,i) = cg%coord(RIGHT, xdim)%r(:)
            cg%gc_xdim(GC3,:,i) = cg%coord(LEFT, xdim)%r(:)
         enddo

         do i = lbound(flind%all_fluids,1), ubound(flind%all_fluids,1)
            cg%gc_xdim(GC1, :, flind%all_fluids(i)%fl%imy) = cg%gc_xdim(GC1, :, flind%all_fluids(i)%fl%imy) * cg%coord(INV_CENTER, xdim)%r(:)
            cg%gc_xdim(GC2, :, flind%all_fluids(i)%fl%imy) = cg%gc_xdim(GC2, :, flind%all_fluids(i)%fl%imy) * cg%coord(RIGHT, xdim)%r(:)
            cg%gc_xdim(GC3, :, flind%all_fluids(i)%fl%imy) = cg%gc_xdim(GC3, :, flind%all_fluids(i)%fl%imy) * cg%coord(LEFT, xdim)%r(:)
         enddo
         cg%gc_ydim(GC2:GC3,:,:) = 1.0     ! [ 1/r , 1 , 1]

         cg%gc_zdim(:,:,:) = 1.0           ! [ 1, 1, 1]

      endif

      select case (sweep)
         case (xdim)
            gc => cg%gc_xdim
         case (ydim)
            cg%gc_ydim(GC1,:,:)   = cg%inv_x(i2)
            gc => cg%gc_ydim
         case (zdim)
            gc => cg%gc_zdim
         case default
            write(msg,'(a,i2)') "[gridgeometry:set_cyl_coeffs] Unknown sweep : ",sweep
            call die(msg)
            if (.false.)  write(0,*) i1,i2 ! suppress compiler warnings
      end select

   end subroutine set_cyl_coeffs
!>
!! \brief routine calculating geometrical source term for cartesian grid
!<
   function cart_geometry_source_terms(u, bb, sweep, i1, i2, cg) result(res)

      use fluidindex, only: flind
      use grid_cont,  only: grid_container

      implicit none

      integer(kind=4),               intent(in) :: sweep
      integer,                       intent(in) :: i1, i2
      real, dimension(:,:), target,  intent(in) :: u, bb
      type(grid_container), pointer, intent(in) :: cg

      real, dimension(size(u,1),flind%fluids)   :: res

      res = 0.0
      return
      if (.false.) write(0,*) sweep, i1, i2, u, bb, cg%inv_x(1)

   end function cart_geometry_source_terms
!>
!! \brief routine calculating geometrical source term for cylindrical grid
!<
   function cyl_geometry_source_terms(u, bb, sweep, i1, i2, cg) result(res)

      use constants,        only: cs_i2_n, xdim
      use fluidindex,       only: flind, iarr_all_dn, iarr_all_my
      use fluidtypes,       only: component_fluid
      use grid_cont,        only: grid_container
      use named_array_list, only: qna

      implicit none

      integer(kind=4),               intent(in) :: sweep
      integer,                       intent(in) :: i1, i2
      real, dimension(:,:), target,  intent(in) :: u, bb
      type(grid_container), pointer, intent(in) :: cg

      real, dimension(size(u,1),flind%fluids)         :: res
      real, dimension(size(u,1),flind%fluids), target :: pp       !< array storing pressure in current sweep (reused later)
      class(component_fluid), pointer                 :: pfl
      real, dimension(:,:),   pointer                 :: puu, pbb
      real, dimension(:),     pointer                 :: ppp, cs2
      integer :: ifl

      if (sweep == xdim) then
         pbb   =>   bb(:,:)
         if (qna%exists(cs_i2_n)) then
            cs2 => cg%q(qna%ind(cs_i2_n))%get_sweep(sweep,i1,i2)
         else
            cs2 => null()
         endif
         do ifl = 1, flind%fluids
            pfl   => flind%all_fluids(ifl)%fl
            puu   =>    u(:, pfl%beg:pfl%end)
            ppp   =>   pp(:, pfl%pos)

            call pfl%compute_pres(cg%n_(sweep), puu, pbb, cs2, ppp)
            res(:, ifl) = cg%inv_x(:) * (u(:, iarr_all_my(ifl))**2 / u(:, iarr_all_dn(ifl)) + ppp(:))
         enddo
      else
         ! Note that there's no additional source term for angular momentum since we're using
         ! modified equation of motion following  Mignone et al. (2007), ApJS 170:228- and
         ! Skinner & Ostriker (2010), ApJSS 188:290-311
         ! The conservative implementation uses the gc(:,:,:) array to modify the azimuthal momentum flux to mimic angular momentum flux.
         res = 0.0
      endif

   end function cyl_geometry_source_terms

!>
!! \brief Detailed execution of geometry source terms designed to use in sources module
!<
   subroutine geometry_source_terms_exec(u, bb, sweep, i1, i2, cg, usrc)

      use fluidindex, only: iarr_all_mx
      use grid_cont,  only: grid_container

      implicit none

      integer(kind=4),               intent(in)    :: sweep   !< direction (x, y or z) we are doing calculations for
      integer,                       intent(in)    :: i1, i2
      real, dimension(:, :),         intent(in)    :: u       !< vector of conservative variables
      real, dimension(:, :),         intent(in)    :: bb      !< local copy of magnetic field
      type(grid_container), pointer, intent(in)    :: cg      !< current grid piece
      real, dimension(:, :),         intent(inout) :: usrc    !< u array update component for sources

      usrc = 0.0
      usrc(:, iarr_all_mx) = geometry_source_terms(u, bb, sweep, i1, i2, cg)

   end subroutine geometry_source_terms_exec

end module gridgeometry
