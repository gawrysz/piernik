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
!! \brief Module containing a routine to compute upper limit of %timestep due to fluids %interactions.
!<
module timestepinteractions
! pulled by ANY
   implicit none
   private
   public :: timestep_interactions

contains
!>
!! \brief Routine that computes upper limit of %timestep due to fluids %interactions.
!! \warning works only with neutrals and dust case !!!!
!! \deprecated BEWARE: works only with neu+dust!!!!
!! \todo check if subtraction of momenta is really the case and rewrite for all fluids
!<
   real function timestep_interactions() result(dt_interact)

      use allreduce,    only: piernik_MPI_Allreduce
      use cg_leaves,    only: leaves
      use cg_list,      only: cg_list_element
      use constants,    only: pMIN, small
      use fluidindex,   only: flind
      use func,         only: L2norm
      use grid_cont,    only: grid_container
      use interactions, only: collfaq, cfl_interact, has_interactions

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      real                           :: val        !< variable used to store the maximum value of relative momentum

      dt_interact = huge(1.)
      if (.not.has_interactions) return
      !    dt_interact_proc = 1.0 / (maxval(collfaq)+small) / maxval(cg%u(iarr_all_dn,:,:,:))

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         val = maxval( L2norm(cg%u(flind%dst%imx,:,:,:),cg%u(flind%dst%imy,:,:,:),cg%u(flind%dst%imz,:,:,:), &
                            &  cg%u(flind%neu%imx,:,:,:),cg%u(flind%neu%imy,:,:,:),cg%u(flind%neu%imz,:,:,:) ) * cg%u(flind%dst%idn,:,:,:) )
         dt_interact = cfl_interact * flind%neu%cs / (maxval(collfaq) * val + small)


         cgl => cgl%nxt
      enddo

      call piernik_MPI_Allreduce(dt_interact, pMIN)

   end function timestep_interactions
!-------------------------------------------------------------------------------
end module timestepinteractions
