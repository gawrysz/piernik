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

!> \brief This module provides a structure for accumulating computational costs of a cg

module cg_cost

   use constants, only: fnamelen

   implicit none

   private
   public :: cg_cost_data_t, cg_cost_t, cost_labels, I_MHD, I_MULTIGRID, I_MULTIPOLE, I_DIFFUSE, I_PARTICLE, I_REFINE, I_IC, I_OTHER

   enum, bind(C)
      enumerator :: I_MHD, I_MULTIGRID, I_MULTIPOLE, I_DIFFUSE, I_PARTICLE, I_REFINE, I_IC, I_OTHER
   end enum

   character(len=*), dimension(I_MHD:I_OTHER), parameter :: cost_labels = &
        [ "MHD       ", &  ! Riemann, RTVD, CT, divB cleaning
        & "multigrid ", &  ! self-gravity multigrid relaxation, residuals
        & "multipole ", &  ! multipole moments <=> potential conversions, not the costs of Q array manipulation
        & "diffusion ", &  ! explicit and multigrid diffussion costs
        & "particles ", &  ! particles in cg
        & "refines   ", &  ! prolongation, restriction, marking criteria
        & "init.cond.", &  ! for use only in the initproblems
        & "other     " ]   ! everything else that is related to cg but not tied to particular algorithm

   type :: cg_cost_data_t
      real, dimension(I_MHD:I_OTHER) :: wtime  ! walltime costs split into different categories
   end type cg_cost_data_t

   type, extends(cg_cost_data_t) :: cg_cost_t
      real, private :: wstart                    ! start value of the timer
      character(len=fnamelen), private :: label  ! used for debugging
   contains
      procedure :: reset  !< Set all counters to 0.
      procedure :: total  !< Return accumulated cost
      procedure :: start  !< Remember start time
      procedure :: stop   !< Stop measuring the time and add to specified timer
      procedure :: get    !< Read specified timer
   end type cg_cost_t

   real, parameter :: T_INVALID = -huge(1.)
   character(len=*), parameter :: invalid_label = "__INVALID__"

contains

!> \brief Set all counters to 0.

   subroutine reset(this)

      implicit none

      class(cg_cost_t), intent(out) :: this

      this%wtime = 0.
      this%wstart = T_INVALID
      this%label = invalid_label

   end subroutine reset

!> \brief Return accumulated cost

   real function total(this)

      implicit none

      class(cg_cost_t), intent(in) :: this

      total = sum(this%wtime(:))

   end function total

!>
!! \brief Remember start time
!!
!! Labels can be set up in the code with a shell script like:
!!     for i in `git grep -l costs%start ` ; do for j in `grep -n costs%start $i | sed 's/:.*//' ` ; do sed -i $j's/costs%start/costs%start("'`basename $i`":$j"'")/' $i ; done ; done
!!
!! I'm afraid that static inclusion of the labels into an often-called routine such as this one may impact the performance.
!! Since it is a debug feature, not required for normal operation, one may consider including it through a precompiler, where necessary.
!! Or create a branch taht will auto-generate the labels but won't nclude them into master.
!!
!<

   subroutine start(this, label)

      use dataio_pub, only: msg, warn
      use MPIF,       only: MPI_Wtime

      implicit none

      class(cg_cost_t),           intent(inout) :: this
      character(len=*), optional, intent(in)    :: label

      real :: t

      if (present(label)) this%label = trim(label)
      t = MPI_Wtime()
      if (this%wstart > T_INVALID) then
         write(msg, '(2(a,f18.6))')"[cg_cost:start] Some counting has already begin by '" // trim(this%label) // "' at ", &
              this%wstart, ". Resetting to ", t
         call warn(msg)
      endif
      this%wstart = t

   end subroutine start

!> \brief Stop measuring the time and add to specified timer

   subroutine stop(this, ind)

      use dataio_pub, only: msg, warn, die
      use MPIF,       only: MPI_Wtime

      implicit none

      class(cg_cost_t), intent(inout) :: this
      integer(kind=4),  intent(in)    :: ind

      real :: t

      t = MPI_Wtime()
      if (this%wstart <= T_INVALID) then
         write(msg, '(a,i3,a,f18.6)')"[cg_cost:start] Counting hasn't been started! index: ", ind, " time: ", t
         call warn(msg)
      endif

      if (ind >= lbound(this%wtime, 1) .and. ind <= ubound(this%wtime, 1)) then
         this%wtime(ind) = this%wtime(ind) + (t - this%wstart)
      else
         call die("[cg_cost:start] invalid index")
      endif

      this%wstart = T_INVALID
      this%label = invalid_label

   end subroutine stop

!> \brief Read specified timer

   real function get(this, ind)

      use dataio_pub, only: die

      implicit none

      class(cg_cost_t), intent(in) :: this
      integer(kind=4),  intent(in) :: ind

      if (ind >= lbound(this%wtime, 1) .and. ind <= ubound(this%wtime, 1)) then
         get = this%wtime(ind)
      else
         call die("[cg_cost:start] invalid index")
         get = 0.
      endif

   end function get

end module cg_cost
