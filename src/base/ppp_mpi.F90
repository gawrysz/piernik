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

!> \brief Module for PPP-guarded MPI calls

module ppp_mpi

   implicit none

   private
   public :: piernik_Waitall

contains

!>
!! \brief a PPP wrapper for MPI_Waitall
!!
!! This routine may sometimes be called only by a subset of MPI ranks.
!! Do not insert MPI_Barrier and do not assume equal call count over all ranks.
!! This may result in differing count of associated PPP timer on different processes.
!<

#define DEBUG_MPI

   subroutine piernik_Waitall(nr, ppp_label, x_mask, use_req2)

      use constants, only: PPP_MPI
      use mpisetup,  only: err_mpi, req, req2, piernik_MPI_Barrier, extra_barriers
      use MPIF,      only: MPI_STATUSES_IGNORE
      use MPIFUN,    only: MPI_Waitall
      use ppp,       only: ppp_main
#ifdef DEBUG_MPI
      use constants, only: INVALID
      use MPIF,      only: MPI_Wtime, MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE, MPI_COMM_WORLD, MPI_Abort, MPI_Request_get_status
      use mpisetup,  only: proc
#endif /* DEBUG_MPI */

      implicit none

      integer(kind=4),           intent(in) :: nr         !< number of requests in req(:) or req2(:)
      character(len=*),          intent(in) :: ppp_label  !< identifier for PPP entry
      integer(kind=4), optional, intent(in) :: x_mask     !< extra mask, if necessary
      logical,         optional, intent(in) :: use_req2   !< use req2 when .true.

      character(len=*), parameter :: mpiw = "MPI_Waitall:"
      integer(kind=4) :: mask
      logical :: r2
#ifdef DEBUG_MPI
      real :: wt0
      logical :: flag
      integer :: i
      integer(kind=4) :: tcnt, cnt_prev
      integer(kind=4), allocatable, dimension(:) :: flags
      real, parameter :: timeout = 10., indecent_time = 0.1 * timeout
      integer(kind=4), parameter :: timeout_code = 17
      logical, parameter :: crash_on_timeout = .false., use_request_get_status = .false.
#endif /* DEBUG_MPI */
      if (nr > 0) then

         mask = PPP_MPI
         if (present(x_mask)) mask = mask + x_mask

#ifdef DEBUG_MPI
         wt0 = MPI_Wtime()
         allocate(flags(nr))
         flags(:) = INVALID
         tcnt = 0
         do while (tcnt < nr)
            cnt_prev = tcnt

            tcnt = 0
            if (use_request_get_status) then
               do i = 1, nr
                  call MPI_Request_get_status(req(i), flag, MPI_STATUS_IGNORE, err_mpi)
                  ! For unknown reasons for OpenMPI 2.1.1-8 (Ubuntu 18.04) the flag is always .false.
                  ! In Ubuntu 20.04 (OpenMPI 4.0.3) it doesn't work too.
                  if (flag) tcnt = tcnt + 1
               enddo
            else
               ! This is less sterile than use of MPI_Request_get_status but may release some internal buffers while waiting for late requests to complete
               call MPI_Testsome(nr, req(:nr), tcnt, flags(cnt_prev+1:), MPI_STATUSES_IGNORE, err_mpi)
               tcnt = cnt_prev + tcnt
            endif
            if (tcnt /= cnt_prev .and. MPI_Wtime() - wt0 > indecent_time) &
                 write(*,*)".@", proc, ppp_label, " : ", tcnt, " out of ", nr, " requests completed in ", MPI_Wtime() - wt0, "s (and counting)"  ! QA_WARN debug

            if (MPI_Wtime() - wt0 > timeout) then
               write(*,*)"-@", proc, ppp_label, " : only ", tcnt, " out of ", nr, " requests completed in ", MPI_Wtime() - wt0, "s (timeout)"  ! QA_WARN debug
               if (crash_on_timeout) then
                  call MPI_Abort(MPI_COMM_WORLD, timeout_code, err_mpi)
               else
                  exit
               endif
            endif

         enddo

         ! write(*,*)"+@", proc, ppp_label, " : ", nr, " requests completed in ", MPI_Wtime() - wt0, "s"  ! QA_WARN debug
         deallocate(flags)
#endif /* DEBUG_MPI */

         call ppp_main%start(mpiw // ppp_label, mask)

         r2 = .false.
         if (present(use_req2)) r2 = use_req2
         if (r2) then
            call MPI_Waitall(nr, req2(:nr), MPI_STATUSES_IGNORE, err_mpi)
         else
            call MPI_Waitall(nr, req(:nr), MPI_STATUSES_IGNORE, err_mpi)
         endif

         call ppp_main%stop(mpiw // ppp_label, mask)
      endif

      if (extra_barriers) call piernik_MPI_Barrier

   end subroutine piernik_Waitall

end module ppp_mpi
