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
!! Main program
!<
program piernik
! pulled by ANY

   use all_boundaries,    only: all_bnd
   use cg_leaves,         only: leaves
   use cg_list_global,    only: all_cg
   use constants,         only: PIERNIK_START, PIERNIK_INITIALIZED, PIERNIK_FINISHED, PIERNIK_CLEANUP, fplen, stdout, I_ONE, CHK, FINAL_DUMP, cbuff_len, PPP_IO, PPP_MPI
   use dataio,            only: write_data, user_msg_handler, check_log, check_tsl, dump, cleanup_dataio
   use dataio_pub,        only: nend, tend, msg, printinfo, warn, die, code_progress, nstep_start
   use div_B,             only: print_divB_norm
   use finalizepiernik,   only: cleanup_piernik
   use fluidindex,        only: flind
   use fluidupdate,       only: fluid_update
   use global,            only: t, nstep, dt, dtm, print_divB, tstep_attempt
   use initpiernik,       only: init_piernik
   use lb_helpers,        only: costs_maintenance
   use list_of_cg_lists,  only: all_lists
   use mpisetup,          only: master, piernik_MPI_Barrier, piernik_MPI_Bcast, cleanup_mpi
   use named_array_list,  only: qna, wna
   use ppp,               only: cleanup_profiling, update_profiling, ppp_main
   use refinement,        only: emergency_fix, updAMR_after
   use refinement_update, only: update_refinement
   use repeatstep,        only: repeat_step
   use timer,             only: walltime_end
   use timestep,          only: time_step, check_cfl_violation
   use user_hooks,        only: finalize_problem, problem_domain_update
#ifdef PERFMON
   use domain,            only: dom
   use timer,             only: timer_start, timer_stop
#endif /* PERFMON */
#if defined DEBUG && defined GRAV && defined NBODY
   use particle_diag,     only: print_all_particles
#endif /* DEBUG && GRAV && NBODY */
#ifdef MAGNETIC
   use all_boundaries,    only: all_mag_boundaries
#endif /* MAGNETIC */

   implicit none

   logical              :: end_sim             !< Used in main loop, to test whether to stop simulation or not
   logical, save        :: tleft = .true.      !< Used in main loop, to test whether to stop simulation or not
   character(len=fplen) :: nstr, tstr
   logical, save        :: first_step = .true., just_expanded
   real                 :: tlast
   logical              :: try_rebalance, rs

   ! ppp-related
   character(len=cbuff_len)    :: label, buf
   character(len=*), parameter :: f_label = "finalize"

   try_rebalance = .false.
   tlast = 0.0

   code_progress = PIERNIK_START

   call init_piernik

   call all_cg%check_na
   !call all_cg%check_for_dirt
#if defined DEBUG && defined GRAV && defined NBODY
   call print_all_particles
#endif /* DEBUG && GRAV && NBODY */

   call piernik_MPI_Barrier
!-------------------------------- MAIN LOOP ----------------------------------
#ifdef PERFMON
   call timer_start
#endif /* PERFMON */

   code_progress = PIERNIK_INITIALIZED

   end_sim = .false.

   if (master) then
      call printinfo("======================================================================================================", .false.)
      call printinfo("###############     Simulation     ###############", .false.)
      call printinfo("Named arrays present at start:", to_stdout=.false.)
      call qna%print_vars(to_stdout=.false.)
      call wna%print_vars(to_stdout=.false.)
      call printinfo("Grid lists present at start:", to_stdout=.false.)
   endif
   call all_lists%print(to_stdout = .false.)  ! needs all procs to participate
   if (master) then
      call printinfo("======================================================================================================", .false.)
   endif

   call print_progress(nstep)
   if (print_divB > 0) call print_divB_norm

   rs = repeat_step()  ! enforce function call
   do while (t < tend .and. nstep < nend .and. .not.(end_sim) .or. rs) ! main loop

      write(buf, '(i10)') nstep
      label = "step " // adjustl(trim(buf))
      if (tstep_attempt /= 0) then
         write(buf, '(i3)') tstep_attempt
         label = "repeated_" // trim(label) // "." // adjustl(trim(buf))
      endif
      call ppp_main%start(label)

      dump(:) = .false.
      if (associated(problem_domain_update)) then
         call problem_domain_update
         if (emergency_fix) try_rebalance = .true.
         just_expanded = emergency_fix
         call update_refinement(refinement_fixup_only=.true.)
#ifdef MAGNETIC
         ! Cheap and dirty fix: an extra update of exteral magnetic boundaries was sometimes needed after expanding magnetized domain.
         ! This is intended to cure insane div B appearing ar fine-coarse boundary touching the external boundary.
         ! The real cause is perhaps due to some data dependencies not fully met.
         if (just_expanded) call all_mag_boundaries
#endif /* MAGNETIC */
      endif

      call all_cg%check_na
      !call all_cg%check_for_dirt

      if (.not. repeat_step()) then
         call time_step(dt, flind, .true.)
         call grace_period

         dtm = dt

         call check_log
         call check_tsl

         tlast = t
      endif

      call fluid_update
      nstep = nstep + I_ONE
      call print_progress(nstep)
      call check_cfl_violation(flind)

      rs = repeat_step()  ! enforce function call
      if ((t - tlast < tiny(1.0)) .and. .not. first_step .and. .not. rs) call die("[piernik] timestep is too small: t == t + 2 * dt")

      call piernik_MPI_Barrier
      call costs_maintenance

      if (.not. repeat_step()) then
         call ppp_main%start('write_data', PPP_IO)
         call write_data(output=CHK)
         call ppp_main%stop('write_data', PPP_IO)

         call user_msg_handler(end_sim)
         call update_refinement
         ! A second call update_refinement here can be used to detect if there are refinement oscillations:
         ! * new "refine" or "correcting" events should not occur
         ! * some "derefine" events are allowed
         ! It can be used for diagnostic purposes. In production runs it may cost too much.
         if (try_rebalance) then
            !> \todo try to rewrite this ugly chain of flags passed through global variables into something more fool-proof
            call leaves%balance_and_update(" (re-balance) ")
            call all_bnd ! For some strange reasons this call prevents MPI-deadlock
            try_rebalance = .false.
         endif

         if (print_divB > 0) then
            if (mod(nstep, print_divB) == 0) call print_divB_norm
         endif
      endif

      if (master) tleft = walltime_end%time_left()
      call ppp_main%start('MPI_Bcast', PPP_MPI)
      call piernik_MPI_Bcast(tleft)
      call ppp_main%stop('MPI_Bcast', PPP_MPI)
      if (.not.tleft) end_sim = .true.

      first_step = .false.

      call ppp_main%stop(label)
      call update_profiling
      rs = repeat_step()  ! enforce function call

      emergency_fix = any(nstep_start + updAMR_after == nstep)
   enddo ! main loop

   if (print_divB > 0) then
      if (mod(nstep, print_divB) /= 0) call print_divB_norm ! print the norm at the end, if it wasn't printed inside the loop above
   endif

   code_progress = PIERNIK_FINISHED

   call ppp_main%start(f_label)
   if (master) then
      write(tstr, '(g14.6)') t
      write(nstr, '(i7)') nstep
      tstr = adjustl(tstr)
      nstr = adjustl(nstr)
      call printinfo("======================================================================================================", .false.)
      call printinfo("###############     Finishing     ###############", .false.)
      if (t >= tend) then
         write(msg, '(2a)') "Simulation has reached final time t = ",trim(tstr)
         call printinfo(msg)
      endif
      if (nstep >= nend) then
         write(msg, '(4a)') "Maximum step count exceeded (",trim(nstr),") at t = ",trim(tstr)
         call warn(msg)
      endif
      if (end_sim) then
         write(msg, '(4a)') "Enforced stop at step ",trim(nstr),", t = ", trim(tstr)
         call warn(msg)
      endif
      if (.not.tleft) then
         write(msg, '(4a)') "Wall time limit exceeded at step ",trim(nstr),", t = ", trim(tstr)
         call warn(msg)
      endif

      call printinfo("Named arrays present at finish:", to_stdout=.false.)
      call qna%print_vars(to_stdout=.false.)
      call wna%print_vars(to_stdout=.false.)
      call printinfo("Grid lists present at finish:", to_stdout=.false.)
   endif
   call all_lists%print(to_stdout = .false.)  ! needs all procs to participate
   if (master) then
      call printinfo("======================================================================================================", .false.)
   endif

   if (associated(finalize_problem)) call finalize_problem

#ifdef PERFMON
   call timer_stop(nstep, dom%total_ncells)
#endif /* PERFMON */
   call write_data(output=FINAL_DUMP)
!---------------------------- END OF MAIN LOOP ----------------------------------

   call piernik_MPI_Barrier

   code_progress = PIERNIK_CLEANUP

   if (master) write(stdout, '(a)', advance='no') "Finishing #"

   call cleanup_piernik

   call ppp_main%stop(f_label)
   call ppp_main%publish  ! we can use HDF5 here because we don't rely on anything that is affected by cleanup_hdf5
   call cleanup_profiling
   call cleanup_mpi
   call cleanup_dataio

   if (master) write(stdout,'(a)')"#"

contains

   subroutine print_progress(nstep)

      use constants,  only: tmr_fu
      use dataio_pub, only: printinfo, msg
      use global,     only: dt, t
      use timer,      only: set_timer, get_timestamp

      implicit none

      integer(kind=4), intent(in) :: nstep
      character(len=*), parameter :: fmt900 = "('   nstep = ',i7,'   dt = ',es23.16,'   t = ',es23.16,'   dWallClock = ',f10.2,' s ',a)"
      logical, save :: first_run = .true.


      if (master) then
         if (first_run) then
            first_run = (set_timer(tmr_fu) < 0.0)
         endif
         write(msg, fmt900) nstep, dt, t, set_timer(tmr_fu), get_timestamp()
         call printinfo(msg, .true.)
      endif

   end subroutine print_progress
!>
!! Meta subroutine responsible for setting proper pointers or doing other magic
!! after relaxation/grace period has passed
!<
   subroutine grace_period

      use all_boundaries, only: all_fluid_boundaries
      use dataio_pub,     only: printinfo
      use global,         only: grace_period_passed, relax_time
      use interactions,   only: interactions_grace_passed
      use mpisetup,       only: master
      use user_hooks,     only: problem_grace_passed

      implicit none

      logical, save     :: runned = .false.

      if (runned) return
      if (grace_period_passed()) then
         if (relax_time > 0.0) then
            ! write info message only if relax_time was set
            write(msg,'(A,ES10.3)') "[piernik:grace_period] grace period has passed.", t
            if (master) call printinfo(msg)
         endif
         call interactions_grace_passed
         if (associated(problem_grace_passed)) call problem_grace_passed
         call all_fluid_boundaries
         runned = .true.
      endif
   end subroutine grace_period

end program piernik
