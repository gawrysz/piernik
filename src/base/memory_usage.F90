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
!! \brief This module provides tools to check memory usage and abort the simulation before it runs into heavy swap usage
!!
!! In this module following namelists of parameters are specified:
!! \copydetails memory_usage::init_memory
!<
module memory_usage

   implicit none

   private
   public :: init_memory, system_mem_usage, check_mem_usage

   integer(kind=4) :: max_mem  !< Maximum allowed RSS memory per thread (in kiB)

   namelist /MEMORY/ max_mem

contains

!-----------------------------------------------------------------------------
!>
!! \brief Routine to set up allowed memory usage for the simulation
!!
!! \n \n
!! @b MEMORY
!! \n \n
!! <table border="+1">
!!   <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!!   <tr><td>max_mem          </td><td>huge(1)</td><td>integer value                        </td><td>\copydoc global::max_mem          </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_memory

      use dataio_pub, only: nh  ! QA_WARN required for diff_nml
      use mpisetup,   only: ibuff, master, slave, piernik_MPI_Bcast

      implicit none

      max_mem     = huge(1_4)

      if (master) then
         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=MEMORY)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=MEMORY, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "MEMORY")
         read(nh%cmdl_nml,nml=MEMORY, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "MEMORY", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=MEMORY)
         close(nh%lun)
         call nh%compare_namelist()

         ibuff(1) = max_mem

      endif

      call piernik_MPI_Bcast(ibuff)

      if (slave) then

         max_mem = ibuff(1)

      endif

   end subroutine init_memory

!>
!! \brief tell the current estimate of memory used by current thread.
!!
!! \todo Add optional intent(in) argument with other interesting labels such as VmPeak or VmSize
!>

   integer(kind=4) function system_mem_usage()

      use constants,  only: INVALID, fnamelen
      use dataio_pub, only: warn
#if defined(__INTEL_COMPILER)
      use ifport,     only: getpid
#endif /* __INTEL_COMPILER */

      implicit none

      integer, parameter :: pidlen = 8
      character(len=fnamelen) :: filename, line
      character(len=pidlen)   :: pid_char
      integer :: stat_lun, io
      logical :: io_exists

      system_mem_usage = INVALID

      write(pid_char, '(i8)') getpid()
      filename = '/proc/' // trim(adjustl(pid_char)) // '/status'

      inquire (file=filename, exist=io_exists)
      if (.not. io_exists) then
         call warn("[memory_usage:system_mem_usage] Cannot find '" // filename // "'")
         return
      endif

      open(newunit=stat_lun, file=filename, status='old')
      do
         read (stat_lun,'(a)', iostat=io) line
         if (io /= 0) exit
         if (line(1:6) == 'VmRSS:') then
            read (line(7:), *) system_mem_usage
            exit
         endif
      enddo

      close(stat_lun)

   end function system_mem_usage

!>
!! \brief Check if current memory usage doex not exceed the limit.
!!
!! \detailed This routine should be called after each potentially large allocaltion.
!! This may be used to crash Piernik before running into swapped memory.
!!
!! Cannot make a clean exit and call a restart dump from here because
!! it is not guaranteed that each process visits this routine.
!<

   subroutine check_mem_usage

      use dataio_pub, only: warn, die, msg
      use mpisetup,   only: proc

      implicit none

      integer :: rss
      real, parameter :: warnlevel = 0.95
      integer, save :: nextwarn = 1, warncnt = 0

      rss = system_mem_usage()

      if (rss > max_mem) then
         write(msg, '(a,2(i11,a),i5)')"[memory_usage:check_mem_usage] RSS exceeded max mem (", rss, " kiB > ", max_mem, " kiB) on process ", proc
         call die(msg)
      else if (rss > warnlevel * max_mem) then
         warncnt = warncnt + 1
         if (warncnt >= nextwarn) then
            write(msg, '(a,2(i11,a),i5)')"[memory_usage:check_mem_usage] RSS is approaching max mem (", rss, " kiB <= ", max_mem, " kiB) on process ", proc
            call warn(msg)
            nextwarn = 2 * nextwarn
         endif
      endif

   end subroutine check_mem_usage

end module memory_usage