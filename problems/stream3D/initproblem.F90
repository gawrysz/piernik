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

! Initial condition for Keplerian disk
! Written by: M. Hanasz, March 2006

   use constants,    only: cbuff_len

   implicit none

   private
   public  :: read_problem_par, problem_initial_conditions, problem_pointers

   real    :: sigma0, Rin, R0, HtoR, eps, amp
   character(len=cbuff_len) :: sigma_model

   namelist /PROBLEM_CONTROL/  sigma0, amp, Rin, R0, HtoR, sigma_model, eps

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use constants,  only: cbuff_len
      use dataio_pub, only: nh
      use mpisetup,   only: cbuff, rbuff, master, slave

      implicit none

      Sigma0  = 1.0
      Rin     = 1.0e-4
      R0      = 1.0
      HtoR    = 1.0
      sigma_model = 'hayashi'
      eps     = 1.0
      amp     = 1.e-5

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

         cbuff(1) =  sigma_model

         rbuff(1) = sigma0
         rbuff(2) = Rin
         rbuff(3) = R0
         rbuff(4) = HtoR
         rbuff(5) = eps
         rbuff(6) = amp

      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         sigma_model  = cbuff(1)

         sigma0       = rbuff(1)
         Rin          = rbuff(2)
         R0           = rbuff(3)
         HtoR         = rbuff(4)
         eps          = rbuff(5)
         amp          = rbuff(6)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   real function dens_Rdistr(R,Rin,n)

      implicit none

      real, intent(in) :: R,Rin,n
      real :: ninv

      ninv = 1./n
      dens_Rdistr = max((R - Rin),Rin)**ninv / R**(2.+ninv)

   end function dens_Rdistr

   subroutine problem_initial_conditions

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: pi, dpi, xdim, ydim, zdim, LO, HI
      use domain,     only: dom
      use fluidindex, only: flind
      use global,     only: smalld
      use gravity,    only: ptmass
      use grid_cont,  only: grid_container
      use units,      only: newtong
#ifndef ISO
      use global,     only: smallei
#endif /* !ISO */

      implicit none

      integer                           :: i, j, k, mnz
      real                              :: xi, yj, zk, rc, H0, sqr_gm, rho0
      real                              :: n,norm, H, ninv
      real                              :: gradP, iOmega, ilook, gradgp
      real, dimension(:,:), allocatable :: noise
      real, dimension(:),   allocatable :: omega,omegad
      type(cg_list_element), pointer    :: cgl
      type(grid_container),  pointer    :: cg
#ifndef ISO
      real                              :: vx, vy, vz
#endif /* !ISO */

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         mnz = max(cg%lhn(zdim,LO)+cg%n_(zdim)/2-1,1)

         allocate(omega(cg%lhn(xdim,LO):cg%lhn(xdim,HI)),omegad(cg%lhn(xdim,LO):cg%lhn(xdim,HI)), noise(3, cg%lhn(zdim,LO):cg%lhn(zdim,HI)))
!   Secondary parameters

         sqr_gm = sqrt(newtong*ptmass)

         n = 0.5*Rin / (R0 - Rin)
         ninv = 1./n
         if (sigma_model == 'hayashi') then
            sigma0 = 0.2* R0**(-1.5)
         endif
         H0 = R0 * HtoR
         flind%neu%cs2 = H0 * (pi*newtong) * sigma0
         flind%neu%cs = sqrt(flind%neu%cs2)

         rho0 = sigma0 / (sqrt(dpi)*H0)

         norm = 1. / dens_Rdistr(R0,Rin,n)

         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            zk = cg%z(k)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               yj = cg%y(j)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  xi = cg%x(i)
                  rc = sqrt(xi**2+yj**2)
                  H = HtoR * rc
                  cg%u(flind%neu%idn,i,j,k) = max(rho0 * norm * dens_RdistR(rc,Rin,n) * exp(- 0.25 * zk**2 / H**2 ), smalld)
                  cg%u(flind%dst%idn,i,j,k) = eps*cg%u(flind%neu%idn,i,j,k)
               enddo
            enddo
         enddo

         do i = cg%lhn(xdim,LO)+1, cg%lhn(xdim,HI)-1   ! 2d
            rc= cg%x(i)*sqrt(2.0)
            gradgp=  0.5*(cg%gp(i+1,i+1,mnz)-cg%gp(i-1,i-1,mnz))/cg%dx/sqrt(2.)
            gradp = -0.5*(cg%u(flind%neu%idn,i+1,i+1,mnz)-cg%u(flind%neu%idn,i-1,i-1,mnz))/cg%dx /sqrt(2.)*flind%neu%cs2
            omega(i)  = sqrt( abs( (gradgp-gradp)/rc ) )
            omegad(i) = sqrt( abs(    gradgp/rc      ) )
         enddo
         omega( cg%lhn(xdim,LO)) = omega( cg%lhn(xdim,LO)+1); omega( cg%lhn(xdim,HI))  = omega(cg%lhn(xdim,HI)-1)
         omegad(cg%lhn(xdim,LO)) = omegad(cg%lhn(xdim,LO)+1); omegad(cg%lhn(xdim,HI)) = omegad(cg%lhn(xdim,HI)-1)

         call random_seed()

         do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
            yj = cg%y(j)
            do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
               xi = cg%x(i)
               rc = sqrt(xi**2+yj**2)
               call random_number(noise)

               ilook = (rc-dom%edge(xdim, LO))/cg%dx/sqrt(2.) + 0.5 + dom%nb + cg%lhn(xdim,LO)-1
               iOmega = omega(int(ilook))+(rc-cg%x(int(ilook))*sqrt(2.))*(omega(int(ilook)+1)-omega(int(ilook))) &
                    &   / (cg%x(int(ilook)+1)-cg%x(int(ilook)))/sqrt(2.)
!
!
               cg%u(flind%neu%imx,i,j,:) = -yj*iOmega*cg%u(flind%neu%idn,i,j,:)
               cg%u(flind%neu%imy,i,j,:) =  xi*iOmega*cg%u(flind%neu%idn,i,j,:)
               cg%u(flind%neu%imz,i,j,:) = 0.0
#ifndef ISO
               cg%u(flind%neu%ien,i,j,:) = flind%neu%cs2/(flind%neu%gam-1.0)*cg%u(flind%neu%idn,i,j,:)
               cg%u(flind%neu%ien,i,j,:) = max(cg%u(flind%neu%ien,i,j,:), smallei)
               cg%u(flind%neu%ien,i,j,:) = cg%u(flind%neu%ien,i,j,:) +0.5*(vx**2+vy**2+vz**2)*cg%u(flind%neu%idn,i,j,:)
#endif /* !ISO */

               iOmega = omegad(int(ilook))+(rc-cg%x(int(ilook))*sqrt(2.))*(omegad(int(ilook)+1)-omegad(int(ilook))) &
                    &   /(cg%x(int(ilook)+1)-cg%x(int(ilook)))/sqrt(2.)
               cg%u(flind%dst%imx,i,j,:) = -yj*iOmega*cg%u(flind%dst%idn,i,j,:) + amp*(noise(1,:)-0.5)
               cg%u(flind%dst%imy,i,j,:) =  xi*iOmega*cg%u(flind%dst%idn,i,j,:) + amp*(noise(2,:)-0.5)
               cg%u(flind%dst%imz,i,j,:) = 0.0 + amp*(noise(3,:)-0.5)

            enddo
         enddo

         deallocate(omega, omegad, noise)

         cgl => cgl%nxt
      enddo

   end subroutine problem_initial_conditions

end module initproblem
