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
!--------------------------------------------------------------------------------------------------------------

#include "piernik.def"

module interpolations
! pulled by RIEMANN

  implicit none
  private
  public :: set_interpolations, interpol

  interface
     subroutine interpolation(prim_var,prim_var_l,prim_var_r)
       
       implicit none
       
       real, dimension(:,:), intent(in)  :: prim_var
       real, dimension(:,:), intent(out) :: prim_var_l
       real, dimension(:,:), intent(out) :: prim_var_r

     end subroutine interpolation

  end interface

  procedure(interpolation), pointer :: interp => null()

contains

   function utoq(u,b_cc) result(q)

     use constants,  only: half, xdim, zdim
     use fluidindex, only: flind
     use fluidtypes, only: component_fluid
     use func,       only: ekin

     implicit none

     real, dimension(:,:),   intent(in)    :: u , b_cc

     real, dimension(size(u,1),size(u,2))  :: q
     integer  :: p

     class(component_fluid), pointer       :: fl

     do p = 1, flind%fluids
        fl => flind%all_fluids(p)%fl

        q(fl%idn,:) =  u(fl%idn,:)
        q(fl%imx,:) =  u(fl%imx,:)/u(fl%idn,:)
        q(fl%imy,:) =  u(fl%imy,:)/u(fl%idn,:)
        q(fl%imz,:) =  u(fl%imz,:)/u(fl%idn,:)
        ! J.CoPhy 208 (2005),Pg 317, Eq. 2. Gas pressure: p = (gamma-1)*(e-half*rho*v^2-half*B^2) and Total pressure: p_T = p + half*B^2. (1) and (2) are markers for HD and MHD.
        if (fl%has_energy) then
            q(fl%ien,:) =  fl%gam_1*(u(fl%ien,:) - ekin(u(fl%imx,:), u(fl%imy,:), u(fl%imz,:), u(fl%idn,:))) ! Primitive variable for gas pressure (p) without magnetic fields. (1)
            if (fl%is_magnetized) then
               q(fl%ien,:) =  q(fl%ien,:) - half*fl%gam_1*sum(b_cc(xdim:zdim,:)**2, dim=1) ! Primitive variable for gas pressure (p) with magnetic fields. The requirement of total pressure is dealt in the fluxes and hlld routines. (2)
            endif
        endif

     enddo

   end function utoq

  subroutine interpol(u, bcc, psi, ql, qr, bccl, bccr, psil, psir)

    implicit none

    real, dimension(:,:), intent(in)     :: u
    real, dimension(:,:), intent(out)    :: ql
    real, dimension(:,:), intent(out)    :: qr

    real, dimension(:,:), intent(in)     :: bcc
    real, dimension(:,:), intent(out)    :: bccl
    real, dimension(:,:), intent(out)    :: bccr

    real, dimension(:,:), intent(in)     :: psi
    real, dimension(:,:), intent(out)    :: psil
    real, dimension(:,:), intent(out)    :: psir

    real, dimension(size(u, 1), size(u, 2)) :: q

    q = utoq(u, bcc)
    call interp(q, ql, qr)
    call interp(bcc, bccl, bccr)
    call interp(psi, psil, psir)

  end subroutine interpol

  subroutine set_interpolations(interpol_str)

    use dataio_pub, only: msg, die

    implicit none

    character(len=*), intent(in) :: interpol_str

    if (associated(interp)) call die("[interpolations:set_interpolations] interp already associated")
    
    select case(interpol_str)
    case('linear', 'LINEAR')
       interp => linear
    case default
       write(msg,'(2a)') "[interpolations:set_interpolations] unknown interpolation ", interpol_str 
       call die(msg)
       interp => null()
    end select
    
  end subroutine set_interpolations

  subroutine linear(q,ql,qr)
    
    use fluxlimiters, only: set_limiters,flimiter,blimiter
    use domain,       only: dom
    use constants,    only: half, GEO_XYZ
    use dataio_pub,   only: die
      
    implicit none

    real, dimension(:,:), intent(in)     :: q
    real, dimension(:,:), intent(out)    :: ql
    real, dimension(:,:), intent(out)    :: qr

    real, dimension(size(q,1),size(q,2)) :: dq_lim, dq_interp

    integer                              :: im1
    integer                              :: n 
    n = size(q,2)
      
    dq_lim = flimiter(q)

    if (dom%geometry_type /= GEO_XYZ) call die("[interpolations:linear] non-cartesian geometry not implemented yet.")

    dq_interp = half*dq_lim
    
    ql = q + dq_interp
    qr = q - dq_interp

    im1 = max(1,n-1) ! neighbouring indices
    !associate(im1 => i - Dom%D_x)
    ! shfit right state
    qr(:,1:im1) = qr(:,2:n)


    ! interpolation for the first and last points
    
    ql(:,1) = q(:,1)
    qr(:,n) = q(:,n)  
      
  end subroutine linear
    
end module interpolations