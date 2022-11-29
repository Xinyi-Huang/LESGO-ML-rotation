!!
!!  Copyright (C) 2009-2013  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
module sim_param
use types, only : rprec
use param, only : ld, ny, nz, lbz
use param, only : zo, dz, vonk, u1flt, v1flt, nx !xiang
implicit none

save
public

logical :: sim_param_initialized = .false.

real (rprec), dimension (:,:), allocatable :: up_swm, vp_swm !xiang dirty swm for dudt
real (rprec), dimension (:, :, :), allocatable :: u, v, w
real (rprec), dimension (:, :, :), allocatable :: dudx, dudy, dudz,  &
                                                  dvdx, dvdy, dvdz,  &
                                                  dwdx, dwdy, dwdz,  &
                                                  RHSx, RHSy, RHSz,  &
                                                  RHSx_f, RHSy_f, RHSz_f

real (rprec), dimension (:, :, :), allocatable :: dpdx, dpdy, dpdz

real (rprec), dimension (:, :, :), allocatable :: txx, txy, tyy
real (rprec), dimension (:, :, :), allocatable :: txz, tyz, tzz

real (rprec), target, dimension (:, :, :), allocatable :: p

real (rprec), dimension (:, :, :), allocatable :: divtx, divty, divtz

real (rprec), dimension (:, :, :), allocatable :: fx, fy, fz, &
                                                  fxa, fya, fza
contains

!
! This subroutine initilizes all global arrays defined in the sim_param
! module. Here they are allocated and initialized to zero.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine sim_param_init ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
allocate ( u(ld, ny, lbz:nz) ); u = 0.0_rprec
allocate ( v(ld, ny, lbz:nz) ); v = 0.0_rprec
allocate ( w(ld, ny, lbz:nz) ); w = 0.0_rprec
allocate( dudx(ld, ny, lbz:nz) ); dudx = 0.0_rprec
allocate( dudy(ld, ny, lbz:nz) ); dudy = 0.0_rprec
allocate( dudz(ld, ny, lbz:nz) ); dudz = 0.0_rprec
allocate( dvdx(ld, ny, lbz:nz) ); dvdx = 0.0_rprec
allocate( dvdy(ld, ny, lbz:nz) ); dvdy = 0.0_rprec
allocate( dvdz(ld, ny, lbz:nz) ); dvdz = 0.0_rprec
allocate( dwdx(ld, ny, lbz:nz) ); dwdx = 0.0_rprec
allocate( dwdy(ld, ny, lbz:nz) ); dwdy = 0.0_rprec
allocate( dwdz(ld, ny, lbz:nz) ); dwdz = 0.0_rprec
allocate( RHSx(ld, ny, lbz:nz) ); RHSx = 0.0_rprec
allocate( RHSy(ld, ny, lbz:nz) ); RHSy = 0.0_rprec
allocate( RHSz(ld, ny, lbz:nz) ); RHSz = 0.0_rprec
allocate( RHSx_f(ld, ny, lbz:nz) ); RHSx_f = 0.0_rprec
allocate( RHSy_f(ld, ny, lbz:nz) ); RHSy_f = 0.0_rprec
allocate( RHSz_f(ld, ny, lbz:nz) ); RHSz_f = 0.0_rprec
allocate ( dpdx(ld, ny, nz) ); dpdx = 0.0_rprec
allocate ( dpdy(ld, ny, nz) ); dpdy = 0.0_rprec
allocate ( dpdz(ld, ny, nz) ); dpdz = 0.0_rprec
allocate ( txx(ld, ny, lbz:nz) ); txx = 0.0_rprec
allocate ( txy(ld, ny, lbz:nz) ); txy = 0.0_rprec
allocate ( tyy(ld, ny, lbz:nz) ); tyy = 0.0_rprec
allocate ( txz(ld, ny, lbz:nz) ); txz = 0.0_rprec
allocate ( tyz(ld, ny, lbz:nz) ); tyz = 0.0_rprec
allocate ( tzz(ld, ny, lbz:nz) ); tzz = 0.0_rprec
allocate ( p(ld, ny, 0:nz) ); p = 0.0_rprec
allocate ( divtx(ld, ny, lbz:nz) ); divtx = 0.0_rprec
allocate ( divty(ld, ny, lbz:nz) ); divty = 0.0_rprec
allocate ( divtz(ld, ny, lbz:nz) ); divtz = 0.0_rprec

#ifdef PPTURBINES
allocate ( fxa(ld, ny, nz) ); fxa = 0.0_rprec
#endif

#if defined(PPLVLSET) || defined(PPATM)
allocate ( fx(ld, ny, nz) ); fx = 0.0_rprec
allocate ( fy(ld, ny, nz) ); fy = 0.0_rprec
allocate ( fz(ld, ny, nz) ); fz = 0.0_rprec
! May already be allocated if using TURBINES
if( .not. allocated(fxa) ) allocate ( fxa(ld, ny, nz) ); fxa = 0.0_rprec
allocate ( fya(ld, ny, nz) ); fya = 0.0_rprec
allocate ( fza(ld, ny, nz) ); fza = 0.0_rprec
#endif 

!xiang
allocate( u1flt(nx,ny) ); u1flt=1/vonk*log(0.5*dz/zo)
allocate( v1flt(nx,ny) ); v1flt=0.0

sim_param_initialized = .true.

return
end subroutine sim_param_init
end module sim_param
