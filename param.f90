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
!!

module param
  use types, only : rprec, point3D_t
#ifdef PPMPI
  use mpi
#endif
  implicit none

  save

  private rprec
  public

!*******************************************************************************
!***  ALL NON-PARAMETER DEFINITIONS READ BY INPUT FILE MUST BE INITIALIZED   ***
!*******************************************************************************

!---------------------------------------------------
! GLOBAL PARAMETERS
!---------------------------------------------------  
  integer, parameter :: CHAR_BUFF_LENGTH = 1024 ! Default size of string buffers with unknown length
  character(*), parameter :: PATH = './'
  character(*), parameter :: checkpoint_file = path // 'vel.out'
  character(*), parameter :: checkpoint_tavg_file = path // 'tavg.out'
#ifdef PPOUTPUT_EXTRA 
  character(*), parameter :: checkpoint_tavg_sgs_file = path // 'tavg_sgs.out'
#endif
  character(*), parameter :: checkpoint_spectra_file = path // 'spectra.out'

!---------------------------------------------------
! MPI PARAMETERS
!---------------------------------------------------

#ifdef PPMPI
  integer :: status(MPI_STATUS_SIZE)
  logical, parameter :: USE_MPI = .true.
  integer, parameter :: lbz = 0  ! overlap level for MPI transfer
#else
  logical, parameter :: USE_MPI = .false.
  integer, parameter :: lbz = 1  ! no overlap level necessary
#endif

  !--this stuff must be defined, even if not using MPI
  ! Setting defaults for ones that can be used even with no MPI
  integer :: nproc = 1 !--this must be 1 if no MPI
  integer :: rank = 0   !--init to 0 (so its defined, even if no MPI)
  integer :: coord = 0  !--same here

  character (8) :: chcoord  !--holds character representation of coord
  integer :: ierr
  integer :: comm
  integer :: up, down
  integer :: global_rank
  integer :: MPI_RPREC, MPI_CPREC
  integer, allocatable, dimension(:) ::  rank_of_coord, coord_of_rank
  integer :: jzmin, jzmax  ! levels that "belong" to this processor, set w/ grid

!---------------------------------------------------
! COMPUTATIONAL DOMAIN PARAMETERS
!---------------------------------------------------  

  integer, parameter :: iBOGUS = -1234567890  !--NOT a new Apple product
  real (rprec), parameter :: BOGUS = -1234567890._rprec
  real(rprec),parameter::pi=3.1415926535897932384626433_rprec

  integer :: Nx=64, Ny=64, Nz=64
  integer :: nz_tot = 64
  integer :: nx2, ny2
  integer :: lh, ld, lh_big, ld_big

  ! this value is dimensional [m]:
  real(rprec) :: z_i = 1000.0_rprec
    
  ! these values should be non-dimensionalized by z_i: 
  ! set as multiple of BL height (z_i) then non-dimensionalized by z_i
  logical :: uniform_spacing = .false.
  real(rprec) :: L_x = 2.0*pi, L_y=2.0*pi, L_z=1.0_rprec

  ! these values are also non-dimensionalized by z_i:
  real(rprec) :: dx, dy, dz
  
!---------------------------------------------------
! MODEL PARAMETERS
!---------------------------------------------------   

  ! Model type: 1->Smagorinsky; 2->Dynamic; 3->Scale dependent
  !             4->Lagrangian scale-sim   5-> Lagragian scale-dep
  integer :: sgs_model=5, wall_damp_exp=2

  ! timesteps between dynamic Cs updates           
  integer :: cs_count = 5

  ! When to start dynamic Cs calculations
  integer :: DYN_init = 100
  
  ! Cs is the Smagorinsky Constant
  ! Co and wall_damp_exp are used in the mason model for smagorisky coeff
  real(rprec) :: Co = 0.16_rprec
  
  ! test filter type: 1->cut off 2->Gaussian 3->Top-hat
  integer :: ifilter = 1

  ! u_star=0.45 m/s if coriolis_forcing=.FALSE. and =ug if coriolis_forcing=.TRUE.
  real(rprec) :: u_star = 0.45_rprec

  ! von Karman constant     
  real(rprec) :: vonk = 0.4_rprec
  
  ! Coriolis stuff
  ! coriol=non-dim coriolis parameter,
  ! ug=horiz geostrophic vel, vg=transverse geostrophic vel
  logical :: coriolis_forcing = .true. 
  real(rprec), dimension(3) :: coriol(3) = (/0.0_rprec, 0.0_rprec, -0.04_rprec/)
  real(rprec) :: ug=1.0_rprec, vg=0.0_rprec

  ! nu_molec is dimensional m^2/s
  real(rprec) :: nu_molec = 1.14e-5_rprec
    
  logical :: molec=.false., sgs=.true., dns_bc=.false.
  
!---------------------------------------------------
! TIMESTEP PARAMETERS
!---------------------------------------------------   

  integer :: nsteps = 50000
  ! -- Maximum runtime in seconds. Simulation will exit if exceeded. (disabled by default)
  integer :: runtime = -1 

  logical :: use_cfl_dt = .false.  
  real(rprec) :: cfl = 0.05
  real(rprec) :: dt_f=2.0e-4, cfl_f=0.05

  real(rprec) :: dt = 2.0e-4_rprec
  real(rprec) :: dt_dim
  
  ! time advance parameters (Adams-Bashforth, 2nd order accurate)
  real(rprec) :: tadv1, tadv2
  
  logical :: cumulative_time = .true.
  character (*), parameter :: fcumulative_time = path // 'total_time.dat'
  
  integer :: jt=0                 ! Local time counter
  integer :: jt_total=0           ! Global time counter
  real(rprec) :: total_time, total_time_dim
  
!---------------------------------------------------
! BOUNDARY/INITIAL CONDITION PARAMETERS
!---------------------------------------------------  

  ! initu = true to read from a file; false to create with random noise
  logical :: initu = .false.
  ! initlag = true to initialize cs, FLM & FMM; false to read from vel.out
  logical :: inilag = .true.

  ! lbc: lower boundary condition:  0 - stress free, 1 - wall 
  integer :: lbc_mom = 1
  ! Xinyi - ubc: upper boundary condition:  0 - stress free, 1 - wall 
  ! integral wall model is not added till now for upper boundary condition 
  integer :: ubc_mom = 1
  
  ! lower boundary condition, roughness length
  real(rprec) :: zo = 0.0001_rprec ! nondimensional
  !xiang 
  real(rprec), allocatable, dimension(:,:) :: u1flt, v1flt

  ! prescribed inflow:   
  logical :: inflow = .false.
  ! Xinyi - intial : Only need initial condition for inflow
  logical :: initial_uniform = .false.
  ! if inflow is true the following should be set:
    ! position of right end of fringe region, as a fraction of L_x
    real(rprec) :: fringe_region_end  = 1.0_rprec
    ! length of fringe region as a fraction of L_x
    real(rprec) :: fringe_region_len = 0.125_rprec

    ! Use uniform inflow instead of concurrent precursor inflow
    logical :: uniform_inflow = .false.
      real(rprec) :: inflow_velocity = 1.0_rprec

  ! if true, imposes a pressure gradient in the x-direction to force the flow
  logical :: use_mean_p_force = .true.
  ! Specify whether mean_p_force should be evaluated as 1/L_z
  logical :: eval_mean_p_force = .false. 
  real(rprec) :: mean_p_force = 1.0_rprec
  
  ! Xinyi : initially scaling constant for mean velocity
  real(rprec) :: initial_velocity_const = 1.0_rprec

!---------------------------------------------------
! DATA OUTPUT PARAMETERS
!---------------------------------------------------

  ! how often to display stdout
  integer :: wbase = 100

  ! how often to write ke to check_ke.out
  integer :: nenergy = 100

  ! how often to display Lagrangian CFL condition of 
  ! dynamic SGS models
  integer :: lag_cfl_count = 1000

  ! Flags for controling checkpointing data
  logical :: checkpoint_data = .false.
  integer :: checkpoint_nskip = 10000

  ! records time-averaged data to files ./output/ *_avg.dat
  logical :: tavg_calc = .false.
  integer :: tavg_nstart = 1, tavg_nend = 50000, tavg_nskip = 100

  ! turns instantaneous velocity recording on or off
  logical :: point_calc = .false.
  integer :: point_nstart=1, point_nend=50000, point_nskip=10
  integer :: point_nloc=1
  type(point3D_t), allocatable, dimension(:) :: point_loc

  ! domain instantaneous output
  logical :: domain_calc=.false.
  integer :: domain_nstart=10000, domain_nend=50000, domain_nskip=10000
  
  ! x-plane instantaneous output
  logical :: xplane_calc=.false.
  integer :: xplane_nstart=10000, xplane_nend=50000, xplane_nskip=10000
  integer :: xplane_nloc=1
  real(rprec), allocatable, dimension(:) :: xplane_loc

  ! y-plane instantaneous output
  logical :: yplane_calc=.false.
  integer :: yplane_nstart=10000, yplane_nend=50000, yplane_nskip=10000
  integer :: yplane_nloc=1
  real(rprec), allocatable, dimension(:) :: yplane_loc

  ! z-plane instantaneous output
  logical :: zplane_calc=.false.
  integer :: zplane_nstart=10000, zplane_nend=50000, zplane_nskip=10000
  integer :: zplane_nloc=1
  real(rprec), allocatable, dimension(:) :: zplane_loc

  logical :: spectra_calc=.false.
  integer :: spectra_nstart=10000, spectra_nend=50000, spectra_nskip=100
  integer :: spectra_nloc=1
  real(rprec), allocatable, dimension(:) :: spectra_loc

  ! Outputs histograms of {Cs^2, Tn, Nu_t, ee} for z-plane locations given below
  logical :: sgs_hist_calc = .false.
  logical :: sgs_hist_cumulative = .false.
  integer :: sgs_hist_nstart = 5000, sgs_hist_nskip = 2
  integer :: sgs_hist_nloc = 3
  real(rprec), allocatable, dimension(:) :: sgs_hist_loc   ! size=sgs_hist_nloc

  ! Limits for Cs^2 (square of Smagorinsky coefficient)
  real(rprec) :: cs2_bmin = 0.0_rprec, cs2_bmax = 0.15_rprec
  integer :: cs2_nbins = 1000

  ! Limits for Tn (Lagrangian time-averaging timescale, models 4,5 only)
  real(rprec) :: tn_bmin = 0.0_rprec, tn_bmax = 1.0_rprec
  integer :: tn_nbins = 1000

  ! Limits for Nu_t (turbulent eddy viscosity)
  real(rprec) :: nu_bmin = 0.0_rprec, nu_bmax = 0.03_rprec
  integer :: nu_nbins = 1000

  ! Limits for ee (error for Germano identity)
  real(rprec) :: ee_bmin = 0.0_rprec, ee_bmax = 100.0_rprec
  integer :: ee_nbins = 10000

end module param

!module for iwmles, each variable is explained!
module iwmles

  use types,only:rprec

  implicit none
  
  integer :: iwm_on
  !                                                  u_tau,x  u_tau,y  tau_wall,x tau_wall,y
  real(kind=rprec), dimension(:,:),   allocatable :: iwm_utx, iwm_uty, iwm_tauwx, iwm_tauwy
  !                                                  filtered tangential velocity, current and previous
  real(kind=rprec), dimension(:,:,:), allocatable :: iwm_flt_tagvel, iwm_flt_tagvel_m
  !                                                  filtered pressure
  real(kind=rprec), dimension(:,:),   allocatable :: iwm_flt_p
  integer :: iwm_dirx=1   !direction x
  integer :: iwm_diry=2   !direction y
  integer :: iwm_DN  =2   !dimension of a surface, wall model always deal with 2D surfaces (because the world is 3D)
  !                                                  integrated profiles, current and previous
  real(kind=rprec), dimension(:,:,:), allocatable :: iwm_inte, iwm_inte_m
  integer :: iwm_Lu  = 1  ! index for integral of u
  integer :: iwm_Luu = 2  ! index for integral of uu
  integer :: iwm_Lv  = 3  ! etc.
  integer :: iwm_Lvv = 4
  integer :: iwm_Luv = 5
  integer :: iwm_LN  = 5  ! the total number of integrals that need to be calculated
  !                                                  unsteady term, convective term, pressure gradient term, turbulent diffusion term, LHS
  real(kind=rprec), dimension(:,:,:), allocatable :: iwm_unsdy, iwm_conv, iwm_PrsGrad, iwm_diff, iwm_LHS
  !                                                  dudz at z=dz/2, dudz at z=zo
  real(kind=rprec), dimension(:,:,:), allocatable :: iwm_dudzT, iwm_dudzB
  !                                                  filtered friction velocity, filtering time scale
  real(kind=rprec), dimension(:,:),   allocatable :: iwm_flt_us, iwm_tR
  !                                                  HALF cell height, zo, linear correction in x, y directions
  real(kind=rprec), dimension(:,:),   allocatable :: iwm_Dz, iwm_z0, iwm_Ax, iwm_Ay

  integer :: iwm_ntime_skip=5   ! number of time steps to skip in between wall stress calculation
  real(kind=rprec) :: iwm_dt    ! time step size seen by the wall model
  integer :: iwm_debug    = 2014  !the file I use to track the performance, this is the year this model is developed
  integer :: iwm_status   = 2016  !the file for check point use, this is the year I graduate...
  
end module iwmles


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Xinyi - WM - rotate : to implement different wall models
!                       for rotating channel
!                       1. Loppi's model - 2018/10
! Here we initialize all parameters used for wall modelling itself
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module wmrotate_param

use types, only: rprec
use param, only: coriolis_forcing, coriol, dz

implicit none

save

public

! If we use a geometric series for points in the wall
  integer       :: wmro_z_n     = 40
  integer       :: wmro_z_ni    = 81
  integer       :: wmro_nstep   = 5
  real(rprec)   :: wmro_z_r     = 1.05_rprec
  real(rprec)   :: Omega        = 0.0_rprec
  real(rprec), allocatable, dimension(:)        :: wmro_zdis
  real(rprec), allocatable, dimension(:)        :: wmro_z

! parameters for equilibrium model
  real(rprec), allocatable, dimension(:)        :: wmro_alpha, wmro_alpha_stable

! parameters for Loppi's model
  real(rprec)   :: wmro_beta    = 3.6_rprec
  real(rprec)   :: wmro_S_n
  real(rprec)   :: A_pressure   = 11.0_rprec, A_stable  = 80.0_rprec
  real(rprec)   :: A_equil      = 17.0_rprec

contains

  subroutine wmrotate_init

  use param, only : vonk
  use param, only : z_i, nu_molec, u_star

  implicit none
  real(rprec)   :: nu
  real(rprec)   :: tmp
  integer       :: i

  nu = nu_molec/(z_i*u_star)
  if(coriolis_forcing) Omega = -coriol(2)/2.0_rprec
  if(abs(Omega) < 1.0e-6_rprec) then
     A_pressure = A_equil
     A_stable = A_equil
  end if
  if(Omega < -1.0e-6_rprec) then
     tmp = A_pressure
     A_pressure = A_stable
     A_stable = tmp
  end if
  wmro_S_n = 2.0_rprec * Omega

  allocate(wmro_zdis(wmro_z_n))
  allocate(wmro_z(wmro_z_ni))
  allocate(wmro_alpha(wmro_z_ni))
  allocate(wmro_alpha_stable(wmro_z_ni))

  wmro_z(1:wmro_z_ni:2) = (/( wmro_z_r**i, i=0,wmro_z_n,1 )/)
  wmro_z(1:wmro_z_ni:2) = 0.5_rprec*dz * &
     (1.0_rprec - wmro_z(1:wmro_z_ni:2))/(1.0_rprec - wmro_z_r**wmro_z_n)
  wmro_z(2:(wmro_z_ni-1):2) = (wmro_z(1:(wmro_z_ni-2):2) + wmro_z(3:wmro_z_ni:2))/2.0_rprec
  wmro_zdis = wmro_z(3:(wmro_z_ni):2) - wmro_z(1:(wmro_z_ni-2):2)

  end subroutine wmrotate_init

  subroutine wmrotate_finalize

  implicit none

  deallocate(wmro_zdis)
  deallocate(wmro_z)
  deallocate(wmro_alpha)
  deallocate(wmro_alpha_stable)

  end subroutine wmrotate_finalize

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Xinyi - WM - ML     : The machine learning result of 
!                       wall model parameters
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Xinyi : This subroutine is only for the subroutine wallstress_ML
function wallstress_ML_cal(ML_in)
use types,only:rprec
implicit none
real(rprec),dimension(2,1),intent(in)   :: ML_in
real(rprec),dimension(1)                :: ML_out
real(rprec)                             :: wallstress_ML_cal
real(rprec):: W1(4,2),W2(4,4),W3(2,4),W4(1,2),B1(4,1),B2(4,1),B3(2,1),B4(1,1)

! Define all the parameters
W1(:,1) = (/-0.000096530509586, -0.005251070687483, -0.002346407669290, -0.006090122176649/)
W1(:,2) = (/0.027263245421098, 2.320018703628134, -0.710667103184319, 0.202751136149669/)
W2(1,:) = (/-0.614244666114551, -0.56875828110715, -768.123583340675, -2.88137027368678/)
W2(2,:) = (/-29.2325471159253, 0.61874311769997, 606.862314031602, -0.979546178847371/)
W2(3,:) = (/-800.686036202809, 540.281718692259, -542.160832895298, -823.615836005128/)
W2(4,:) = (/8.87452448422246, -1.57586699125869, -254.733947116813, 2.8570684215825/)
W3(1,:) = (/-8.61032891969321, -52.5081877227373, 238.702964339122, -378.726549169011/)
W3(2,:) = (/5.69137173878176, 6.43463663134031, 15.2708054748494, 5.6704967141009/)
W4(1,:) = (/0.243608988604323, 3.64392478869163/)

B1(:,1) = (/2.46643851383815, -0.556583366474677, -3.03634911516902, -0.624446440839603/)
B2(:,1) = (/-766.225624615761, 632.606691471833, 544.387077019739, -262.657497454692/)
B3(:,1) = (/70.3836160446785, -9.05649341076006/)
B4(1,1) = 1.39916470522558

! Calculate the result
ML_out = reshape( matmul(W4, tanh(matmul(W3, tanh(matmul(W2, &
       & tanh(matmul(W1,ML_in)+B1) )+B2) )+B3) )+B4 , (/1/))
wallstress_ML_cal = ML_out(1)

end function wallstress_ML_cal

end module  wmrotate_param
