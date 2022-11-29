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

! For use with staggered grid LES
! JDA, 23 Jan 96
!--provides txz, tyz (w-nodes) and dudz, dvdz (w-nodes) at jz=1

!xiang: when this function is called, it already ensures coord=0
subroutine wallstress ()
use grid_m, only : grid
use types,only:rprec
use param,only:dz,ld,lh,nx,ny,nz,vonk,lbc_mom,zo, domain_nskip, jt_total
use param,only: dt, u1flt, v1flt
use sim_param,only:u,v,dudz,dvdz,txz,tyz
use test_filtermodule
use iwmles, only : iwm_on !xiang for the option of integral wall model
implicit none
integer::jx,jy
real(rprec),dimension(nx,ny)::denom,u_avg,ustar
real(rprec),dimension(ld,ny)::u1,v1
! No need to define phi_m or psi_m as a matrix as only a constant value is used
real(rprec)::const,phi_m,psi_m
integer :: iswm_on, digits_jt
CHARACTER*50 :: fname
integer, pointer, dimension(:) :: autowrap_i, autowrap_j !useful array for autowraped index
real(rprec) :: Cd
integer :: jjx, jjy, ctjj
real, dimension(:,:,:), allocatable :: randnum
real, dimension(:,:), allocatable :: tfluc
real :: Tflt, eps_flt, sumtmp
iswm_on=0
Cd=0.0


select case (lbc_mom)
  case (0) ! Stress free
    txz(:, :, 1) = 0._rprec
    tyz(:, :, 1) = 0._rprec
    dudz(:, :, 1) = 0._rprec
    dvdz(:, :, 1) = 0._rprec
  case (1) ! Wall
   if(iwm_on == 0)then !if not using integral wall model...
     if(iswm_on == 0)then
      ! See John D. Albertson's dissertation, eqns (2.46)-(2.52)
      ! For dudz and dvdz at wall, we should use derivwall.f90
      ! Also, see:
      ! E. Bou-Zeid, C. Meneveau & M.B. Parlange, "A scale-dependent Lagrangian dynamic model
      !   for large eddy simulation of complex turbulent flows" (2005) -- Appendix    
     
      !TS Remove the following line when obukhov.f is used
      psi_m=0._rprec
      phi_m=1._rprec
  
      u1=u(:,:,1)
      v1=v(:,:,1)
      
      !test the 2nd point==================
      !u1=u(:,:,2)
      !v1=v(:,:,2)
      !test temporal filtering=============
      !Tflt=dz/2.0/1.0/vonk
      !do jx=1,nx
      !do jy=1,ny
      !    eps_flt=dt/Tflt
      !    if(eps_flt>1.0_rprec)then
      !        eps_flt=1.0_rprec
      !    end if
      !    u1flt(jx,jy)=u1flt(jx,jy)*(1.0_rprec-eps_flt)+u1(jx,jy)*eps_flt
      !    v1flt(jx,jy)=v1flt(jx,jy)*(1.0_rprec-eps_flt)+v1(jx,jy)*eps_flt
      !end do
      !end do
      !do jx=1,nx
      !do jy=1,ny
      !    u1(jx,jy)=u1flt(jx,jy)
      !    v1(jx,jy)=v1flt(jx,jy)
      !end do
      !end do
      !original filtering
      call test_filter ( u1 )
      call test_filter ( v1 )
      !test spatial filtering===========================
      !nullify(autowrap_i, autowrap_j)
      !autowrap_i => grid % autowrap_i
      !autowrap_j => grid % autowrap_j
      !do jy=1,ny
      !do jx=1,nx
      !   u1(jx,jy)=0.0_rprec
      !   v1(jx,jy)=0.0_rprec
      !   do jjx=-1,1
      !   do jjy=-1,1
      !       u1(jx,jy)=u1(jx,jy)+u(autowrap_i(jx+jjx),autowrap_j(jy+jjy),1)
      !       v1(jx,jy)=v1(jx,jy)+v(autowrap_i(jx+jjx),autowrap_j(jy+jjy),1)
      !   end do
      !   end do
      !   u1(jx,jy)=u1(jx,jy)/9.0_rprec
      !   v1(jx,jy)=v1(jx,jy)/9.0_rprec
      !end do
      !end do
      
      denom=log(0.5_rprec*dz/zo)-psi_m
      !test the 2nd point ==================
      !denom=log(1.5_rprec*dz/zo)-psi_m
      u_avg=sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)
      ustar=u_avg*vonk/denom  

      !generate random number=============
      !allocate(randnum(nx,ny,2))
      !call init_random_seed
      !call random_number(randnum)
      !do jy=1,ny
      !do jx=1,nx
      !  randnum(jx,jy,1)=(randnum(jx,jy,1)-0.5)*2.0*0.37
      !  randnum(jx,jy,2)=(randnum(jx,jy,2)-0.5)*2.0*0.37
      !end do
      !end do     

      !for artificial fluctuations
      !allocate(tfluc(nx,ny))
      !sumtmp=0.0_rprec
      !do jy=1,ny
      !do jx=1,nx
      !    sumtmp=sumtmp+u1(jx,jy)
      !end do
      !end do    
      !sumtmp=sumtmp/nx/ny
      !do jy=1,ny
      !do jx=1,nx
      !    tfluc(jx,jy)=-0.15_rprec*(u1(jx,jy)-1.0/vonk*log(dz/2.0/zo))
      !end do
      !end do
      do jy=1,ny
      do jx=1,nx
         const=-(ustar(jx,jy)**2) /u_avg(jx,jy)
         txz(jx,jy,1)=const *u1(jx,jy)
         tyz(jx,jy,1)=const *v1(jx,jy)
         !test random number wall-stress
         !txz(jx,jy,1)=-1.00!+tfluc(jx,jy)!+randnum(jx,jy,1)
         !tyz(jx,jy,1)=0.00!randnum(jx,jy,2)
      !TS REMOVE derivwall.f90 and add it here
      !do finite difference here============================
      !   dudz(jx,jy,1)=(u(jx,jy,2)-u(jx,jy,1))/dz
      !   dvdz(jx,jy,1)=(v(jx,jy,2)-v(jx,jy,1))/dz
      !this is as in Moeng 84
         dudz(jx,jy,1)=ustar(jx,jy)/(0.5_rprec*dz*vonK)*u(jx,jy,1)/u_avg(jx,jy)&
      !TS ADD for non-neutral case
             *phi_m
         dvdz(jx,jy,1)=ustar(jx,jy)/(0.5_rprec*dz*vonK)*v(jx,jy,1)/u_avg(jx,jy)&
      !TS ADD for non-neutral case
             *phi_m
      !   dudz(jx,jy,1)=merge(0._rprec,dudz(jx,jy,1),u(jx,jy,1).eq.0._rprec)
      !   dvdz(jx,jy,1)=merge(0._rprec,dvdz(jx,jy,1),v(jx,jy,1).eq.0._rprec)
      end do
      end do
      if( mod(jt_total,domain_nskip)==1)then
        digits_jt=int(log(real(jt_total))/log(10.0))+1
        write(fname,'(A,I<digits_jt>,A)') 'equil_output.',jt_total,'.dat' 
        open(375,file=fname)  
        do jy=1,ny 
        do jx=1,nx
          write(375,*) txz(jx,jy,1), tyz(jx,jy,1), u1(jx,jy), v1(jx,jy)
        end do
        end do
        close(375)
      endif
      !for filtering=============
      nullify(autowrap_i, autowrap_j)
      !for random number========
      !deallocate(randnum)
      !for filtering
      else
        call swm_wallstress() !xiang wallstress from the slip wall model
     end if
   else
    call iwm_wallstress() !xiang: calculate stress using the integral wall model...
   endif
end select
end subroutine wallstress


! Xinyi : This subroutine is for upper boundary layer, and almost exactly the same as
! the one for the lower boundary layer. Take care with nz and nz-1
subroutine wallstress_upper ()
use types,only:rprec
use param,only:dz,ld,lh,nx,ny,nz,vonk,ubc_mom,zo, domain_nskip, jt_total
use sim_param,only:u,v,dudz,dvdz,txz,tyz
use test_filtermodule
implicit none
integer::jx,jy
real(rprec),dimension(nx,ny)::denom,u_avg,ustar
real(rprec),dimension(ld,ny)::u1,v1
! No need to define phi_m or psi_m as a matrix as only a constant value is used
real(rprec)::const
integer :: digits_jt
CHARACTER*50 :: fname

select case (ubc_mom)
  case (0) ! Stress free
    txz(:, :, nz) = 0._rprec
    tyz(:, :, nz) = 0._rprec
    dudz(:, :, nz) = 0._rprec
    dvdz(:, :, nz) = 0._rprec
  case (1) ! Wall
  
      u1=u(:,:,nz-1)
      v1=v(:,:,nz-1)
      
      call test_filter ( u1 )
      call test_filter ( v1 )
      
      denom=log(0.5_rprec*dz/zo)
      u_avg=sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)
      ustar=u_avg*vonk/denom  

      do jy=1,ny
      do jx=1,nx
         const=(ustar(jx,jy)**2) /u_avg(jx,jy)
         txz(jx,jy,nz)=const *u1(jx,jy)
         tyz(jx,jy,nz)=const *v1(jx,jy)
         dudz(jx,jy,nz)=-ustar(jx,jy)/(0.5_rprec*dz*vonK)*u(jx,jy,nz-1)/u_avg(jx,jy)
         dvdz(jx,jy,nz)=-ustar(jx,jy)/(0.5_rprec*dz*vonK)*v(jx,jy,nz-1)/u_avg(jx,jy)
      end do
      end do
      if( mod(jt_total,domain_nskip)==1)then
        digits_jt=int(log(real(jt_total))/log(10.0))+1
        write(fname,'(A,I<digits_jt>,A)') 'equil_output_upper.',jt_total,'.dat' 
        open(375,file=fname)  
        do jy=1,ny 
        do jx=1,nx
          write(375,*) txz(jx,jy,nz), tyz(jx,jy,nz), u1(jx,jy), v1(jx,jy)
        end do
        end do
        close(375)
      endif
end select
end subroutine wallstress_upper

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Xinyi - WM - rotate : to implement different wall models
!                       for rotating channel
!                       1. Loppi's model - 2018/10
!                       2. Rotation min model - 2018/12
!                       3. Harmonic model - 2018/12
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!xiang: when this function is called, it already ensures coord=0
subroutine wallstress_wmro ()
use types,only:rprec
use param,only:dz,ld,lh,nx,ny,nz,vonk,lbc_mom,zo, domain_nskip, jt_total, jt
use param,only: nu_molec,z_i,u_star
use sim_param,only:u,v,dudz,dvdz,txz,tyz
use wmrotate_param
use test_filtermodule
implicit none
integer     ::jx,jy
real(rprec),dimension(ld,ny)    :: u1,v1
real(rprec),dimension(wmro_z_ni):: wmro_a, wmro_b, wmro_c ! quadratic parameter
real(rprec),dimension(wmro_z_ni):: wmro_mdudz, wmro_dudz, wmro_dudzf, wmro_coeff_i, wmro_Ri
real(rprec) :: const, wmro_judge, wmro_judge2, wmro_usol, wmro_coeff, u_avg
real(rprec) :: wmro_sign, wmro_nu
integer     :: digits_jt
CHARACTER*50 :: fname

wmro_nu = nu_molec/(z_i*u_star)

select case (lbc_mom)
  case (0) ! Stress free
    txz(:, :, 1) = 0._rprec
    tyz(:, :, 1) = 0._rprec
    dudz(:, :, 1) = 0._rprec
    dvdz(:, :, 1) = 0._rprec
  case (1) ! Wall
  if(modulo(jt,wmro_nstep) == 1) then
      u1=u(:,:,1)
      v1=v(:,:,1)

      !original filtering
      !call test_filter ( u1 )
      !call test_filter ( v1 )

         ! Xinyi - trick : initialization done only for the very first point
         u_avg = sqrt(u1(1,1)**2 + v1(1,1)**2)
         const = wmro_nu * (u_avg / (0.5_rprec*dz))
         !wmro_mdudz = - u_avg / (0.5_rprec*dz)
         wmro_dudz = u_avg / (0.5_rprec*dz)

      do jy=1,ny
      do jx=1,nx
         u_avg = sqrt(u1(jx,jy)**2 + v1(jx,jy)**2)
         wmro_judge = 1.0_rprec

         do while (abs(wmro_judge) > 1.0e-6_rprec)
            ! To generate damping parameters, we need to turn to the friction scale
            ! Model Loppi :
            !wmro_alpha = vonk*wmro_z*sqrt(const)*(1.0_rprec - exp(-wmro_z*sqrt(const)/wmro_nu/A_pressure))**2

            ! Model harmonic:
            wmro_judge2 = 1.0_rprec
            do while (abs(wmro_judge2) > 1.0e-6_rprec)
               wmro_dudzf = wmro_dudz
               if (Omega > 1.0e-6_rprec) then
                  wmro_alpha = 1.0_rprec / (1.0_rprec / (vonk*wmro_z*(1.0_rprec + 7.0_rprec*(2.0_rprec*Omega)/wmro_dudz)) + 1.0_rprec / (sqrt(const)/(2.0_rprec*Omega)))
               else
                  wmro_alpha = vonk*wmro_z
               end if
               wmro_alpha = 1.0_rprec / (wmro_nu + wmro_alpha * sqrt(const) * (1.0_rprec - exp(-wmro_z*sqrt(const)/wmro_nu/A_equil))**2)
               wmro_dudz = const * wmro_alpha
               wmro_judge2 = maxval(abs(wmro_dudz - wmro_dudzf))
            end do

            ! Model min:
            !if (Omega > 1.0e-6_rprec) then
            !    wmro_alpha = min(vonk*wmro_z,const / (2.0_rprec*abs(Omega)))
            !else
            !    wmro_alpha = vonk*wmro_z
            !end if
            !wmro_alpha = wmro_nu + wmro_alpha * sqrt(const) * (1.0_rprec - exp(-wmro_z*sqrt(const)/wmro_nu/A_equil))**2
            !wmro_alpha = 1.0_rprec / wmro_alpha

            ! For solving quadratic equations
            ! Model Loppi : 
            !wmro_a = wmro_nu + wmro_alpha
            !wmro_b = -(-const + wmro_beta*wmro_alpha*wmro_S_n)
            !wmro_c = -wmro_beta*wmro_alpha*wmro_S_n**2
            !wmro_Ri = wmro_S_n/wmro_mdudz
            !wmro_Ri = wmro_Ri*(wmro_Ri+1.0_rprec)
            !where(wmro_Ri*wmro_beta > 1.0_rprec) 
            !   wmro_mdudz = -const/wmro_nu
            !elsewhere(abs(wmro_Ri) < 1.0e-10_rprec)
            !   wmro_mdudz = -const/wmro_a
            !elsewhere
            !   wmro_mdudz = (-wmro_b - &
            !      sqrt(wmro_b**2-4.0_rprec*wmro_a*wmro_c))/ &
            !      (2.0_rprec*wmro_a)
            !end where
            ! Model min, harmonic : 
            !wmro_dudz = const * wmro_alpha

            ! converge for shooting method
            ! Model Loppi : 
            !wmro_usol = -sum(wmro_zdis/6.0_rprec * &
            !   (wmro_mdudz(1:(size(wmro_mdudz)-2):2) + wmro_mdudz(3::2) + &
            !   wmro_mdudz(2::2)*4.0_rprec))
            ! Model harmonic : 
            wmro_usol = sum(wmro_zdis/6.0_rprec * &
               (wmro_dudz(1:(size(wmro_dudz)-2):2) + wmro_dudz(3::2) + &
               wmro_dudz(2::2)*4.0_rprec))
            ! Model Loppi : 
            !wmro_coeff_i = wmro_mdudz**2/(wmro_mdudz**2*wmro_a - wmro_c)
            ! Model harmonic:
            wmro_coeff_i = wmro_alpha
            wmro_coeff = sum(wmro_zdis/6.0_rprec * &
               (wmro_coeff_i(1:(size(wmro_coeff_i)-2):2) + wmro_coeff_i(3::2) + &
               wmro_coeff_i(2::2)*4.0_rprec))
            ! Model Loppi : 
            !wmro_judge = (u_avg - wmro_usol)/wmro_coeff
            ! Model harmonic: 
            wmro_judge = (u_avg - wmro_usol)/wmro_coeff
            const = const + wmro_judge
         end do

         txz(jx,jy,1) = -const*u1(jx,jy)/u_avg
         tyz(jx,jy,1) = -const*v1(jx,jy)/u_avg
         !dudz(jx,jy,1) = -wmro_mdudz(wmro_z_ni)*u1(jx,jy)/u_avg
         !dvdz(jx,jy,1) = -wmro_mdudz(wmro_z_ni)*v1(jx,jy)/u_avg
         dudz(jx,jy,1) = wmro_dudz(wmro_z_ni)*u1(jx,jy)/u_avg
         dvdz(jx,jy,1) = wmro_dudz(wmro_z_ni)*v1(jx,jy)/u_avg
      end do
      end do
  end if

      ! Xinyi - Suggestion : Might want to use some matrices to store several layers as the initial
      ! condition of the next time step

      ! Output for equil_output.dat
      if( mod(jt_total,domain_nskip)==1)then
        digits_jt=int(log(real(jt_total))/log(10.0))+1
        write(fname,'(A,I<digits_jt>,A)') 'equil_output.',jt_total,'.dat' 
        open(375,file=fname)  
        do jy=1,ny 
        do jx=1,nx
          write(375,*) txz(jx,jy,1), tyz(jx,jy,1), u1(jx,jy), v1(jx,jy), wmro_nu*dudz(jx,jy,1), wmro_nu*dvdz(jx,jy,1)
        end do
        end do
        close(375)
      end if
end select
end subroutine wallstress_wmro


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Xinyi - WM - rotate : to implement different wall models
!                       from machine learning
!                       1. Developed from ML - 2019/04
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!xiang: when this function is called, it already ensures coord=0
subroutine wallstress_ML ()
use types,only:rprec
use param,only:dz,ld,lh,nx,ny,nz,vonk,lbc_mom,zo, domain_nskip, jt_total, jt
use param,only: nu_molec,z_i,u_star
use sim_param,only:u,v,dudz,dvdz,txz,tyz
use wmrotate_param
use test_filtermodule
implicit none
integer     ::jx,jy
real(rprec),dimension(ld,ny)    :: u1,v1
real(rprec) :: utaup_f, utaup, utaup_d, utaup_tmp, wmro_judge, wmro_y, wmro_dudz, u_avg
real(rprec) :: wmro_usol_f, wmro_mlin_f(2,1), wmro_usol, wmro_usol_d, wmro_mlin(2,1), wmro_flag
real(rprec) :: wmro_nu
integer     :: digits_jt
integer     :: i, j                     ! Xinyi - debug : Where is it stuck
CHARACTER*50 :: fname

wmro_nu = nu_molec/(z_i*u_star)

select case (lbc_mom)
  case (0) ! Stress free
    txz(:, :, 1) = 0._rprec
    tyz(:, :, 1) = 0._rprec
    dudz(:, :, 1) = 0._rprec
    dvdz(:, :, 1) = 0._rprec
  case (1) ! Wall
  if(modulo(jt,wmro_nstep) == 1) then
      u1=u(:,:,1)
      v1=v(:,:,1)
      wmro_y=0.5_rprec*dz

      !original filtering
      call test_filter ( u1 )
      call test_filter ( v1 )

         ! Xinyi - trick : initialization done only for the very first point
         u_avg = sqrt(u1(1,1)**2 + v1(1,1)**2)
         utaup_f = sqrt(wmro_nu * (u_avg / (0.5_rprec*dz)))

      do jy=1,ny
      do jx=1,nx
         u_avg = sqrt(u1(jx,jy)**2 + v1(jx,jy)**2)
         wmro_judge = 1.0_rprec

         ! Xinyi - debug :
         i = 0

         wmro_mlin_f = reshape((/wmro_y*utaup_f/wmro_nu, &
                     & wmro_y/utaup_f*wmro_S_n/), (/2,1/))
         if (abs(wmro_S_n) .LE. 1e-6) then
             wmro_usol_f = log(wmro_y*utaup_f/wmro_nu)/vonk
         else 
             wmro_usol_f = (sign(0.5_rprec,wmro_y-utaup_f/wmro_S_n/vonk)+0.5_rprec) &
                 & *(wmro_S_n*wmro_y/utaup_f+log(utaup_f**2/wmro_S_n/wmro_nu)/vonk-(log(vonk)+1.0_rprec)/vonk) &
                 &       + (sign(0.5_rprec,-wmro_y+utaup_f/wmro_S_n/vonk)+0.5_rprec) &
                 & *log(wmro_y*utaup_f/wmro_nu)/vonk
         end if
         wmro_usol_f = (wallstress_ML_cal(wmro_mlin_f) & 
                   & + wmro_usol_f) * utaup_f

         utaup       = 1.0_rprec

         do while (abs(wmro_judge) > 1.0e-6_rprec)
            ! To generate damping parameters, we need to turn to the friction scale

            ! Xinyi - debug :
            i = i + 1
            if(i > 1000) print *, 'Stuck for u_avg ', jx, jy, u_avg, utaup, wmro_mlin, wmro_usol, & 
                                & '//y,nu,Sn,MLcal', wmro_y, wmro_nu, wmro_S_n, wallstress_ML_cal(wmro_mlin)

            ! Model ML :
             wmro_mlin = reshape((/wmro_y*utaup/wmro_nu, &
                         & wmro_y/utaup*wmro_S_n/), (/2,1/))
             if (abs(wmro_S_n) .LE. 1e-6) then
                 wmro_usol = log(wmro_y*utaup/wmro_nu)/vonk
             else 
                 wmro_usol = (sign(0.5_rprec,wmro_y-utaup/wmro_S_n/vonk)+0.5_rprec) &
                     & *(wmro_S_n*wmro_y/utaup+log(utaup**2/wmro_S_n/wmro_nu)/vonk-(log(vonk)+1.0_rprec)/vonk) &
                     &       + (sign(0.5_rprec,-wmro_y+utaup/wmro_S_n/vonk)+0.5_rprec) &
                     & *log(wmro_y*utaup/wmro_nu)/vonk
             end if
             wmro_usol = (wallstress_ML_cal(wmro_mlin) & 
                   & + wmro_usol) * utaup

             ! In case we are stuck, use the bisection method
             if (i > 100 .and. (((wmro_usol_f - u_avg) > 0.0_rprec .and. (wmro_usol - u_avg) < 0.0_rprec) .or. ((wmro_usol_f - u_avg) < 0.0_rprec .and. (wmro_usol - u_avg) > 0.0_rprec))) then
                 exit
             end if

             ! Normal loop
             utaup_d = utaup+1.0e-10
             wmro_mlin = reshape((/wmro_y*utaup_d/wmro_nu, &
                         & wmro_y/utaup_d*wmro_S_n/), (/2,1/))
             if (abs(wmro_S_n) .LE. 1e-6) then
                 wmro_usol_d = log(wmro_y*utaup_d/wmro_nu)/vonk
             else 
                 wmro_usol_d = (sign(0.5_rprec,wmro_y-utaup_d/wmro_S_n/vonk)+0.5_rprec) &
                     & *(wmro_S_n*wmro_y/utaup_d+log(utaup_d**2/wmro_S_n/wmro_nu)/vonk-(log(vonk)+1.0_rprec)/vonk) &
                     &       + (sign(0.5_rprec,-wmro_y+utaup_d/wmro_S_n/vonk)+0.5_rprec) &
                     & *log(wmro_y*utaup_d/wmro_nu)/vonk
             end if
             wmro_usol_d = (wallstress_ML_cal(wmro_mlin) & 
                   & + wmro_usol_d) * utaup_d

            utaup_tmp   = utaup
            wmro_judge  = (u_avg - wmro_usol)
            utaup       = utaup + sign(1.0_rprec,wmro_judge)*min(abs(wmro_judge)* &
                        & abs(utaup-utaup_d)/abs(wmro_usol-wmro_usol_d), 0.1_rprec*utaup, &
                        & abs(utaup-utaup_f*3.0_rprec))*0.5_rprec

            wmro_usol_f = wmro_usol
            utaup_f     = utaup_tmp

         end do

         do while (abs(wmro_judge) > 1.0e-2_rprec)

            ! Xinyi - debug :
            i = i + 1
            if(i > 1000) print *, 'Stuck for u_avg ', jx, jy, u_avg, utaup, wmro_mlin, wmro_usol, & 
                                & '//y,nu,Sn,MLcal', wmro_y, wmro_nu, wmro_S_n, wallstress_ML_cal(wmro_mlin)

            utaup_d = 0.5_rprec * (utaup_f + utaup)
             ! In case we are stuck, use the bisection method
             wmro_mlin = reshape((/wmro_y*utaup_d/wmro_nu, &
                         & wmro_y/utaup_d*wmro_S_n/), (/2,1/))
             if (abs(wmro_S_n) .LE. 1e-6) then
                 wmro_usol_d = log(wmro_y*utaup_d/wmro_nu)/vonk
             else 
                 wmro_usol_d = (sign(0.5_rprec,wmro_y-utaup_d/wmro_S_n/vonk)+0.5_rprec) &
                     & *(wmro_S_n*wmro_y/utaup_d+log(utaup_d**2/wmro_S_n/wmro_nu)/vonk-(log(vonk)+1.0_rprec)/vonk) &
                     &       + (sign(0.5_rprec,-wmro_y+utaup_d/wmro_S_n/vonk)+0.5_rprec) &
                     & *log(wmro_y*utaup_d/wmro_nu)/vonk
             end if
             wmro_usol_d = (wallstress_ML_cal(wmro_mlin) & 
                   & + wmro_usol_d) * utaup_d
             wmro_judge      = (u_avg - wmro_usol_d)

             print *, 'Bisection method', utaup_f, wmro_usol_f, utaup, wmro_usol, utaup_d, wmro_usol_d
             if (((wmro_usol_d - u_avg) > 0.0_rprec .and. (wmro_usol - u_avg) < 0.0_rprec) .or. ((wmro_usol_d - u_avg) < 0.0_rprec .and. (wmro_usol - u_avg) > 0.0_rprec)) then
                 utaup_f     = utaup
                 wmro_usol_f = wmro_usol
             end if
             utaup           = utaup_d
             wmro_usol       = wmro_usol_d
         end do

         txz(jx,jy,1) = -utaup**2*u1(jx,jy)/u_avg
         tyz(jx,jy,1) = -utaup**2*v1(jx,jy)/u_avg
         wmro_y = wmro_y +0.01_rprec*dz
         wmro_mlin = reshape((/wmro_y*utaup/wmro_nu, &
                     & wmro_y/utaup*wmro_S_n/), (/2,1/))
         if (abs(wmro_S_n) .LE. 1e-6) then
             wmro_usol = log(wmro_y*utaup/wmro_nu)/vonk
         else 
             wmro_usol = (sign(0.5_rprec,wmro_y-utaup/wmro_S_n/vonk)+0.5_rprec) &
                 & *(wmro_S_n*wmro_y/utaup+log(utaup**2/wmro_S_n/wmro_nu)/vonk-(log(vonk)+1.0_rprec)/vonk) &
                 &       + (sign(0.5_rprec,-wmro_y+utaup/wmro_S_n/vonk)+0.5_rprec) &
                 & *log(wmro_y*utaup/wmro_nu)/vonk
         end if
         wmro_dudz = (wallstress_ML_cal(wmro_mlin) & 
               & + wmro_usol) * utaup

         wmro_y = wmro_y - 0.02_rprec*dz
         wmro_mlin = reshape((/wmro_y*utaup/wmro_nu, &
                     & wmro_y/utaup*wmro_S_n/), (/2,1/))
         if (abs(wmro_S_n) .LE. 1e-6) then
             wmro_usol = log(wmro_y*utaup/wmro_nu)/vonk
         else 
             wmro_usol = (sign(0.5_rprec,wmro_y-utaup/wmro_S_n/vonk)+0.5_rprec) &
                 & *(wmro_S_n*wmro_y/utaup+log(utaup**2/wmro_S_n/wmro_nu)/vonk-(log(vonk)+1.0_rprec)/vonk) &
                 &       + (sign(0.5_rprec,-wmro_y+utaup/wmro_S_n/vonk)+0.5_rprec) &
                 & *log(wmro_y*utaup/wmro_nu)/vonk
         end if
         wmro_dudz = ( - (wallstress_ML_cal(wmro_mlin)+wmro_usol) * &
                     & utaup + wmro_dudz)/(0.02_rprec*dz)
         dudz(jx,jy,1) = wmro_dudz*u1(jx,jy)/u_avg
         dvdz(jx,jy,1) = wmro_dudz*v1(jx,jy)/u_avg

         wmro_y = wmro_y + 0.01_rprec*dz

      end do
      end do
  end if

      ! Xinyi - Suggestion : Might want to use some matrices to store several layers as the initial
      ! condition of the next time step

      ! Output for equil_output.dat
      if( mod(jt_total,domain_nskip)==1)then
        digits_jt=int(log(real(jt_total))/log(10.0))+1
        write(fname,'(A,I<digits_jt>,A)') 'equil_output.',jt_total,'.dat' 
        open(375,file=fname)  
        do jy=1,ny 
        do jx=1,nx
          write(375,*) txz(jx,jy,1), tyz(jx,jy,1), u1(jx,jy), v1(jx,jy), wmro_nu*dudz(jx,jy,1), wmro_nu*dvdz(jx,jy,1)
        end do
        end do
        close(375)
      end if
end select
end subroutine wallstress_ML

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Xinyi - WM - rotate : to implement different wall models
!                       for rotating channel
!                       1. Loppi's model - 2018/10
!                       2. Non-harmonic model - 2018/11
!                         (Assuming positive spanwise rotation)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Xinyi : This subroutine is for upper boundary layer, and almost exactly the same as
! the one for the lower boundary layer. Take care with nz and nz-1
subroutine wallstress_upper_wmro ()
use types,only:rprec
use param,only:dz,ld,lh,nx,ny,nz,vonk,ubc_mom,zo, domain_nskip, jt_total, jt
use param,only: nu_molec,z_i,u_star
use sim_param,only:u,v,dudz,dvdz,txz,tyz
use wmrotate_param
use test_filtermodule
implicit none
integer     ::jx,jy
real(rprec),dimension(ld,ny)    :: u1,v1
real(rprec),dimension(wmro_z_ni):: wmro_a, wmro_b, wmro_c ! quadratic parameter
real(rprec),dimension(wmro_z_ni):: wmro_mdudz, wmro_dudz, wmro_dudzf, wmro_coeff_i, wmro_Ri
real(rprec) :: const, wmro_judge, wmro_usol, wmro_coeff, u_avg
real(rprec) :: wmro_sign, wmro_nu
integer     :: digits_jt
CHARACTER*50 :: fname

wmro_nu = nu_molec/(z_i*u_star)

select case (ubc_mom)
  case (0) ! Stress free
    txz(:, :, nz) = 0._rprec
    tyz(:, :, nz) = 0._rprec
    dudz(:, :, nz) = 0._rprec
    dvdz(:, :, nz) = 0._rprec
  case (1) ! Wall
  
  if(modulo(jt,wmro_nstep) == 1) then
      u1=u(:,:,nz-1)
      v1=v(:,:,nz-1)

      !original filtering
      call test_filter ( u1 )
      call test_filter ( v1 )

      do jy=1,ny
      do jx=1,nx

         ! Model No-slip upper boundary layer :
!         txz(jx,jy,nz) = wmro_nu*u1(jx,jy)/(0.5_rprec*dz)
!         tyz(jx,jy,nz) = wmro_nu*v1(jx,jy)/(0.5_rprec*dz)
!         dudz(jx,jy,nz) = -u1(jx,jy)/(0.5_rprec*dz)
!         dvdz(jx,jy,nz) = -v1(jx,jy)/(0.5_rprec*dz)
         
         u_avg = sqrt(u1(jx,jy)**2 + v1(jx,jy)**2)
         const = wmro_nu * (u_avg / (0.5_rprec*dz))
         wmro_judge = 1.0_rprec
        !wmro_mdudz = - u_avg / (0.5_rprec*dz)
         wmro_dudz = u_avg / (0.5_rprec*dz)

         do while (abs(wmro_judge) > 1.0e-6_rprec)

            ! To generate damping parameters, we need to turn to the friction scale
            ! Model Loppi : 
            !wmro_alpha_stable = vonk*wmro_z*sqrt(const)*(1.0_rprec - exp(-wmro_z*sqrt(const)/wmro_nu/A_stable))**2

            ! Model min:
            if (Omega > 1.0e-6_rprec) then
                wmro_alpha = min(vonk*wmro_z,const / (2.0_rprec*abs(Omega)))
            else
                wmro_alpha = vonk*wmro_z
            end if
            wmro_alpha = wmro_nu + wmro_alpha * sqrt(const) * (1.0_rprec - exp(-wmro_z*sqrt(const)/wmro_nu/A_equil))**2
            wmro_alpha = 1.0_rprec / wmro_alpha

            ! For solving quadratic equations
            ! Model Loppi :
            !wmro_a = wmro_nu + wmro_alpha_stable
            !wmro_b = -(-const + wmro_beta*wmro_alpha_stable*wmro_S_n)
            !wmro_c = -wmro_beta*wmro_alpha_stable*wmro_S_n**2
            !wmro_Ri = wmro_S_n/wmro_mdudz
            !wmro_Ri = wmro_Ri*(wmro_Ri+1.0_rprec)
            !where(wmro_Ri*wmro_beta > 1.0_rprec) 
            !   wmro_mdudz = -const/wmro_nu
            !elsewhere(abs(wmro_Ri) < 1.0e-10_rprec)
            !   wmro_mdudz = -const/wmro_a
            !elsewhere
            !   wmro_mdudz = (-wmro_b - &
            !      sqrt(wmro_b**2-4.0_rprec*wmro_a*wmro_c))/ &
            !      (2.0_rprec*wmro_a)
            !end where
            ! Model min, harmonic:
            wmro_dudz = const * wmro_alpha

            ! converge for shooting method
            !wmro_usol = -sum(wmro_zdis/6.0_rprec * &
            !   (wmro_mdudz(1:(size(wmro_mdudz)-2):2) + wmro_mdudz(3::2) + &
            !   wmro_mdudz(2::2)*4.0_rprec))
            wmro_usol = sum(wmro_zdis/6.0_rprec * &
               (wmro_dudz(1:(size(wmro_dudz)-2):2) + wmro_dudz(3::2) + &
               wmro_dudz(2::2)*4.0_rprec))
            !wmro_coeff_i = wmro_mdudz**2/(wmro_mdudz**2*wmro_a - wmro_c)
            wmro_coeff_i = wmro_alpha
            wmro_coeff = sum(wmro_zdis/6.0_rprec * &
               (wmro_coeff_i(1:(size(wmro_coeff_i)-2):2) + wmro_coeff_i(3::2) + &
               wmro_coeff_i(2::2)*4.0_rprec))
            !wmro_judge = (u_avg - wmro_usol)/wmro_coeff
            wmro_judge = (u_avg - wmro_usol)/wmro_coeff
            const = const + wmro_judge
         end do

         txz(jx,jy,nz) = const*u1(jx,jy)/u_avg
         tyz(jx,jy,nz) = const*v1(jx,jy)/u_avg
         !dudz(jx,jy,nz) = wmro_mdudz(wmro_z_ni)*u1(jx,jy)/u_avg
         !dvdz(jx,jy,nz) = wmro_mdudz(wmro_z_ni)*v1(jx,jy)/u_avg
         dudz(jx,jy,nz) = -wmro_dudz(wmro_z_ni)*u1(jx,jy)/u_avg
         dvdz(jx,jy,nz) = -wmro_dudz(wmro_z_ni)*v1(jx,jy)/u_avg
      end do
      end do
  end if 

      ! Xinyi - Suggestion : Might want to use some matrices to store several layers as the initial
      ! condition of the next time step

      ! Output for equil_output
      if( mod(jt_total,domain_nskip)==1)then
        digits_jt=int(log(real(jt_total))/log(10.0))+1
        write(fname,'(A,I<digits_jt>,A)') 'equil_output_upper.',jt_total,'.dat' 
        open(375,file=fname)  
        do jy=1,ny 
        do jx=1,nx
          write(375,*) txz(jx,jy,nz), tyz(jx,jy,nz), u1(jx,jy), v1(jx,jy), wmro_nu*dudz(jx,jy,nz), wmro_nu*dvdz(jx,jy,nz)
        end do
        end do
        close(375)
      endif
end select
end subroutine wallstress_upper_wmro
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Xinyi - WM - rotate : to implement different wall models
!                       1. Loppi's model - 2018/10
!                       2. Non-harmonic model - 2018/11
!                         (Assuming positive spanwise rotation)
! The end of both boundaries
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!===========================the following subroutine & modules are slip wall model and further developements============!
subroutine  swm_wallstress()
  use grid_m, only : grid
  use types,only:rprec
  use param,only:dz,ld,lh,nx,ny,nz,vonk,zo,jt_total,domain_nskip, coord, dx, dy, dt
  use sim_param,only:u,v,w,p,dudz,dvdz,txz,tyz, up_swm, vp_swm
  use sgs_param,only:Nu_t
  use test_filtermodule
  implicit none
  integer::jx,jy
  real(rprec) :: usMean, lp, nu_v, B
  real(rprec),dimension(nx,ny)::u_slip, v_slip, us_WM, uSlp, Neqx, Neqy
  real(rprec),dimension(ld,ny)::u1,v1,dudz1,dvdz1, Nut1, p1
  CHARACTER*50 :: fname
  integer :: swm_file
  integer :: digits_jt
  integer, pointer, dimension(:) :: autowrap_i, autowrap_j !useful array for autowraped index
  real(rprec) :: Cd
  Cd=0.0
  
    u1=u(:,:,1)
    v1=v(:,:,1)   
    !p1 is used for non equilibrium effects---
    p1=p(:,:,1)
    do jx=1,nx
    do jy=1,ny
        p1(jx,jy)=p1(jx,jy)-0.5*(u1(jx,jy)**2+v1(jx,jy)**2+(w(jx,jy,2)*0.25)**2)
    end do
    end do
    !-----
    call test_filter ( u1 )
    call test_filter ( v1 )    
    call test_filter ( p1 ) !p1 is used for non equilibrium effects---
    Nut1=Nu_t(:,:,2)
    usMean=1.0_rprec
    
    lp=dz/2.0*log(dz/2.0/zo)
    do jy=1,ny
    do jx=1,nx
        u_slip(jx,jy)=lp/(lp+dz/2.0)*u1(jx,jy)
        v_slip(jx,jy)=lp/(lp+dz/2.0)*v1(jx,jy)
    end do
    end do
    
    !=====this is to include non-equilibrium effects
    nullify(autowrap_i, autowrap_j)
    autowrap_i => grid % autowrap_i
    autowrap_j => grid % autowrap_j
    do jy=1,ny
    do jx=1,nx
      Neqx(jx,jy)=-( u1(autowrap_i(jx+1),jy)*u1(autowrap_i(jx+1),jy) &
                    -u1(autowrap_i(jx-1),jy)*u1(autowrap_i(jx-1),jy) )/(2.0*dx)*(dz/2.0) &
                  -( u1(jx,autowrap_j(jy+1))*v1(jx,autowrap_j(jy+1)) &
                    -u1(jx,autowrap_j(jy-1))*v1(jx,autowrap_j(jy-1)) )/(2.0*dy)*(dz/2.0) &
                  -( p1(autowrap_i(jx+1),jy) &
                    -p1(autowrap_i(jx-1),jy) )/(2.0*dx)*(dz/2.0)
      Neqy(jx,jy)=-( u1(autowrap_i(jx+1),jy)*v1(autowrap_i(jx+1),jy) &
                    -u1(autowrap_i(jx-1),jy)*v1(autowrap_i(jx-1),jy) )/(2.0*dx)*(dz/2.0) &
                  -( v1(jx,autowrap_j(jy+1))*v1(jx,autowrap_j(jy+1)) &
                    -v1(jx,autowrap_j(jy-1))*v1(jx,autowrap_j(jy-1)) )/(2.0*dy)*(dz/2.0) &
                  -( p1(jx,autowrap_j(jy+1)) &
                    -p1(jx,autowrap_j(jy-1)) )/(2.0*dy)*(dz/2.0)
    end do
    end do
      !this is to include dudt effects
    if(abs(up_swm(1,1))>0.000001)then
      do jy=1,ny
      do jx=1,nx
         Neqx(jx,jy)=Neqx(jx,jy)-(u1(jx,jy)-up_swm(jx,jy))/dt*(dz/2.0)
         Neqy(jx,jy)=Neqy(jx,jy)-(v1(jx,jy)-vp_swm(jx,jy))/dt*(dz/2.0) 
      end do
      end do
    endif
    !-----------------
    
    nu_v=zo*usMean*exp(vonk*B)
    do jy=1,ny
    do jx=1,nx
       txz(jx,jy,1)=-(usMean*vonk*dz/2.0+nu_v)*(u1(jx,jy)-u_slip(jx,jy))/(dz/2.0)
       tyz(jx,jy,1)=-(usMean*vonk*dz/2.0+nu_v)*(v1(jx,jy)-v_slip(jx,jy))/(dz/2.0)
       txz(jx,jy,1)=txz(jx,jy,1)+Neqx(jx,jy)/log(dz/2.0/zo) !for non-equilibrium effects
       tyz(jx,jy,1)=tyz(jx,jy,1)+Neqy(jx,jy)/log(dz/2.0/zo) !for non-equilibrium effects
       dudz(jx,jy,1)=u(jx,jy,1)/(dz/2.0*log(dz/zo/2.0))
       dvdz(jx,jy,1)=v(jx,jy,1)/(dz/2.0*log(dz/zo/2.0))
    end do
    end do
    
    if( mod(jt_total,domain_nskip)==1)then
        digits_jt=int(log(real(jt_total))/log(10.0))+1
        write(fname,'(A,I<digits_jt>,A)') 'swm_output.',jt_total,'.dat' 
        swm_file=125
        open(swm_file,file=fname)             
        do jy=1,ny 
        do jx=1,nx
            write(swm_file,*) u_slip(jx,jy), v_slip(jx,jy), txz(jx,jy,1), tyz(jx,jy,1), u1(jx,jy), v1(jx,jy)
        end do
        end do
        close(swm_file)
    endif
	
    !update the last time step velocity use for non-equilibrium effects
    do jy=1,ny
	do jx=1,nx
		up_swm(jx,jy)=u1(jx,jy) !xiang dirty swm for dudt
		vp_swm(jx,jy)=v1(jx,jy) !xiang dirty swm for dudt
	end do
	end do
    nullify(autowrap_i, autowrap_j) !nullify always after using
end subroutine swm_wallstress

!===========================the following subroutines & modules are integral wall models================================!

subroutine iwm_wallstress
  use param, only : jt, nx,ny, dt
  use iwmles, only : iwm_ntime_skip, iwm_tauwx, iwm_tauwy, iwm_dudzT, iwm_dirx, iwm_diry, iwm_dt
  use sim_param,only:dudz,dvdz,txz,tyz
  
  implicit none
  
  integer :: iwm_i, iwm_j
  
  !xiang: calculate the time step used in the integral wall model 
  !! DO NOT USE iwm_ntime_skip=1 !! !! this number is hard coded to prevent any mis-use...
  if(mod(jt,iwm_ntime_skip)==1)then
    iwm_dt=dt
  else
    iwm_dt=iwm_dt+dt
  end if
  
  !xiang: compute the wall stress
  if(mod(jt,iwm_ntime_skip)==0)then
    call iwm_calc_lhs() !gather flow status, update the integrated unsteady term, convective term, turbulent diffusion term etc.
    call iwm_calc_wallstress() !the subroutine to calculate wall stress 
  endif
  call iwm_monitor() !this is to monitor any quantity from the iwm, useful debugging tool
  
  !xiang: imposing txz, tyz, dudz, dvdz every time step even iwm_* are not computed every time step. 
  do iwm_i=1,nx
  do iwm_j=1,ny
    txz(iwm_i,iwm_j,1)=-iwm_tauwx(iwm_i,iwm_j)   !wall stress, use the value calculated in iwm, note the negative sign
    tyz(iwm_i,iwm_j,1)=-iwm_tauwy(iwm_i,iwm_j)   !wall stress, use the value calculated in iwm, note the negative sign
    
    !xiang: I think those quantities should be consistent with the iwm.
    !       Users could switch back to equilibrium values, which I have tested and is fine.
    dudz(iwm_i,iwm_j,1)=iwm_dudzT(iwm_i,iwm_j,iwm_dirx) !dudz, use the value calculated in iwm, note the positive sign
    dvdz(iwm_i,iwm_j,1)=iwm_dudzT(iwm_i,iwm_j,iwm_diry) !dudz, use the value calculated in iwm, note the positive sign
  end do
  end do
  
end subroutine iwm_wallstress

!xiang: memory allocation and initialize everything with plug flow conditions
subroutine iwm_malloc
  use types,only : rprec
  use param,only : nx, ny, dz, vonk, zo, cfl, L_x
  use iwmles

  implicit none

  real(rprec) :: usinit, uinit, vinit, Dzp
  
  usinit= 1.0_rprec  !initial value for us (the friction velocity)
  uinit = usinit/vonk*log(dz/2.0_rprec/zo) !initial value for the x-velocity at first grid point
  vinit = 0.0_rprec !initial value for the y-velocity at first grid point
  Dzp=dz/2.0_rprec  !at the height of the first grid point.
  
  !us in x, y directions
  allocate(iwm_utx(nx,ny))
  allocate(iwm_uty(nx,ny))
  iwm_utx = usinit
  iwm_uty = 0._rprec

  !wall stress in x, y directions
  allocate(iwm_tauwx(nx,ny))
  allocate(iwm_tauwy(nx,ny))
  iwm_tauwx = usinit**2.0_rprec
  iwm_tauwy = 0._rprec

  !flitered velocity at the first grid point in x, y directions
  allocate(iwm_flt_tagvel  (nx,ny,iwm_DN))
  allocate(iwm_flt_tagvel_m(nx,ny,iwm_DN))
  iwm_flt_tagvel  (:,:,iwm_dirx) = uinit
  iwm_flt_tagvel  (:,:,iwm_diry) = vinit
  iwm_flt_tagvel_m(:,:,iwm_dirx) = uinit
  iwm_flt_tagvel_m(:,:,iwm_diry) = vinit
  
  !pressure at first grid point
  allocate(iwm_flt_p(nx,ny)) 
  iwm_flt_p = 0._rprec
  
  !integrals of Lu, Lv, etc.
  allocate(iwm_inte  (nx,ny,iwm_LN))
  allocate(iwm_inte_m(nx,ny,iwm_LN))
  iwm_inte   (:,:,iwm_Lu) =uinit*Dzp
  iwm_inte   (:,:,iwm_Lv) =0._rprec
  iwm_inte   (:,:,iwm_Luu)=uinit**2.0_rprec*Dzp
  iwm_inte   (:,:,iwm_Lvv)=0._rprec
  iwm_inte   (:,:,iwm_Luv)=0._rprec
   iwm_inte_m(:,:,iwm_Lu) =uinit*Dzp
   iwm_inte_m(:,:,iwm_Lv) =0._rprec
   iwm_inte_m(:,:,iwm_Luu)=uinit**2.0_rprec*Dzp
   iwm_inte_m(:,:,iwm_Lvv)=0._rprec
   iwm_inte_m(:,:,iwm_Luv)=0._rprec

  !each term in the integral equation and top/bottom derivatives
  allocate(iwm_unsdy  (nx,ny,iWM_DN))
  allocate(iwm_conv   (nx,ny,iWM_DN))
  allocate(iwm_PrsGrad(nx,ny,iWM_DN))
  allocate(iwm_diff   (nx,ny,iwm_DN))
  allocate(iwm_LHS    (nx,ny,iWM_DN))
  allocate(iwm_dudzT  (nx,ny,iwm_DN))
  allocate(iwm_dudzB  (nx,ny,iwm_DN))
  iWM_unsdy   = 0._rprec
  iWM_conv    = 0._rprec
  iWM_PrsGrad = 0._rprec
  iwm_diff    = 0._rprec
  iWM_LHS     = -uinit*Dzp
  iwm_dudzT(:,:,iwm_dirx) = usinit/vonk/Dzp
  iwm_dudzT(:,:,iwm_diry) = 0._rprec
  iwm_dudzB(:,:,iwm_dirx) = usinit/vonk/zo
  iwm_dudzB(:,:,iwm_diry) = 0._rprec
  
  !filtered friction velocity and the filtering time scale, tR<1
  allocate(iwm_flt_us(nx,ny))
  allocate(iwm_tR    (nx,ny))
  iwm_flt_us = usinit
  iwm_tR = (cfl*L_x/nx/uinit)/(dz/2._rprec/vonk/usinit)  

  !cell height and imposed roughness length
  allocate(iwm_Dz(nx,ny))
  allocate(iwm_z0(nx,ny))
  iWM_Dz = dz/2._rprec
  iWM_z0 = zo !we leave the possibility of zo as a function of x-y
  
  !linear correction to the log profile
  allocate(iwm_Ax(nx,ny))
  allocate(iwm_Ay(nx,ny))
  iwm_Ax=0._rprec
  iwm_Ay=0._rprec
  
  !time step seen by the iwm
  iwm_dt=iwm_ntime_skip*cfl*L_x/nx/uinit
  
end subroutine iwm_malloc

!xiang: this subroutine deallocate memory used for iwm
subroutine iwm_demalloc

  use iwmles

  implicit none

  deallocate(iwm_utx)
  deallocate(iwm_uty)

  deallocate(iwm_tauwx)
  deallocate(iwm_tauwy)

  deallocate(iwm_flt_tagvel  )
  deallocate(iwm_flt_tagvel_m)

  deallocate(iwm_inte  )
  deallocate(iwm_inte_m)

  deallocate(iwm_unsdy  )
  deallocate(iwm_conv   )
  deallocate(iwm_PrsGrad)
  deallocate(iwm_diff   )
  deallocate(iwm_LHS    )
  deallocate(iwm_dudzT  )
  deallocate(iwm_dudzB  )

  deallocate(iwm_flt_us)
  deallocate(iwm_tR    )

  deallocate(iwm_Dz)
  deallocate(iwm_z0)
  deallocate(iwm_Ax)
  deallocate(iwm_Ay)
  
end subroutine iwm_demalloc

!xiang: calculate the left hand side of the iwm system.
subroutine iwm_calc_lhs()
  use iwmles
  use grid_m, only : grid
  use types,only : rprec
  use param,only : nx,ny,dx,dy,ld
  use sim_param,only : u,v,w,p
  use test_filtermodule
  
  implicit none
  
  integer, pointer, dimension(:) :: autowrap_i, autowrap_j !useful array for autowraped index
  integer :: iwm_i,iwm_j 
  real(kind=rprec), dimension(ld,ny) :: u_inst, v_inst, w_inst, p_inst !the instantaneous field
  real(kind=rprec) :: p_bar !the mean pressure at first grid point 
  real(kind=rprec) :: Luux, Luvx, Luvy, Lvvy, Lux, Lvy !for temporary storage of derivativex of integrals like dLudx, dLvdx...
  real(kind=rprec) :: phip, phim 
  

  nullify(autowrap_i, autowrap_j)
  autowrap_i => grid % autowrap_i
  autowrap_j => grid % autowrap_j
 
  !updat the u, v for previous time step
  do iwm_i=1,nx
  do iwm_j=1,ny
    iwm_flt_tagvel_m(iwm_i,iwm_j,iwm_dirx) = iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)
    iwm_flt_tagvel_m(iwm_i,iwm_j,iwm_diry) = iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)
  end do
  end do

  !get the instantaneous field
  u_inst=u(:,:,1)
  v_inst=v(:,:,1) !let us do not worry about the half cell displacement in v
  w_inst=w(:,:,2)*0.25_rprec !w is quadrantic near the wall
  p_inst=p(:,:,1)
  do iwm_i=1,nx
  do iwm_j=1,ny
    p_inst(iwm_i,iwm_j)= p_inst(iwm_i,iwm_j) &
                        -0.5_rprec*( (u_inst(iwm_i,iwm_j))**2.0_rprec &
                                    +(v_inst(iwm_i,iwm_j))**2.0_rprec &
                                    +(w_inst(iwm_i,iwm_j))**2.0_rprec   )  ! the real pressure is needed, this step is CODE SPECIFIC!
  end do
  end do
  !obtain the pressure fluctuations
  
  p_bar=0._rprec
  do iwm_i=1,nx
  do iwm_j=1,ny
    p_bar=p_bar+p_inst(iwm_i,iwm_j)
  end do
  end do
  p_bar=p_bar/nx/ny
  do iwm_i=1,nx
  do iwm_j=1,ny
    p_inst(iwm_i,iwm_j)=p_inst(iwm_i,iwm_j)-p_bar
  end do
  end do
   
  !all the data enters must be filtered (see Anderson&Meneveau 2011 JFM)
  call test_filter ( u_inst )
  call test_filter ( v_inst )
  call test_filter ( w_inst )
  call test_filter ( p_inst )

  !temporal filtering 
  do iwm_i=1,nx
  do iwm_j=1,ny
    iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx) = iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)*(1._rprec-iwm_tR(iwm_i,iwm_j)) &
                                          +u_inst        (iwm_i,iwm_j)         *iwm_tR(iwm_i,iwm_j)
    iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry) = iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)*(1._rprec-iwm_tR(iwm_i,iwm_j)) &
                                          +v_inst        (iwm_i,iwm_j)         *iwm_tR(iwm_i,iwm_j)
    iwm_flt_p     (iwm_i,iwm_j)          = iwm_flt_p     (iwm_i,iwm_j)         *(1._rprec-iwm_tR(iwm_i,iwm_j)) &
                                          +p_inst        (iwm_i,iwm_j)         *iwm_tR(iwm_i,iwm_j)
  end do
  end do

  !Calculate LHS, calculation of the integrals is done from the last time step in the subroutine iwm_calc_wallstress, so is iwm_diff
  do iwm_i=1,nx
  do iwm_j=1,ny
    !the unsteady term
    iwm_unsdy(iwm_i,iwm_j,iwm_dirx)=(iwm_inte(iwm_i,iwm_j,iwm_Lu)-iwm_inte_m(iwm_i,iwm_j,iwm_Lu))/iwm_dt
    iwm_unsdy(iwm_i,iwm_j,iwm_diry)=(iwm_inte(iwm_i,iwm_j,iwm_Lv)-iwm_inte_m(iwm_i,iwm_j,iwm_Lv))/iwm_dt  
    !the convective term
    phip=iwm_inte(autowrap_i(iwm_i+1),iwm_j,iwm_Luu)
    phim=iwm_inte(autowrap_i(iwm_i-1),iwm_j,iwm_Luu)
    Luux=(phip-phim)/dx/2._rprec
    phip=iwm_inte(iwm_i,autowrap_j(iwm_j+1),iwm_Luv)
    phim=iwm_inte(iwm_i,autowrap_j(iwm_j-1),iwm_Luv)
    Luvy=(phip-phim)/dy/2._rprec
    phip=iwm_inte(autowrap_i(iwm_i+1),iwm_j,iwm_Luv)
    phim=iwm_inte(autowrap_i(iwm_i-1),iwm_j,iwm_Luv)
    Luvx=(phip-phim)/dx/2._rprec
    phip=iwm_inte(iwm_i,autowrap_j(iwm_j+1),iwm_Lvv)
    phim=iwm_inte(iwm_i,autowrap_j(iwm_j-1),iwm_Lvv)
    Lvvy=(phip-phim)/dy/2._rprec
    phip=iwm_inte(autowrap_i(iwm_i+1),iwm_j,iwm_Lu )
    phim=iwm_inte(autowrap_i(iwm_i-1),iwm_j,iwm_Lu )
    Lux =(phip-phim)/dx/2._rprec
    phip=iwm_inte(iwm_i,autowrap_j(iwm_j+1),iwm_Lv )
    phim=iwm_inte(iwm_i,autowrap_j(iwm_j-1),iwm_Lv )
    Lvy=(phip-phim)/dy/2._rprec
    iwm_conv(iwm_i,iwm_j,iwm_dirx)=Luux+Luvy-iwm_flt_tagvel_m(iwm_i,iwm_j,iWM_dirx)*(Lux+Lvy)
    iwm_conv(iwm_i,iwm_j,iwm_diry)=Luvx+Lvvy-iwm_flt_tagvel_m(iwm_i,iwm_j,iWM_diry)*(Lux+Lvy)    
    !the pressure gradient term
    phip=iwm_flt_p(autowrap_i(iwm_i+1),iwm_j)
    phim=iwm_flt_p(autowrap_i(iwm_i-1),iwm_j)
    iwm_PrsGrad(iwm_i,iwm_j,iwm_dirx)=(phip-phim)/dx/2._rprec*iwm_Dz(iwm_i,iwm_j)-1.0_rprec*iwm_Dz(iwm_i,iwm_j) !including the mean unit pressure gradient
    phip=iwm_flt_p(iwm_i,autowrap_j(iwm_j+1))
    phim=iwm_flt_p(iwm_i,autowrap_j(iwm_j-1))
    iwm_PrsGrad(iwm_i,iwm_j,iwm_diry)=(phip-phim)/dy/2._rprec*iwm_Dz(iwm_i,iwm_j)     
    !the left hand side.
    iwm_lhs(iwm_i,iwm_j,iwm_dirx)=  -iwm_inte(iwm_i,iwm_j,iwm_Lu) &
                                  +iwm_dt*( iwm_conv(iwm_i,iwm_j,iwm_dirx) &
                                           +iwm_PrsGrad(iwm_i,iwm_j,iwm_dirx) &
                                           -iwm_diff(iwm_i,iwm_j,iwm_dirx)     )  ! this is the integrated momentum equation, except for the Lu term
    iwm_lhs(iwm_i,iwm_j,iwm_diry)=  -iwm_inte(iwm_i,iwm_j,iwm_Lv) &
                                  +iwm_dt*( iwm_conv(iwm_i,iwm_j,iwm_diry) &
                                           +iwm_PrsGrad(iwm_i,iwm_j,iwm_diry) &
                                           -iwm_diff(iwm_i,iwm_j,iwm_diry)     ) ! this is the integrated momentum equation, except for the Lv term
  end do
  end do

  nullify(autowrap_i, autowrap_j)
end subroutine iwm_calc_lhs  

subroutine iwm_slv(lhsx,lhsy,Ux,Uy,Dz,z0,utx,uty,fx,fy)
  use types,only:rprec
  use param,only:vonk
  implicit none
  real(kind=rprec), intent(in)  :: lhsx,lhsy,Ux,Uy,Dz,z0,utx,uty
  real(kind=rprec), intent(out) :: fx,fy
  real(kind=rprec) :: Ax, Ay, Vel, inteLu, inteLv

  Ax=(Ux-utx/vonk*log(Dz/z0))/((1.0-z0/Dz))
  Ay=(Uy-uty/vonk*log(Dz/z0))/((1.0-z0/Dz))
  Vel=sqrt(Ux**2.0+Uy**2.0)
  inteLu = 1.0/2.0*Dz*Ax*(1.0-z0/Dz)**2.0 &
          +1.0/vonk*utx*Dz*(z0/Dz-1.0+log(Dz/z0))
  inteLv = 1.0/2.0*Dz*Ay*(1.0-z0/Dz)**2.0 &
          +1.0/vonk*uty*Dz*(z0/Dz-1.0+log(Dz/z0))
  fx=inteLu+lhsx
  fy=inteLv+lhsy
end subroutine iwm_slv

subroutine iwm_calc_wallstress
  use iwmles
  use types,only:rprec
  use param,only:vonk,nx,ny,coord,ld
  use sim_param,only:u,v
  use test_filtermodule
  
  implicit none

  integer :: iwm_i, iwm_j
  real(kind=rprec) :: fx,fy,fxp,fyp
  real(kind=rprec) :: iwm_tol, iwm_eps
  real(kind=rprec) :: a11, a12, a21, a22
  real(kind=rprec) :: iwmutxP,iwmutyP
  integer          :: iter, MaxIter, equil_flag, div_flag
  real(kind=rprec) :: equilWMpara,equilutx,equiluty
  CHARACTER*50     :: fname
  real(kind=rprec) :: iwmpAx, iwmpAy, iwmputx,iwmputy,iwmpz0,iwmpDz
  real(kind=rprec) :: utaup
  real(kind=rprec) :: dVelzT, dVelzB, Vel
  real(kind=rprec) :: Cd
  Cd=0.0
  
  MaxIter=1500

  iwm_tol = 0.000001_rprec
  iwm_eps = 0.000000001_rprec
  
  Do iwm_i=1,nx
  Do iwm_j=1,ny

    !use Newton method to solve the system
    iwm_utx(iwm_i,iwm_j)=1.0_rprec*sign(1.0_rprec,iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx))
    iwm_uty(iwm_i,iwm_j)=0.1_rprec*sign(1.0_rprec,iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry))

    Call iwm_slv(iwm_lhs       (iwm_i,iwm_j,iwm_dirx),iwm_lhs       (iwm_i,iwm_j,iwm_diry), &
                 iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx),iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry), &
                 iwm_Dz        (iwm_i,iwm_j),         iwm_z0        (iwm_i,iwm_j),          &
                 iwm_utx       (iwm_i,iwm_j),         iwm_uty       (iwm_i,iwm_j),          &
                 fx,fy                                                                     )
    
    iter=0
    equil_flag=0
    div_flag=0
    do while (max(abs(fx),abs(fy))>iwm_tol)      
      iwmutxP=iwm_utx(iwm_i,iwm_j)+iWM_eps
      iwmutyP=iwm_uty(iwm_i,iwm_j)
      Call iwm_slv(iwm_lhs       (iwm_i,iwm_j,iwm_dirx),iwm_lhs       (iwm_i,iwm_j,iwm_diry), &
                   iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx),iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry), &
                   iwm_Dz        (iwm_i,iwm_j),         iwm_z0        (iwm_i,iwm_j),          &
                   iwmutxP,                             iwmutyP,                              &
                   fxp,fyp                                                                  )
      a11=(fxp-fx)/iwm_eps
      a21=(fyp-fy)/iwm_eps
      iwmutxP=iwm_utx(iwm_i,iwm_j)
      iwmutyP=iwm_uty(iwm_i,iwm_j)+iwm_eps
      Call iwm_slv(iwm_lhs       (iwm_i,iwm_j,iwm_dirx),iwm_lhs       (iwm_i,iwm_j,iwm_diry), &
                   iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx),iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry), &
                   iwm_Dz        (iwm_i,iwm_j),         iwm_z0        (iwm_i,iwm_j),          &
                   iwmutxP,                             iwmutyP,                              &
                   fxp,fyp                                                                  )
      a12=(fxp-fx)/iwm_eps
      a22=(fyp-fy)/iwm_eps
      iwm_utx(iwm_i,iwm_j)=iwm_utx(iwm_i,iwm_j)-0.50*( a22*fx-a12*fy)/(a11*a22-a12*a21)
      iwm_uty(iwm_i,iwm_j)=iwm_uty(iwm_i,iwm_j)-0.50*(-a21*fx+a11*fy)/(a11*a22-a12*a21)
      Call iwm_slv(iwm_lhs       (iwm_i,iwm_j,iwm_dirx),iwm_lhs       (iwm_i,iwm_j,iwm_diry), &
                   iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx),iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry), &
                   iwm_Dz        (iwm_i,iwm_j),         iwm_z0        (iwm_i,iwm_j),          &
                   iwm_utx       (iwm_i,iwm_j),         iwm_uty       (iwm_i,iwm_j),          &
                   fx,fy                                                                     )
      iter=iter+1
      if(iter>MaxIter)then !maximum iteration reached, according to me, this check is never invoked
        equil_flag=1
        div_flag=1;
        exit
      endif
    end do  
    
    !infinity check  !!this check is never invoked, but it is useful to keep!!
    if(iwm_utx(iwm_i,iwm_j)-1.0==iwm_utx(iwm_i,iwm_j) .or. iwm_uty(iwm_i,iwm_j)-1.0==iwm_uty(iwm_i,iwm_j))then  
        equil_flag=1
        div_flag  =1
    endif
    
    !calculate equilibrium us for equil_flag=1 use
    equilutx=vonk*iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)/log(iwm_Dz(iwm_i,iwm_j)/iwm_z0(iwm_i,iwm_j))  
    equiluty=vonk*iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)/log(iwm_Dz(iwm_i,iwm_j)/iwm_z0(iwm_i,iwm_j))    
    if(equil_flag==1)then
      iwm_utx(iwm_i,iwm_j)=equilutx
      iwm_uty(iwm_i,iwm_j)=equiluty
    endif
    
    !calculate Ax, Ay
    if(equil_flag==1)then
      iwmpAx=0.0_rprec
      iwmpAy=0.0_rprec
    else
      iwmpAx= ( iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)  &
               -iwm_utx(iwm_i,iwm_j)/vonk*log(iwm_Dz(iwm_i,iwm_j)/iwm_z0(iwm_i,iwm_j)))   &
             /((1.0_rprec-iwm_z0(iwm_i,iwm_j)/iwm_Dz(iwm_i,iwm_j)))                              !eq. D2 in Yang et al 2015
      iwmpAy= ( iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)  &
               -iwm_uty(iwm_i,iwm_j)/vonk*log(iwm_Dz(iwm_i,iwm_j)/iwm_z0(iwm_i,iwm_j)))   &
             /((1.0_rprec-iwm_z0(iwm_i,iwm_j)/iwm_Dz(iwm_i,iwm_j)))                              !eq. D2 in Yang et al 2015
    endif
    
    !check for excessive linear term correction !!after first 100 step this check is rarely invoked!!
    if(abs(iwmpAx)>1.0_rprec .or. abs(iwmpAy)>1.0_rprec)then
        equil_flag=1
        iwm_utx(iwm_i,iwm_j)=equilutx
        iwm_uty(iwm_i,iwm_j)=equiluty
        iwmpAx=0.0_rprec
        iwmpAy=0.0_rprec
    endif
    
    !store the linear correction
    iwm_Ax(iwm_i,iwm_j)=iwmpAx
    iwm_Ay(iwm_i,iwm_j)=iwmpAy
    
    !update integral for last time step
    iwm_inte_m(iwm_i,iwm_j,iwm_Lu ) = iwm_inte(iwm_i,iwm_j,iwm_Lu )
    iwm_inte_m(iwm_i,iwm_j,iwm_Lv ) = iwm_inte(iwm_i,iwm_j,iwm_Lv )
    iwm_inte_m(iwm_i,iwm_j,iwm_Luv) = iwm_inte(iwm_i,iwm_j,iwm_Luv)
    iwm_inte_m(iwm_i,iwm_j,iwm_Luu) = iwm_inte(iwm_i,iwm_j,iwm_Luu)
    iwm_inte_m(iwm_i,iwm_j,iwm_Lvv) = iwm_inte(iwm_i,iwm_j,iwm_Lvv)
    
    !those temporary variables are used for convenient reference
    iwmputx = iwm_utx(iwm_i,iwm_j)
    iwmputy = iwm_uty(iwm_i,iwm_j)
    iwmpDz  = iwm_Dz (iwm_i,iwm_j)
    iwmpz0  = iwm_z0 (iwm_i,iwm_j)
    
    !calculate the needed integrals
    iwm_inte(iwm_i,iwm_j,iwm_Lu) = 1.0/2.0*iwmpDz*iwmpAx*(1.0-iwmpz0/iwmpDz)**2.0 &
           +1.0/vonk*iwmputx*iwmpDz*(iwmpz0/iwmpDz-1.0+log(iwmpDz/iwmpz0))     ! Eq. D7 in Yang et al 2015
    iwm_inte(iwm_i,iwm_j,iwm_Lv) = 1.0/2.0*iwmpDz*iwmpAy*(1.0-iwmpz0/iwmpDz)**2.0 &
           +1.0/vonk*iwmputy*iwmpDz*(iwmpz0/iwmpDz-1.0+log(iwmpDz/iwmpz0))     ! Eq. D7 in Yang et al 2015
    iwm_inte(iwm_i,iwm_j,iwm_Luv) = 1.0/vonk**2.0*iwmputx*iwmputy*iwmpDz*(1.0-2*iwmpz0/iwmpDz+(1.0-log(iwmpDz/iwmpz0))**2.0) &
           +1.0/3.0*iwmpAx*iwmpAy*iwmpDz*(1.0 - iwmpz0/iwmpDz)**3.0 &
           -1.0/4.0/vonk*(iwmpAx*iwmputy+iwmpAy*iwmputx)*iwmpDz &
            *(1.0-4.0*iwmpz0/iwmpDz+3.0*iwmpz0**2.0/iwmpDz**2.0-2.0*log(iwmpDz/iwmpz0)+4.0*iwmpz0/iwmpDz*log(iwmpDz/iwmpz0))   ! Eq. D8 in Yang et al 2015
    iwm_inte(iwm_i,iwm_j,iwm_Luu) = 1.0/vonk**2.0*iwmputx**2.0*iwmpDz*((log(iwmpDz/iwmpz0)-1.0)**2.0-2.0*iwmpz0/iwmpDz+1.0) &
           +1.0/3.0*iwmpAx**2.0*iwmpDz*(1.0-iwmpz0/iwmpDz)**3.0 &
           -1.0/2.0/vonk*iwmputx*iwmpAx*iwmpDz &
            *(1.0-4.0*iwmpz0/iwmpDz+3.0*iwmpz0**2.0/iwmpDz**2.0-2.0*log(iwmpDz/iwmpz0)+4.0*iwmpz0/iwmpDz*log(iwmpDz/iwmpz0))   ! Eq. D8 in Yang et al 2015           
    iwm_inte(iwm_i,iwm_j,iwm_Lvv) = 1.0/vonk**2.0*iwmputy**2.0*iwmpDz*((log(iwmpDz/iwmpz0)-1.0)**2.0-2.0*iwmpz0/iwmpDz+1.0) &
           +1.0/3.0*iwmpAy**2.0*iwmpDz*(1.0-iwmpz0/iwmpDz)**3.0 &
           -1.0/2.0/vonk*iwmputy*iwmpAy*iwmpDz&
           *(1.0-4.0*iwmpz0/iwmpDz-3.0*iwmpz0**2.0/iwmpDz**2.0-2.0*log(iwmpDz/iwmpz0)+4.0*iwmpz0/iwmpDz*log(iwmpDz/iwmpz0))    ! Eq. D9 in Yang et al 2015
    
    !calculate top and bottom derivatives
    iwm_dudzT(iwm_i,iwm_j,iwm_dirx)=1.0/iwmpDz*(iwmpAx+iwmputx/vonk)  ! Eq. D5 (a)
    iwm_dudzT(iwm_i,iwm_j,iwm_diry)=1.0/iwmpDz*(iwmpAy+iwmputy/vonk)
    iwm_dudzB(iwm_i,iwm_j,iwm_dirx)=1.0/iwmpDz*iwmpAx+iwmputx/vonk/iwmpz0  ! Eq. D5 (b)
    iwm_dudzB(iwm_i,iwm_j,iwm_diry)=1.0/iwmpDz*iwmpAy+iwmputy/vonk/iwmpz0
    
    !calculte the turbulent diffusion term
    Vel=sqrt(iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)**2.0_rprec+iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)**2.0_rprec) !total velocity
    dVelzT=abs(iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)/Vel*iwm_dudzT(iwm_i,iwm_j,iwm_dirx) &
              +iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)/Vel*iwm_dudzT(iwm_i,iwm_j,iwm_diry))    ! Eq. D6
    dvelzB=sqrt(iwm_dudzB(iwm_i,iwm_j,iwm_dirx)**2.0_rprec+iwm_dudzB(iwm_i,iwm_j,iwm_diry)**2.0_rprec)  !Eq. D6
    iwm_diff(iwm_i,iwm_j,iwm_dirx)= (vonk*iwmpDz)**2.0_rprec*dVelzT*iwm_dudzT(iwm_i,iwm_j,iwm_dirx) &
                                   -(vonk*iwmpz0)**2.0_rprec*dVelzB*iwm_dudzB(iwm_i,iwm_j,iwm_dirx)    !Eq. D4, the eddy viscosity is nu_T=(vonk*y)^2*dudy, hence the formule
    iwm_diff(iwm_i,iwm_j,iwm_diry)= (vonk*iwmpDz)**2.0_rprec*dVelzT*iwm_dudzT(iwm_i,iwm_j,iwm_diry) &
                                   -(vonk*iwmpz0)**2.0_rprec*dVelzB*iwm_dudzB(iwm_i,iwm_j,iwm_diry)    !Eq. D4
  
    !calculate the wall stress
    if(equil_flag==1)then
      equilWMpara= sqrt(iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)**2.0_rprec+iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)**2.0_rprec) &
                  *vonk**2.0_rprec/(log(iwm_Dz(iwm_i,iwm_j)/iwm_z0(iwm_i,iwm_j)))**2.0_rprec
      iwm_tauwx(iwm_i,iwm_j)=equilWMpara*iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)+Cd*iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)**2.0
      iwm_tauwy(iwm_i,iwm_j)=equilWMpara*iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)+Cd*iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)**2.0
    else
      iwm_tauwx(iwm_i,iwm_j)=(vonk*iwmpz0)**2.0_rprec*dVelzB*iwm_dudzB(iwm_i,iwm_j,iwm_dirx)+Cd*iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)**2.0  !Eq. D4
      iwm_tauwy(iwm_i,iwm_j)=(vonk*iwmpz0)**2.0_rprec*dVelzB*iwm_dudzB(iwm_i,iwm_j,iwm_diry)+Cd*iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)**2.0
    endif
    
    !calculate the friciton velocity 
    utaup=(iwm_tauwx(iwm_i,iwm_j)**2.0_rprec +iwm_tauwy(iwm_i,iwm_j)**2.0_rprec )**0.25_rprec  !definition of friction velocity
    iwm_flt_us(iwm_i,iwm_i)=iwm_flt_us(iwm_i,iwm_j)*(1.0_rprec -iwm_tR(iwm_i,iwm_j))+utaup*iwm_tR(iwm_i,iwm_j) !the filtered friction velocity used for filtering time scale

    !update the filtering time scale
    iwm_tR(iwm_i,iwm_j)=iwm_dt/(iwm_Dz(iwm_i,iwm_j)/iwm_flt_us(iwm_i,iwm_j)/vonk)  !Eq. 26
    !filtering time scale can only be larger than the time step, if not, then just use the instantaneous flow field to do the model
    if(iwm_tR(iwm_i,iwm_j)>1.0_rprec)then
      iwm_tR(iwm_i,iwm_j)=1.0_rprec
    endif
    
  end do
  end do
  
end subroutine iwm_calc_wallstress

!this subroutine is to monitor the parameters at one point, do not call this subroutine if you are not interested in how the model works
subroutine iwm_monitor
  use iwmles
  use param,only:coord,nx,ny, domain_nskip
  use io,only:jt_total
  
  implicit none
  
  integer :: iwm_i,iwm_j,dmpPrd, digits_jt
  CHARACTER*50 :: fname
  
  dmpPrd=domain_nskip
  if( mod(jt_total,dmpPrd)==1)then
    digits_jt=int(log(real(jt_total))/log(10.0))+1
    write(fname,'(A,I<digits_jt>,A)') 'iwm_output.',jt_total,'.dat' 
    open(iwm_debug+coord,file=fname)  
    do iwm_i=1,nx 
    do iwm_j=1,ny
      write(iwm_debug+coord,*) iwm_tauwx(iwm_i,iwm_j), iwm_tauwy(iwm_i,iwm_j), iwm_Ax(iwm_i,iwm_j), iwm_Ay(iwm_i,iwm_j), iwm_flt_tagvel(iwm_i,iwm_j,:)
    end do
    end do
    close(iwm_debug+coord)
  endif
  
end subroutine iwm_monitor

!put down a check point for the integral wall model, this subroutine is called after maksing sure iwm_on=1
subroutine iwm_checkPoint()
    use iwmles !all need to be used for next time step, al all need to be recorded
    implicit none
    open(iwm_status,file='iwm_checkPoint.dat', &
                    form='unformatted', &
                    status='unknown', &
                    position='rewind')
    write(iwm_status) iwm_utx(:,:), iwm_uty(:,:), iwm_tauwx(:,:), iwm_tauwy(:,:), &
                      iwm_flt_tagvel(:,:,1:iwm_DN), iwm_flt_tagvel_m(:,:,1:iwm_DN), &
                      iwm_flt_p(:,:), &
                      iwm_inte(:,:,1:iwm_LN), iWM_inte_m(:,:,1:iwm_LN),    &
                      iwm_unsdy(:,:,1:iwm_DN), iwm_conv(:,:,1:iwm_DN), iwm_PrsGrad(:,:,1:iwm_DN), &
                      iwm_diff(:,:,1:iwm_DN), iwm_LHS(:,:,1:iwm_DN), &
                      iwm_dudzT(:,:,1:iwm_DN), iwm_dudzB(:,:,1:iwm_DN), &
                      iwm_flt_us(:,:), iwm_tR(:,:), &
                      iwm_Dz(:,:), iwm_z0(:,:), iwm_Ax(:,:), iwm_Ay(:,:), &
                      iwm_dt
    close(iwm_status)
end subroutine iwm_checkPoint

!read a check point for the integral wall model, this subroutine is called after maksing sure iwm_on=1
subroutine iwm_read_checkPoint()
    use iwmles !all need to be used for next time step, al all need to be read
    implicit none
    open(iwm_status+1,file='iwm_checkPoint.dat', &
                    form='unformatted')
    read(iwm_status+1)iwm_utx(:,:), iwm_uty(:,:), iwm_tauwx(:,:), iwm_tauwy(:,:), &
                      iwm_flt_tagvel(:,:,1:iwm_DN), iwm_flt_tagvel_m(:,:,1:iwm_DN), &
                      iwm_flt_p(:,:), &
                      iwm_inte(:,:,1:iwm_LN), iWM_inte_m(:,:,1:iwm_LN),    &
                      iwm_unsdy(:,:,1:iwm_DN), iwm_conv(:,:,1:iwm_DN), iwm_PrsGrad(:,:,1:iwm_DN), &
                      iwm_diff(:,:,1:iwm_DN), iwm_LHS(:,:,1:iwm_DN), &
                      iwm_dudzT(:,:,1:iwm_DN), iwm_dudzB(:,:,1:iwm_DN), &
                      iwm_flt_us(:,:), iwm_tR(:,:), &
                      iwm_Dz(:,:), iwm_z0(:,:), iwm_Ax(:,:), iwm_Ay(:,:), &
                      iwm_dt
    close(iwm_status+1)
end subroutine iwm_read_checkPoint
