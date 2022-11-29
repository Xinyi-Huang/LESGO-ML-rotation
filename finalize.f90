!!
!!  Copyright (C) 2012-2013  Johns Hopkins University
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

!*******************************************************************************
subroutine finalize()
!*******************************************************************************
! 
! This subroutine is called by the main program. It is a driver subroutine for
! calling all the finalize routines of the various lesgo modules.
!
use param, only : coord
use iwmles, only : iwm_on !xiang for iwm finalize
use param, only : lbc_mom !xiang: always ensure iwm_on=1 only when lbc_mom=1
use wmrotate_param, only : wmrotate_finalize ! Xinyi - WM - rotate : finalization
#ifdef PPMPI
use param, only : MPI_COMM_WORLD, ierr
#endif
#ifdef PPTURBINES
use turbines, only : turbines_finalize
#endif
#ifdef PPATM
use atm_lesgo_interface, only : atm_lesgo_finalize
#endif

implicit none

! Turbines:
#ifdef PPTURBINES
call turbines_finalize ()   ! must come before MPI finalize
#endif   

!finalize for integral wall model xiang
if(lbc_mom == 1)then
if(iwm_on ==  1)then
	if(coord==0) call iwm_demalloc()
endif
endif 

! Xinyi - WM - rotate : Finalization for wall model of rotating channel
call wmrotate_finalize()

! Actuator Turbine Model:
#ifdef PPATM
call atm_lesgo_finalize ()   ! write the restart files
#endif   

! SGS variable histograms
!if (sgs_hist_calc) then
!  call sgs_hist_finalize()
!endif

! MPI:
#ifdef PPMPI
! First make sure everyone in has finished
call mpi_barrier( MPI_COMM_WORLD, ierr )
call mpi_finalize (ierr)
#endif

return
end subroutine finalize
