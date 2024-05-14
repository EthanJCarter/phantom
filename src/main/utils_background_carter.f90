!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module background_carter
!
! This module contains utility routines for disc setup when using
! radiative transfer approximation (ieos=21,icooling=8) but NOT
! using stellar heating (via L_star) to set the background.
!
! :References:
!   Young et. al (2024), Submitted.
!
! :Owner: Ethan Carter
!
! :Runtime parameters:
!   - G           : *in code units*
!   - M_disc      : *disc mass*
!   - M_star      : *mass of central star*
!   - Qmin        : *minimum Toomre Q parameter*
!   - R_c         : *characteristic radius of the exponential taper*
!   - R_in        : *inner disc boundary*
!   - R_out       : *outer disc boundary*
!   - R_ref       : *reference radius*
!   - R_warp      : *position of warp*
!   - T_in        : *temperature (K) at R=R_in*
!   - T_out       : *temperature (K) at R=R_out*
!   - T_ref       : *temperature (K) at R=R_ref*
!   - alphaSS_max : *maximum Shakura-Sunyaev alpha viscosity in disc*
!   - alphaSS_min : *minimum Shakura-Sunyaev alpha viscosity in disc*
!   - c           : *in code units*
!   - cs0         : *sound speed at R=1*
!   - n           : *number of particles in the disc*
!   - p_index     : *power law index of surface density profile*
!   - psi_max     : *maximum warp amplitude*
!   - q_index     : *power law index of sound speed profile*
!   - sig_in      : *surface density (g/cm^2) at R=R_in*
!   - sig_max     : *maximum surface density (g/cm^2)*
!   - sig_out     : *surface density (g/cm^2) at R=R_out*
!   - sig_ref     : *surface density (g/cm^2) at R=R_ref*
!   - udist       : *distance units (cgs)*
!   - umass       : *mass units (cgs)*
!   - utime       : *time units (cgs)*
!
! :Dependencies: centreofmass, dim, eos, externalforces, infile_utils, io,
!   mpidomain, mpiutils, options, part, physcon, random, units, vectorutils
!
 use part,    only:igas,labeltype
 use units,   only:umass,udist,utime
 use physcon, only:gg,kb_on_mh,kboltz,solarm

 implicit none

 public calc_csr, calc_Tr, read_setupfile, get_Tr
 real, public :: R_in,R_ref,R_out,p_index,q_index,M_star,H_R

 private

 contains

 !----------------------------------------------------------------
 !
 ! This is a wrapper subroutine for the module.
 !
 !----------------------------------------------------------------
 subroutine get_Tr(i,Tr)
   use infile_utils, only: inopts
   integer, intent(in) :: i
   real, intent(out) :: Tr
   character(len=100)  :: filename
   real              :: csr
   type(inopts), allocatable :: db(:)
   integer           :: size_db, ierr

   filename = 'disc.setup'

   !defaults
   Tr = 0
   H_R = 0.0556
   M_star = 0.8
   q_index = 0.25
   R_ref = 10

   call read_setupfile(filename,ierr)

   call calc_csr(i, csr)
   call calc_Tr(csr, Tr)

end subroutine get_Tr

 !----------------------------------------------------------------
 !
 ! This is a subroutine for calculating csr, the sound speed for
 ! particle i at radial distance ri from the sink. This is set by
 ! H/r, M_star, and q_index.
 !
 !----------------------------------------------------------------
subroutine calc_csr(i,csr,ri_2)
 use io,       only:warning, fatal
 use eos,      only:gmw
 use units,    only:umass,udist,unit_density,unit_ergg,utime,unit_velocity
 use part,       only:eos_vars,igasP,xyzh,xyzmh_ptmass,igamma

 real,intent(in), optional     :: ri_2
 integer, intent(in)           :: i
 real,intent(out)              :: csr
 real                          :: T1AU
 real                          :: ri,ri2,cs0,cs_sq,G
 integer                       :: iline, iparams=10, ierr
 character (len=120)           :: discprefix

 if (present(ri_2)) then
    ! If ri2 is passed as an argument, we can skip calculating it.
    !print*, 'ri_2 is present, so skipping calculation of ri2.'
    ri = ri_2**0.5
 else
    ! If ri2 is not passed as an argumnt, calculate from
    ! xyzh and position of sink particle (star).
    !print*, 'ri_2 is not present, calculating ri2 from position.'
    ri2 = ((xyzh(1,i)-xyzmh_ptmass(1,1))**2d0 &
    + (xyzh(2,i)-xyzmh_ptmass(2,1))**2d0 & 
    + (xyzh(3,i)-xyzmh_ptmass(3,1))**2d0)
    ri=ri2**0.5
 endif

 !Work in code units until sound speed, then convert using unit_velocity
 cs0 = H_R*((M_star/R_ref)**0.5)*(R_ref**(q_index))
 !Convert here
 cs0 = cs0
 csr = cs0*ri**(-q_index)

end subroutine calc_csr

 !----------------------------------------------------------------
 !
 ! This is a subroutine for calculating Tr from the sound speed.
 !
 !----------------------------------------------------------------
subroutine calc_Tr(csr,Tr)
 use io,       only:warning, fatal
 use eos,      only:gmw
 use physcon,  only:kb_on_mh
 use units, only: unit_velocity
 real, intent(in) :: csr
 real, intent(out) :: Tr
 real             :: cs_sq

 cs_sq = (csr**2d0)*unit_velocity*unit_velocity
 !print*, 'cs_sq: ', cs_sq
 Tr = cs_sq*(gmw/kb_on_mh)
 !print*, gmw/kb_on_mh
 !print*, Tr

end subroutine calc_Tr

!------------------------------------------------------------------------
!
! read setup file - repurposed from setup_unifdis.f90
!
!------------------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
   use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
   use setunits,     only:read_options_and_set_units
   character(len=*), intent(in)  :: filename
   integer,          intent(out) :: ierr
   integer, parameter :: iunit = 21
   integer :: nerr
   integer          :: size_db
   type(inopts), allocatable :: db(:)
  
   print "(a)",' reading setup options from '//trim(filename)
   nerr = 0
   ierr = 0
   call open_db_from_file(db,filename,iunit,ierr)
   !
   ! units
   !
   call read_options_and_set_units(db,nerr,.false.)
   !
   ! reference radii
   !
   call read_inopt(R_in,'R_in',db,errcount=nerr)
   call read_inopt(R_ref,'R_ref',db,errcount=nerr)
   !
   ! other parameters
   !
   call read_inopt(M_star,'m1',db,errcount=nerr)
   call read_inopt(H_R,'H_R',db,errcount=nerr)
   call read_inopt(q_index,'qindex',db,errcount=nerr)

   size_db = size(db)
   !print*, db

   if (nerr > 0) then
      print "(1x,i2,a)",nerr,' error(s) during read of setup file...'
      ierr = nerr
      return
   endif
  
  end subroutine read_setupfile

end module background_carter
