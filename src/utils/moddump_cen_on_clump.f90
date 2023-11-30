!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module moddump
    !
    ! Centre on a given particle based on user-input
    !
    ! :References: None
    !
    ! :Owner: Ethan Carter
    !
    ! :Runtime parameters: None
    !
    ! :Dependencies: part
    !
     implicit none
    
    contains
    
    subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
     use part, only: xyzmh_ptmass,vxyz_ptmass,nptmass
     integer, intent(inout) :: npart
     integer, intent(inout) :: npartoftype(:)
     real,    intent(inout) :: massoftype(:)
     real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
     real,dimension(3)   :: pid_xyz,pid_vxyz
     integer :: i, pid

     !Ask for particle id
     print*, 'Total number of particles: ', npart
     print*, 'Total number of sinks: ', nptmass

     write(*,*) 'Please enter a suitable integer particle ID: '
     read(*,*) pid

     print*, 'Chosen PID is', pid, '. centering particles relative to chosen centre...'
     
     !Obtain a copy of xyz and vxyz of target centre
     pid_xyz  = xyzh(1:3, pid)
     pid_vxyz = vxyzu(1:3,pid)

          !if (nptmass >= 3) then
     print*, 'Centering sinks relative to chosen centre...'
     do i=1,nptmass
         xyzmh_ptmass(1:3,i) = (xyzmh_ptmass(1:3,i) - pid_xyz)
         vxyz_ptmass(1:3,i) = (vxyz_ptmass(1:3,i) - pid_vxyz)
         print*, xyzmh_ptmass(1:3,i)
         !print*, i
     enddo
     !endif

     do i=1, npart
      if (pid >= 1 .and. pid <= npart) then
         !print*, 'Before re-centering - Particle ', i, ':', xyzh(1:3,i), vxyzu(1:3,i)
         xyzh(1:3,i) = (xyzh(1:3,i) - pid_xyz)
         vxyzu(1:3,i) = (vxyzu(1:3,i) - pid_vxyz)
      else
          ! Handle invalid particle ID
          print*, 'Invalid particle ID'
          return
      endif
     enddo

     print*, 'Particles centered on ', pid

     return
     end subroutine modify_dump
end module moddump
