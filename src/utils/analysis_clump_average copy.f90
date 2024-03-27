!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
    !
    ! Analysis routine for finding and outputting clump information and averages wrt time
    !
    !
    ! :References: None
    !
    ! :Owner: Ethan Carter | Adam Fenton (see specific_output.f90)
    !
    ! :Runtime parameters: None
    !
    ! :Dependencies: part
    !
     implicit none
     character(len=20), parameter, public :: analysistype = 'clump_history'
     public :: do_analysis
    
     private
    
    contains
    
    subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)

        use io,         only:fatal
        use units,      only:unit_density
        use part,       only:xyzmh_ptmass,vxyz_ptmass,ihacc,rhoh,massoftype,igas, itemp, nptmass,xyzmh_ptmass,isdead_or_accreted, eos_vars
        use physcon,    only:pi
        use io,         only:fatal

        !types for OOP
        !Define a derived data type for each set of values
        type :: clump
          real(8), dimension(1000) :: num, rho, temp, x, y, z
          integer, dimension(1000) :: pid
        end type clump

        !Generic declarations for analysis subroutine
        character(len=*), intent(in)    :: dumpfile
        real,dimension(:,:), intent(in) :: xyzh,vxyzu
        real, intent(in)                :: pmass,time
        integer, intent(in)             :: npart,iunit,numfile

        !For tracking
        integer :: i,j,nlines, n_clumps_in_restart, clump_id
        integer :: k = 1
        integer :: n_clumps = 0
        logical :: file_exists
        real    :: rhoi, dist_sq
        real, dimension(1000) :: sink_flag_debug, clump_flag_debug
        real, dimension(1000) :: distance2, distance2_sinks
        type(clump) :: clump_hist

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

        print*, 'Acquiring neighbour particles to clump centre...'

        do i=1,npart
          distance2 = ((xyzh(1,i))**2) &
          + ((xyzh(2,i))**2) &
          + (xyzh(3,i)**2)

          if (distance2 .LE. 100) then
            clump_hist(j)%pid = i
            clump_hist(j)%x = xyzh(1,i)
            clump_hist(j)%y = xyzh(2,i)
            clump_hist(j)%z = xyzh(3,i)

            j = j + 1
            print*, 'Neighbour added with pid ', i

          endif

        enddo


    end subroutine do_analysis

end module analysis
