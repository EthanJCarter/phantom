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
     character(len=20), parameter, public :: analysistype = 'clump_avereages'
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
          real(8), dimension(40) :: num, rho, temp, x, y, z
          integer, dimension(40) :: pid
        end type clump

        !Define a derived data type for each set of values
        type :: sink
          real(8) :: num, rho, temp, x, y, z
        end type sink

        !Generic declarations for analysis subroutine
        character(len=*), intent(in)    :: dumpfile
        real,dimension(:,:), intent(in) :: xyzh,vxyzu
        real, intent(in)                :: pmass,time
        integer, intent(in)             :: npart,iunit,numfile

        !For tracking
        integer :: i,j,k,l,nlines, n_clumps_in_restart, clump_id
        integer :: n_clumps = 0
        real    :: selected_rho, new_exponent, averho, tempi, avetemp
        real,    dimension(50) ::  clump_output_density, clump_dens, clump_temp
        integer, dimension(50) :: clump_pid, rhotemp
        integer :: restart_file_read_counter = 0
        integer :: count,density_check,a,density_new, w
        real    :: exp_min = 1e-9
        logical :: file_exists
        character(len=13) :: clump_info_file
        real    :: chosen, iexp, density_specifier
        real    :: particle_radius, keplerian_velocity, particle_velocity
        integer :: flag1,flag2, new_clump, away_from_sinks
        real    :: rhoi, dist_sq
        real, dimension(1000) :: sink_flag_debug, clump_flag_debug
        real, dimension(1000) :: distance2, distance2_sinks
        ! Find the index of the minimum value in clumpi%rho
        integer :: min_index, max_index
        type(clump) :: clumpi


        if (n_clumps == 0) then
            print*, exp_min
            do i=1, npart
                rhoi = rhoh(xyzh(4,i),massoftype(igas))
                  if ((rhoi*unit_density) > exp_min) then
                    clump_dens(1)= rhoi
                    clump_pid(1) = i
                    n_clumps = 1
                    clump_output_density(1) = exp_min
                    clump_temp(1) = eos_vars(itemp,i)
                    print '(A29,E10.1)', "Candidate particle with rho >",clump_output_density(1)

                    do k=1,1
                        !Ugly print statement as cocacenating
                        !the format causes a fault for some reason
                        print '(A9)', "Location:"
                        print '(A2,(F10.3),1x)', 'x=',xyzh(1,clump_pid(k))
                        print '(A2(F10.3),1x)', 'y=',xyzh(2,clump_pid(k))
                        print '(A2(F10.3))', 'z=',xyzh(3,clump_pid(k))
                        print '(A, E10.3, E10.3)', 'rhoi=', (clump_dens(k) *unit_density)
                        print '(A, E10.3)', 'temp=', clump_temp(k)
                        !print '(A9, F10.3, F10.3, F10.3)', "Location:", &
                        !xyzh(1,clump_pid(k)), ',', &
                        !xyzh(2,clump_pid(k)), ',', &
                        !xyzh(3,clump_pid(k))
                    enddo

                exit
                  endif
            enddo
        endif

        if (n_clumps > 0) then
            iexp = exp_min
            do i=1, npart
              if (.not. isdead_or_accreted(xyzh(4,i))) then ! i.e. if the particle is alive and hasn't been accreted by any sink
                rhoi = rhoh(xyzh(4,i),massoftype(igas))
              if ((rhoi *unit_density) > exp_min) then !If above dens_min, check if close to other clumps or sinks
                do k=1, n_clumps
                  distance2(k) = ((xyzh(1,i) - xyzh(1,clump_pid(k)))**2 &
                                  + (xyzh(2,i) - xyzh(2,clump_pid(k)))**2 &
                                  + (xyzh(3,i) - xyzh(3,clump_pid(k)))**2)
                enddo

                do j=1, nptmass
                  distance2_sinks(j) = (xyzh(1,i) - xyzmh_ptmass(1,j))**2 &
                                        + (xyzh(2,i) - xyzmh_ptmass(2,j))**2 &
                                        + (xyzh(3,i) - xyzmh_ptmass(3,j))**2
                enddo

                new_clump = 1       ! Flag to determine if new clump is away from other clumps
                away_from_sinks = 1 ! Flag to determine if new clump is away from sinks
                flag1 = 1 !Default flags to true
                flag2 = 1
                do k=1,n_clumps
                  !if within a certain distance (100 code units?)
                  if (distance2(k) < 100) then
                    flag1 = 0 !Set flag to false, don't do anything
                  endif
                enddo
  
                new_clump=new_clump*flag1

                if (nptmass > 1) then
                    do j=1,nptmass
                      if (distance2_sinks(j) < 100) then
                        flag2 = 0 !Set flag to false, don't do anything
                      endif
                    enddo
    
                    away_from_sinks=away_from_sinks*flag2 
    
                    !!! DEBUGGING
                    do j=1, nptmass
                      sink_flag_debug(j) = away_from_sinks
    
                    enddo
                  endif
    
                  !if it is new and far away enough from existing sinks
                  if (new_clump==1 .and. away_from_sinks==1) then
                    !Add to clump count
                    n_clumps = n_clumps + 1
                    !redefine output density - is this needed?
                    clump_output_density(k) = exp_min
    
                    !set clump density and pid to current particle values
                    clump_dens(k)= rhoi
                    clump_pid(k) = i
                    clump_temp(k) = eos_vars(itemp,i)

                    print '(A29, E10.1)', "Candidate particle with rho >",exp_min
                    print '(A9)', "Location:"
                    print '(A2, F10.3,1x)', 'x=',xyzh(1,clump_pid(k))
                    print '(A2, F10.3,1x)', 'y=',xyzh(2,clump_pid(k))
                    print '(A2, F10.3)', 'z=',xyzh(3,clump_pid(k))
                    print '(A, E10.3, E10.3)', 'rhoi=', (clump_dens(k) *unit_density)
                    print '(A, E10.3)', 'temp=', clump_temp(k)
    
                  endif
    
                  !If it is NOT a new clump and is close enough to an existing clump
                  if (new_clump == 0) then
                    do k=1,n_clumps
                      if (rhoi > clump_dens(k) .and. (distance2(k) < 1)) then ! check if i is alive
                        clump_dens(k)= rhoi !Define new clump dens as particle dens
                        clump_pid(k) = i !Define clump_pid as current particle id
                      endif
                    enddo
                  endif
                endif
              endif
            enddo
          endif

          do k=1,n_clumps
            ! Set all values in temp clump to zero
            clumpi%num = 0.0d0
            clumpi%rho = 0.0d0
            clumpi%temp = 0.0d0
            clumpi%x = 0.0d0
            clumpi%y = 0.0d0
            clumpi%z = 0.0d0
            clumpi%pid = 0
            do i=1,npart
              rhoi = rhoh(xyzh(4,i),massoftype(igas)) * unit_density
              tempi = eos_vars(itemp,i)
              dist_sq = ((xyzh(1,i) - xyzh(1,clump_pid(k)))**2 &
              + (xyzh(2,i) - xyzh(2,clump_pid(k)))**2 &
              + (xyzh(3,i) - xyzh(3,clump_pid(k)))**2)
              min_index = minloc(clumpi%rho, dim=1)
              !print '(A, I3.2)', 'Index: ', min_index
              if (dist_sq < 100 .and. rhoi > minval(clumpi%rho)) then
                clumpi%rho(min_index) = rhoi ! Replace smallest value with the desired value
                clumpi%temp(min_index) = tempi
                !clumpi%x(min_index) = xyzh(1,i)
                !clumpi%y(min_index) = xyzh(2,i)
                !clumpi%z(min_index) = xyzh(3,i)
                !clumpi%pid(min_index) = i
              endif
            enddo
            max_index = maxloc(clumpi%rho, dim=1)
            averho  = (sum(clumpi%rho))/40 !Size of clumpi is hardcoded for ease for now.
            avetemp = (sum(clumpi%temp))/40
            print '(A, I, E10.3)', 'Average density for clump', k, averho
            print '(A,I3.3,A10)', 'Writing to: ', k, '.small.dat'
            write(clump_info_file, '(I3.3,A10)') k, ".small.dat"
            open(499,file=clump_info_file,position='append')
            write(499, '(I10, 4E10.3, 6F12.3)') clump_pid(k), (clump_dens(k)*unit_density), &
            clump_temp(k), averho, avetemp,  &
            xyzh(1,clump_pid(k)), xyzh(2,clump_pid(k)), xyzh(3,clump_pid(k)), &
            vxyzu(1,clump_pid(k)), vxyzu(2,clump_pid(k)), vxyzu(3,clump_pid(k))
            enddo

    end subroutine do_analysis

end module analysis
