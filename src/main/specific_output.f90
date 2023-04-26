module density
    implicit none

     INTEGER, public :: i,j,k,l,den_min,den_max,nlines, n_clumps_in_restart, clump_id, n_clumps
     REAL ,public :: selected_rho, new_exponent
     real, dimension(50) ::  clump_output_density, clump_dens
     integer, dimension(50) :: clump_pid
     integer :: restart_file_read_counter = 0
     integer, private :: w

   contains


     subroutine specific_output(den_min,den_max,N)
       use timestep,  only:time
       use units,         only:unit_density
       use part,             only:igas,massoftype,xyzh,rhoh, npart, vxyzu, nptmass, xyzmh_ptmass,isdead_or_accreted
       use readwrite_dumps,  only:write_fulldump
       use io,   only:fatal

       INTEGER :: den_min, den_max, count,density_check,a,density_new
       integer, intent (in) :: N
       logical :: file_exists

       CHARACTER(LEN=30) :: Format, clumps_and_sinks
       character(len=7) :: clump_info_file
       character(len=200) ::dumpfile_extension,clump_id, dumpfile_prefix,runid,particle_id, sub_file
       character(len=500) :: dumpfile, dumpfile_check,dumpfile_check_start

       real  :: chosen, den,density_specifier, new_exponent
       real , dimension(N) :: values
       real :: dyn_time_inner_disc, particle_radius, keplerian_velocity, particle_velocity
       integer :: flag1,flag2, new_clump, away_from_sinks
       ! real, dimension(50) :: clump_dens
       real :: rhoi
       real, dimension(50) :: sink_flag_debug, clump_flag_debug
       real, dimension(50) :: distance2, distance2_sinks

       dyn_time_inner_disc =(10.0**1.5) * (3.15E7/5.023E6)
         IF (n_clumps == 0 .and. time > dyn_time_inner_disc) then
           DO i=1, npart
             rhoi = rhoh(xyzh(4,i),massoftype(igas))
               IF ((rhoi *unit_density) > 1E-9) then
                 clump_dens(1)= rhoi
                 clump_pid(1) = i
                 n_clumps = 1
                 den = den_min
                 clump_output_density(1) = 10.0**den
                 write(clumps_and_sinks,'(A16,I3.3,A4)')"clumps_and_sinks",n_clumps, ".dat"
                 open(7728,file=clumps_and_sinks,position='append')
                 write(7728,*) "Number of clumps:", n_clumps
                 write(7728,*) "Number of sinks:", nptmass
                 do k = 1, n_clumps
                   do j = 1, nptmass
                     write(7728,*)'Clumps: ' // NEW_LINE('A'),&
                     k,',', &
                     clump_pid(k), &
                     clump_dens(k) * unit_density,',', &
                     xyzh(1,clump_pid(k)),',', &
                     xyzh(2,clump_pid(k)),',', &
                     xyzh(3,clump_pid(k)), NEW_LINE('A'), &
                     'Sinks: ' // NEW_LINE('A'),&
                     j,',', &
                     xyzmh_ptmass(1,j),',',&
                     xyzmh_ptmass(2,j),',',&
                     xyzmh_ptmass(3,j)
   
                   enddo
                 enddo
                 close(7728)
   
                 exit
               END IF
           END DO
         END IF
   
   
   
   
         IF (n_clumps > 0 .and. time > dyn_time_inner_disc) then
           den = -9
           DO i=1, npart
             IF (.not. isdead_or_accreted(xyzh(4,i))) then ! i.e. if the particle is alive and hasn't been accreted by any sink
               rhoi = rhoh(xyzh(4,i),massoftype(igas))
                 IF ((rhoi *unit_density) > 1E-9) then
                   DO k=1, n_clumps
   
                     distance2(k) = ((xyzh(1,i) - xyzh(1,clump_pid(k)))**2 &
                                   + (xyzh(2,i) - xyzh(2,clump_pid(k)))**2 &
                                   + (xyzh(3,i) - xyzh(3,clump_pid(k)))**2)
                   enddo
                   DO j=1, nptmass
                     distance2_sinks(j) = (xyzh(1,i) - xyzmh_ptmass(1,j))**2 &
                                        + (xyzh(2,i) - xyzmh_ptmass(2,j))**2 &
                                        + (xyzh(3,i) - xyzmh_ptmass(3,j))**2
                   enddo
   
   
   
                   new_clump = 1       ! Flag to determine if new clump is away from other clumps
                   away_from_sinks = 1 ! Flag to determine if new clump is away from sinks
                   flag1 = 1
                   flag2 = 1
                   do k=1,n_clumps
                       if (distance2(k) < 100) then
                         flag1 = 0
                         exit
                       endif
                   enddo
                   new_clump=new_clump*flag1
   
                   !!! DEBUGGING
                   do k=1,n_clumps
                     clump_flag_debug(k) = new_clump
                   enddo
   
                   if (nptmass > 1) then
                     do j=1,nptmass
                         if (distance2_sinks(j) < 100) then
                           flag2 = 0
                           exit
                         endif
                     enddo
                     away_from_sinks=away_from_sinks*flag2
   
                   !!! DEBUGGING
                     do j=1, nptmass
                       sink_flag_debug(j) = away_from_sinks
                     enddo
                   endif
   
   
                   if (new_clump==1 .and. away_from_sinks==1) then
                     n_clumps = n_clumps + 1
                     clump_output_density(n_clumps) = 10.0**den
                     clump_dens(n_clumps)= rhoi
                     clump_pid(n_clumps) = i
                     write(clumps_and_sinks,'(A16,I3.3,A4)')"clumps_and_sinks",n_clumps, ".dat"
                     open(7728,file=clumps_and_sinks,position='append')
                     write(7728,*) "Number of clumps:", n_clumps
                     write(7728,*) "Number of sinks:", nptmass
                     do k = 1, n_clumps
                       write(7728,*) k,',', &
                       clump_pid(k), &
                       clump_dens(k) * unit_density,',', &
                       xyzh(1,clump_pid(k)),',', &
                       xyzh(2,clump_pid(k)),',', &
                       xyzh(3,clump_pid(k))
                     enddo
                       do j = 1, nptmass
                         write(7728,*)j,',', &
                         xyzmh_ptmass(1,j),',',&
                         xyzmh_ptmass(2,j),',',&
                         xyzmh_ptmass(3,j)
                     enddo
                     close(7728)
                   endif
   
   
   
                   if (new_clump == 0) then
                     do k=1,n_clumps
                       if (rhoi > clump_dens(k) .and. (distance2(k) < 1)) then ! check is i is alive
                         clump_dens(k)= rhoi
                         clump_pid(k) = i
   
   
                       end if
   
                     enddo
                   endif
                 end if
             END IF
           end DO
         ENDIF
   
         do w = 1, n_clumps
   
           if ((clump_dens(w) * unit_density) .GE. clump_output_density(w) .and. (clump_output_density(w) .LE. 1e-3 )) then
             runid = 'run1'
             write(dumpfile_extension, ' (I2)')int(abs(log10(clump_output_density(w))) * 10)
             write(clump_id, ' (I5)')w
             write(clump_info_file, '(I3.3,A4)')w, ".dat"
             write(particle_id,' (I7)')clump_pid(w)
   
             Format = "(A4,A1,I3.3,A1,I7.7,A1,I3.3)"
             write(dumpfile,Format)runid,".",w, ".",clump_pid(w),".",(int(abs(log10(clump_output_density(w))) * 10))
             call write_fulldump(time,dumpfile)
   !          call write_restart_file()
   
   
   
             open(499,file=clump_info_file,position='append')
             write(499,*) (clump_dens(w) * unit_density),",",&
                        xyzh(1,clump_pid(w)),",",xyzh(2,clump_pid(w)),",",xyzh(3,clump_pid(w)),",",&
                       vxyzu(1,clump_pid(w)),",",vxyzu(2,clump_pid(w)),",",vxyzu(3,clump_pid(w))
             close(499)
             new_exponent = (log10(clump_output_density(w)) + 1)
             clump_output_density(w) = 10.**new_exponent
             ! if (new_exponent == -7) then
             !   write(*,*)
             !   stop 0
             ! endif
   
           end if
   
         end do
   
   
   
   
     END subroutine specific_output
   
   
   
     subroutine write_restart_file()
       ! Write a file containing the number of clumps and time of creation with clumpID, clumpPID and next output density
       ! for each clump. This file will be read on a restart to continue clump tracking where a previous run left off.
       use timestep,  only:time
       character(len=12) :: clump_restart_file
   
       real :: new_exponent
   
   
   
       write(*,*) "Writing restart file"
       write(clump_restart_file, '(A7,A1,A4)') "restart","_","file" ! Create a file with name restart_file
       open(2,file=clump_restart_file,status='replace')
       ! This is where the file is written
       write(2,*) n_clumps, time
       do k=1, n_clumps ! Loop over all clumps
         new_exponent = (log10(clump_output_density(k)))
         write(2,*) k, clump_pid(k), new_exponent
       enddo
       close(2)
     end subroutine write_restart_file
   
   
   
     subroutine read_restart_file()
       ! We need to read the restart file to initialise a resumed run. The following subroutine reads the
       ! restart_file and assigns the stored values to new arrays.
       ! TO:DO
         ! I need to integrate the read_restart_file subroutine with specific output routine by
         ! allocating the read values from the file to the arrays in the specific output subroutine
   
       integer :: temp,n, clump_ID,clump_particle_ID,io
   
       ! integer :: temp,n, clump_ID,clump_part_ID,nlines,io
       real :: time,next_density
       integer, dimension(20) :: clump_ids,clump_pids
       real, dimension(20) :: clump_densities,restart_read_test
   
   
       ! Count number of lines in file, this is required for correct file reading.
       nlines = 0
       open(1,file='restart_file')
       DO
         READ(1,*,iostat=io)
         IF (io/=0) EXIT
         nlines = nlines + 1
       END DO
   
   
   
      ! Rewind to the start of file to re-read it now that we know how many lines there are.
      REWIND(1)
      !
      read(1,*) n_clumps_in_restart,time ! get number of clumps
      !
      close(1)
   
     end subroutine read_restart_file
   
   
   
     subroutine assign_values_from_restart()
       integer :: temp,n, clump_ID,clump_particle_ID,io
       real :: time,next_density
       integer, dimension(20) :: clump_ids,clump_pids
       real, dimension(20) :: clump_densities,restart_read_test
       open(1,file='restart_file')
       do n=1,nlines
         if (n==1) then               ! This will read just the first line of the file containing
           read(1,*) n_clumps, time   ! the number of clumps and the time the file was created
         else
           read(1,*) clump_id,clump_particle_id,next_density   ! For all other lines in the file, from 2
           clump_pid(clump_id) = clump_particle_id             ! to n_clumps, read in clumpid, clumpPid
           clump_output_density(clump_id) = 10**next_density   ! and next density and assign the latter 2
         endif                                                 ! to the arrays
       end do
       close(1)
   
     end subroutine assign_values_from_restart
   
   end module density