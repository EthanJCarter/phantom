!
!
!
!
!
module clump_tracking
  !
  ! Routine for in-situ output of clump information
  ! Produces fulldumps every order of magnitude of density per clump
  ! & produces a smaller 'average' output file every 0.1 magnitudes of density per clump
  !
  ! :References: None
  !
  ! :Owner: Adam Fenton | Further adapted by Ethan Carter (current maintainer)
  !
  ! :Dependencies: timestep, units, part, readwrite_dumps, io
  !
    implicit none

    integer, public :: i,j,k,l,w,nlines, n_clumps_in_restart, clump_id, n_clumps
    real ,public :: selected_rho, new_exponent
    real,public, dimension(50) ::  full_output_density, clump_dens, clump_temp, small_output_density
    integer,public, dimension(50) :: clump_pid
    integer, public :: restart_file_read_counter = 0

    public specific_output, write_restart_file, read_restart_file, assign_values_from_restart, read_params_file

    private

    contains


     subroutine specific_output(exp_min,exp_max)
       use timestep,  only:time
       use units,         only:unit_density
       use part,             only:igas,massoftype,xyzh,rhoh, npart, vxyzu, nptmass, xyzmh_ptmass,isdead_or_accreted,eos_vars,itemp
       use readwrite_dumps,  only:write_fulldump
       use io,   only:fatal

       !types for OOP
       !Define a derived data type for each set of values
       type :: clump
          real(8), dimension(40) :: num, rho, temp, x, y, z
          integer, dimension(40) :: pid
       end type clump

       integer :: count,density_check,a,density_new, min_index, max_index
       real, intent (in) :: exp_min,exp_max
       logical :: file_exists

       character(len=30) :: format, clump_info
       character(len=13) :: clump_info_file
       character(len=200) ::dumpfile_extension, dumpfile_prefix,runid,particle_id, sub_file
       character(len=500) :: dumpfile, dumpfile_check,dumpfile_check_start

       real  :: chosen, iexp, density_specifier
       !real , dimension(N) :: values !depreciated
       real :: dyn_time_inner_disc, particle_radius, keplerian_velocity, particle_velocity
       integer :: flag1,flag2, new_clump, away_from_sinks
       real :: rhoi, tempi, averho, avetemp, dist2_clumps
       real, dimension(1000) :: sink_flag_debug, clump_flag_debug
       real, dimension(1000) :: distance2, distance2_sinks

       type(clump) :: clumpi


       dyn_time_inner_disc =(10.0**1.5) * (3.15E7/5.023E6) 
         if (n_clumps == 0 .and. time > dyn_time_inner_disc) then
           do i=1, npart
             rhoi = rhoh(xyzh(4,i),massoftype(igas))
               if ((rhoi *unit_density) > exp_min) then
                 clump_dens(1)= rhoi
                 clump_pid(1) = i
                 clump_temp(1) = eos_vars(itemp,i)
                 n_clumps = 1
                 full_output_density(1) = exp_min
                 small_output_density(1) = exp_min
                 write(clump_info,'(A10,I3.3,A4)')"info", n_clumps, ".dat"
                 open(7228,file=clump_info,position='append')
                 write(7228,*) "Number of clumps:", n_clumps
                 write(7228,*) "Number of sinks:", nptmass
                 do k = 1, n_clumps
                   write(7228,*) k,',', &
                   clump_pid(k), &
                   clump_dens(k) * unit_density,',', &
                   xyzh(1,clump_pid(k)),',', &
                   xyzh(2,clump_pid(k)),',', &
                   xyzh(3,clump_pid(k))
                 enddo

                 do j = 1, nptmass
                   write(7228,*)j,',', &
                   xyzmh_ptmass(1,j),',',&
                   xyzmh_ptmass(2,j),',',&
                   xyzmh_ptmass(3,j)

                 enddo
                 close(7228)

                    exit
               endif
             enddo
           endif

           if (n_clumps > 0 .and. time > dyn_time_inner_disc) then
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
                       flag2 = 0 !Set flag to false, don't do anything
                       exit

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
                   !Define output densities for newest clump
                   full_output_density(n_clumps) = exp_min
                   small_output_density(n_clumps) = exp_min

                   !Set clump density and pid to current particle values
                   clump_dens(n_clumps)= rhoi
                   clump_pid(n_clumps) = i
                   clump_temp(n_clumps) = eos_vars(itemp,i)
                   !Write-out information to clump_info.dat
                   !Writes a new file everytime there is a new clump
                   write(clump_info,'(A10,I3.3,A4)') "info",n_clumps,".dat"
                   open(7228,file=clump_info,position='append')
                   write(7228,*) "Number of clumps:", n_clumps
                   write(7228,*) "Number of sinks:", nptmass
                   do k = 1, n_clumps
                     write(7228,*) k,',', &
                     clump_pid(k), &
                     clump_dens(k) * unit_density,',', &
                     xyzh(1,clump_pid(k)),',', &
                     xyzh(2,clump_pid(k)),',', &
                     xyzh(3,clump_pid(k))

                   enddo

                   do j = 1, nptmass
                     write(7228,*)j,',', &
                     xyzmh_ptmass(1,j),',',&
                     xyzmh_ptmass(2,j),',',&
                     xyzmh_ptmass(3,j)
                   enddo
                   close(7228)
                 endif

                 !If it is NOT a new clump and is close enough to an existing clump
                 if (new_clump == 0) then
                   do k=1,n_clumps
                     if (rhoi > clump_dens(k) .and. (distance2(k) < 1)) then ! check if i is alive
                       clump_dens(k)= rhoi !Define new clump dens as particle dens
                       clump_pid(k) = i !Define clump_pid as current particle id
                       clump_temp(k) = eos_vars(itemp,i)
                     endif
                   enddo
                 endif
               endif
             endif
           enddo
         endif

         !Iterate over all clumps found by above checks
         do w = 1, n_clumps

           !print*, "clump_density: ", (clump_dens(w)*unit_density), "clump output: ", full_output_density(w)
           if ((clump_dens(w) * unit_density) .GE. (small_output_density(w)) .and. (small_output_density(w) .LE. exp_max )) then
            !Set all values in temp clump to zero
            clumpi%num = 0.0d0
            clumpi%rho = 0.0d0
            clumpi%temp = 0.0d0
            clumpi%x = 0.0d0
            clumpi%y = 0.0d0
            clumpi%z = 0.0d0
            clumpi%pid = 0

            do i=1,npart
              rhoi = rhoh(xyzh(4,i),massoftype(igas))
              tempi = eos_vars(itemp,i)
              dist2_clumps = ((xyzh(1,i) - xyzh(1,clump_pid(w)))**2 &
              + (xyzh(2,i) - xyzh(2,clump_pid(w)))**2 &
              + (xyzh(3,i) - xyzh(3,clump_pid(w)))**2)
              min_index = minloc(clumpi%rho, dim=1)
              !print '(A, I3.2)', 'Index: ', min_index
              if (dist2_clumps < 100 .and. (rhoi * unit_density) > (minval(clumpi%rho) * unit_density)) then
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
            print '(A, I, E10.3)', 'Average density for clump', w, (averho*unit_density)
            print '(A,I3.3,A8)', 'Writing to: ', w, '.ave.dat'
            write(clump_info_file, '(I3.3,A8)') w, ".ave.dat"
            open(499,file=clump_info_file,position='append')
            write(499, '(I10, 5E10.3, 6F12.3)') clump_pid(w), time,  &
             (clump_dens(w)*unit_density), clump_temp(w), (averho*unit_density), avetemp,  &
            xyzh(1,clump_pid(w)), xyzh(2,clump_pid(w)), xyzh(3,clump_pid(w)), &
            vxyzu(1,clump_pid(w)), vxyzu(2,clump_pid(w)), vxyzu(3,clump_pid(w))

            !Update small_output_dens for that clump
            new_exponent = (log10(small_output_density(w)) + 0.1)
            small_output_density(w) = 10.**new_exponent

            print*, new_exponent
          endif

           !If the clump density is between dens_min and dens_max and is above the current full_output_density, write out a fulldump
           if ((clump_dens(w) * unit_density) .GE. full_output_density(w) .and. (full_output_density(w) .LE. exp_max )) then
             print*, "Trying to write up clump of density: ", full_output_density(w), "and PID: ", clump_pid(w)
             runid = 'lom'

             write(dumpfile_extension, ' (I2)')int(abs(log10(full_output_density(w))) * 10)
             write(clump_id, ' (I5)')w
             write(particle_id,' (I7)')clump_pid(w)
             format = "(A3,A1,I3.3,A1,I3.3,A1,I7.7)"
             write(dumpfile,format)runid,".",w, ".",(int(abs(log10(full_output_density(w))) * 10)),".",clump_pid(w)
             call write_fulldump(time,dumpfile)
             !call write_restart_file()
             !Changed to a logfile-style format
             write(clump_info,'(A4,A1,I3.3,A4)') "info",".",(int(abs(log10(full_output_density(w))) * 10)),".log"
             open(7228,file=clump_info,position='append')
             write(7228,*) "Number of clumps:", n_clumps
             write(7228,*) "Number of sinks:", nptmass

             do k = 1, n_clumps
               write(7228,*) k,',', &
               clump_pid(k), &
               clump_dens(k) * unit_density,',', &
               xyzh(1,clump_pid(k)),',', &
               xyzh(2,clump_pid(k)),',', &
               xyzh(3,clump_pid(k))
             enddo

             do j = 1, nptmass
               write(7228,*)j,',', &
               xyzmh_ptmass(1,j),',',&
               xyzmh_ptmass(2,j),',',&
               xyzmh_ptmass(3,j)
             enddo
             close(7228)

             new_exponent = (log10(full_output_density(w)) + 1)
             full_output_density(w) = 10.**new_exponent

           endif

         enddo


       end subroutine specific_output


       subroutine write_restart_file()
         ! Write a file containing the number of clumps and time of creation with clumpID, clumpPID and next output density
         ! for each clump. This file will be read on a restart to continue clump tracking where a previous run left off.
         use timestep,  only:time
         character(len=12) :: clump_restart_file

         write(*,*) "Writing restart file"
         write(clump_restart_file, '(A7,A1,A4)') "restart","_","file" ! Create a file with name restart_file
         open(2,file=clump_restart_file,status='replace')
         ! This is where the file is written
         write(2,*) n_clumps, time
         do k=1, n_clumps ! Loop over all clumps
           new_exponent = (log10(full_output_density(w)))
           write(2,*) k, clump_pid(w), new_exponent

         enddo

         close(2)
       end subroutine write_restart_file

       !-----DEPRECIATED-----
       subroutine read_restart_file()
         ! We need to read the restart file to initialise a resumed run. The following subroutine reads the
         ! restart_file and assigns the stored values to new arrays.
         ! TO:do
         ! I need to integrate the read_restart_file subroutine with specific output routine by
         ! allocating the read values from the file to the arrays in the specific output subroutine

         integer :: temp,n,clump_particle_ID,io

         ! integer :: temp,n, clump_ID,clump_part_ID,nlines,io
         real :: time,next_density
         integer, dimension(20) :: clump_ids,clump_pids
         real, dimension(20) :: clump_densities,restart_read_test

         ! Count number of lines in file, this is required for correct file reading.
         nlines = 0
         open(1,file='restart_file')

         do
           read(1,*,iostat=io)
           if (io/=0) EXIT
           nlines = nlines + 1

         enddo

         ! Rewind to the start of file to re-read it now that we know how many lines there are.
         rewind(1)
         read(1,*) n_clumps_in_restart,time ! get number of clumps
         close(1)

       end subroutine read_restart_file

       !-----DEPRECIATED-----
       subroutine assign_values_from_restart()
         integer :: temp,n,clump_particle_ID,io
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
             full_output_density(clump_id) = 10**next_density   ! and next density and assign the latter 2
   
           endif                                                 ! to the arrays
         enddo
   
         close(1)
   
       end subroutine assign_values_from_restart
   
       function read_params_file(line1,line2) result(out_vals)
         !Obtain lower and upper limits for tracking from params file
     
           integer, intent(in) :: line1, line2 !Currently 3 and 5 for lower and upper exps
           character(len=128) :: line !current line as string
           integer :: fid, ierr !Unit ID and input error flag
           real(kind=8) :: out1, out2 !Output args so we only call this once
           real, dimension(2) :: out_vals !Output args array
     
           open(newunit=fid, file=trim('tracking.params.dat'), status='old', action='read', iostat=ierr)
           if (ierr/=0) then
             write(*,*) "STOP: error opening file"
             stop
           endif
     
           do i=1, line2
     
             read(fid, '(A)', iostat=ierr) line
             if (ierr/=0) then
               write(*,*) "STOP: error reading line"
               stop
             elseif (i == line1) then !Read expected line for low den
               read(line, *) out1
             elseif (i==line2) then !Read expected line for high den
               read(line, *) out2
             endif
     
           enddo
     
           ! Convert the read string to a real number
           !read(line, *) value
     
           !write(*, *) "Line ", "     :  ", trim(line)
     
           !Display the value converted to real
           !write(*, *) "Read value:", value
     
           close(fid)
     
           out_vals(1) = out1
           out_vals(2) = out2
     
           return
     
         end function read_params_file
   
     end module clump_tracking
