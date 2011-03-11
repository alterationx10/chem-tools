        program sos

! This program will calculate the SOS from a ".cddat" type
! output file (from turbomole).       

!Decalartions go here        
        implicit none
        integer :: i, j, k ! These are for dummy variables
        integer :: narg, nl ! These are more specific
        double precision :: sovers, mr
        double precision, parameter :: rcm2au=8065.4804d0*27.21138d0, &
        au2ev = 27.21138d0, Na = 0.07732d0, conv = 2d0*235.726372d0, &
        bau2cgs = 1.342219d0*0.0001d0
        character(LEN=50) :: datafile
        double precision, allocatable :: data(:,:)

!Help the user
        narg = iargc()

        if (narg .eq. 0) then
                write (*,*) "This program uses a '.cddat'-type file &
                (from a Turbomole output). This file should have FIVE &
               lines of header (!) and two columns of data: RCM, R(cgs)"
                write (*,*) "Usage:"
                write (*,*) "sos file.cddat"
        end if

!Find out about the data file
        call getarg(1,datafile)
        call count(nl,datafile,15)
        ! The first 5 lines are header here
        nl = nl - 5

!Allocate the array and load in the data
        allocate(data(nl,2)) 
        read (15,*)
        read (15,*)
        read (15,*)
        read (15,*)
        read (15,*)
        do i = 1, nl
                read(15,*) data(i,1:2)
        end do


!Lets convert to au here, and help the user out some more

        if (data(nl,1) < 1000) then
                stop 'Are you sure your .cddat file is in RCM?'
        end if

        do i = 1, nl
                data(i,1) = data(i,1)/rcm2au
        end do 


        do i = 1, nl-1
                if (data(i,1) > data(i+1,1)) then
                   stop 'Your data is not sorted from lowest to highest'
                end if
        end do


!Lets do the math!
        do j = 1, nl
                sovers=0d0
                mr=0d0
                call subsos(j, sovers)
                mr=bau2cgs*Na*Na*rcm2au*rcm2au*sovers/100d0

                write (*,*) j, data(j,1)*au2ev, sovers, mr
        end do                

! Clean up
        deallocate(data)
        stop 'Normal Termination'
 6666   stop 'Error reading cdspectrum file'

        contains        

        subroutine count(lines, file, fn)
        integer :: lines, fn
        character(LEN=50) :: file

        lines = 0
        open (fn,FILE=file, STATUS="OLD",ERR=6660)
 100    read (fn,*,end=200,ERR=6660)
        nl = nl + 1
        go to 100
 6660   stop 'Error opening/reading file'
 200    continue
        rewind(fn)
        return

        end subroutine

        subroutine subsos(nexcit, sum)
        integer :: nexcit, l
        double precision :: sum
        sum = 0d0
        do l = 1, nexcit
                sum = sum + (data(l,2)/conv) / &
                (data(l,1)*data(l,1) - Na*Na)
        end do
        sum = (2d0/3d0)*sum

        return

        end subroutine

        end program sos
