      program rot2zero

      implicit none
      double precision, parameter :: pi = 3.1415926535897932385d0
      integer :: i, j, k, narg, nl
      character(LEN=50) :: filename
      double precision :: x, y, z, dummy, rsum
      double precision, allocatable :: rotsum(:,:)

      narg = iargc()

      !A little help for the user
      if (narg .eq. 0) then
      write (*,*) "This program uses a '.cddat' format.
     + Simply use a file that you would read into plotspec"
      write (*,*) "Example:"
      write (*,*) "rot2zero fenchone.cddat"
      write (*,*) "Note:"
      stop
      end if

      !Get the filename
      call getarg(1,filename)

      !Count the lines
      nl=0
      open (15,FILE=filename, STATUS="OLD",ERR=6666)
 100  read (15,*,end=200,ERR=6001)
      nl = nl + 1
      go to 100
 6666 stop 'Error opening cddat file'
 200  continue
      rewind(15)

      !We want our array to be 5 lines shorter than the .cddat file (header is useless)
      allocate(rotsum(nl-5,2))
     


       !We will read past the first 5 lines
       read (15,*,ERR=6002)
       read (15,*,ERR=6002)
       read (15,*,ERR=6002)
       read (15,*,ERR=6002)
       read (15,*,ERR=6002)


       !Read in Rot. Strengths to the array
       do i=1,nl-5
         read (15,*,ERR=6003) dummy , rotsum(i,1)
       end do

      !Keep a running sum
      rsum=0d0
      do j=1,nl-5
         rsum=rsum+rotsum(j,1)
         rotsum(j,2)=rsum
      end do

      !Write to std. out
      write (*,'(A, A)') "	   Rot. Str. (x10^-40 cgs)	Running Sum"
      do k=1,nl-5
      write (*,*) k, rotsum(k,1:2)
      end do 

 5999 stop 'Normal Termination'
 6001 stop 'Error 1' 
 6002 stop 'Error 2'
 6003 stop 'Error 3'
 6004 stop 'Error 4'



      end program rot2zero
