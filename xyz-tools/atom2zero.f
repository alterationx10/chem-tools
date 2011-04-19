      program a20

      implicit none
      double precision, parameter :: pi = 3.1415926535897932385d0
      integer :: i, j, k, narg, nl, natoms, i_a2c
      character(LEN=50) :: filename
      character(LEN=3) :: a2c
      character(LEN=3), allocatable :: atoms(:,:)
      double precision :: x, y, z, dummy, csum, xsum, ysum, zsum
     +, xav, yav, zav
      double precision, allocatable :: cent(:,:)

      narg = iargc()

      if (narg .eq. 0) then
      write (*,*) "This program uses the XYZ format.
     + Simply use the XYZ filname as the argument"
      write (*,*) "Example:"
      write (*,*) "a20 filename.xyz atom#"
      write (*,*) "Note:"
      stop
      end if

      ! Get arguments
      call getarg(1,filename)
      call getarg(2,a2c)
      read (a2c,*) i_a2c

      !Read in XYZ file (count lines)
      nl=0
      open (15,FILE=filename, STATUS="OLD",ERR=6666)
 100  read (15,*,end=200,ERR=6001)
      nl = nl + 1
      go to 100
 6666 stop 'Error opening xyz file'
 200  continue
      rewind(15)

      ! Allocate arrays
      allocate(cent(nl-2,4),atoms(nl-2,1))

      ! Read XYZ into array
      read (15,*,ERR=6002) natoms
      read (15,*,ERR=6003)
      if (natoms .ne. nl-2) then
        stop 'Number atoms and actual number of 
     + atoms in XYZ file do no agree'
      end if
      do i=1,nl-2
        read (15,*,ERR=6004) atoms(i,1), cent(i,1:3)
      end do

      ! Dump to std. out
      write (*,*) natoms
      write (*,*) 
      do i = 1, natoms
      write (*,*) atoms(i,1),
     +            cent(i,1) - cent(i_a2c,1),
     +            cent(i,2) - cent(i_a2c,2), 
     +            cent(i,3) - cent(i_a2c,3)
      end do


 5999 stop 'Normal Termination'
 6001 stop 'Error 1' 
 6002 stop 'Error 2'
 6003 stop 'Error 3'
 6004 stop 'Error 4'
 6005 stop 'Error 5'


      end program a20
