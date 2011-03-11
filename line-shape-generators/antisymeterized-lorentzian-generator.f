      program cdgenerator

      implicit none
      double precision, parameter :: pi = 3.1415926535897932385d0
      integer :: i, j, k, m, nl, narg
      double precision :: n, ex, kappa, gammma, nu, space, cd, ord, 
     + sum, temp
      character(LEN=50) :: pf, fakespace, fakegamma
      double precision, allocatable :: values(:,:), cdout(:,:),  
     + ordout(:,:), hilb(:,:)

! My pateneted help protocol!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      narg = iargc()
      if (narg .eq. 0) then
         write (*,*) "This program generates a (linear combination of)
     + symmeterized lorentzian(s) and the corresponding antisymmeterized
     + pair"
         write (*,*)
         write (*,*) "Use with a filename and a desired spacing of
     + points followed by the half-width-at-half-height."
         write (*,*)
         write (*,*) "The file you call should have the position of the
     + desired peak in the first column and the desired peak height in
     + the second."
         write (*,*)
         write (*,*) "Example:"
         write (*,*) "antisymmeterized-lorentzian-generator filename
     + spacing hwhh"
         write (*,*)
         write (*,*) "Note: Data will be generated from +/- 20 w.r.t
     + the last position-value"
      stop
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the file names, count lines and all that jazz!!!!!!!!!!!!!!!!!!!
      call getarg(1,pf)
      call getarg(2,fakespace)
      call getarg(3,fakegamma)

      open (16,STATUS="SCRATCH")
      write (16,*) fakespace, fakegamma
      rewind(16)
      read (16,*) space, gammma
 
      nl = 0
 300  open (11,FILE=pf, STATUS="OLD", ERR=666)
      read (11,*,end=200,ERR=666)
      nl = nl + 1
      go to 300
 666  stop 'Error opening file'
 200  continue
      rewind(11)

      allocate(values(nl*2,2), hilb(nl*2,2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do m = 1, nl
         read (11,*,ERR=6000) values(nl+m,1:2)
      end do
      do m = 1, nl
         values(m,1)=-values(nl+1-m,1)
         values(m,2)=-values(nl+1-m,2)
      end do

      open(14,FILE="symmeterized-lorentzian.dat",STATUS="replace")
      call cdgenerate
      open(15,FILE="antisymmeterized-lorentzian.dat",STATUS="replace")
      call ordgenerate

      stop 'Normal termination'
 6000 stop 'Error reading file'
      contains




      subroutine cdgenerate

      ex = -values(nl*2,1)-20
      do while (ex < values(nl*2,1)+20)
      sum = 0
      do k = 1, nl*2
      kappa = values(k,2)
      nu = values(k,1)
      temp = (kappa*(gammma)**2)/( (gammma)**2 + (ex-nu)**2) -
     +    (kappa*(gammma)**2)/( (gammma)**2 + (ex+nu)**2)
      sum = sum + temp
      end do
      write (14,*) ex, sum
      ex = ex + space
      end do

      end subroutine cdgenerate



      subroutine ordgenerate
      
      ex = -values(nl*2,1)-20
      do while (ex < values(nl*2,1)+20)
      sum = 0
      do k = 1, nl*2
      kappa = values(k,2)
      nu = values(k,1)
      temp = (-kappa*(gammma)*(ex-nu))/( (gammma)**2 + (ex-nu)**2) +
     +       (kappa*(gammma)*(ex+nu))/( (gammma)**2 + (ex+nu)**2)

      sum = sum + temp
      end do
      write (15,*) ex, sum
      ex = ex + space
      end do

      end subroutine ordgenerate


      end program cdgenerator
