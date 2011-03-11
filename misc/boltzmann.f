      program boltzmann

      implicit none
      double precision, parameter :: pi=3.1415926535897932385d0,
     + r=8.314472d0, na=6.02214179d23 
      integer :: i, j, k, narg, nl, nconf
      real :: printtemp
      character(LEN=50) :: filename, atoms
      character(LEN=50), allocatable:: conformer(:,:)
      double precision :: x, y, z, dummy, distsum, temp, kb, norm
      double precision, allocatable :: energy(:,:),dist(:,:),unnorm(:,:)

!Lets get the basics out of the way!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      narg = iargc()
      if (narg .eq. 0) then
         write (*,*) "call this program with a properly formated input
     + file"
         write (*,*) "boltzmann filename"
         write (*,*) "Example file:"
         write (*,*) "3 298 #Number of coformers and the Temperature in
     +  Kelvin"
         write (*,*) "ConfomerName1 ConformerEnergy1 #(in kJ/mol)"
         write (*,*) "ConfomerName1 ConformerEnergy1 #(in kJ/mol)"
         write (*,*) "ConfomerName1 ConformerEnergy1 #(in kJ/mol)"
         go to 6002
      end if

      call getarg(1,filename)

      nl=0
      open (15,FILE=filename, STATUS="OLD",ERR=6666)
 100  read (15,*,end=200,ERR=6000)
      nl = nl + 1
      go to 100
 6666 stop 'Error opening file'
 200  continue
      rewind(15)

      allocate(conformer(nl-1,1),energy(nl-1,1),dist(nl-1,1),
     + unnorm(nl-1,1))

      read (15,*) nconf, temp
      printtemp = temp
      if (nconf .ne. nl-1) then
         go to 6003
      end if

      do i=1,nl-1
         read (15,*,ERR=6000) conformer(i,1), energy(i,1)
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Energy Check!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      k=0
      do i=1,nl-1
         if (energy(i,1) .eq. 0d0) then
            k = k + 1
         end if
      end do

      if ( k .ne. 1) then
         go to 6001
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate the boltzmann factors for the confomers
      do i=1,nl-1
          dist(i,1) = Exp(-(energy(i,1)*1000d0)/(r*temp))
      end do
      distsum=0d0
      do i=1,nl-1
         distsum = distsum + dist(i,1)
      end do
      do i=1,nl-1
         unnorm(i,1) = (dist(i,1) / distsum)
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Make sure it is normalized!!!!!!!!!!!!!!!!!!!!!!!!!
      norm=0d0
      do i=1,nl-1
         norm = norm + unnorm(i,1)
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Print the results!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write (*,*) "The results are:"
      do i=1,nl-1
      write (*,409) conformer(i,1), ((unnorm(i,1))/norm)*100
     + , "% at ", temp, "K"
 409  format (A,F5.2,A,F6.2,A)
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!FINISH IT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 5999 stop 'Normal Termination'
 6000 stop 'Error with File'
 6001 stop 'Error: Lowest conformation not set to 0' 
 6002 stop
 6003 stop 'The number of conformers that you think are in your file
     + and the number there actually are
     + are not in agreement'
      end program boltzmann
