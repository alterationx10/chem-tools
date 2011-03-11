      program cnc

      implicit none
      double precision, parameter :: pi = 3.1415926535897932385d0
      integer :: i, j, k, narg, nl, natoms
      character(LEN=50) :: filename
      character(LEN=3), allocatable :: atoms(:,:)
      double precision :: x, y, z, dummy, csum, xsum, ysum, zsum
     +, xav, yav, zav
      double precision, allocatable :: cent(:,:)

      narg = iargc()

      if (narg .eq. 0) then
      write (*,*) "This program uses the XYZ format.
     + Simply use the XYZ filname as the argument"
      write (*,*) "Example:"
      write (*,*) "cnc filename.xyz"
      write (*,*) "Note:"
      write (*,*) "Currently supports elements 1-99
     + (and the elements are CASE SENSITIVE!)"
      stop
      end if

      call getarg(1,filename)

      nl=0
      open (15,FILE=filename, STATUS="OLD",ERR=6666)
 100  read (15,*,end=200,ERR=6001)
      nl = nl + 1
      go to 100
 6666 stop 'Error opening xyz file'
 200  continue
      rewind(15)

      allocate(cent(nl-2,4),atoms(nl-2,1))

       read (15,*,ERR=6002) natoms
       read (15,*,ERR=6003)

       if (natoms .ne. nl-2) then
          stop 'Number atoms and actual number of 
     + atoms in XYZ file do no agree'
       end if
       

       do i=1,nl-2
         read (15,*,ERR=6004) atoms(i,1), cent(i,1:3)
         call atomcharge
       end do

      csum=0d0
      do i = 1, nl-2
         csum = csum + cent(i,4)
      end do

      xsum=0d0
      do i = 1, nl-2
         xsum = xsum + cent(i,1)*cent(i,4)
      end do

      ysum=0d0
      do i = 1, nl-2
         ysum = ysum + cent(i,2)*cent(i,4)
      end do

      zsum=0d0
      do i = 1, nl-2
         zsum = zsum + cent(i,3)*cent(i,4)
      end do

      xav = xsum/csum
      yav = ysum/csum
      zav = zsum/csum


      write (*,*) "The center of Nuclear Charge is:"
      write (*,*) xav, yav, zav
!      write (*,*) xsum, ysum, zsum, csum





 5999 stop 'Normal Termination'
 6001 stop 'Error 1' 
 6002 stop 'Error 2'
 6003 stop 'Error 3'
 6004 stop 'Error 4'
 6005 stop 'Error 5'

      contains

      subroutine atomcharge

      if (atoms(i,1) .eq. "H") then
      cent(i,4) = 1
      elseif (atoms(i,1) .eq. "He" ) then 
      cent(i,4) = 2
      elseif (atoms(i,1) .eq. "Li" ) then
      cent(i,4) = 3
      elseif (atoms(i,1) .eq. "Be" ) then 
      cent(i,4) = 4
      elseif (atoms(i,1) .eq. "B" ) then 
      cent(i,4) = 5
      elseif (atoms(i,1) .eq. "C" ) then 
      cent(i,4) = 6
      elseif (atoms(i,1) .eq. "N" ) then 
      cent(i,4) = 7
      elseif (atoms(i,1) .eq. "O" ) then 
      cent(i,4) = 8
      elseif (atoms(i,1) .eq. "F" ) then 
      cent(i,4) = 9
      elseif (atoms(i,1) .eq. "Ne" ) then 
      cent(i,4) = 10
      elseif (atoms(i,1) .eq. "Na" ) then 
      cent(i,4) = 11
      elseif (atoms(i,1) .eq. "Mg" ) then 
      cent(i,4) = 12
      elseif (atoms(i,1) .eq. "Al" ) then 
      cent(i,4) = 13
      elseif (atoms(i,1) .eq. "Si" ) then 
      cent(i,4) = 14
      elseif (atoms(i,1) .eq. "P" ) then 
      cent(i,4) = 15
      elseif (atoms(i,1) .eq. "S" ) then 
      cent(i,4) = 16
      elseif (atoms(i,1) .eq. "Cl" ) then 
      cent(i,4) = 17
      elseif (atoms(i,1) .eq. "Ar" ) then 
      cent(i,4) = 18
      elseif (atoms(i,1) .eq. "K" ) then 
      cent(i,4) = 19
      elseif (atoms(i,1) .eq. "Ca" ) then 
      cent(i,4) = 20
      elseif (atoms(i,1) .eq. "Sc" ) then 
      cent(i,4) = 21
      elseif (atoms(i,1) .eq. "Ti" ) then 
      cent(i,4) = 22
      elseif (atoms(i,1) .eq. "V" ) then 
      cent(i,4) = 23
      elseif (atoms(i,1) .eq. "Cr" ) then 
      cent(i,4) = 24
      elseif (atoms(i,1) .eq. "Mn" ) then 
      cent(i,4) = 25
      elseif (atoms(i,1) .eq. "Fe" ) then 
      cent(i,4) = 26
      elseif (atoms(i,1) .eq. "Co" ) then 
      cent(i,4) = 27
      elseif (atoms(i,1) .eq. "Ni" ) then 
      cent(i,4) = 28
      elseif (atoms(i,1) .eq. "Cu" ) then 
      cent(i,4) = 29
      elseif (atoms(i,1) .eq. "Zn" ) then 
      cent(i,4) = 30
      elseif (atoms(i,1) .eq. "Ga" ) then 
      cent(i,4) = 31
      elseif (atoms(i,1) .eq. "Ge" ) then 
      cent(i,4) = 32
      elseif (atoms(i,1) .eq. "As" ) then
       cent(i,4) = 33
      elseif (atoms(i,1) .eq. "Se" ) then
       cent(i,4) = 34
      elseif (atoms(i,1) .eq. "Br" ) then 
      cent(i,4) = 35
      elseif (atoms(i,1) .eq. "Kr" ) then 
      cent(i,4) = 36
      elseif (atoms(i,1) .eq. "Rb" ) then 
      cent(i,4) = 37
      elseif (atoms(i,1) .eq. "Sr" ) then 
      cent(i,4) = 38
      elseif (atoms(i,1) .eq. "Y" ) then 
      cent(i,4) = 39
      elseif (atoms(i,1) .eq. "Zr" ) then 
      cent(i,4) = 40
      elseif (atoms(i,1) .eq. "Nb" ) then 
      cent(i,4) = 41
      elseif (atoms(i,1) .eq. "Mo" ) then 
      cent(i,4) = 42
      elseif (atoms(i,1) .eq. "Tc" ) then 
      cent(i,4) = 43
      elseif (atoms(i,1) .eq. "Ru" ) then 
      cent(i,4) = 44
      elseif (atoms(i,1) .eq. "Rh" ) then 
      cent(i,4) = 45
      elseif (atoms(i,1) .eq. "Pd" ) then 
      cent(i,4) = 46
      elseif (atoms(i,1) .eq. "Ag" ) then 
      cent(i,4) = 47
      elseif (atoms(i,1) .eq. "Cd" ) then 
      cent(i,4) = 48
      elseif (atoms(i,1) .eq. "In" ) then 
      cent(i,4) = 49
      elseif (atoms(i,1) .eq. "Sn" ) then 
      cent(i,4) = 50
      elseif (atoms(i,1) .eq. "Sb" ) then 
      cent(i,4) = 51
      elseif (atoms(i,1) .eq. "Te" ) then 
      cent(i,4) = 52
      elseif (atoms(i,1) .eq. "I" ) then 
      cent(i,4) = 53
      elseif (atoms(i,1) .eq. "Xe" ) then 
      cent(i,4) = 54
      elseif (atoms(i,1) .eq. "Cs" ) then 
      cent(i,4) = 55
      elseif (atoms(i,1) .eq. "Ba" ) then 
      cent(i,4) = 56
      elseif (atoms(i,1) .eq. "La" ) then 
      cent(i,4) = 57
      elseif (atoms(i,1) .eq. "Ce" ) then 
      cent(i,4) = 58
      elseif (atoms(i,1) .eq. "Pr" ) then 
      cent(i,4) = 59
      elseif (atoms(i,1) .eq. "Nd" ) then 
      cent(i,4) = 60
      elseif (atoms(i,1) .eq. "Pm" ) then 
      cent(i,4) = 61
      elseif (atoms(i,1) .eq. "Sm" ) then 
      cent(i,4) = 62
      elseif (atoms(i,1) .eq. "Eu" ) then 
      cent(i,4) = 63
      elseif (atoms(i,1) .eq. "Gd" ) then 
      cent(i,4) = 64
      elseif (atoms(i,1) .eq. "Tb" ) then 
      cent(i,4) = 65
      elseif (atoms(i,1) .eq. "Dy" ) then 
      cent(i,4) = 66
      elseif (atoms(i,1) .eq. "Ho" ) then 
      cent(i,4) = 67
      elseif (atoms(i,1) .eq. "Er" ) then 
      cent(i,4) = 68
      elseif (atoms(i,1) .eq. "Tm" ) then 
      cent(i,4) = 69
      elseif (atoms(i,1) .eq. "Yb" ) then 
      cent(i,4) = 70
      elseif (atoms(i,1) .eq. "Lu" ) then 
      cent(i,4) = 71
      elseif (atoms(i,1) .eq. "Hf" ) then 
      cent(i,4) = 72
      elseif (atoms(i,1) .eq. "Ta" ) then 
      cent(i,4) = 73
      elseif (atoms(i,1) .eq. "W" ) then 
      cent(i,4) = 74
      elseif (atoms(i,1) .eq. "Re" ) then 
      cent(i,4) = 75
      elseif (atoms(i,1) .eq. "Os" ) then 
      cent(i,4) = 76
      elseif (atoms(i,1) .eq. "Ir" ) then 
      cent(i,4) = 77
      elseif (atoms(i,1) .eq. "Pt" ) then 
      cent(i,4) = 78
      elseif (atoms(i,1) .eq. "Au" ) then 
      cent(i,4) = 79
      elseif (atoms(i,1) .eq. "Hg" ) then 
      cent(i,4) = 80
      elseif (atoms(i,1) .eq. "Tl" ) then 
      cent(i,4) = 81
      elseif (atoms(i,1) .eq. "Pb" ) then
       cent(i,4) = 82
      elseif (atoms(i,1) .eq. "Bi" ) then 
      cent(i,4) = 83
      elseif (atoms(i,1) .eq. "Po" ) then 
      cent(i,4) = 84
      elseif (atoms(i,1) .eq. "At" ) then 
      cent(i,4) = 85
      elseif (atoms(i,1) .eq. "Rn" ) then 
      cent(i,4) = 86
      elseif (atoms(i,1) .eq. "Fr" ) then 
      cent(i,4) = 87
      elseif (atoms(i,1) .eq. "Ra" ) then 
      cent(i,4) = 88
      elseif (atoms(i,1) .eq. "Ac" ) then 
      cent(i,4) = 89
      elseif (atoms(i,1) .eq. "Th" ) then 
      cent(i,4) = 90
      elseif (atoms(i,1) .eq. "Pa" ) then 
      cent(i,4) = 91
      elseif (atoms(i,1) .eq. "U" ) then 
      cent(i,4) = 92
      elseif (atoms(i,1) .eq. "Np" ) then 
      cent(i,4) = 93
      elseif (atoms(i,1) .eq. "Pu" ) then 
      cent(i,4) = 94
      elseif (atoms(i,1) .eq. "Am" ) then
       cent(i,4) = 95
      elseif (atoms(i,1) .eq. "Cm" ) then 
      cent(i,4) = 96
      elseif (atoms(i,1) .eq. "Bk" ) then 
      cent(i,4) = 97
      elseif (atoms(i,1) .eq. "Cf" ) then 
      cent(i,4) = 98
      elseif (atoms(i,1) .eq. "Es" ) then
       cent(i,4) = 99
      else 
      stop 'You have invalid elements in your XYZ file'
      end if

      return
      end subroutine atomcharge


      end program cnc
