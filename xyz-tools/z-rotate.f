      program zrotate

      implicit none
      double precision, parameter :: pi = 3.1415926535897932385d0
      character(LEN=1) :: lettern="N"
      character(LEN=2), allocatable :: atoms(:,:)
      character(LEN=50) :: ang, xyz
      integer :: i, j, k, narg, natoms
      double precision :: param(1,2), distance, angle, Degree
      double precision :: x, y, z, nvector(1,3), norm, rotmat(3,3)
      double precision, allocatable :: original(:,:), expanded(:,:)
     +, rotated(:,:), final(:,:), vector(:,:)

      Degree = pi/180.0d0

      narg = iargc()
      if (narg .eq. 0) then
         write (*,*) 'This program takes an XYZ and rotates it around 
     + the Z-axis' 
         write (*,*)
         write (*,*) "Example:"
         write (*,*) "z-rotate angle filename.xyz"
         stop
      end if

      call getarg(1,ang)
      call getarg(2,xyz)

!     I don't like reading from standard in
      open (14,STATUS="SCRATCH")
      write (14,*) ang
      rewind(14)
      read (14,*) angle
      angle = angle*Degree
!     ---

!     Get the number of atoms, coordinates, and the "N" vector
      open(15,FILE=xyz,STATUS="OLD")
      read(15,*) natoms
      read(15,*)

      allocate(original(natoms,3),expanded(natoms,3), rotated(natoms,3)
     +, final(natoms,3), atoms(natoms,1), vector(natoms,3))

      i = 1
      do while (i .le. natoms)
      read(15,*) atoms(i,1), original(i,1), original(i,2), original(i,3)
      i = i + 1
      end do
!     ---


      nvector(1,1) = 0d0
      nvector(1,2) = 0d0
      nvector(1,3) = 1d0
!     ---

      rotmat(1,1) = (1-COS(angle))*nvector(1,1)*nvector(1,1)+COS(angle)
      rotmat(1,2) = (1-COS(angle))*nvector(1,1)*nvector(1,2)
     + - SIN(angle)*nvector(1,3)
      rotmat(1,3) = (1-COS(angle))*nvector(1,1)*nvector(1,3)
     + + SIN(angle)*nvector(1,2)
      rotmat(2,1) = (1-COS(angle))*nvector(1,1)*nvector(1,2)
     + + SIN(angle)*nvector(1,3)
      rotmat(2,2) = (1-COS(angle))*nvector(1,2)*nvector(1,2)+COS(angle)
      rotmat(2,3) = (1-COS(angle))*nvector(1,2)*nvector(1,3)
     + - SIN(angle)*nvector(1,1)
      rotmat(3,1) = (1-COS(angle))*nvector(1,1)*nvector(1,3)
     + - SIN(angle)*nvector(1,2)
      rotmat(3,2) = (1-COS(angle))*nvector(1,2)*nvector(1,3)
     + + SIN(angle)*nvector(1,1)
      rotmat(3,3) = (1-COS(angle))*nvector(1,3)*nvector(1,3)+COS(angle)

      !write (*,*) rotmat(1,1:3)
      !write (*,*) rotmat(2,1:3)
      !write (*,*) rotmat(3,1:3)
 

      i = 1
      do while (i .le. natoms)
           rotated(i,1) = original(i,1)*rotmat(1,1)
     ++ original(i,2)*rotmat(2,1) +  original(i,3)*rotmat(3,1)

           rotated(i,2) = original(i,1)*rotmat(1,2)
     ++ original(i,2)*rotmat(2,2) +  original(i,3)*rotmat(3,2)

           rotated(i,3) = original(i,1)*rotmat(1,3)
     ++ original(i,2)*rotmat(2,3) +  original(i,3)*rotmat(3,3)

      i = i + 1
      end do



!     i = 1
!     do while (i .le. natoms)
!     rotated(i,1) = original(i,1)
!     rotated(i,2) = original(i,2)
!     rotated(i,3) = original(i,3)
!     i = i + 1
!     end do
!     ---

!     Add amount of "N" vector to the rotated coordinates
      i = 1
      do while (i .le. natoms)
      final(i,1) = rotated(i,1) !+ (distance*nvector(1,1))
      final(i,2) = rotated(i,2) !+ (distance*nvector(1,2))
      final(i,3) = rotated(i,3) !+ (distance*nvector(1,3))
      i = i + 1
      end do
!     ---  

!     Dump back an XYZ file
!      open(16,FILE="mod-ligand.xyz",STATUS="REPLACE")
      write(*,*) natoms
      write (*,*) "XYZ spacer" !"nvectro =", nvector(1,1:3)
      i = 1
      do while (i .le. natoms)
      write(*,*) atoms(i,1), final(i,1:3)
      i = i + 1
      end do
!     ---

 5999 stop 'Normal Termination'
 6000 stop 'Error' 

      end program zrotate
