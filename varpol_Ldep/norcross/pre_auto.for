      program as_error

      implicit none
      character*20 filename

      write(*,*) "das file name (max 20 characters):"
      read(*,*) filename
      call open_files(filename)
      call inp_das()
      call inp_exact()
      call system("~/autovarlambda/asdeck25.x < "//filename)
      call inp_computed()
      call ener_out()
      call print_ener()

      stop
      end

!***********************************************************************
      subroutine open_files(filename)
!***********************************************************************
      implicit none
      character*20 filename

      open(unit=20,file=filename,status='unknown')
      open(unit=30,file='exactvalues.dat',status='unknown')
      open(unit=25,file='relat_error.dat',status="unknown")
!... later      open(unit=40,file='TERMS',status='unknown')
!... later      open(unit=25,file='tmp',status='unknown')

      return
      end

