      program as_error

      implicit none
      character*20 filename
      character*1 abener
      integer ieabs
      common/enercomp/ieabs

c      iebas=1 ! compute er% of absolute energies 
      ieabs=0 ! compute er% of energies relative to the ground state
      write(*,*) " das file name (max 20 characters):"
      read(*,*) filename
      write(*,*) " compute er% of absolute energies? (y/n)"
      read(*,*) abener
      if(abener.eq."Y".or.abener.eq."y") ieabs=1
      call open_files(filename)
      call inp_das()
      call inp_obs()
      call system("~/autovarlambda/asdeck25.x < "//filename)
      call inp_comp()
      call compare_ei()
      call compare_aki()
      call print_ener()
      call print_aki()

      stop
      end

!***********************************************************************
      subroutine open_files(filename)
!***********************************************************************
      implicit none
      character*20 filename

      open(unit=20,file=filename,status='unknown')
      open(unit=25,file='relat_error.dat',status="unknown")

      return
      end

