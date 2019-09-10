!
!  set of subroutines for varlambda
!
!  Alejandra Mendez - 6/03/19
!  v1.3
!

!***********************************************************************
      subroutine open_files()
!***********************************************************************
      implicit none

      open(unit=20,file='das',status='unknown')
      open(unit=30,file='exactvalues.dat',status='unknown')
!... later      open(unit=40,file='TERMS',status='unknown')
!... later      open(unit=25,file='tmp',status='unknown')

      return
      end

!***********************************************************************
      subroutine run_as(lam0,NVMX)
!***********************************************************************
      integer, intent (in) :: NVMX
      real*8, intent(in) :: lam0(NVMX)
cf2py integer, intent (in) :: NVMX
cf2py real*8, intent(in) :: lam0(NVMX)

      call lambdas(lam0,NVMX)
      call modify_das
      call system("~/autovarlambda/asdeck25.x < tmp")
      call inp_computed
      call ener_out

      return
      end

!***********************************************************************
      subroutine lambdas(lam0,NVMX)
!***********************************************************************
!  this subroutine arranges the values of lambda in the array lam
!    nvlam : number of lambda parameters that will be varied
!    lam   : array of lambda parameters
!    ips   : flags and gives the number of pseudostates
!    ntlam  : total number of lambdas
      implicit none
      integer, intent (in) :: NVMX
      real*8, intent (in) :: lam0(NVMX)
      integer NMX
      parameter (NMX=45)
      integer DMX
      parameter (DMX=100)

      integer i
!      integer invi,invf
      real*8 linp

      real*8 lam
      integer nvlam,ntlam,ips,NLAM
      common/lambck/lam(NMX),nvlam,ips,ntlam,NLAM
      common/laminpbck/linp(DMX)

      if(NLAM.gt.NMX) then
         write(*,*) "Increase NMX! stop. NMX=",NLAM
         stop
      endif

!     Initialize lam array
      do 100 i=1,NLAM
         lam(i) = 1.0
100   continue

!     write varying lambda values
      nvlam = NVMX
      do 200 i=1,nvlam
         lam(i) = lam0(i)
!         print*,ips,i,lam0(i)
         if(ips.gt.0.and.i.ge.ips) lam(i) = -1.0*lam0(i)
!         print*,lam(i)                   !test
200   continue
      do i=nvlam+1,NLAM
         lam(i) = linp(i)
!         print*,i,lam(i)
      enddo

      return
!999     return
      end

!***********************************************************************
      subroutine modify_das()
!***********************************************************************
      implicit none
      integer DMX,NMX
      parameter (DMX=100,NMX=45)

      integer i
      integer ndas,nvlam,ips,ntlam,NLAM
      integer MXCONF,NZION,MCFMX
      character*250 das
      real*8 lam
      common/dasbck/das(DMX),ndas
      common/lambck/lam(NMX),nvlam,ips,ntlam,NLAM
      common/salgebbck/MXCONF,NZION,MCFMX

      open(unit=25,file='tmp',status='unknown')
      call write_das()
      write(25,1000) (lam(i),i=1,NLAM)
      !print*,ntlam,(lam(i),i=1,NLAM)
      if(MCFMX.ne.0) then
         write(25,1001) (i,i=1,MCFMX)
      endif
1000  format(5(f12.8,1x))
1001  format(10(i2,1x))
      close(unit=30)
      close(unit=25)
      return
      end

!***********************************************************************
      subroutine write_das()
!***********************************************************************
      implicit none
      integer DMX
      parameter (DMX=100)

      integer i,j
      integer ndas
      character*250 das
      character*250 line
      common/dasbck/das(DMX),ndas
      
      do 100 i=1,ndas
         write(25,'(A)') das(i)
         line = das(i)
         do j=1,4
            if(line(j:j+6).eq.'SMINIM') goto 101
         enddo
100   continue
101   continue

      return
      end


