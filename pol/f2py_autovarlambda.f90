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
!...later       open(unit=15,file='tmp',status="unknown")
      open(unit=25,file='relat_error.dat',status="unknown")
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
      call write_das
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
      integer nvlam,ntlam,ips
      common/lambck/lam(NMX),nvlam,ips,ntlam
      common/laminpbck/linp(DMX)

      if(ntlam.gt.NMX) then
         write(*,*) "Increase NMX! stop. NMX=",ntlam
         stop
      endif

!     Initialize lam array
      do 100 i=1,ntlam
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
      do i=nvlam+1,ntlam
         lam(i) = linp(i)
!         print*,i,lam(i)
      enddo

      return
      end

!***********************************************************************
      subroutine run_as_pol(pol)
!***********************************************************************
      implicit none
      integer PDIM 
      parameter (PDIM=2)
      real*8, intent(in) :: pol(PDIM)
cf2py real*8, intent(in) :: pol(PDIM)

      call varpol(pol)
      call write_das()
      call system("~/autovarlambda/asdeck25.x < tmp")
      call inp_computed
      call ener_out

      return
      end

!***********************************************************************
      subroutine varpol(pol)
!***********************************************************************
      implicit none
      integer PDIM 
      parameter (PDIM=2)
      real*8, intent (in) :: pol(PDIM)
      integer ipol
      real*8 alfd,rcut
      common/polbck/alfd,rcut,ipol

      alfd=pol(1)
      rcut=pol(2)

      return
      end

!***********************************************************************
      subroutine write_das()
!***********************************************************************
      implicit none

      integer DMX
      parameter (DMX=100)
      integer NMX
      parameter (NMX=45)

      integer i
      character fhead*50,fsalgebini*34,fmxvorb*16,fmxconf*16
      character fsalgebfin*15,fsalgeb*131
      character arad*4,acup*5
      character anorb*2,forb*13,fcfg*13
      character fsminim*14,fnzion*15,fend*7,fmaxe*14
      character forthog*16,fprint*15,fradout*16
      character fnlam*14,fsto*15,fpol*49
      character aorthog*5,aform*6,ayes*5
      character fformat*200,fformat0*200
      integer removeblanks,ib
      character*16 falfd,frcut

      character*250 das
      integer ndas
      CHARACTER(LEN=5) RAD
      CHARACTER(LEN=4) CUP
      CHARACTER(LEN=3) ORTHOG
      integer MXVORB,MXCONF,KUTSO,NZION,MAXE,NLAM,MCFMX
      real*8 lam
      integer nvlam,ntlam,ips
      integer ipol
      real*8 alfd,rcut

      common/dasbck/das(DMX),ndas
      common/salgebbck/MXVORB,MXCONF,KUTSO,RAD,CUP
      common/sminimbck/NZION,MAXE,NLAM,MCFMX,ORTHOG
      common/lambck/lam(NMX),nvlam,ips,ntlam
      common/polbck/alfd,rcut,ipol

      open(unit=15,file='tmp',status="unknown")

      do i=1,ntlam
         lam(i) = 1.000
      enddo

!... write header
      fhead="('A.S. Automatically generated with oompaloompa',/"
!... write namelist SALGEB
      fsalgebini="'&SALGEB RAD=',a4,1x,'CUP=',a5,1x,"
      if(MXVORB.lt.10) fmxvorb="'MXVORB=',i1,1x,"
      if(MXVORB.ge.10) fmxvorb="'MXVORB=',i2,1x,"
      if(MXVORB.ge.100) fmxvorb="'MXVORB=',i3,1x,"
      if(MXCONF.lt.10) fmxconf="'MXCONF=',i1,1x,"
      if(MXCONF.ge.10) fmxconf="'MXCONF=',i2,1x,"
      if(MXCONF.ge.100) fmxconf="'MXCONF=',i3,1x,"
      fsalgebfin="'KUTSO=0 &END')"
      fsalgeb=fhead//fsalgebini//fmxvorb//fmxconf//fsalgebfin
      arad="'NO'"
      acup="'ICM'"
      write(15,fsalgeb) arad,acup,MXVORB,MXCONF
      if(MXVORB.lt.10) write(anorb,'(i1)') MXVORB
      if(MXVORB.ge.10) write(anorb,'(i2)') MXVORB
      forb="("//anorb//"(i3,i2))"
      fcfg="(i4,"//anorb//"(i5))"
!... write orbitals + electronic configurations
      do 100 i=1,MXCONF+1
        write(15,'(a)') das(i+2)
100   continue
!... write namelist SMINIM
      fsminim="('&SMINIM',1x,"
      if(abs(NZION).lt.10) fnzion="'NZION=',i2,1x,"
      if(abs(NZION).ge.10) fnzion="'NZION=',i3,1x,"
      if(MAXE.lt.10) fmaxe="'MAXE=',i1,1x,"
      if(MAXE.gt.10) fmaxe="'MAXE=',i2,1x,"
      if(NLAM.lt.10) fnlam="'NLAM=',i1,1x,"
      if(NLAM.ge.10) fnlam="'NLAM=',i2,1x,"
      fsto=' '
      if(NZION.lt.0.and.MXCONF.lt.10) fsto="'MCFMX=',i1,1x,"
      if(NZION.lt.0.and.MXCONF.ge.10) fsto="'MCFMX=',i2,1x,"
      forthog="'ORTHOG=',a5,1x,"
      fprint="'PRINT=',a6,1x,"
      fradout="'RADOUT=',a5,1x,"
      if(NZION.gt.0) fformat=fsminim//fnzion//fmaxe//fnlam//forthog//
     +      fprint//fradout
      if(NZION.lt.0) fformat=fsminim//fnzion//fmaxe//fnlam//fsto//
     +      forthog//fprint//fradout
      ib=removeblanks(fformat)
      if(ipol.eq.1) then
         falfd="'ALFD=',f5.3,1x,"
         frcut="'RCUT=',f5.3,1x,"
!         fpol=falfd//frcut//"'IPOLFN=9999',1x,"
         fpol=falfd//frcut//"1x,"
         fformat=fformat(1:ib)//fpol
         ib=removeblanks(fformat)
      endif
      fend="'&END')"
      fformat0=fformat(1:ib)//fend
      aorthog="'YES'"
      aform="'FORM'"
      ayes="'YES'"
      if(ipol.eq.0) then
         if(NZION.gt.0) 
     +      write(15,fformat0) NZION,MAXE,NLAM,aorthog,aform,ayes
         if(NZION.lt.0) 
     +      write(15,fformat0) NZION,MAXE,NLAM,MXCONF,aorthog,aform,
     +                          ayes
      elseif(ipol.eq.1) then
         if(NZION.gt.0) 
     +      write(15,fformat0) NZION,MAXE,NLAM,aorthog,aform,ayes,
     +                          alfd,rcut
         if(NZION.lt.0) 
     +      write(15,fformat0) NZION,MAXE,NLAM,MXCONF,aorthog,aform,
     +                          ayes,alfd,rcut
            print*,alfd,rcut
      endif
!      write(15,1008) (lam(i),i=1,nvlam)
      write(15,1008) (lam(i),i=1,ntlam)
      if(NZION.lt.0) write(15,1009) (i,i=1,MXVORB)

1008  format(5(f12.9,1x))
1009  format(10(i3))

      close(unit=15)
      return
      end

c***********************************************************************
      integer function removeblanks(string)
c***********************************************************************
      character*200 string
      integer i
      do i=1,200
         if(string(i:i+1).eq.'  ') then
            removeblanks=i-1
            return
         endif
      enddo
      return
      end
