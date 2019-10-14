      program testprint

      implicit none
      integer ipol,lrcut
      real*8 valfd,vrcut
      common/polbck/valfd,vrcut,ipol,lrcut

      ipol=1
      lrcut=1
      valfd=0.006
      vrcut=0.410

      call open_files
      call inp_das
      call inp_exact
      call write_das

      stop
      end


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
      character fnlam*14,fsto*15,fpol*56
      character aorthog*5,aform*6,ayes*5
      character fformat*200,fformat0*200
      integer removeblanks,ib
      character*28 falfd,frcut

      character*250 das
      integer ndas
      CHARACTER(LEN=5) RAD
      CHARACTER(LEN=4) CUP
      CHARACTER(LEN=3) ORTHOG
      integer MXVORB,MXCONF,KUTSO,NZION,MAXE,NLAM,MCFMX
      real*8 ALFD(3),RCUT(3)
      real*8 lam
      integer nvlam,ntlam,ips
      integer ipol,lrcut
      real*8 valfd,vrcut

      common/dasbck/das(DMX),ndas
      common/salgebbck/MXVORB,MXCONF,KUTSO,RAD,CUP
      common/sminimbck/ALFD,RCUT,NZION,MAXE,NLAM,MCFMX,ORTHOG
      common/lambck/lam(NMX),nvlam,ips,ntlam
      common/polbck/valfd,vrcut,ipol,lrcut

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
         alfd(lrcut)=valfd
         rcut(lrcut)=vrcut
         falfd="'ALFD=',2(f6.4,','),f6.4,1x,"
         frcut="'RCUT=',2(f6.4,','),f6.4,1x,"
!         fpol=falfd//frcut//"'IPOLFN=9999',1x,"
!         print*,'*',falfd,'*'
!         print*,'*',frcut,'*'
         fpol=falfd//frcut
!         print*,'*',fpol,'*'
         fformat=fformat(1:ib)//fpol
         ib=removeblanks(fformat)
      endif
      fend="'&END')"
      fformat0=fformat(1:ib)//fend
      print*,'*',fformat(1:ib),'*'
      print*,'*',fformat0,'*'
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
      implicit none
      character*200 string
      integer i
      do i=1,200
         if(string(i:i+1).eq.'  ') then
            removeblanks=i-1
            go to 999
         endif
      enddo

999   return
      end
