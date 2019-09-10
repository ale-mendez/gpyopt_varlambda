      program test

      call open_files()
      call inp_das()
      call inp_exact()
      call write_das_2()

      stop
      end
!***********************************************************************
      subroutine write_das_2()
!***********************************************************************
      implicit none

      integer DMX
      parameter (DMX=100)

      integer i,j
      integer ipol,ismin
      integer ncfg,MXCONF,NZION,MCFMX
      character fhead*50,fsalgebini*34,fmxvorb*16,fmxconf*16
      character fsalgebfin*15,fsalgeb*131
      character arad*4,acup*5
      character anorb*2,forb*13,fcfg*13
      integer ishift
      character fsminim*14,fnzion*15,fend*7,fmaxe*14
      character forthog*16,fprint*15,fradout*16
      character fnlam*14,fsto*15,fpol*46
      character aorthog*5,aform*6,ayes*5
      character fformat*200,fformat0*200
      integer removeblanks,ib
      character dum*250,dumm5*5,dumm4*4
      real*8 alfd,rcut
      integer ialfd,ircut,iend
      character*5 chalfd,chrcut

      character*250 das
      integer ndas
      common/dasbck/das(DMX),ndas
      common/salgebbck/MXCONF,NZION,MCFMX

      ipol=0
      ncfg=mxconf
      do 100 i=1,ncfg+3
         write(666,'(a)') das(i+ishift)
100   continue
      ismin=ncfg+4
      if(ipol.eq.1) then
         alfd=1.000
         rcut=1.000
         dum=das(ismin)
         do 110 i=1,245
            dumm5=dum(i:i+5)
            dumm4=dum(i:i+4)
            if (dumm5.eq.'ALFD=') ialfd=i
            if (dumm5.eq.'RCUT=') ircut=i
            if (dumm4.eq.'&END') iend=i
110      continue
         write(*,*) iend
         write(chalfd,'(f5.3)') alfd
         write(chrcut,'(f5.3)') rcut
         write(666,*)dum(1:ialfd-1)
         write(666,*)dum(ialfd:ialfd+4)//chalfd
         write(666,*)dum(ircut:ircut+4)//chrcut
         write(666,*)dum(iend:iend+3)
      else
         write(666,*) das(ismin)
      endif
      do 120 i=ismin+1,ndas
         write(666,*) das(i)
120   continue

1008  format(5(f12.9,1x))
1009  format(10(i3))

      return
      end

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

