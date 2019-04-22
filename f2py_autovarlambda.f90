!
!  set of subroutines for varlam_hyperopt.py
!
!  Alejandra Mendez - 6/03/19
!  v1.3
!

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
      parameter (NMX=25)
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
!     if
!         ips.eq.0 => no pseudostates are used
!         ips.gt.0 => pseudostates are used
!         ips.gt.nvlam+1 => lambda values are not varied and pseudostates are used
!                          >> print lam=1. for i=nvlam+1 to ips
!                          >> print lam=-1. for ips
!      ntlam=nvlam+1
!      if(ips.eq.0) then
!         lam(ntlam)=1.0
!         goto 999
!      elseif(ips.gt.0) then
!         if(ips.eq.ntlam) then
!            lam(ntlam)=-1.0*lam0(ntlam)
!         elseif(ips.gt.ntlam) then
!            invi=ntlam
!            invf=ips-1
!            do 300 i=invi,invf
!               lam(i)=1.0
!300         continue
!            ntlam=ips
!            lam(ips)=-1.0
!          endif
!       endif

      return
!999     return
      end
!***********************************************************************
      subroutine inp_das()
!***********************************************************************
      IMPLICIT REAL*8 (A-H,O-P,R-Z)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT INTEGER*4 (Q)

      integer DMX
      parameter (DMX=100)
      integer NMX
      parameter (NMX=25)

      integer i,ic
      real*8 linp
      character*250 line
      CHARACTER(LEN=8) TITLE,RTITLE
      CHARACTER(LEN=5) RAD
      CHARACTER(LEN=4) PHASE,RUN
      CHARACTER(LEN=4) CUP,BORN,BASIS
      CHARACTER(LEN=3) CPU,TARGET
      CHARACTER(LEN=1) XDR
      CHARACTER(LEN=3) RADOUT,POTOUT,POTL,ORTHOG,STONLZ,FAC,MIXBV
      CHARACTER(LEN=4) POTIN,PPOT             !TERMINATOR    ,MYRGE
      CHARACTER(LEN=6) PRINT,TCC

      character*250 das
      integer ndas
      real*8 lam
      integer nvlam,ips,ntlam

      common/dasbck/das(DMX),ndas
      common/salgebbck/MXCONF,NZION,MCFMX
      common/lambck/lam(NMX),nvlam,ips,ntlam,NLAM
      common/laminpbck/linp(DMX)

      namelist/SALGEB/
     X BASIS,BDISK,BFOT,BORN,BPASS,
     X CPU,CUP,
     X ECNTRB,
     X FACTJ,FACTL,
     X ICFG,IDIAG,IDW,IFILL,INAST,INASTJ,IPAR,ISCALR,ITANAL,
     X J2MIN,J2MAX,JNAST,IJNAST,                     !RM STGJK VARIABLES
     X K1,K2,KCOR1,KCOR2,KCUT,KCUT0,KCUTCC,KCUTI,KCUTP,KORB1,KORB2,
     X KPOL0,KPOLE,KPOLM,KSUBCF,KUTDSK,KUTM1,KUTOO,KUTOOX,KUTSO,
     X KUTSS,KUTSSX,
     X LCMIN,LCON,LCON0,LCONI,LCONDW,LCONDWJ,LMAX
     X,LRGLAM,LVMIN,LVMAX,LXTRA,
     X MBP1MX,MBP2MX,MEKVMX,             !NOT FOR USER JOE TO MESS WITH
     X MAXJFS,MAXLAM,MXLAMX,MAXLX,MAXLOO,MCFSS,MDEL,MENGB,
     X MINJT,MAXJT,MINJTP,MAXJTP,
     X MINLT,MAXLT,MINLTP,MAXLTP,MINLTS,MAXLTS,
     X MINST,MAXST,MINSTP,MAXSTP,MINSTS,MAXSTS,
     X MODE,MPRINT,MSTART,MSTRT0,MXCCF,MXCONF,MXVORB,
     X NAST,NASTJ,NASTJP,NASTP,NASTS,
     X NMETA,NMETAJ,NMETAG,NMETGJ,NXTRA,
     X PHASE,
     X QCL0,QCS0,QCUT,QQCUT,
     X RAD,RTITLE,RUN,
     X TARGET,TITLE,
     X XDR
      NAMELIST/SMINIM/
     X AJUSTX,ALAV,ALFD,ATM,
     X BALAN,
     X CMXLSA,CMXLSR,CMXICA,CMXICR,
     X DHNS0,
     X ECNTRB,ECORR,EIMXLS,EIMXIC,ESKPL,ESKPH,
     X FAC,
     X ibreit,IBWRM,ICAV,IDIAG,IFIX,IGAGR,IGAUGE,IMAXIT,INCLUD,INUKE,
     X IOPTIM,IORT,IPOLFN,IREL,irtard,ISHFTLS,ISHFTIC,ITANAL,ITOL,
     X IWGHT,IXTRA,IZESP,
     X JPRINT,JRAD,JLOWMN,JLOWMX,JUPMN,JUPMX,
     X KCUT,KTCC,KUTCA,
     X LLOWMN,LLOWMX,LUPMN,LUPMX,
     X M,MAXE,MAXLAM,MCFMX,MCFSTO,MDELE,MDEN,MEXPOT,MEXTRE,MGRP,MHF,
     X MIXBV,MPNCH,MPRINT,MPSEUD,MRAD,MRED,MSTEP,MULTS,
     X NDEN,NFIX,NLAM,NLAMD,NLAMQ,NMETAP,NMETAPJ,NMETAR,NMETARJ,NOCC,
     X NPE,NRSLMX,NVAR,NVARD,NVARQ,NZION,
     X ORTHOG,
     X POTIN,POTL,POTOUT,PPOT,PRINT,
     X QED,
     X RAD,RADOUT,RCAV,RCUT,RMIN1,RMIN2,RNUKE,RZERO,
     X SCALER,SKIN,STOLB,STONLZ,
     X TCC,TINORB,TK0,TOLB,TOLE,TOLTCC,TVARY,
     X WLG1,WLG2,
     X xmax,
     X ZESP

      MCFMX=0
      ic=0
      do 100 i=1,DMX
         read(20,'(A)',end=101) line
         ic=ic+1
         das(i)=line
!         write(111,'(A)') das(i)         !test
         linp(i)=1.d0
100   continue
101   ndas=ic
      rewind(20)
      NZION=0
      ips=0
      read(20,SALGEB)
      read(20,SMINIM)
      read(20,*) (linp(i),i=1,NLAM)
      do 200 i=1,NLAM
!        print*,i,linp(i)
        if(linp(i).lt.0) then
            ips=i
            goto 201
        endif
200   continue
201   if(ips.ne.0) write(*,*) " >> Pseudostates are considered << "

      close(unit=20)
      return
      end

!***********************************************************************
      subroutine inp_exact()
!***********************************************************************
      implicit none
      integer LMX
      parameter (LMX=300)

      integer i,ir
      integer sr,lr,pr,cfr,ni
      real*8 er

      integer se,le,pe,cfe,ne
      real*8 ee,eegr
      common/exactbck/se(LMX),le(LMX),pe(LMX),cfe(LMX),ee(LMX),eegr,ne

      read(30,*)
      ir=0
      do 100 i=1,LMX
         read(30,*,end=101) sr,lr,pr,cfr,ni,er
         if(sr.eq.0) goto 101
         if(sr.ne.0.and.i.eq.LMX) then
            print*,'Increase LMX'
            stop
         endif
         ir=ir+1
         se(i)=sr
         le(i)=lr
         cfe(i)=cfr
         ee(i)=er
100   continue
101   ne=ir
      eegr=er
!      print*,ne,eegr                     !test

      return
      end

!***********************************************************************
      subroutine modify_das()
!***********************************************************************
      implicit none
      integer DMX,NMX
      parameter (DMX=100,NMX=25)

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

!***********************************************************************
      subroutine inp_computed()
!***********************************************************************
      implicit none
      integer LMX
      parameter (LMX=300)

      integer i,j,ir
      integer sr,lr,pr,cfr,ni
      real*8 er
      integer ss,ll,jmin,jmax,gi

      integer sc,lc,pc,cfc,nc,gic
      real*8 ec,ecgr
      common/compbck/sc(LMX),lc(LMX),pc(LMX),cfc(LMX),ec(LMX),ecgr,nc,
     |                  gic(LMX)

      open(unit=40,file='TERMS',status='unknown')
      read(40,*)
      ir=0
      do 100 i=1,LMX
         read(40,*) sr,lr,pr,cfr,ni,er
         if(sr.eq.0) goto 101
         if(sr.ne.0.and.i.eq.LMX) then
            print*,'Increase LMX'
            stop
         endif
         ir=ir+1
         sc(i)=sr
         lc(i)=lr
         cfc(i)=cfr
         ec(i)=er
100   continue
101   nc=ir
      ecgr=er
!      print*,nc,ecgr                     !test

!     compute statistic weigth gi
      do 200 i=1,nc
         ss=(sc(i)-1)/2
         ll=lc(i)
         jmin=abs(ss-ll)
         jmax=ss+ll
         gi=0
         do 210 j=jmin,jmax,1
            gi=gi+2*j+1
210      continue
         gic(i)=gi
200   continue

      close(unit=40)

      return
      end

!***********************************************************************
      subroutine ener_out()
!***********************************************************************
      implicit none
      integer LMX
      parameter (LMX=300)

      real*8 enere,enerc
      integer neex,ncex

!     details of output:
!
!       enere -- exact energy values
!       enerc -- computed energy values
!       enere(1),enerc(1) -- ground state energies
!       enere(>1),enerc(>1) -- excited energies respect to ground
!       neex -- number of exact excited energies
!       ncex -- number of computed excited energies

      integer i,j

      integer se,le,pe,cfe,ne
      integer sc,lc,pc,cfc,nc,gic
      real*8 ee,eegr
      real*8 ec,ecgr
      common/exactbck/se(LMX),le(LMX),pe(LMX),cfe(LMX),ee(LMX),eegr,ne
      common/compbck/sc(LMX),lc(LMX),pc(LMX),cfc(LMX),ec(LMX),ecgr,nc,
     |                  gic(LMX)
      common/eneroutbck/enere(LMX),enerc(LMX),neex,ncex

      neex=ne
      ncex=nc

!.... write exact and computed ground state energies in first row of enert
      enere(1) = eegr
      enerc(1) = ecgr

!.... check config and quantum numbers and write excited energies
      do 100 i=2,ne
         do 110 j=2,nc
            if(se(i).eq.sc(j).and.le(i).eq.lc(j).and.
     |         pe(i).eq.pc(j).and.cfe(i).eq.cfc(j)) then
               enere(i) = ee(i)
               enerc(i) = ec(j)
            endif
110      continue
100   continue

      return
      end
