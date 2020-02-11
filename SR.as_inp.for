
!***********************************************************************
      subroutine inp_das()
!***********************************************************************
      IMPLICIT REAL*8 (A-H,O-P,R-Z)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT INTEGER*4 (Q)

      integer DMX
      parameter (DMX=100)
      integer NMX
      parameter (NMX=45)

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
      CHARACTER(LEN=4) POTIN,PPOT
      CHARACTER(LEN=6) PRINT,TCC

      real*8 ALFD(3),RCUT(3)

      character*250 das
      integer ndas
      real*8 lam
      integer nvlam,ips,ntlam

      common/dasbck/das(DMX),ndas
      common/salgebbck/MXVORB,MXCONF,KUTSO,RAD,CUP
      common/sminimbck/ALFD,RCUT,NZION,MAXE,NLAM,MCFMX,ORTHOG
      common/lambck/lam(NMX),nvlam,ips,ntlam
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
      ntlam=NLAM
      read(20,*) (linp(i),i=1,ntlam)
      do 200 i=1,ntlam
!        print*,i,linp(i)
        if(linp(i).lt.0) then
            ips=i
            goto 201
        endif
200   continue
201   close(unit=20)
      
      return
      end

!***********************************************************************
      subroutine inp_comp
!***********************************************************************

c...later      open(unit=41,file='ols',status='unknown')
c...later      open(unit=42,file='oic',status='unknown')

      call inp_ols_as()
      call inp_oic_as()

      return
      end

!***********************************************************************
      subroutine inp_ols_as()
!***********************************************************************
      implicit none
      integer MXT,MXTT
      parameter (MXT=150,MXTT=200)

      character line*80,dum*2,dum2*6,aterm*6,aegr*15
      integer i,j,itr,ss,ll,jmin,jmax,gi
      integer nconfigs
      real*8 aak
      integer k,t(MXT)
      integer upcf,upw,upt,lwt,lwcf,lww

      integer lwcfc(MXTT),lwsc(MXTT),lwlc(MXTT),lwpc(MXTT),!lwgic(MXTT),
     |        upcfc(MXTT),upsc(MXTT),uplc(MXTT),uppc(MXTT),!upgic(MXTT),
     |        ntranc
      real*8 akic(MXTT)

      integer ierr
      integer sc,lc,pc,cfc,gic,nc
      real*8 ec,egrc
      common/errlogbck/ierr
      common/cei_ls/sc(MXT),lc(MXT),pc(MXT),cfc(MXT),gic(MXT),ec(MXT),
     |              egrc,nc

      if (ierr.eq.1) goto 999

      open(unit=41,file='ols',status='unknown')

!... read configurations
      read(41,*)
      read(41,'(A)') line
      do 100 i=1,80
         dum=line(i:i+1)
         if(dum.eq.'CF') then
            dum=line(i-2:i-1)
            goto 101
         endif
100   continue
101   read(dum,*) nconfigs
      do 110 i=1,nconfigs+3
         read(41,*)
110   continue

!... read ground energy, energy levels and quantum numbers
      read(41,'(A)') line
      do 200 i=1,80
         dum2=line(i:i+6)
         if(dum2.eq.'NTERM=') then
            aterm=line(i+7:i+12)
         elseif(dum2.eq.'E1/RY=') then
            aegr=line(i+7:i+23)
            goto 201
         endif
200   continue
201   read(aterm,*) nc
      read(aegr,*) egrc
      read(41,*)
      if (nc.ge.MXT) then
         write(*,*) 'Increase MXT. Stop.'
         stop
      endif
      do 210 i=1,nc
         read(41,*) k,t(i),ss,lc(i),cfc(i),ec(i)
         sc(i)=abs(ss)
         if(ss.lt.0) pc(i)=1
         if(ss.gt.0) pc(i)=0
210   continue
      read(41,*)
      read(41,*)

!... compute statistic weigth gi
      do 220 i=1,nc
         ss=(sc(i)-1)/2
         ll=lc(i)
         jmin=abs(ss-ll)
         jmax=ss+ll
         gi=0
         do 230 j=jmin,jmax,1
            gi=gi+2*j+1
230      continue
         gic(i)=gi
220   continue

!... read computed transition data
      itr=0
      do 300 i=1,MXTT
         if(i.gt.MXTT) then
            write(*,*) 'Increase MXTT. Stop.'
            stop
         endif
         read(41,*,end=301) upcf,upt,upw,lwcf,lwt,lww,aak
         if(aak.eq.0.d0) goto 300
         itr=itr+1
         akic(itr)=aak
         do 310 j=1,nc
            if(upt.eq.t(j)) then
               upcfc(itr)=cfc(j)
               upsc(itr)=sc(j)
               uplc(itr)=lc(j)
               uppc(itr)=pc(j)
            endif
            if(lwt.eq.t(j)) then
               lwcfc(itr)=cfc(j)
               lwsc(itr)=sc(j)
               lwlc(itr)=lc(j)
               lwpc(itr)=pc(j)
            endif
310      continue
300   continue
301   ntranc=itr

      close(unit=41)
      return
999    return
      end

!***********************************************************************
      subroutine inp_oic_as()
!***********************************************************************
      implicit none
      integer MXL,MXLT
      parameter(MXL=200,MXLT=3500)

      character line*80,dum*2,dum2*7,alev*6,aegr*15
      integer i,j,ss,jj,itr
      integer nconfigs,nlevels
      real*8 egr,aak
      integer k,lv(MXL),tt,is(MXL),il(MXL),igi(MXL),ip(MXL),icf(MXL)
      real*8 dek(MXL)
      integer upcf,upw,uplv,lwlv,lwcf,lww

      integer ierr
      integer lwcfc,lwsc,lwlc,lwpc,lwgic,upcfc,upsc,uplc,uppc,upgic
      integer ntranc
      real*8 akic
      common/errlogbck/ierr
      common/caki_ic/lwcfc(MXLT),lwsc(MXLT),lwlc(MXLT),lwpc(MXLT),
     |               lwgic(MXLT),upcfc(MXLT),upsc(MXLT),uplc(MXLT),
     |               uppc(MXLT),upgic(MXLT),akic(MXLT),ntranc

      if (ierr.eq.1) goto 999

      open(unit=42,file='oic',status='unknown')

c... read configurations
      read(42,*)
      read(42,'(A)') line
      do 100 i=1,80
         dum=line(i:i+1)
         if(dum.eq.'CF') then
            dum=line(i-2:i-1)
            goto 101
         endif
100   continue
101   read(dum,*) nconfigs
      do 110 i=1,nconfigs+3
         read(42,*)
110   continue

c... read ground energy, energy levels and quantum numbers
      read(42,'(A)') line
      do 200 i=1,80
         dum2=line(i:i+6)
         if(dum2.eq.'NLEVEL=') then
            alev=line(i+7:i+12)
         elseif(dum2.eq.'E1/RY= ') then
            aegr=line(i+7:i+23)
            goto 201
         endif
200   continue
201   read(alev,*) nlevels
      read(aegr,*) egr
      read(42,*)
      if (nlevels.ge.MXL) then
         write(*,*) 'Increase MXL. Stop.'
         stop
      endif
      do 210 i=1,nlevels
         read(42,*) k,lv(i),tt,ss,il(i),jj,icf(i),dek(i)
         is(i)=abs(ss)
         igi(i)=jj+1
         if(ss.lt.0) ip(i)=1
         if(ss.gt.0) ip(i)=0
210   continue
      read(42,*)
      read(42,*)

c... read computed transition data
      itr=0
      do 300 i=1,MXLT
         if(i.gt.MXLT) then
            write(*,*) 'Increase MXLT. Stop.'
            stop
         endif
         read(42,*,end=301) upcf,uplv,upw,lwcf,lwlv,lww,aak
         itr=itr+1
         akic(itr)=dabs(aak)
         do 310 j=1,nlevels
            if(uplv.eq.lv(j)) then
               upcfc(itr)=icf(j)
               upsc(itr)=is(j)
               uplc(itr)=il(j)
               uppc(itr)=ip(j)
               upgic(itr)=igi(j)
            endif
            if(lwlv.eq.lv(j)) then
               lwcfc(itr)=icf(j)
               lwsc(itr)=is(j)
               lwlc(itr)=il(j)
               lwpc(itr)=ip(j)
               lwgic(itr)=igi(j)
            endif
310      continue
!         if (lwcfc(itr).eq.1) then
!         write(*,*) itr,akic(itr),
!     |             lwcfc(itr),lwsc(itr),lwlc(itr),lwpc(itr),lwgic(itr),
!     |             upcfc(itr),upsc(itr),uplc(itr),uppc(itr),upgic(itr)
!         endif
300   continue
301   ntranc=itr
!      print*,"read ",ntranc," transitions from AS"

      close(unit=42)
      return
999    return
      end

