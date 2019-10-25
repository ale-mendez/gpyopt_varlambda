!***********************************************************************
      subroutine compare_ei()
!***********************************************************************
!
!     details of output:
!
!       enere -- exact energy values
!       enerc -- computed energy values
!       enere(1),enerc(1) -- ground state energies
!       enere(>1),enerc(>1) -- excited energies respect to ground
!       neex -- number of exact excited energies
!       ncex -- number of computed excited energies
!
      implicit none
      integer MXTN,MXT,LMX
      parameter (MXTN=30,MXT=100,LMX=50)

      integer i,j,iprint
      real*8 diff,errp

      integer ierr
      integer se,le,pe,cfe,ne
      real*8 ee,egre
      integer sc,lc,pc,cfc,gic,nc
      real*8 ec,egrc
      real*8 enere,enerc
      integer neex,ncex

!      common/exactbck/se(LMX),le(LMX),pe(LMX),cfe(LMX),ee(LMX),eegr,ne
!      common/compbck/sc(LMX),lc(LMX),pc(LMX),cfc(LMX),ec(LMX),ecgr,nc,
!     |                  gic(LMX)
      common/errlogbck/ierr
      common/eei_ls/se(MXTN),le(MXTN),pe(MXTN),cfe(MXTN),ee(MXTN),
     |              egre,ne
      common/cei_ls/sc(MXT),lc(MXT),pc(MXT),cfc(MXT),gic(MXT),ec(MXT),
     |              egrc,nc
      common/eicompare/enere(LMX),enerc(LMX),neex,ncex

      neex=ne
      ncex=nc

      if (ierr.eq.1) then
         ncex=ne
         do i=1,ne
             enere(i) = 0.0d0
             enerc(i) = 0.0d0
         enddo
         goto 999
      endif

!... write exact and computed ground state energies in first row of enert
      enere(1) = egre
      enerc(1) = egrc

!... check config and quantum numbers and write absolute energies of excited levels
      do 100 i=2,neex
         do 110 j=2,nc
            if(se(i).eq.sc(j).and.le(i).eq.lc(j).and.
     |         pe(i).eq.pc(j).and.cfe(i).eq.cfc(j)) then
               enere(i) = ee(i)+egre
               enerc(i) = ec(j)+egrc
            endif
110      continue
100   continue

      iprint=0
      if (iprint.eq.1) then
      write(6,1000) 
      do 200 i=1,neex
         diff = enere(i)-enerc(i)
         if(enere(i).eq.0.d0) goto 200
         errp = dabs(diff/enere(i))*100.d0
         write(6,1001) i,enere(i),enerc(i),errp
200   continue
1000  format(3x,'i',6x,'observed',4x,'computed',5x,'erp')
1001  format(1x,i3,1x,3(f12.4))
      endif

      return
999    return
      end

!***********************************************************************
      subroutine print_ener()
!***********************************************************************
      implicit none
      integer LMX
      parameter (LMX=50)

      integer i
      real*8 errp,diff
      real*8 enere,enerc
      integer neex,ncex

      common/eicompare/enere(LMX),enerc(LMX),neex,ncex

      write(25,1001)
      do 100 i=1,neex
         if(enere(i).eq.0.d0) goto 100
         diff = enere(i)-enerc(i)
         errp = dabs(diff/enere(i))*100.d0
         write(25,1000) i,errp,enere(i),enerc(i)
100   continue

1000  format(i4,3(f12.6,2x))
1001  format(3x,"i",6x,"er%",10x,"observed",5x,"computed")
      close(unit=25)
      return
      end

!***********************************************************************
      subroutine compare_aki
!***********************************************************************
      implicit none
      integer MXLTN,MXLT
      parameter (MXLTN=10,MXLT=250)

      integer i,j,iprint
      real*8 errp
      character*3 ires

      integer lwcfe,lwse,lwle,lwpe,lwgie,upcfe,upse,uple,uppe,upgie
      integer ntrane
      real*8 akie,facce
      integer lwcfc,lwsc,lwlc,lwpc,lwgic,upcfc,upsc,uplc,uppc,upgic
      integer ntranc
      real*8 akic
      real*8 vakie,vakic
      integer ntrtot

      common/eaki_ic/lwcfe(MXLTN),lwse(MXLTN),lwle(MXLTN),lwpe(MXLTN),
     |               lwgie(MXLTN),upcfe(MXLTN),upse(MXLTN),uple(MXLTN),
     |               uppe(MXLTN),upgie(MXLTN),akie(MXLTN),facce(MXLTN),
     |               ntrane
      common/caki_ic/lwcfc(MXLT),lwsc(MXLT),lwlc(MXLT),lwpc(MXLT),
     |               lwgic(MXLT),upcfc(MXLT),upsc(MXLT),uplc(MXLT),
     |               uppc(MXLT),upgic(MXLT),akic(MXLT),ntranc
      common/akicompare/vakie(MXLTN),vakic(MXLTN),ntrtot

      ntrtot=0
      do 100 i=1,ntrane
         do 110 j=1,ntranc
            if (lwcfe(i).eq.lwcfc(j).and.upcfe(i).eq.upcfc(j).and.
     |          lwse(i).eq.lwsc(j).and.upse(i).eq.upsc(j).and.
     |          lwle(i).eq.lwlc(j).and.uple(i).eq.uplc(j).and.
     |          lwgie(i).eq.lwgic(j).and.upgie(i).eq.upgic(j)) then
               vakie(i)=akie(i)
               vakic(i)=akic(j)
               ntrtot=ntrtot+1
            endif
110      continue
100   continue

      iprint=0
      if (iprint.eq.1) then
      write(6,1000)
      do 200 i=1,ntrtot
         errp=dabs(vakie(i)-vakic(i))/vakie(i)
         if(errp.gt.facce(i)) ires='out'
         if(errp.le.facce(i)) ires='in'
         write(6,1001) i,vakie(i),vakic(i),errp,facce(i),ires
200   continue
1000  format(3x,'i',6x,'observed',4x,'computed',5x,'erp',7x,'acc')
1001  format(1x,i3,3x,1p,2(e12.3),f10.4,0p,f10.4,a5)
      endif

      return
      end
