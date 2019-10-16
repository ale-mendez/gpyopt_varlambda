!***********************************************************************
      subroutine ener_out()
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

      integer i,j

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
      common/eneroutbck/enere(LMX),enerc(LMX),neex,ncex

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
      do 100 i=2,ne
         do 110 j=2,nc
            if(se(i).eq.sc(j).and.le(i).eq.lc(j).and.
     |         pe(i).eq.pc(j).and.cfe(i).eq.cfc(j)) then
               enere(i) = ee(i)+egre
               enerc(i) = ec(j)+egrc
            endif
110      continue
100   continue

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

      common/eneroutbck/enere(LMX),enerc(LMX),neex,ncex

      write(25,1001)
      do 100 i=1,neex
         if(enere(i).eq.0) goto 100
         diff = enere(i)-enerc(i)
         errp = dabs(diff/enere(i))*100.d0
         write(25,1000) i,errp,enere(i),enerc(i)
100   continue

1000  format(i4,3(f12.6,2x))
1001  format(3x,"i",6x,"er%",10x,"exact",8x,"computed")
      close(unit=25)
      return
      end
