
!***********************************************************************
      subroutine inp_exact()
!***********************************************************************
!
!  Input of NIST values with autostructure's TERM file format
!
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
            print*,'Increase LMX. Stop.'
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
!      print*,ne,eegr
      close(30)
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
      integer ierr
      common/compbck/sc(LMX),lc(LMX),pc(LMX),cfc(LMX),ec(LMX),ecgr,nc,
     |                  gic(LMX)
      common/errlogbck/ierr

      if (ierr.eq.1) goto 999

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
!      print*,nc,ecgr

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
999    return
      end

