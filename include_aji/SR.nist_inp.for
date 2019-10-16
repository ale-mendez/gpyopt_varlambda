      subroutine inp_exper

      open(unit=31,file='NIST_cfgs.dat',status='unknown')
      open(unit=32,file='NIST_terms.dat',status='unknown')
      open(unit=33,file='NIST_energies.dat',status='unknown')
      open(unit=34,file='NIST_lines.dat',status='unknown')

      call inp_cfgs_nist
      call inp_terms_nist
      call inp_lines_nist

      return
      end

!***********************************************************************
      subroutine inp_cfgs_nist
!***********************************************************************
!
!  Input of configuration list with NIST format
!
      implicit none
      integer MXCF
      parameter (MXCF=18)

      integer i,ic

      integer icf,ncf
      character*10 cf
      common/cfbck/icf(MXCF),cf(MXCF),ncf

!... read header
      do 100 i=1,3
         read(31,*)
100   continue

      ic=0
!... read configurations
      do 200 i=1,MXCF
         if(i.gt.MXCF) then
            write(*,*) 'Increase MXCFN. Stop.'
            stop
         endif
         read(31,*,end=201) icf(i),cf(i)
!         print*,icf(i),cf(i)
         ic=ic+1
200   continue
201   ncf=ic
!      print*,ncf

      close(unit=31)
      return
      end

!***********************************************************************
      subroutine inp_terms_nist
!***********************************************************************
!
!  Input of NIST terms with NIST format
!
      implicit none
      integer MXCF,MXTN
      parameter (MXCF=18,MXTN=30)

      integer i,j,iss,ill,ipp,it
      character icfterms*10,islterms*3
      real*8 enernist

      character*200 line1,line2
      character dum*1,abin*30,aion*30

      integer icf,ncf
      character*10 cf
      integer cfe,se,le,pe,ne
      real*8 ee,egre,eion
      common/cfbck/icf(MXCF),cf(MXCF),ncf
      common/eei_ls/se(MXTN),le(MXTN),pe(MXTN),cfe(MXTN),ee(MXTN),
     |              egre,ne

!... read headers
      do 100 i=1,3
         read(32,*)
         read(33,*)
100   continue

!... read list of terms
      it=0
      do 200 i=1,MXTN
         if(i.gt.MXTN) then
            write(*,*) 'Increase MXTN. Stop.'
            stop
         endif
         read(32,*,end=201) icfterms,islterms,enernist
         do 210 j=1,ncf
           if(icfterms.eq.cf(j)) then
              call SLP(islterms,iss,ill,ipp)
              it=it+1
              cfe(it)=icf(j)
              se(it)=iss
              le(it)=ill
              pe(it)=ipp
              ee(it)=enernist
           endif
210      continue
!         print*,icfterms,islterms,enernist
!         print*,cfe(it),se(it),le(it),pe(it),ee(it)
200   continue
201   ne=it

c... read binding and ionization energies 
      read(33,'(A)') line1
      read(33,'(A)') line2
      do 300 i=1,200
         dum=line1(i:i)
         if(dum.eq.'>') then
            abin=line1(i+1:i+31)
            aion=line2(i+1:i+31)
         endif
300   continue
      read(abin,*) egre
      read(aion,*) eion
!      print*,'binener=',egre
!      print*,'ionener=',eion

      close(unit=32)
      close(unit=33)
      return
      end

!***********************************************************************
      subroutine inp_lines_nist
!***********************************************************************
!
!  Input of NIST lines with NIST format
!
      implicit none
      integer MXCF,MXTN,MXLTN
      parameter (MXCF=18,MXTN=30,MXLTN=10)

      integer i,j,itr
      real*8 aji,ei,ek,fvalue
      character acc*3,dum*1,lcf*10,ucf*10,lt*3,ut*3
      integer lj,uj,ntrane
      integer iss,ill,ipp

      integer icf,ncf
      character*10 cf
      integer cfe,se,le,pe,ne
      real*8 ee,egre
      integer lwcfe,lwse,lwle,lwpe,lwgie,upcfe,upse,uple,uppe,upgie
      real*8 akie,facce
      common/cfbck/icf(MXCF),cf(MXCF),ncf
      common/eei_ls/se(MXTN),le(MXTN),pe(MXTN),cfe(MXTN),ee(MXTN),
     |              egre,ne
      common/eaki_ic/lwcfe(MXLTN),lwse(MXLTN),lwle(MXLTN),lwpe(MXLTN),
     |               lwgie(MXLTN),upcfe(MXLTN),upse(MXLTN),uple(MXLTN),
     |               uppe(MXLTN),upgie(MXLTN),akie(MXLTN),facce(MXLTN),
     |               ntrane

!... read header
      do 100 i=1,5
         read(34,*)
100   continue

      itr=0
      do 200 i=1,MXLTN
         if(i.gt.MXLTN) then
            write(*,*) 'Increase MXLTN. Stop.'
            stop
         endif
         read(34,*,end=201) aji,acc,ei,dum,ek,lcf,lt,lj,ucf,ut,uj
!         print*,aki,acc,lcf,lt,lj,ucf,ut,uj
         itr=itr+1
         akie(itr)=aji
         lwgie(itr)=2*lj+1
         upgie(itr)=2*uj+1
         call accuracy_nist(acc,fvalue)
         facce=fvalue
         do 210 j=1,ncf
           if(lcf.eq.cf(j)) then
              call SLP(lt,iss,ill,ipp)
              lwcfe(itr)=icf(j)
              lwse(itr)=iss
              lwle(itr)=ill
              lwpe(itr)=ipp
           endif
           if(ucf.eq.cf(j)) then
              call SLP(ut,iss,ill,ipp)
              upcfe(itr)=icf(j)
              upse(itr)=iss
              uple(itr)=ill
              uppe(itr)=ipp
           endif
210      continue
!         print*,lcf,lt,lwcf(i),lws(i),lwl(i),lwp(i),lww(i)
!         print*,ucf,ut,upcf(i),ups(i),upl(i),upp(i),upw(i)
200   continue
201   ntrane=itr

      close(unit=34)
      return
      end

!***********************************************************************
      subroutine SLP(islterms,iss,ill,ipp)
!***********************************************************************
      implicit none
      character islterms*3
      character*1 as,al,ap
      integer iss,ill,ipp

      as=islterms(1:1)
      al=islterms(2:2)
      ap=islterms(3:3)
      read(as,*) iss
      if(al.eq.'S') ill=0
      if(al.eq.'P') ill=1
      if(al.eq.'D') ill=2
      if(al.eq.'F') ill=3
      if(al.eq.'G') ill=4
      if(al.eq.'H') ill=5
      if(ap.eq.' ') ipp=0
      if(ap.eq.'*') ipp=1

      return
      end

!***********************************************************************
      subroutine accuracy_nist(acc,fvalue)
!***********************************************************************
      implicit none
      character*3 acc
      real*8 fvalue

!... accuracy is best or equal than fvalue%
      if(acc.eq.'AAA') fvalue=0.3d0
      if(acc.eq.'AA') fvalue=1.d0
      if(acc.eq.'A+') fvalue=2.d0
      if(acc.eq.'A') fvalue=3.d0
      if(acc.eq.'B+') fvalue=7.d0
      if(acc.eq.'B') fvalue=10.d0
      if(acc.eq.'C+') fvalue=18.d0
      if(acc.eq.'C') fvalue=25.d0
      if(acc.eq.'D+') fvalue=40.d0
      if(acc.eq.'D') fvalue=50.d0
!... accuracy is worst or equal than fvalue%
      if(acc.eq.'E') fvalue=50.d0

      return
      end
