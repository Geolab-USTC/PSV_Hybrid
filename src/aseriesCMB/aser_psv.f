      subroutine GRT(isourcec,iflat,ispc,mtdc,nbc,jo0,
     *      cc0,ss0,dd0,tth0,
     *      lfinal0,nen0,na0,nray0,ncoun0,linec,iwhc,
     *      dpc,nnc,xxx,tsc,nstressc,mc,so,greenr,greenz)
c   isourcec The source type =0: isotropic source
c                            =1: dislocation source
c   iflat = 0       flat layer 1: flatting
c   nbc      The layer where the source is
c   jo0
c   cc0
c   ss0
c   dd0
c   tth0
c   dpc      The time interval
c   nnc      The number of time steps
c   xxx      The epicentral distance
c   tsc      The starting time
c   so       The radiation pattern
c   nstressc The control of output Greens function
c            For dislocation source:
c               =0   w    velocity
c               =12  Txz  stress
c               =22  Tzz  stress
c               =11  Txx  stress
c            For isotropic source:
c               <=1  The Green's function
c               > 1  The z-deravitive of Green's function
c   greenr   The output of the r- Green's function
c   greenz   The output of the z- Green's function

      include 'psv.p'

      integer nbc,jo0,lfinal0,nnc,nstressc,isourcec
      integer ispc, mtdc
      integer linec
      real dpc,xxx,tsc
      real mc
      real so(*), greenr(*)
      real greenz(*)
      real cc0(*),ss0(*),dd0(*),tth0(*)
      integer nen0(*),na0(*),nray0(*),ncoun0(*)

      common/trace/tracing
      common/down/ndown, nstress
      common/pathc/po,to,k
      common/flt/phh(iflt,inz),ns,beta,nz,nup,nst
      common/str/ph(itt,inz)
      common/orst/cc(icc),ss(icc),dd(icc),tth(icc),xx
      common/rays/na(ina1,ina2),nray(ina1,ina2),nend,lmax,nsp,love,
     *nzr,iwh
      common/tfix/tn1,tn2,tn3,tn4,jn1,jn2,jn3,jn4,det,ncont
      common/term/iordr
      common/coeff/itrnm
      common/travel/rlp(icc),rls(icc),ll,fdp,map,nplnw
      common/thy/tt(itt),pp(itt),ff(itt),l
      common/plotc/con
      common/stuff/c(icc),s(icc),d(icc),th(icc),rcsq(icc),rssq(icc),x
      common/cagcon/tstart,dp,nn,ndirt,nnend,dltp,nnn,dtim,mtd,dltm,nfm
     +,mkp
      common/soray/nsorce,ntype,ndiret,tqc,tqs,tqd,jpt,nconjt,dconjt,nb
      common/dscfil/far1
      common/tmpwr/iwr
      common/sourcetype/isource
      common/scase/ispecial
      common/linepoint/line

      dimension nen(ina2),ncoun(ina2)
      real moment

      common/lprint/prnt,prnts,prntc,flat,prcon,prampl
      logical prnt,prnts,prntc,flat,prcon,prampl

c DLTM and mtd used in the contor
c dltm max seperation between time points
c mtd=3 controls rate of approaching this limit
c Dltp the last time point before TO must be less than
c Dtim the seperation between TO and the first time point after TO
c muse be less than dtim
      iwr=0
      nfm=0

      prnt=.false.
      prnts=.false.
      prntc=.false.
      prcon=.false.
      prampl=.false.
      tracing = -1

c        write(*,*)"nnc=",nnc
      do j=1,nnc
        greenr(j)  = 0.0
        greenz(j) = 0.0
      enddo

c    isource=1 dislocation source
c           =0 explosition source
      isource=isourcec
      nplnw=2
      itrnm=2
      iordr=1
      ndirt=2

      iwh=iwhc
      line = linec
      nb =nbc
      dp =dpc
      nstress=nstressc
      tstart=tsc
      jo=jo0
      dt=dp
      ispecial = ispc
      moment = mc

      nn=nnc
      depth=0.0
      do jj=1,nbc
        depth=depth+tth0(jj)
      enddo
      sdepth=depth
      rhos=dd0(nbc)

      do j=1,jo0
        cc(j)=cc0(j)
        ss(j)=ss0(j)
        dd(j)=dd0(j)
        tth(j)=tth0(j)
      enddo

      ns=nb
      ncont = 1

      flat = .true.
      if(iflat .eq. 0)flat=.false.

      lfinal=lfinal0
      kk=1
      do  nsp=1,lfinal0
        nen(nsp)=nen0(nsp)
        ncoun(nsp)=1
        do j=1,nen(nsp)
          na(j,nsp)=na0(kk)
          nray(j,nsp)=nray0(kk)
          kk=kk+1
        enddo
      enddo

      do i=1,jo
        cc(i)=cc0(i)
        ss(i)=ss0(i)
        dd(i)=dd0(i)
        tth(i)=tth0(i)
      enddo

      if(tracing .gt. 0) then
        write(*,*)'nb=',nb,' dp=',dp,' nstress=',nstress
        write(*,*)'tstart=',tstart,' jo=',jo,' nn=',nnc
        write(*,*)"jo=",jo
        do j=1,jo
          write(*,*)cc(j), ss(j), dd(j), tth(j)
        enddo
        write(*,*)"lfinal=",lfinal
        do nsp=1,lfinal
          write(*,*)nen(nsp),(na(j,nsp),j=1,nen(nsp))
          write(*,*)ncoun(nsp), (nray(j,nsp),j=1,nen(nsp))
        enddo
      endif

      tracing = -1
      jpt = 0
      nz=8
      nconjt = 0
      lmax = jo+1
      call curay(jo)
      ll=lmax
      dtim = 0.001
      dtim = dp/2.
      dltp = dtim
c      dtim = dp/5.
      dltm =5.0
c Changed temporatorly
      mtd=mtdc
      nnn=15
      ns=nz
      beta=s(nb)

      x=xxx
      x=abs(x)
      far1=0.0
      xx=x
c      con=1.
c      ra = (xx/111.195) * 0.01745
c      con=sqrt(ra/sin(ra))*sqrt(2./xx)*moment
c      if(line .eq. 1) con=sqrt(ra/sin(ra))*moment
      con=sqrt(2./xx)*moment
      if(line .eq. 1) con=moment
      hs =0.
      l=1
      nnsp=1
      do  nsp=nnsp,lfinal
        ndown=ncoun0(nsp)
c comments, the nrec is the layer where the receivers are,
c it needs caution, when an interface is present.
        nrec=na(nen(nsp)-1,nsp)
c          /* parameters at receiver */
        tqc=cc(nrec)
        tqs=ss(nrec)
        tqd=dd(nrec)
c          write(*,*)tqc,tqs,tqd

        love=0
        if(nray(1,nsp).eq.4) love=2

        nst=nray(1,nsp)

        nsorce=nray(1,nsp)
        if(nsorce.eq.4) nsorce=3

        nup=0
        if(na(1,nsp).gt.nb) nup=2

        ndiret=0
        ncount=ncoun(nsp)
        nend=nen(nsp)-1
c       nend=nen(nsp)

        nx=1
        do nm=1,nend
           nx=max0(nx,na(nm,nsp))
        enddo

        lmax=nx
        tq=0
        call trav(ncount)
        if (love.eq.2) nray(1,nsp)=4
      enddo

      if(love .eq. 2)then
c  jz = 1 vss (if isource == 0, v)
c  jz = 2 vds (if isource == 0, dv/dz)
        do jk=1,nn
          pp(jk)=so(5)*ph(jk,1)
        enddo
      else
          jz=5
c     write(*,*)"jz=",jz
          do jk=1,nn
              pp(jk)=so(1)*ph(jk,jz)+so(2)*ph(jk,jz-1)+so(3)*ph(jk,jz-2)
          enddo
      endif
      if(line .ne. 1)then
          con = 2./3.1415926
          nnz = nn
          call step(nnz,1,0,2,0,dp,2)
          dpc  = 2.0*dp
          nnc  = nn/2
      endif

      do jk=1,nnc
          greenr(jk)=pp(jk)
      enddo
      ymax=0.0
      do jk=1,nn
          if(abs(pp(jk)).gt.ymax)ymax=abs(pp(jk))
      enddo
c      write(*,*)'ymax=',ymax
      if(love .ne. 2)then
          jz=8
          do jk=1,nn
              pp(jk)=so(1)*ph(jk,jz)+so(2)*ph(jk,jz-1)+so(3)*ph(jk,jz-2)
          enddo
      else
          do jk=1,nn
              pp(jk)=so(4)*ph(jk,2)
          enddo
      endif
      ymax=0.0
      do jk=1,nn
         if(abs(pp(jk)).gt.ymax)ymax=abs(pp(jk))
      enddo
c       write(*,*)'ymax=',ymax
      if(line .ne. 1)then
        call step(nnz,1,0,2,0,dp,2)
      endif
      if(isource .eq. 0)then
c          /* this is d/dz */
        call diff(nn,pp,dp)
      endif
      do jk=1,nnc
         greenz(jk)=pp(jk)
      enddo
      return
      end
