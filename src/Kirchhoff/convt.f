
c      call step(n,1,dt,p1,p2,p)
c      dt=dt*2.
c      n=n/2

      subroutine step(nfa,nfad,dp,p1,p2,p)

          dimension p1(nfa), p2(nfa), p(nfa)

          nff=nfa/2

	  l=0
          do n=1,nff,nfad
              l=l+1
              p(l)=convs(dp,nfa-1,n-1,p1,p2)
          enddo

          return
          end

      function convs(del,nf,n,fa,fp)

      dimension fa(*), fp(*)

      nn=n
      dn=del
      if(nn.lt.1) go to 2
      ndo=min0(nn,(nf-1)/2)
      ip=2
      np=2*nn
      even=fp(ip)*fa(np)
      odd=0
      if(ndo.lt.2) go to 11
      do 10 i=2,ndo
      ip=ip+1
      np=np-1
      odd=odd+fp(ip)*fa(np)
      ip=ip+1
      np=np-1
      even=even+fp(ip)*fa(np)
  10      continue
 11     continue
      ends=fp(1)*fa(2*nn+1)+fp(ip+1)*fa(np-1)
      convs=dn*(ends+4.*even+2.*odd)/3.
	go to 5
 2      convs=0.
 5	continue
      return
      end 

      subroutine convt1(x,nx,y,ny,z,nz,dt)
      dimension x(*),y(*),z(*),c(160200),d(160200)
      nz=max0(nx,ny)
      call log2fd(nz,n,l2n)
      l2n=l2n+1
      n=2*n
      nx1=nx+1
      ny1=ny+1
      do 1 i=1,nx
      j2=2*i
      j1=j2-1
      c(j1)=x(i)
1     c(j2)=0.
      do 2 i=nx1,n
      j2=2*i
      j1=j2-1
      c(j1)=0.
2     c(j2)=0.
      do 3 i=1,ny
      j2=2*i
      j1=j2-1
      d(j1)=y(i)
3     d(j2)=0.
      do 4 i=ny1,n
      j2=2*i
      j1=j2-1
      d(j1)=0.
4     d(j2)=0.
      call coolb(l2n,c,-1.)
      call coolb(l2n,d,-1.)
      nhalf=n/2
      ncent=nhalf+1
      do 6 i=1,ncent
      j2=2*i
      j1=j2-1
      e1=c(j1)*d(j1)-c(j2)*d(j2)
      e2=c(j1)*d(j2)+c(j2)*d(j1)
      c(j1)=e1
6     c(j2)=e2
      call conj(c,ncent)
      call coolb(l2n,c,1.)
      fscl=dt/float(n)
      do 7 i=1,nz
      j1=2*i-1
7     z(i)=c(j1)*fscl
      return
      end

      subroutine coolb(nn,datai,signi)
      dimension datai(*)
      n=2**(nn+1)
      j=1
      do 5 i=1,n,2
      if(i-j)1,2,2
    1 tempr=datai(j)
      tempi=datai(j+1)
      datai(j)=datai(i)
      datai(j+1)=datai(i+1)
      datai(i)=tempr
      datai(i+1)=tempi
    2 m=n/2
    3 if(j-m)5,5,4
    4 j=j-m
      m=m/2
      if(m-2)5,3,3
    5 j=j+m
      mmax=2
    6 if(mmax-n)7,10,10
    7 istep=2*mmax
      theta=signi*6.28318531/float(mmax)
      sinth=sin(theta/2.)
      wstpr=-2.0  *sinth*sinth
      wstpi= sin(theta)
      wr=1.
      wi=0.
      do 9 m=1,mmax,2
      do 8 i=m,n,istep
      j=i+mmax
      tempr=wr*datai(j)-wi*datai(j+1)
      tempi=wr*datai(j+1)+wi*datai(j)
      datai(j)=datai(i)-tempr
      datai(j+1)=datai(i+1)-tempi
      datai(i)=datai(i)+tempr
    8 datai(i+1)=datai(i+1)+tempi
      tempr=wr
      wr=wr*wstpr-wi*wstpi+wr
    9 wi=wi*wstpr+tempr*wstpi+wi
      mmax=istep
      go to 6
   10 return
      end
