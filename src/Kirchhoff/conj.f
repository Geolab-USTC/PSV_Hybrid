      subroutine conj(c,ncent)
      dimension c(*)
      nhm1=ncent-2
      ncent2=ncent*2
      do 1 i=1,nhm1
      l2=2*i
      l1=l2-1
      l3=l2+1
      k1=ncent2+l1
      k2=ncent2+l2
      j1=ncent2-l3
      j2=ncent2-l2
      c(k1)=c(j1)
1     c(k2)=-c(j2)
      return
      end
      subroutine log2fd(np,n,l2n)
      n1=0
      n=1
      l2n=0
1     continue
      if((np.gt.n1).and.(np.le.n)) go to 2
      n1=n1*2
      n=n*2
      l2n=l2n+1
      go to 1
2     return
      end
