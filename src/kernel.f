c     
c     $Header: kernel.f,v 1.2 89/07/04 18:43:07 keith Exp $
c     
c     $Log:	kernel.f,v $
c     Revision 1.2  89/07/04  18:43:07  keith
c     Fixed error in kernel and force which led to sites being allocated the
c     wrong potential parameters.  Needed extra parameter to kernel.
c     
c     Revision 1.1  89/04/20  16:01:50  keith
c     Initial revision
c     
      subroutine kernel(j0, nnab, forceij, pe, rsqr, nabchg, chg,
     +     norm, alpha, ptype, n, potpar)
      implicit real*8 (A-H,O-Z)
      integer      j0, nnab, ptype, n
      real*8 pe, norm, alpha, chg
      real*8 forceij(*), rsqr(*), nabchg(*), potpar(n,*)
      
      integer jsite
      parameter(a1=  0.254829592)
      parameter(a2= -0.284496736)
      parameter(a3=  1.421413741)
      parameter(a4= -1.453152027)
      parameter(a5=  1.061405429)
      parameter(PP=  0.3275911)
      
      POLY(t) = (t*(a1+t*(a2+t*(a3+t*(a4+t*(a5))))))
      
      goto (1000, 2000, 3000) ptype+1
C     Unknown potential type
      write(*,*) ' Unknown potential type ',ptype
      stop
C     Lennard-Jones case
 1000 continue
      if(alpha .ge. 0.0) then
         do 100 jsite = j0+1, nnab
            r     = sqrt(rsqr(jsite))
            rr    = 1.0 / r
            rsqrr = rr**2
            r6r = (potpar(jsite,2)**2 * rsqrr)**3
            r12r  = (r6r)**2
            expa2r2 = nabchg(jsite)* chg * exp(-(alpha*r)**2)
            t = 1.0/(1.0+PP*alpha*r)
            erfcterm = POLY(t) * expa2r2 * rr
            pe = pe + erfcterm + potpar(jsite,1)*(r12r - r6r)
            forceij(jsite) = rsqrr*(6.0*potpar(jsite,1)
     +           *(2*r12r- r6r)
     +           + erfcterm + norm * expa2r2)
 100     continue
      else
         do 101 jsite = j0+1, nnab
            r     = sqrt(rsqr(jsite))
            rr    = 1.0 / r
            rsqrr = rr**2
            r6r = (potpar(jsite,2)**2 * rsqrr)**3
            r12r  = (r6r)**2
            pe = pe + potpar(jsite,1)*(r12r - r6r)
            forceij(jsite) = 6.0*rsqrr*potpar(jsite,1)*
     +           (2*r12r - r6r)
 101     continue
      end if   
      return
C     6-exp potential
 2000 continue
      if(alpha .ge. 0.0) then
         do 200 jsite = j0+1, nnab
            r     = sqrt(rsqr(jsite))
            rr    = 1.0 / r
            rsqrr = rr**2
            r6r   = potpar(jsite,1) * rsqrr**3
            expa2r2 = nabchg(jsite)* chg * exp(-(alpha*r)**2)
            t = 1.0/(1.0+PP*alpha*r)
            erfcterm = POLY(t) * expa2r2 * rr
            expf1 = potpar(jsite,2) * exp(-potpar(jsite,3) * r)
            pe = pe + erfcterm - r6r + expf1
            forceij(jsite) = (-6.0*r6r + erfcterm + norm * expa2r2)
     +           *rsqrr + potpar(jsite,3) *expf1 * rr
 200     continue
      else
         do 201 jsite = j0+1, nnab
            r     = sqrt(rsqr(jsite))
            rr    = 1.0 / r
            rsqrr = rr**2
            r6r   = potpar(jsite,1) * rsqrr**3
            expf1 = potpar(jsite,2) * exp(-potpar(jsite,3) * r)
            pe = pe - r6r + expf1
            forceij(jsite) = -6.0*r6r*rsqrr + potpar(jsite,3)*expf1 * rr
 201     continue
      end if
      return
C     MCY potential
 3000 continue
      if(alpha .gt. 0) then
         do 300 jsite = j0+1, nnab
            r     = sqrt(rsqr(jsite))
            rr    = 1.0 / r
            rsqrr = rr**2
            expa2r2 = nabchg(jsite)* chg * exp(-(alpha*r)**2)
            t = 1.0/(1.0+PP*alpha*r)
            erfcterm = POLY(t) * expa2r2 * rr
            expf1 =  potpar(jsite,1) * exp(-potpar(jsite,2) * r)
            expf2 = -potpar(jsite,3) * exp(-potpar(jsite,4) * r)
            pe = pe + erfcterm + expf1 + expf2
            forceij(jsite) =
     +           (potpar(jsite,2)*expf1 + potpar(jsite,4)*expf2) * rr
     +           + (erfcterm + norm * expa2r2) * rsqrr
 300     continue
      else
         do 301 jsite = j0+1, nnab
            r     = sqrt(rsqr(jsite))
            rr    = 1.0 / r
            rsqrr = rr**2
            expf1 =  potpar(jsite,1) * exp(-potpar(jsite,2) * r)
            expf2 = -potpar(jsite,3) * exp(-potpar(jsite,4) * r)
            pe = pe + expf1 + expf2
            forceij(jsite) =
     +           (potpar(jsite,2)*expf1 + potpar(jsite,4)*expf2) * rr
 301     continue
      end if
      return
      
      end
