      program test
C this is just perform a simpe test
C omit main program from implementation in c(++) code  
      implicit none 
      real*8 ecms,x1,x2
      real*8 pjet1(0:3),pjet2(0:3),phiggs(0:3)
      integer iopthc,i,j,flavour1,flavour2
      real*8 rsmi,di,dbi,dti,dtbi,lambdai
      real*8 a1hwwi,a2hwwi,a3hwwi,a1haai,a2haai,a3haai
      real*8 a1hazi,a2hazi,a3hazi,a1hzzi,a2hzzi,a3hzzi
      real*8 oo,oo2,weight
C
      iopthc = 1
      rsmi    = 1.0d0
      di      = 0.0d0
      dbi     = 0.0d0
      dti     = 0.0d0
      dtbi    = 0.0d0
      lambdai = 0.0d0
C
      ecms = 8000d0
      x1 = 1.0d0/8.0d0   !   500 GeV
      x2 = 1.0d0/16.0d0   ! -250 GeV
      pjet1(3)   =   1400d0
      pjet2(3)   =  -1200d0
      phiggs(3)  =     50d0 
      pjet1(1)   =     30d0
      pjet2(1)   =    -50d0    
      phiggs(1)  =     20d0 
      pjet1(2)   =     30d0
      pjet2(2)   =     20d0
      phiggs(2)  =    -50d0 
      pjet1(0) = sqrt(pjet1(1)*pjet1(1)+pjet1(2)*pjet1(2)
     +              +pjet1(3)*pjet1(3))
      pjet2(0) = sqrt(pjet2(1)*pjet2(1)+pjet2(2)*pjet2(2)
     +               +pjet2(3)*pjet2(3))
      phiggs(0) = sqrt(phiggs(1)*phiggs(1)+phiggs(2)*phiggs(2)
     +                +phiggs(3)*phiggs(3)+120d0*120d0)
      flavour1 = 1
      flavour2 = -1

      do i =0,3
         write(*,*)'i ',i,pjet1(i)+pjet2(i)+phiggs(i)
      enddo
C
      call optobs(ecms,x1,x2,pjet1,pjet2,phiggs,
C      call optobs(ecms,x1,x2,pjet1,flavour1,pjet2,flavour2,phiggs,
C     +     iopthc,rsmi,di,dbi,dti,dtbi,lambdai,
C     +     a1hwwi,a2hwwi,a3hwwi,a1haai,a2haai,a3haai,
C     +     a1hazi,a2hazi,a3hazi,a1hzzi,a2hzzi,a3hzzi,
     +     oo,oo2)
C    +      ,weight) 
C
      write (*,*)' optimal observble ',oo
      write (*,*)' reweighting weight ',weight
      return
      end
********************************************************************************
      subroutine optobs(ecms,x1,x2,pjet1,pjet2,phiggs,oo,oo2) 
C      subroutine optobs(ecms,x1,x2,pjet1,flav1,pjet2,flav2,phiggs,oo) 
********************************************************************************
C      input: ecms proton-pton center-of mass energy in GeV 
C             x1, x2: Bjorken x of incoming partons
C             pjet1(0:3):  E,px,py,pz of 1st jet
C             pjet1(0:3):  E,px,py,pz of 2nd jet
C             phiggs(0:3): E,px,py,pz of Higgs boson
C      -----------------------------------------------------------------------------
C      output: oo - optimal observable
C              oo2 - optimal observable 2nd order (testing)
*********************************************************************************
      implicit none
C input and output variables
      real*8 ecms,m2sm,m2full,m2cpo,x1,x2
      real*8, intent(out) :: oo
      !f2py intent(out) :: oo
      real*8, intent(out) :: oo2
      !f2py intent(out) :: oo2
      real*8 pjet1(0:3),pjet2(0:3),phiggs(0:3)
      integer i1,i2,i,j,flav1,flav2
C
      real*8 xp1,xp2,wmat2
      real*8 p(6,0:3),pdf1(-6:6),pdf2(-6:6)
      real*8 m2i(-5:5,-5:5),m2if(-5:5,-5:5,-5:5,-5:5)
C
      integer qborn,qw,qz,qschan,qtchan,qch2,qchint,
     &     qbini,qbfin,qwidth,qfact,qbos,qferm,qsoft,qhh2,
     &     qqcddiag,qqcdnondiag,qqcdgsplit,qqcdggfus,qcp
      common/rcoptions/qborn,qw,qz,qschan,qtchan,qch2,qchint,
     &     qbini,qbfin,qwidth,qfact,qbos,qferm,qsoft,qhh2,
     &     qqcddiag,qqcdnondiag,qqcdgsplit,qqcdggfus,qcp
C
      real*8 pi,el,alpha,alpha0,alphaz,GF,alphas
      real*8 cw,cw2,sw,sw2,mw,mw2,gw,mz,mz2,gz,mh,mh2
      real*8 ml,ml2,mqp,mqp2,mqm,mqm2

      complex*16 v,cv
      common/param/pi,el,alpha,alpha0,alphaz,GF,alphas,
     &     v(3,3),cv(3,3),cw,cw2,sw,sw2,mw,mw2,gw,mz,mz2,gz,mh,mh2,
     &     ml(3),ml2(3),mqp(3),mqp2(3),mqm(3),mqm2(3)
C     
      real*8 sinthetac,mfact
      real*8 ppecms
      common/ppecms/ppecms
C
      complex*16 ccw,ccw2,csw,csw2,cmw,cmw2,cmz,cmz2
      common/cparam/ccw,ccw2,csw,csw2,cmw,cmw2,cmz,cmz2
C
      complex*16 xcw,xcw2,xsw,xsw2,xmw,xmw2,xmz,xmz2,xmh,xmh2,null
      complex*16 xml,xml2,xmqp,xmqp2,xmqm,xmqm2
C
      common/xparam/xcw,xcw2,xsw,xsw2,xmw,xmw2,xmz,xmz2,xmh,xmh2,null, 
     &     xml(3),xml2(3),xmqp(3),xmqp2(3),xmqm(3),xmqm2(3) 
C
      real*8 qu,qd,ql,qn,qf,mu,mu2,md,md2,mlep,mlep2
      complex*16 guu,gdd,gnn,gll     
      common/qf/qu,qd,ql,qn,qf(4),mu,mu2,md,md2,mlep,mlep2,
     +     guu(-1:1),gdd(-1:1),gnn(-1:1),gll(-1:1)
C
      integer qhvv
      real*8 rsm,d,db,dt,dtb,lambdahvv
      real*8 a1hww,a2hww,a3hww,a1haa,a2haa,a3haa
      real*8 a1haz,a2haz,a3haz,a1hzz,a2hzz,a3hzz

      common/hvv/rsm,d,db,dt,dtb,lambdahvv,
     &     a1hww,a2hww,a3hww,a1haa,a2haa,a3haa,
     &     a1haz,a2haz,a3haz,a1hzz,a2hzz,a3hzz,qhvv
C     
      ppecms = ecms             !  pp-center of mass energy 
c
C P(1,x)  = 1st incoming proton P(2,x) = 2nd incoming proton 
C P(3,x)  = 1st outgoing jet    P(3,x) = 2nd outgoing jet 
C P(5,x)  = Higgs boosn candidate       
C
      xp1=x1                    !  Bjorken x1
      xp2=x2                    !  Bjorken x2
      p(1,0) = xp1*ppecms/2.0d0
      p(1,1) = 0.0d0
      p(1,2) = 0.0d0
      p(1,3) = xp1*ppecms/2.0d0
      p(2,0) = xp2*ppecms/2.0d0
      p(2,1) = 0.0d0
      p(2,2) = 0.0d0
      p(2,3) = -xp2*ppecms/2.0d0
      do I = 0,3
         p(3,i) = pjet1(i)
         p(4,i) = pjet2(i)
         p(5,i) = phiggs(i)
         p(6,i) = 0.0d0
      enddo
C set HAWK "constants" (ony those needed by LO calculation)
      alpha = 1.0d0/137.0d0
      pi = 4d0*datan(1d0)
      el = sqrt(4d0*pi*alpha)   ! alpha, alpha0 oder alphaz
C
      mz   = 91.1876d0
      cmz = mz
      xmz = mz
      mz2 = mz*mz
      xmz2 = mz2 
      cmz2  = mz2
C
      mw   = 80.425d0
      xmw  = mw
      cmw  = mw
      mw2  = mw*mw
      xmw2 = mw2
      cmw2 = mw2
C
      mh   = 125d0
C
      cw  = mw/mz
      xcw = cw
      ccw = cw 
      cw2 = cw*cw
      xcw2 = cw2
      ccw2 = cw2
      sw2  = 1d0-cw2
      csw2 = sw2
      xsw2 = sw2
      sw  = sqrt(sw2)
      xsw  = sw 
      csw  = sw 
C fermion el-mag charges
      qu =  2d0/3d0
      qd =  -1d0/3d0
      qn =   0d0
      ql =   -1d0
C fermion couplings 
      guu(+1) = -qu*xsw/xcw
      guu(-1) = -qu*xsw/xcw+1d0/2d0/xsw/xcw
      gdd(+1) = -qd*xsw/xcw
      gdd(-1) = -qd*xsw/xcw-1d0/2d0/xsw/xcw
      gnn(+1) = -qn*xsw/xcw
      gnn(-1) = -qn*xsw/xcw+1d0/2d0/xsw/xcw
      gll(+1) = -ql*xsw/xcw
      gll(-1) = -ql*xsw/xcw-1d0/2d0/xsw/xcw
C CKM matrix and its conjugated
      sinthetac = 0.226548d0
      v(1,2) = sinthetac
      v(1,1) = sqrt(1d0-v(1,2)*v(1,2)) ! bug v(1,1)
      v(2,1) = -v(1,2)
      v(2,2) = v(1,1)
      v(1,3) = 0d0
      v(2,3) = 0d0
      v(3,1) = 0d0
      v(3,2) = 0d0
      v(3,3) = 1d0
      do I =1,3
         do J = 1,3
            cv(i,j) = dconjg(v(i,j))
         enddo
      enddo
C set HAWK flags (only those needed by LO calculation)
      qw     = 1 
      qz     = 1 
C
      qschan = 0 
      qtchan = 1 
      qch2   = 1
      qchint = 0 
C     
      qbini       = 0 
      qbfin       = 0 
      qwidth      = 0 
      qcp         = 0 
C anomalous couplings
C dt = 1,0 rest=0
      a1hww = mw/sw
      a2hww = 0.0d0
      a3hww = 2d0/sw/mw
      a1haa = 0d0
      a2haa = 0d0
      a3haa = 0d0
      a1haz = 0d0
      a2haz = 0d0
      a3haz = 0d0
      a1hzz = mw/sw/cw2
      a2hzz = 0d0
      a3hzz = 4d0*cw2/2d0/sw/mw
C otherwise use direct in put value in terms if a1 a2 a3
c PDFs
      mfact = mh
      call getpdf(xp1,mfact,pdf1)
      call getpdf(xp2,mfact,pdf2)
C for a test withour pdf weighting
C      do i = -6,6
C         pdf1(i) = 1.0d0
C         pdf2(i) = 1.0d0
C      enddo
C     SM:
      qhvv = 0 
c squared matrix elements
      call Mat2born0(p,m2i,m2if,1)
c weighted |M|^2
      wmat2 = 0d0
      do 100 i1=-4,4
         do 100 i2=-4,4
            if (i1*i2.eq.0) goto 100
            wmat2 = wmat2 + m2i(i1,i2)*pdf1(i1)*pdf2(i2) 
 100          continue
      m2sm = wmat2
C      m2sm = m2i(flav1,flav2)
C     full
      qhvv = 1 
C
c squared matrix elements
      call Mat2born0(p,m2i,m2if,1)
c weighted |M|^2
      wmat2 = 0d0
      do 200 i1=-4,4
         do 200 i2=-4,4
            if (i1*i2.eq.0) goto 200
            wmat2 = wmat2 + m2i(i1,i2)*pdf1(i1)*pdf2(i2) 
 200          continue
      m2full = wmat2
C      m2full = m2i(flav1,flav2)
C pure anom. CP-odd, ohne Formfaktor:
      qhvv  = 1
      a1hww = 0d0
      a1hzz = 0D0
c squared matrix elements
      call Mat2born0(p,m2i,m2if,1)
c weighted |M|^2
      wmat2 = 0d0
      do 300 i1=-4,4
         do 300 i2=-4,4
            if (i1*i2.eq.0) goto 300
            wmat2 = wmat2 + m2i(i1,i2)*pdf1(i1)*pdf2(i2) 
 300          continue
      m2cpo = wmat2
C      m2cpo = m2i(flav1,flav2)
C     
      oo = (m2full-m2sm-m2cpo)/m2sm/2.0d0
      oo2 = m2cpo/m2sm
C      write(*,*)'lambdahvv ',lambdahvv
      return
      end

**********************************************************************************
      subroutine reweight(ecms,ipara,rsmin,din,dbin,dtin,dtbin,
     +     lambdahvvin,
     +     a1hwwin,a2hwwin,a3hwwin,a1haain,a2haain,a3haain,
     +     a1hazin,a2hazin,a3hazin,a1hzzin,a2hzzin,a3hzzin,
     +     x1,x2,pjet1,pjet2,phiggs,weight) 
**********************************************************************************
C      input: ecms proton-pton center-of mass energy in GeV 
C             x1, x2: Bjorken x of incoming partons
C             ipara = 1 use parametrization in terms of d, db dt, dbt etc.
C                     else use parametrization in a1, a2, a3
C             anomalous couplings: rsmin,din,dbin,dtin,dtbin, lambdahvvin,
C             a1hwwin,a2hwwin,a3hwwin,a1haain,a2haain,a3haain,
C             a1hazin,a2hazin,a3hazin,a1hzzin,a2hzzin,a3hzzin
C             pjet1(0:3):  E,px,py,pz of 1st jet
C             pjet1(0:3):  E,px,py,pz of 2nd jet
C             phiggs(0:3): E,px,py,pz of Higgs boson
C      -----------------------------------------------------------------------------
C      output: weight - weight for reweighting from SM to anomalous coupling choice
************************************************************************************
      implicit none
C input and output variables
      integer ipara
      real*8 rsmin,din,dbin,dtin,dtbin,lambdahvvin
      real*8 a1hwwin,a2hwwin,a3hwwin,a1haain,a2haain,a3haain
      real*8 a1hazin,a2hazin,a3hazin,a1hzzin,a2hzzin,a3hzzin
C
      real*8 ecms,m2sm,m2full,x1,x2
      real*8, intent(out) :: weight
      !f2py intent(out) :: weight
      real*8 pjet1(0:3),pjet2(0:3),phiggs(0:3)
      integer i1,i2,i,j
C
      real*8 xp1,xp2,wmat2
      real*8 p(6,0:3),pdf1(-6:6),pdf2(-6:6)
      real*8 m2i(-5:5,-5:5),m2if(-5:5,-5:5,-5:5,-5:5)
C
      integer qborn,qw,qz,qschan,qtchan,qch2,qchint,
     &     qbini,qbfin,qwidth,qfact,qbos,qferm,qsoft,qhh2,
     &     qqcddiag,qqcdnondiag,qqcdgsplit,qqcdggfus,qcp
      common/rcoptions/qborn,qw,qz,qschan,qtchan,qch2,qchint,
     &     qbini,qbfin,qwidth,qfact,qbos,qferm,qsoft,qhh2,
     &     qqcddiag,qqcdnondiag,qqcdgsplit,qqcdggfus,qcp
C
      real*8 pi,el,alpha,alpha0,alphaz,GF,alphas
      real*8 cw,cw2,sw,sw2,mw,mw2,gw,mz,mz2,gz,mh,mh2
      real*8 ml,ml2,mqp,mqp2,mqm,mqm2

      complex*16 v,cv
      common/param/pi,el,alpha,alpha0,alphaz,GF,alphas,
     &     v(3,3),cv(3,3),cw,cw2,sw,sw2,mw,mw2,gw,mz,mz2,gz,mh,mh2,
     &     ml(3),ml2(3),mqp(3),mqp2(3),mqm(3),mqm2(3)
C     
      real*8 sinthetac,mfact
      real*8 ppecms
      common/ppecms/ppecms
C
      complex*16 ccw,ccw2,csw,csw2,cmw,cmw2,cmz,cmz2
      common/cparam/ccw,ccw2,csw,csw2,cmw,cmw2,cmz,cmz2
C
      complex*16 xcw,xcw2,xsw,xsw2,xmw,xmw2,xmz,xmz2,xmh,xmh2,null
      complex*16 xml,xml2,xmqp,xmqp2,xmqm,xmqm2
C
      common/xparam/xcw,xcw2,xsw,xsw2,xmw,xmw2,xmz,xmz2,xmh,xmh2,null, 
     &     xml(3),xml2(3),xmqp(3),xmqp2(3),xmqm(3),xmqm2(3) 
C
      real*8 qu,qd,ql,qn,qf,mu,mu2,md,md2,mlep,mlep2
      complex*16 guu,gdd,gnn,gll     
      common/qf/qu,qd,ql,qn,qf(4),mu,mu2,md,md2,mlep,mlep2,
     +     guu(-1:1),gdd(-1:1),gnn(-1:1),gll(-1:1)
C
      integer qhvv

      real*8 rsm,d,db,dt,dtb,lambdahvv
      real*8 a1hww,a2hww,a3hww,a1haa,a2haa,a3haa
      real*8 a1haz,a2haz,a3haz,a1hzz,a2hzz,a3hzz


      common/hvv/rsm,d,db,dt,dtb,lambdahvv,
     &     a1hww,a2hww,a3hww,a1haa,a2haa,a3haa,
     &     a1haz,a2haz,a3haz,a1hzz,a2hzz,a3hzz,qhvv
C     
      ppecms = ecms             !  pp-center of mass energy 
c
C P(1,x)  = 1st incoming proton P(2,x) = 2nd incoming proton 
C P(3,x)  = 1st outgoing jet    P(3,x) = 2nd outgoing jet 
C P(5,x)  = Higgs boosn candidate       
C
      xp1=x1                    !  Bjorken x1
      xp2=x2                    !  Bjorken x2
      p(1,0) = xp1*ppecms/2.0d0
      p(1,1) = 0.0d0
      p(1,2) = 0.0d0
      p(1,3) = xp1*ppecms/2.0d0
      p(2,0) = xp2*ppecms/2.0d0
      p(2,1) = 0.0d0
      p(2,2) = 0.0d0
      p(2,3) = -xp2*ppecms/2.0d0
      do I = 0,3
         p(3,i) = pjet1(i)
         p(4,i) = pjet2(i)
         p(5,i) = phiggs(i)
         p(6,i) = 0.0d0
      enddo


      rsm       = rsmin
      d         = din
      db        = dbin
      dt        = dtin
      dtb       = dtbin
      lambdahvv = lambdahvvin
      a1hww     = a1hwwin
      a2hww     = a2hwwin
      a3hww     = a3hwwin
      a1haa     = a1haain
      a2haa     = a2haain
      a3haa     = a3haain
      a1haz     = a1hazin
      a2haz     = a2hazin
      a3haz     = a3hazin
      a1hzz     = a1hzzin
      a2hzz     = a2hzzin
      a3hzz     = a3hzzin

C set HAWK "constants" (ony those needed by LO calculation)
      alpha = 1.0d0/137.0d0
      pi = 4d0*datan(1d0)
      el = sqrt(4d0*pi*alpha)   ! alpha, alpha0 oder alphaz
C
      mz   = 91.1876d0
      cmz = mz
      xmz = mz
      mz2 = mz*mz
      xmz2 = mz2 
      cmz2  = mz2
C
      mw   = 80.425d0
      xmw  = mw
      cmw  = mw
      mw2  = mw*mw
      xmw2 = mw2
      cmw2 = mw2
C
      mh   = 125d0
C
      cw  = mw/mz
      xcw = cw
      ccw = cw 
      cw2 = cw*cw
      xcw2 = cw2
      ccw2 = cw2
      sw2  = 1d0-cw2
      csw2 = sw2
      xsw2 = sw2
      sw  = sqrt(sw2)
      xsw  = sw 
      csw  = sw 
C fermion el-mag charges
      qu =  2d0/3d0
      qd =  -1d0/3d0
      qn =   0d0
      ql =   -1d0
C fermion couplings 
      guu(+1) = -qu*xsw/xcw
      guu(-1) = -qu*xsw/xcw+1d0/2d0/xsw/xcw
      gdd(+1) = -qd*xsw/xcw
      gdd(-1) = -qd*xsw/xcw-1d0/2d0/xsw/xcw
      gnn(+1) = -qn*xsw/xcw
      gnn(-1) = -qn*xsw/xcw+1d0/2d0/xsw/xcw
      gll(+1) = -ql*xsw/xcw
      gll(-1) = -ql*xsw/xcw-1d0/2d0/xsw/xcw
C CKM matrix and its conjugated
      sinthetac = 0.226548d0
      v(1,2) = sinthetac
      v(1,1) = sqrt(1d0-v(1,2)*v(1,2)) ! bug v(1,1)
      v(2,1) = -v(1,2)
      v(2,2) = v(1,1)
      v(1,3) = 0d0
      v(2,3) = 0d0
      v(3,1) = 0d0
      v(3,2) = 0d0
      v(3,3) = 1d0
      do I =1,3
         do J = 1,3
            cv(i,j) = dconjg(v(i,j))
         enddo
      enddo
C set HAWK flags (only those needed by LO calculation)
      qw     = 1 
      qz     = 1 
C
      qschan = 0 
      qtchan = 1 
      qch2   = 1
      qchint = 0 
C     
      qbini       = 0 
      qbfin       = 0 
      qwidth      = 0 
      qcp         = 0 
C anomalous couplings
      if (ipara.eq.1) then  
         a1hww = mw/sw*rsm
         a2hww = 2d0*d /sw/mw
         a3hww = 2d0*dt/sw/mw
         a1haa = 0d0
         a2haa = 4d0*(d *sw2+db *cw2)/2d0/sw/mw
         a3haa = 4d0*(dt*sw2+dtb*cw2)/2d0/sw/mw
         a1haz = 0d0
         a2haz = -2d0*cw*(d -db )/mw ! sign change due to BHS convention
         a3haz = -2d0*cw*(dt-dtb)/mw ! sign change due to BHS convention
         a1hzz = mw/sw/cw2*rsm
         a2hzz = 4d0*(d *cw2+db *sw2)/2d0/sw/mw
         a3hzz = 4d0*(dt*cw2+dtb*sw2)/2d0/sw/mw
      endif
C otherwise use direct in put value in terms if a1 a2 a3
c PDFs
      mfact = mh
      call getpdf(xp1,mfact,pdf1)
      call getpdf(xp2,mfact,pdf2)
C for a test withour pdf weighting
C      do i = -6,6
C         pdf1(i) = 1.0d0
C         pdf2(i) = 1.0d0
C      enddo
C     SM:
      qhvv = 0 
c squared matrix elements
      call Mat2born0(p,m2i,m2if,1)
c weighted |M|^2
      wmat2 = 0d0
      do 100 i1=-4,4
         do 100 i2=-4,4
            if (i1*i2.eq.0) goto 100
            wmat2 = wmat2 + m2i(i1,i2)*pdf1(i1)*pdf2(i2) 
 100  continue
      m2sm = wmat2
C     full
      qhvv = 1 
C
c squared matrix elements
      call Mat2born0(p,m2i,m2if,1)
c weighted |M|^2
      wmat2 = 0d0
      do 200 i1=-4,4
         do 200 i2=-4,4
            if (i1*i2.eq.0) goto 200
            wmat2 = wmat2 + m2i(i1,i2)*pdf1(i1)*pdf2(i2) 
 200  continue
      m2full = wmat2
C
      weight = m2full/m2sm
      return
      end
********************************************************************************

************************************************************************
      subroutine Mat2born0(p,m2i0,m2if0,qbcalc)
************************************************************************
*     Born structures needed for subtraction function for
*     parton1(p1) + parton2(p2) --> f(p3) + f'(p4) + H(p5)
*     qbcalc = 0/1: b-quarks ex-/in-cluded in evaluation
*     
*     Note: anomalous HVV couplings for qhvv>0 included
*-----------------------------------------------------------------------
*     13.10.06 Stefan Dittmaier
************************************************************************
      implicit real*8 (a-z)
      complex*16 v,cv
      complex*16 ccw,ccw2,csw,csw2,cmw,cmw2,cmz,cmz2,cmh,cmh2
      complex*16 pws,pws12,pwt11,pwt22,pwt12,pwt21
      complex*16 pzs,pzs12,pzt11,pzt22,pzt12,pzt21
      real*8 p(6,0:3)
      real*8 m2aqq0ss(5,5,5,5),m2qaq0ss(5,5,5,5)
      real*8 m2qq0ss(5,5,5,5),m2aqaq0ss(5,5,5,5)
      real*8 m2aqq0tt(5,5,5,5),m2qaq0tt(5,5,5,5)
      real*8 m2qq0tt(5,5,5,5),m2aqaq0tt(5,5,5,5)
      real*8 m2aqq0st(5,5,5,5),m2qaq0st(5,5,5,5)
      real*8 m2qq0st(5,5,5,5),m2aqaq0st(5,5,5,5)
      real*8 m2i0(-5:5,-5:5)
      real*8 m2if0(-5:5,-5:5,-5:5,-5:5)
      integer i1,i2,i3,i4,i,j,qb,qbcalc
      integer qborn,qw,qz,qschan,qtchan,qch2,qchint,
     &     qbini,qbfin,qwidth,qfact,qbos,qferm,qsoft,qhh2,
     &     qqcddiag,qqcdnondiag,qqcdgsplit,qqcdggfus,qcp
      integer qhvv
      
      common/param/pi,el,alpha,alpha0,alphaz,GF,alphas,
     &     v(3,3),cv(3,3),cw,cw2,sw,sw2,mw,mw2,gw,mz,mz2,gz,mh,mh2,
     &     ml(3),ml2(3),mqp(3),mqp2(3),mqm(3),mqm2(3)
      common/cparam/ccw,ccw2,csw,csw2,cmw,cmw2,cmz,cmz2
      common/rcoptions/qborn,qw,qz,qschan,qtchan,qch2,qchint,
     &     qbini,qbfin,qwidth,qfact,qbos,qferm,qsoft,qhh2,
     &     qqcddiag,qqcdnondiag,qqcdgsplit,qqcdggfus,qcp
      common/hvv/rsm,d,db,dt,dtb,lambdahvv,
     &     a1hww,a2hww,a3hww,a1haa,a2haa,a3haa,
     &     a1haz,a2haz,a3haz,a1hzz,a2hzz,a3hzz,qhvv
      
      data m2aqq0ss/625*0/,m2aqq0st/625*0/,m2aqq0tt/625*0/
      data m2qq0ss/625*0/,m2qq0st/625*0/,m2qq0tt/625*0/
      data m2aqaq0ss/625*0/,m2aqaq0st/625*0/,m2aqaq0tt/625*0/
      data m2qaq0ss/625*0/,m2qaq0st/625*0/,m2qaq0tt/625*0/
      
c energies and angles
      eb1  = p(1,0)
      eb2  = p(2,0)
      e1   = p(3,0)
      e2   = p(4,0)
      cth1 = p(3,3)/p(3,0)
      cth2 = p(4,3)/p(4,0)
      if (p(1,3).lt.0d0) then
         cth1 = -cth1
         cth2 = -cth2
      endif
      sth1 = sqrt(p(3,1)**2+p(3,2)**2)/p(3,0)
      sth2 = sqrt(p(4,1)**2+p(4,2)**2)/p(4,0)
      if (cth1.gt.0d0) then
         c22th1 = (1d0+cth1)*.5d0
         s22th1 = sth1**2/(4d0*c22th1)
      else
         s22th1 = (1d0-cth1)*.5d0
         c22th1 = sth1**2/(4d0*s22th1)
      endif
      if (cth2.gt.0d0) then
         c22th2 = (1d0+cth2)*.5d0
         s22th2 = sth2**2/(4d0*c22th2)
      else
         s22th2 = (1d0-cth2)*.5d0
         c22th2 = sth2**2/(4d0*s22th2)
      endif

c invariants
      s    =  4d0*eb1*eb2
      t11  = -4d0*eb1*e1*s22th1
      t12  = -4d0*eb1*e2*s22th2
      t21  = -4d0*eb2*e1*c22th1
      t22  = -4d0*eb2*e2*c22th2
c     s12  = mh2-s+4d0*e1*eb1*s22th1+4d0*e1*eb2*c22th1
c     &		    +4d0*e2*eb1*s22th2+4d0*e2*eb2*c22th2
      s12  = 2d0*( p(3,0)*p(4,0)-p(3,1)*p(4,1)
     &     -p(3,2)*p(4,2)-p(3,3)*p(4,3))
      
c propagators
      if(qwidth.eq.1) then
         pws12 = 1d0/(s12-cmw2)
         pzs12 = 1d0/(s12-cmz2)
         pws = 1d0/(s-cmw2)
         pzs = 1d0/(s-cmz2)
          pwt11 = 1d0/(t11-cmw2)
          pzt11 = 1d0/(t11-cmz2)
          pwt22 = 1d0/(t22-cmw2)
          pzt22 = 1d0/(t22-cmz2)
          pwt12 = 1d0/(t12-cmw2)
          pzt12 = 1d0/(t12-cmz2)
          pwt21 = 1d0/(t21-cmw2)
          pzt21 = 1d0/(t21-cmz2)
        else

          if ((qwidth.eq.0).or.((qwidth.eq.2).and.(s12.lt.0d0))) then
            pws12 = 1d0/(s12-mw2)
            pzs12 = 1d0/(s12-mz2)
          else
            pws12 = 1d0/(s12-cmw2)
            pzs12 = 1d0/(s12-cmz2)
          endif
          
          if ((qwidth.eq.0).or.((qwidth.eq.2).and.(s.lt.0d0))) then
            pws = 1d0/(s-mw2)
            pzs = 1d0/(s-mz2)
          else
            pws = 1d0/(s-cmw2)
            pzs = 1d0/(s-cmz2)
          endif
          
          if ((qwidth.eq.0).or.((qwidth.eq.2).and.(t11.lt.0d0))) then
            pwt11 = 1d0/(t11-mw2)
            pzt11 = 1d0/(t11-mz2)
          else
            pwt11 = 1d0/(t11-cmw2)
            pzt11 = 1d0/(t11-cmz2)
          endif
          
          if ((qwidth.eq.0).or.((qwidth.eq.2).and.(t22.lt.0d0))) then
            pwt22 = 1d0/(t22-mw2)
            pzt22 = 1d0/(t22-mz2)
          else
            pwt22 = 1d0/(t22-cmw2)
            pzt22 = 1d0/(t22-cmz2)
          endif
          
          if ((qwidth.eq.0).or.((qwidth.eq.2).and.(t12.lt.0d0))) then
            pwt12 = 1d0/(t12-mw2)
            pzt12 = 1d0/(t12-mz2)
          else
            pwt12 = 1d0/(t12-cmw2)
            pzt12 = 1d0/(t12-cmz2)
          endif
          
          if ((qwidth.eq.0).or.((qwidth.eq.2).and.(t21.lt.0d0))) then
            pwt21 = 1d0/(t21-mw2)
            pzt21 = 1d0/(t21-mz2)
          else
            pwt21 = 1d0/(t21-cmw2)
            pzt21 = 1d0/(t21-cmz2)
          endif
        endif
C
	qb = qbcalc
	if ((qbini.eq.0).and.(qbfin.eq.0)) qb = 0
C
	if (qhvv.eq.0) then

c qbar-q --> qbar-q
           call Mat2aqqborn(m2aqq0ss,m2aqq0tt,m2aqq0st,
     &          s,s12,t11,t21,t12,t22,
     &          pws,pws12,pwt11,pwt22,pzs,pzs12,pzt11,pzt22,qb)
c q-q    --> q-q
           call Mat2aqqborn(m2qq0ss,m2qq0tt,m2qq0st,
     &          t21,t12,t11,s,s12,t22,
     &          pwt21,pwt12,pwt11,pwt22,pzt21,pzt12,pzt11,pzt22,qb)

c other channels from CP symmetry for qcp=1
           if (qcp.eq.0) then
c q-qbar --> q-qbar
              call Mat2aqqborn(m2qaq0ss,m2qaq0tt,m2qaq0st,
     &             s,s12,t22,t12,t21,t11,
     &             pws,pws12,pwt22,pwt11,pzs,pzs12,pzt11,pzt22,qb)
c qbar-qbar --> qbar-qbar
              call Mat2aqqborn(m2aqaq0ss,m2aqaq0tt,m2aqaq0st,
     &             t12,t21,t11,s12,s,t22,
     &             pwt12,pwt21,pwt11,pwt22,pzt21,pzt12,pzt11,pzt22,qb)
           else
              do 103 i4=1,4+qb
                 do 103 i3=1,4+qb
                    do 103 i2=1,4+qb
                       do 103 i1=1,4+qb
                          m2qaq0ss(i1,i2,i3,i4)  = m2aqq0ss(i1,i2,i3,i4) 
                          m2qaq0tt(i1,i2,i3,i4)  = m2aqq0tt(i1,i2,i3,i4) 
                          m2qaq0st(i1,i2,i3,i4)  = m2aqq0st(i1,i2,i3,i4) 
                          m2aqaq0ss(i1,i2,i3,i4) = m2qq0ss(i1,i2,i3,i4)  
                          m2aqaq0tt(i1,i2,i3,i4) = m2qq0tt(i1,i2,i3,i4)  
                          m2aqaq0st(i1,i2,i3,i4) = m2qq0st(i1,i2,i3,i4)  
 103          continue
           endif
                    
        else

	call setprods3(p)

c qbar-q --> qbar-q
        call Mat2aqqbornhvv(m2aqq0ss,m2aqq0tt,m2aqq0st,-1,-2,3,4,qb)
c q-q    --> q-q
        call Mat2aqqbornhvv(m2qq0ss,m2qq0tt,m2qq0st,3,-2,-1,4,qb)
c q-qbar --> q-qbar
        call Mat2aqqbornhvv(m2qaq0ss,m2qaq0tt,m2qaq0st,-2,-1,4,3,qb)
c qbar-qbar --> qbar-qbar
        call Mat2aqqbornhvv(m2aqaq0ss,m2aqaq0tt,m2aqaq0st,-1,4,3,-2,qb)

	endif

c add final states for fixed initial states
        do 200 i1=-4-qbini,4+qbini
        do 200 i2=-4-qbini,4+qbini
          m2i0(i1,i2) = 0d0
          if (i1*i2.eq.0) goto 200
        do 300 i3=1,4+qbfin
        do 300 i4=1,4+qbfin
          if(i1.gt.0) then
            if ((i2.gt.0)) then
              m2if0(i1,i2,i3,i4)    = m2qq0ss(i3,i2,i1,i4)*.5d0
     &           + m2qq0tt(i3,i2,i1,i4)*.5d0 + m2qq0st(i3,i2,i1,i4)*.5d0
              m2i0(i1,i2) = m2i0(i1,i2) + m2if0(i1,i2,i3,i4)
            else
              m2if0(i1,i2,i3,-i4)    = m2qaq0ss(-i2,i1,i4,i3)
     &           + m2qaq0tt(-i2,i1,i4,i3) + m2qaq0st(-i2,i1,i4,i3)
              m2i0(i1,i2) = m2i0(i1,i2) + m2if0(i1,i2,i3,-i4)
            endif
          elseif (i2.gt.0) then
            m2if0(i1,i2,-i3,i4)    = m2aqq0ss(-i1,i2,i3,i4)
     &         + m2aqq0tt(-i1,i2,i3,i4) + m2aqq0st(-i1,i2,i3,i4)
            m2i0(i1,i2) = m2i0(i1,i2) + m2if0(i1,i2,-i3,i4)
          else
            m2if0(i1,i2,-i3,-i4)    = m2aqaq0ss(-i1,i4,i3,-i2)*.5d0
     &         + m2aqaq0tt(-i1,i4,i3,-i2)*.5d0 
     &	       + m2aqaq0st(-i1,i4,i3,-i2)*.5d0
            m2i0(i1,i2) = m2i0(i1,i2) + m2if0(i1,i2,-i3,-i4)
	  endif
300     continue
200     continue

	end

************************************************************************
************************************************************************
      subroutine Mat2aqqborn(m2aqq0ss,m2aqq0tt,m2aqq0st,
     &     s,s12,t11,t21,t12,t22,
     &     pws,pws12,pwt11,pwt22,pzs,pzs12,pzt11,pzt22,qb)
************************************************************************
*       generic Born structures needed for subtraction function for
*	    anti-q(p1) + q(p2) --> f(p3) + f'(p4) + H(p5)
*
*       from generic amplitudes for
*	    H(p) --> fa(k1) + anti-fb(k2) + fc(k3) + anti-fd(k4)
*	         CC: f1     + anti-f2     + f3     + anti-f4     
*	         NC: f1     + anti-f1     + f3     + anti-f3     
*	             f2     + anti-f2     + f3     + anti-f3     
*	             f1     + anti-f1     + f4     + anti-f4     
*
*	fermions: f1,f4 = generic   up-type fermions
*	          f2,f3 = generic down-type fermions
*-----------------------------------------------------------------------
*       13.10.06 Stefan Dittmaier
************************************************************************
      implicit real*8 (a-z)
      complex*16 pre,pres,pret,mats0,matt0,wmats0,wmatt0
	complex*16 v,cv
	complex*16 pws,pws12,pwt11,pwt22
	complex*16 pzs,pzs12,pzt11,pzt22
        complex*16 ccw,ccw2,csw,csw2,cmw,cmw2,cmz,cmz2,cmh,cmh2
        complex*16 xcw,xcw2,xsw,xsw2,xmw,xmw2,xmz,xmz2,xmh,xmh2,null
        complex*16 xml,xml2,xmqp,xmqp2,xmqm,xmqm2
	real*8 m1243sq(-1:1,-1:1),m4213sq(-1:1,-1:1)
        real*8 m2aqq0ss(5,5,5,5),m2aqq0tt(5,5,5,5),m2aqq0st(5,5,5,5)
	integer ia,ib,ic,id
	integer i1,i2,i3,i4,j,j1,j2,f1,f2,f3,f4,g1,g2,g3,g4,qb
        integer qborn,qw,qz,qschan,qtchan,qch2,qchint,
     &                   qbini,qbfin,qwidth,qfact,qbos,qferm,qsoft,qhh2,
     &			 qqcddiag,qqcdnondiag,qqcdgsplit,qqcdggfus,qcp
	complex*16 guu,gdd,gnn,gll     
	logical checkss12,checkt11t22

        common/param/pi,el,alpha,alpha0,alphaz,GF,alphas,
     &         v(3,3),cv(3,3),cw,cw2,sw,sw2,mw,mw2,gw,mz,mz2,gz,mh,mh2,
     &         ml(3),ml2(3),mqp(3),mqp2(3),mqm(3),mqm2(3)
        common/rcoptions/qborn,qw,qz,qschan,qtchan,qch2,qchint,
     &                   qbini,qbfin,qwidth,qfact,qbos,qferm,qsoft,qhh2,
     &			 qqcddiag,qqcdnondiag,qqcdgsplit,qqcdggfus,qcp
        common/cparam/ccw,ccw2,csw,csw2,cmw,cmw2,cmz,cmz2
        common/xparam/xcw,xcw2,xsw,xsw2,xmw,xmw2,xmz,xmz2,xmh,xmh2,null,
     &       xml(3),xml2(3),xmqp(3),xmqp2(3),xmqm(3),xmqm2(3)
        common/qf/qu,qd,ql,qn,qf(4),mu,mu2,md,md2,mlep,mlep2,
     +       guu(-1:1),gdd(-1:1),gnn(-1:1),gll(-1:1)

	el3 = el**3

        m1243sq(1,1)   = 4D0*t12*t21
        m1243sq(-1,1)  = 4D0*t11*t22
        m1243sq(1,-1)  = m1243sq(-1,1)
        m1243sq(-1,-1) = m1243sq(1,1)

        M4213sq(1,1)   = m1243sq(1,1)
        M4213sq(-1,1)  = 4D0*s*s12
        M4213sq(1,-1)  = M4213sq(-1,1)
        M4213sq(-1,-1) = m1243sq(-1,-1)

	checkss12   = ((  s*qschan.gt.0d0).or.(  s*qtchan.lt.0d0)).and.
     &	              ((s12*qschan.gt.0d0).or.(s12*qtchan.lt.0d0))
	checkt11t22 = ((t11*qschan.gt.0d0).or.(t11*qtchan.lt.0d0)).and.
     &	              ((t22*qschan.gt.0d0).or.(t22*qtchan.lt.0d0))


c qbar-q scattering
c
c   s-channel:                             t-channel:
c
c   fa(p1) \_          / fd(p3)            fa(p1) --<--*--<-- fb(p3)
c          |\        |/                                >
c            \       /-                                >
c             *vvvvv*                                  >
c           _/       \                                 >
c           /|        \|                               >
c   fb(p2) /          -\ fc(p4)            fd(p2) -->--*-->-- fc(p4)


c Type: ubar-d -> ubar-d
c-----------------------

c*** generic amplitudes
	pres = 0d0
	pret = 0d0

c ubar-d -W-> cbar-s    (k1->-p1,k2->-p2,k3->p4,k4->p3)   
        if (checkss12) pres = qw*el3*PWS12*PWS*XMW/(2D0*XSW*XSW2)

c ubar-s -Z-> ubar-s    (k1->-p1,k2->p3,k3->p4,k4->-p2)   
        if (checkt11t22) pret = qz*el3*PZT11*PZT22/XSW/XCW2*XMW

c*** individual channels
 	do 200 g1=1,2
 	do 200 g2=1,2+qb
 	do 200 g3=1,2
 	do 200 g4=1,2+qb
 	  f1 = 2*g1
 	  f2 = 2*g2-1
 	  f3 = 2*g3
 	  f4 = 2*g4-1
 	  m2aqq0ss(f1,f2,f3,f4) = 0d0
 	  m2aqq0tt(f1,f2,f3,f4) = 0d0
 	  m2aqq0st(f1,f2,f3,f4) = 0d0

	do 200 i1=1,4

	  if (i1.eq.1) then
	    mats0  = 0d0
	    matt0  = pret*guu(-1)*gdd(+1)
	    fac    = m4213sq(-1,+1)
	  elseif (i1.eq.2) then
	    mats0  = -pres
	    matt0  = pret*guu(-1)*gdd(-1)
	    fac    = m4213sq(-1,-1)
	  elseif (i1.eq.3) then
	    mats0  = 0d0
	    matt0  = pret*guu(+1)*gdd(+1)
	    fac    = m4213sq(+1,+1)
	  elseif (i1.eq.4) then
	    mats0  = 0d0
	    matt0  = pret*guu(+1)*gdd(-1)
	    fac    = m4213sq(+1,-1)
	  endif

 	  wmats0   = cv(g1,g2)*v(g3,g4)*mats0

 	  if ((f1.eq.f3).and.(f2.eq.f4)) then
 	    m2aqq0ss(f1,f2,f3,f4) = m2aqq0ss(f1,f2,f3,f4) 
     &		         + fac*9d0*abs(wmats0)**2 *qch2
 	    m2aqq0tt(f1,f2,f3,f4) = m2aqq0tt(f1,f2,f3,f4) 
     &		         + fac*9d0*abs(matt0)**2 *qch2
 	    m2aqq0st(f1,f2,f3,f4) = m2aqq0st(f1,f2,f3,f4) 
     &	          - fac*3d0*2d0*dreal(wmats0*dconjg(matt0)) *qchint 
 	  else
 	    m2aqq0ss(f1,f2,f3,f4) = m2aqq0ss(f1,f2,f3,f4) 
     &			 + fac*9d0*abs(wmats0)**2 *qch2 
 	  endif

200	continue

c Type: dbar-u -> dbar-u
c-----------------------

c*** generic amplitudes
	pres = 0d0
	pret = 0d0

c dbar-u -W-> sbar-c    (k3->-p1,k4->-p2,k1->p4,k2->p3)   
	if (checkss12) pres = qw*el3*PWS12*PWS*XMW/(2D0*XSW*XSW2)

c dbar-c -Z-> dbar-c    (k3->-p1,k4->p3,k1->p4,k2->-p2)   
        if (checkt11t22) pret = qz*el3*PZT11*PZT22*XMW/(XSW*XCW2)

 	do 205 g1=1,2+qb
 	do 205 g2=1,2
 	do 205 g3=1,2+qb
 	do 205 g4=1,2
 	  f1 = 2*g1-1
 	  f2 = 2*g2
 	  f3 = 2*g3-1
 	  f4 = 2*g4
 	  m2aqq0ss(f1,f2,f3,f4) = 0d0
 	  m2aqq0tt(f1,f2,f3,f4) = 0d0
 	  m2aqq0st(f1,f2,f3,f4) = 0d0

	do 205 i1=1,4

	  if (i1.eq.1) then
	    mats0  = 0d0
	    matt0  = pret*guu(+1)*gdd(-1)
	    fac    = m4213sq(-1,+1)
	  elseif (i1.eq.2) then
	    mats0  = -pres
	    matt0  = pret*guu(-1)*gdd(-1)
	    fac    = m4213sq(-1,-1)
	  elseif (i1.eq.3) then
	    mats0  = 0d0
	    matt0  = pret*guu(+1)*gdd(+1)
	    fac    = m4213sq(+1,+1)
	  elseif (i1.eq.4) then
	    mats0  = 0d0
	    matt0  = pret*guu(-1)*gdd(+1)
	    fac    = m4213sq(+1,-1)
	  endif

 	  wmats0 = v(g2,g1)*cv(g4,g3)*mats0

 	  if ((f1.eq.f3).and.(f2.eq.f4)) then
 	    m2aqq0ss(f1,f2,f3,f4) = m2aqq0ss(f1,f2,f3,f4) 
     &		         + fac*9d0*abs(wmats0)**2 *qch2 
 	    m2aqq0tt(f1,f2,f3,f4) = m2aqq0tt(f1,f2,f3,f4) 
     &		         + fac*9d0*abs(matt0)**2 *qch2 
 	    m2aqq0st(f1,f2,f3,f4) = m2aqq0st(f1,f2,f3,f4) 
     &	          - fac*3d0*2d0*dreal(wmats0*dconjg(matt0)) *qchint
 	  else
 	    m2aqq0ss(f1,f2,f3,f4) = m2aqq0ss(f1,f2,f3,f4) 
     &			 + fac*9d0*abs(wmats0)**2 *qch2 
 	  endif
205	continue


c Type: ubar-u -> ubar-u
c----------------------- 

c*** generic amplitudes
	pres = 0d0
	pret = 0d0

c ubar-u -Z-> cbar-c  (k1->-p1,k2->-p2,k3->p4,k4->p3)
	if (checkss12) pres = qz*el3*PZS*PZS12*XMW/(XSW*XCW2)

c ubar-c -Z-> ubar-c  (k1->-p1,k2->p3,k3->p4,k4->-p2)
	if (checkt11t22) pret = qz*el3*PZT11*PZT22*XMW/(XSW*XCW2)

	m2aqq0ss(2,2,2,2) = 0d0
	m2aqq0tt(2,2,2,2) = 0d0
	m2aqq0st(2,2,2,2) = 0d0
	m2aqq0ss(2,2,4,4) = 0d0
	m2aqq0tt(2,2,4,4) = 0d0
	m2aqq0st(2,2,4,4) = 0d0
	m2aqq0ss(2,4,2,4) = 0d0
	m2aqq0tt(2,4,2,4) = 0d0
	m2aqq0st(2,4,2,4) = 0d0

	do 201 i1=1,6

	if (i1.eq.1) then
	  mats0  = 0d0
	  matt0  = pret*guu(+1)*guu(-1)
	  fac    = m4213sq(-1,+1)
	elseif (i1.eq.2) then
	  mats0  = -pres*guu(-1)*guu(-1)
	  matt0  =  pret*guu(-1)*guu(-1)
	  fac    = m4213sq(-1,-1)
	elseif (i1.eq.3) then
	  mats0  = -pres*guu(+1)*guu(+1)
	  matt0  =  pret*guu(+1)*guu(+1)
	  fac    = m4213sq(+1,+1)
	elseif (i1.eq.4) then
	  mats0  = 0d0
	  matt0  = pret*guu(-1)*guu(+1)
	  fac    = m4213sq(+1,-1)
	elseif (i1.eq.5) then
	  mats0  = pres*guu(-1)*guu(+1)
	  matt0  = 0d0
	  fac    = m1243sq(-1,+1)
	elseif (i1.eq.6) then
	  mats0  = pres*guu(+1)*guu(-1)
	  matt0  = 0d0
	  fac    = m1243sq(+1,-1)
	endif

c*** ubar-u -ZZ-> ubar-u    
	m2aqq0ss(2,2,2,2) = m2aqq0ss(2,2,2,2) 
     &			    + fac*9d0*abs(mats0)**2 *qch2
	m2aqq0tt(2,2,2,2) = m2aqq0tt(2,2,2,2) 
     &			    + fac*9d0*abs(matt0)**2 *qch2 
	m2aqq0st(2,2,2,2) = m2aqq0st(2,2,2,2) 
     &	      - fac*3d0*2d0*dreal(mats0*dconjg(matt0)) *qchint
c*** ubar-u -Z-> cbar-c 
	m2aqq0ss(2,2,4,4) = m2aqq0ss(2,2,4,4) 
     &			    + fac*9d0*abs(mats0)**2 *qch2 
c*** ubar-c -Z-> ubar-c 
	m2aqq0tt(2,4,2,4) = m2aqq0tt(2,4,2,4) 
     &			    + fac*9d0*abs(matt0)**2 *qch2 
201	continue
	m2aqq0ss(4,4,4,4) = m2aqq0ss(2,2,2,2) 
	m2aqq0tt(4,4,4,4) = m2aqq0tt(2,2,2,2) 
	m2aqq0st(4,4,4,4) = m2aqq0st(2,2,2,2) 
	m2aqq0ss(4,4,2,2) = m2aqq0ss(2,2,4,4)
	m2aqq0tt(4,4,2,2) = m2aqq0tt(2,2,4,4)
	m2aqq0st(4,4,2,2) = m2aqq0st(2,2,4,4)
	m2aqq0ss(4,2,4,2) = m2aqq0ss(2,4,2,4)
	m2aqq0tt(4,2,4,2) = m2aqq0tt(2,4,2,4)
	m2aqq0st(4,2,4,2) = m2aqq0st(2,4,2,4)


c Type: ubar-u -> dbar-d
c----------------------- 

c*** generic amplitudes
	pres = 0d0
	pret = 0d0

c ubar-u -Z-> sbar-s  (k1->-p1,k2->-p2,k3->p4,k4->p3)
        if (checkss12) pres = qz*el3*PZS*PZS12*XMW/(XSW*XCW2)

c ubar-u -W-> dbar-d  (k1->-p1,k2->p3,k3->p4,k4->-p2)
	if (checkt11t22) pret = qw*el3*PWT11*PWT22*XMW/(2D0*XSW*XSW2)

 	do 202 g1=1,2
 	do 202 g2=1,2
 	do 202 g3=1,2+qb
 	do 202 g4=1,2+qb
 	  f1 = 2*g1
 	  f2 = 2*g2
 	  f3 = 2*g3-1
 	  f4 = 2*g4-1
 	  m2aqq0ss(f1,f2,f3,f4) = 0d0
 	  m2aqq0tt(f1,f2,f3,f4) = 0d0
 	  m2aqq0st(f1,f2,f3,f4) = 0d0

	do 202 i1=1,4

	  if (i1.eq.1) then
	    matt0  = 0d0
	    mats0  = pres*guu(+1)*gdd(-1)
	    fac    = m1243sq(-1,+1)
	  elseif (i1.eq.2) then
	    matt0  = -pret
	    mats0  = pres*guu(-1)*gdd(-1)
	    fac    = m1243sq(-1,-1)
	  elseif (i1.eq.3) then
	    matt0  = 0d0
	    mats0  = pres*guu(+1)*gdd(+1)
	    fac    = m1243sq(+1,+1)
	  elseif (i1.eq.4) then
	    matt0  = 0d0
	    mats0  = pres*guu(-1)*gdd(+1)
	    fac    = m1243sq(+1,-1)
	  endif

 	  wmatt0 = cv(g1,g3)*v(g2,g4)*matt0

 	  if ((f1.eq.f2).and.(f3.eq.f4)) then
 	    m2aqq0ss(f1,f2,f3,f4) = m2aqq0ss(f1,f2,f3,f4) 
     &		          + fac*9d0*abs(mats0)**2 *qch2 
 	    m2aqq0tt(f1,f2,f3,f4) = m2aqq0tt(f1,f2,f3,f4) 
     &		          + fac*9d0*abs(wmatt0)**2 *qch2
 	    m2aqq0st(f1,f2,f3,f4) = m2aqq0st(f1,f2,f3,f4) 
     &	          - fac*3d0*2d0*dreal(wmatt0*dconjg(mats0)) *qchint
 	  else
 	    m2aqq0tt(f1,f2,f3,f4) = m2aqq0tt(f1,f2,f3,f4) 
     &			 + fac*9d0*abs(wmatt0)**2 *qch2 
 	  endif
202	  continue


c Type: dbar-d -> dbar-d
c-----------------------

c*** generic amplitudes
	pres = 0d0
	pret = 0d0

c ubar-u -Z-> cbar-c  (k1->-p1,k2->-p2,k3->p4,k4->p3)
	if (checkss12) pres = qz*el3*PZS*PZS12*XMW/(XSW*XCW2)

c ubar-c -Z-> ubar-c  (k1->-p1,k2->p3,k3->p4,k4->-p2)
	if (checkt11t22) pret = qz*el3*PZT11*PZT22*XMW/(XSW*XCW2)

	m2aqq0ss(1,1,1,1) = 0d0
	m2aqq0tt(1,1,1,1) = 0d0
	m2aqq0st(1,1,1,1) = 0d0
	m2aqq0ss(1,1,3,3) = 0d0
	m2aqq0tt(1,1,3,3) = 0d0
	m2aqq0st(1,1,3,3) = 0d0
	m2aqq0ss(1,3,1,3) = 0d0
	m2aqq0tt(1,3,1,3) = 0d0
	m2aqq0st(1,3,1,3) = 0d0

	do 203 i1=1,6

	if (i1.eq.1) then
	  mats0  = 0d0
	  matt0  = pret*gdd(+1)*gdd(-1)
	  fac    = m4213sq(-1,+1)
	elseif (i1.eq.2) then
	  mats0  = -pres*gdd(-1)*gdd(-1)
	  matt0  =  pret*gdd(-1)*gdd(-1)
	  fac    = m4213sq(-1,-1)
	elseif (i1.eq.3) then
	  mats0  = -pres*gdd(+1)*gdd(+1)
	  matt0  =  pret*gdd(+1)*gdd(+1)
	  fac    = m4213sq(+1,+1)
	elseif (i1.eq.4) then
	  mats0  = 0d0
	  matt0  = pret*gdd(-1)*gdd(+1)
	  fac    = m4213sq(+1,-1)
	elseif (i1.eq.5) then
	  mats0  = pres*gdd(-1)*gdd(+1)
	  matt0  = 0d0
	  fac    = m1243sq(-1,+1)
	elseif (i1.eq.6) then
	  mats0  = pres*gdd(+1)*gdd(-1)
	  matt0  = 0d0
	  fac    = m1243sq(+1,-1)
	endif

c*** dbar-d -ZZ-> dbar-d    
	m2aqq0ss(1,1,1,1) = m2aqq0ss(1,1,1,1) 
     &			    + fac*9d0*abs(mats0)**2 *qch2
	m2aqq0tt(1,1,1,1) = m2aqq0tt(1,1,1,1) 
     &			    + fac*9d0*abs(matt0)**2 *qch2
	m2aqq0st(1,1,1,1) = m2aqq0st(1,1,1,1) 
     &	       - fac*3d0*2d0*dreal(mats0*dconjg(matt0)) *qchint
c*** dbar-d -Z-> sbar-s 
	m2aqq0ss(1,1,3,3) = m2aqq0ss(1,1,3,3) 
     &			    + fac*9d0*abs(mats0)**2 *qch2 
c*** dbar-s -Z-> dbar-s 
	m2aqq0tt(1,3,1,3) = m2aqq0tt(1,3,1,3) 
     &			    + fac*9d0*abs(matt0)**2 *qch2 
203	continue
	m2aqq0ss(3,3,3,3) = m2aqq0ss(1,1,1,1) 
	m2aqq0tt(3,3,3,3) = m2aqq0tt(1,1,1,1) 
	m2aqq0st(3,3,3,3) = m2aqq0st(1,1,1,1) 
	m2aqq0ss(3,3,1,1) = m2aqq0ss(1,1,3,3)
	m2aqq0tt(3,3,1,1) = m2aqq0tt(1,1,3,3)
	m2aqq0st(3,3,1,1) = m2aqq0st(1,1,3,3)
	m2aqq0ss(3,1,3,1) = m2aqq0ss(1,3,1,3)
	m2aqq0tt(3,1,3,1) = m2aqq0tt(1,3,1,3)
	m2aqq0st(3,1,3,1) = m2aqq0st(1,3,1,3)
	if (qb.eq.1) then
        m2aqq0ss(5,5,5,5) = m2aqq0ss(1,1,1,1)
        m2aqq0tt(5,5,5,5) = m2aqq0tt(1,1,1,1)
        m2aqq0st(5,5,5,5) = m2aqq0st(1,1,1,1)
        m2aqq0ss(1,1,5,5) = m2aqq0ss(1,1,3,3)
        m2aqq0tt(1,1,5,5) = m2aqq0tt(1,1,3,3)
        m2aqq0st(1,1,5,5) = m2aqq0st(1,1,3,3)
        m2aqq0ss(3,3,5,5) = m2aqq0ss(1,1,3,3)
        m2aqq0tt(3,3,5,5) = m2aqq0tt(1,1,3,3)
        m2aqq0st(3,3,5,5) = m2aqq0st(1,1,3,3)
        m2aqq0ss(5,5,1,1) = m2aqq0ss(1,1,3,3)
        m2aqq0tt(5,5,1,1) = m2aqq0tt(1,1,3,3)
        m2aqq0st(5,5,1,1) = m2aqq0st(1,1,3,3)
        m2aqq0ss(5,5,3,3) = m2aqq0ss(1,1,3,3)
        m2aqq0tt(5,5,3,3) = m2aqq0tt(1,1,3,3)
        m2aqq0st(5,5,3,3) = m2aqq0st(1,1,3,3)
        m2aqq0ss(1,5,1,5) = m2aqq0ss(1,3,1,3)
        m2aqq0tt(1,5,1,5) = m2aqq0tt(1,3,1,3)
        m2aqq0st(1,5,1,5) = m2aqq0st(1,3,1,3)
        m2aqq0ss(3,5,3,5) = m2aqq0ss(1,3,1,3)
        m2aqq0tt(3,5,3,5) = m2aqq0tt(1,3,1,3)
        m2aqq0st(3,5,3,5) = m2aqq0st(1,3,1,3)
        m2aqq0ss(5,1,5,1) = m2aqq0ss(1,3,1,3)
        m2aqq0tt(5,1,5,1) = m2aqq0tt(1,3,1,3)
        m2aqq0st(5,1,5,1) = m2aqq0st(1,3,1,3)
        m2aqq0ss(5,3,5,3) = m2aqq0ss(1,3,1,3)
        m2aqq0tt(5,3,5,3) = m2aqq0tt(1,3,1,3)
        m2aqq0st(5,3,5,3) = m2aqq0st(1,3,1,3)
	endif

c Type: dbar-d -> ubar-u
c-----------------------
	pres = 0d0
	pret = 0d0

c*** generic amplitudes

c dbar-d -Z-> cbar-c  (k3->-p1,k4->-p2,k1->p4,k2->p3)
        if (checkss12) pres = qz*el3*PZS*PZS12*XMW/(XSW*XCW2)

c dbar-s -W-> ubar-c  (k3->-p1,k4->p3,k1->p4,k2->-p2)
	if (checkt11t22) pret = qw*el3*PWT11*PWT22*XMW/(2D0*XSW*XSW2)

 	do 204 g1=1,2+qb
 	do 204 g2=1,2+qb
 	do 204 g3=1,2
 	do 204 g4=1,2
 	  f1 = 2*g1-1
 	  f2 = 2*g2-1
 	  f3 = 2*g3
 	  f4 = 2*g4
 	  m2aqq0ss(f1,f2,f3,f4) = 0d0
 	  m2aqq0tt(f1,f2,f3,f4) = 0d0
 	  m2aqq0st(f1,f2,f3,f4) = 0d0

	do 204 i1=1,4

	  if (i1.eq.1) then
	    matt0  = 0d0
	    mats0  = pres*guu(-1)*gdd(+1)
	    fac    = m1243sq(-1,+1)
	  elseif (i1.eq.2) then
	    matt0  = -pret
	    mats0  = pres*guu(-1)*gdd(-1)
	    fac    = m1243sq(-1,-1)
	  elseif (i1.eq.3) then
	    matt0  = 0d0
	    mats0  = pres*guu(+1)*gdd(+1)
	    fac    = m1243sq(+1,+1)
	  elseif (i1.eq.4) then
	    matt0  = 0d0
	    mats0  = pres*guu(+1)*gdd(-1)
	    fac    = m1243sq(+1,-1)
	  endif

 	  wmatt0 = v(g3,g1)*cv(g4,g2)*matt0

 	  if ((f1.eq.f2).and.(f3.eq.f4)) then
 	    m2aqq0ss(f1,f2,f3,f4) = m2aqq0ss(f1,f2,f3,f4) 
     &		         + fac*9d0*abs(mats0)**2 *qch2
 	    m2aqq0tt(f1,f2,f3,f4) = m2aqq0tt(f1,f2,f3,f4) 
     &		         + fac*9d0*abs(wmatt0)**2 *qch2
 	    m2aqq0st(f1,f2,f3,f4) = m2aqq0st(f1,f2,f3,f4) 
     &	          - fac*3d0*2d0*dreal(wmatt0*dconjg(mats0)) *qchint
 	  else
 	    m2aqq0tt(f1,f2,f3,f4) = m2aqq0tt(f1,f2,f3,f4) 
     &			 + fac*9d0*abs(wmatt0)**2 *qch2 
 	  endif
204	  continue

	end
************************************************************************
        subroutine Mat2aqqbornhvv(m2aqq0ss,m2aqq0tt,m2aqq0st,
     &	                          q1,q2,q3,q4,qb)
************************************************************************
*       generic Born structures needed for subtraction function for
*	    anti-q(p1) + q(p2) --> f(p3) + f'(p4) + H(p5)
*
*	including anomalous HVV couplings !!!
*
*       from generic amplitudes for
*	    H(p) --> fa(k1) + anti-fb(k2) + fc(k3) + anti-fd(k4)
*	         CC: f1     + anti-f2     + f3     + anti-f4     
*	         NC: f1     + anti-f1     + f3     + anti-f3     
*	             f2     + anti-f2     + f3     + anti-f3     
*	             f1     + anti-f1     + f4     + anti-f4     
*
*	fermions: f1,f4 = generic   up-type fermions
*	          f2,f3 = generic down-type fermions
*-----------------------------------------------------------------------
*       14.6.12 Stefan Dittmaier
************************************************************************
        implicit real*8 (a-z)
	complex*16 msH4f(-1:1,-1:1),mtH4f(-1:1,-1:1)
	complex*16 mats0,matt0,wmats0,wmatt0
	complex*16 sp(-4:4,-4:4),csp(-4:4,-4:4)
	complex*16 v,cv
        complex*16 ccw,ccw2,csw,csw2,cmw,cmw2,cmz,cmz2,cmh,cmh2
        complex*16 xcw,xcw2,xsw,xsw2,xmw,xmw2,xmz,xmz2,xmh,xmh2,null
        complex*16 xml,xml2,xmqp,xmqp2,xmqm,xmqm2
        real*8 m2aqq0ss(5,5,5,5),m2aqq0tt(5,5,5,5),m2aqq0st(5,5,5,5)
        real*8 vp(-4:4,-4:4)
	integer ia,ib,ic,id,q1,q2,q3,q4
	integer i1,i2,i3,i4,j,j1,j2,f1,f2,f3,f4,g1,g2,g3,g4,qb
        integer qborn,qw,qz,qschan,qtchan,qch2,qchint,
     &                   qbini,qbfin,qwidth,qfact,qbos,qferm,qsoft,qhh2,
     &			 qqcddiag,qqcdnondiag,qqcdgsplit,qqcdggfus,qcp
	complex*16 guu,gdd,gnn,gll     
C
        common/param/pi,el,alpha,alpha0,alphaz,GF,alphas,
     &         v(3,3),cv(3,3),cw,cw2,sw,sw2,mw,mw2,gw,mz,mz2,gz,mh,mh2,
     &         ml(3),ml2(3),mqp(3),mqp2(3),mqm(3),mqm2(3)
        common/rcoptions/qborn,qw,qz,qschan,qtchan,qch2,qchint,
     &                   qbini,qbfin,qwidth,qfact,qbos,qferm,qsoft,qhh2,
     &			 qqcddiag,qqcdnondiag,qqcdgsplit,qqcdggfus,qcp
        common/cparam/ccw,ccw2,csw,csw2,cmw,cmw2,cmz,cmz2
        common/xparam/xcw,xcw2,xsw,xsw2,xmw,xmw2,xmz,xmz2,xmh,xmh2,null,
     &       xml(3),xml2(3),xmqp(3),xmqp2(3),xmqm(3),xmqm2(3)
        common/qf/qu,qd,ql,qn,qf(4),mu,mu2,md,md2,mlep,mlep2,
     &       guu(-1:1),gdd(-1:1),gnn(-1:1),gll(-1:1)

	common/prods3/sp,csp,vp

	fac = el**6

c invariants
        s   = vp(q1,q2)
	s12 = vp(q3,q4)
	t11 = vp(q1,q3)
	t22 = vp(q2,q4)

	call fixahvv(s,s12,
     &             a1hwws,a2hwws,a3hwws,a1haas,a2haas,a3haas,
     &             a1hazs,a2hazs,a3hazs,a1hzzs,a2hzzs,a3hzzs)
	call fixahvv(t11,t22,
     &             a1hwwt,a2hwwt,a3hwwt,a1haat,a2haat,a3haat,
     &             a1hazt,a2hazt,a3hazt,a1hzzt,a2hzzt,a3hzzt)

c qbar-q scattering
c
c   s-channel:                             t-channel:
c
c   fa(p1) \_          / fd(p3)            fa(p1) --<--*--<-- fb(p3)
c          |\        |/                                >
c            \       /-                                >
c             *vvvvv*                                  >
c           _/       \                                 >
c           /|        \|                               >
c   fb(p2) /          -\ fc(p4)            fd(p2) -->--*-->-- fc(p4)


c Type: ubar-d -> ubar-d
c-----------------------

c*** generic amplitudes
	do i1=-1,1,2
	do i2=-1,1,2
	  msH4f(i1,i2) = 0d0
	  mtH4f(i1,i2) = 0d0
	enddo
	enddo

c ubar-d -W-> cbar-s    (k1->-p1,k2->-p2,k3->p4,k4->p3)   
	if (qw.eq.1) call M_WW_hvv(msH4f,q1,q2,q4,q3,
     &	                  a1hwws,a2hwws,a3hwws)

c ubar-s -Z-> ubar-s    (k1->-p1,k2->p3,k3->p4,k4->-p2)   
	if (qz.eq.1) call M_ZZ_hvv(mtH4f,guu,gdd,qu,qd,q1,q3,q4,q2,
     &	                  a1haat,a2haat,a3haat,a1hazt,a2hazt,a3hazt,
     &	                  a1hzzt,a2hzzt,a3hzzt)

c*** individual channels
 	do 200 g1=1,2
 	do 200 g2=1,2+qb
 	do 200 g3=1,2
 	do 200 g4=1,2+qb
 	  f1 = 2*g1
 	  f2 = 2*g2-1
 	  f3 = 2*g3
 	  f4 = 2*g4-1
 	  m2aqq0ss(f1,f2,f3,f4) = 0d0
 	  m2aqq0tt(f1,f2,f3,f4) = 0d0
 	  m2aqq0st(f1,f2,f3,f4) = 0d0

	do 200 i1=-1,1,2
	do 200 i2=-1,1,2
	do 200 i3=-1,1,2
	  mats0  = 0d0
	  matt0  = 0d0

          if (i1.eq.-i2) mats0  = msH4f(i2,i3)
          if (i1.eq.-i3) matt0  = mtH4f(i3,i2)

 	  wmats0 = cv(g1,g2)*v(g3,g4)*mats0

 	  if ((f1.eq.f3).and.(f2.eq.f4)) then
 	    m2aqq0ss(f1,f2,f3,f4) = m2aqq0ss(f1,f2,f3,f4) 
     &		         + fac*9d0*abs(wmats0)**2 *qch2
 	    m2aqq0tt(f1,f2,f3,f4) = m2aqq0tt(f1,f2,f3,f4) 
     &		         + fac*9d0*abs(matt0)**2 *qch2
 	    m2aqq0st(f1,f2,f3,f4) = m2aqq0st(f1,f2,f3,f4) 
     &	          - fac*3d0*2d0*dreal(wmats0*dconjg(matt0)) *qchint 
 	  else
 	    m2aqq0ss(f1,f2,f3,f4) = m2aqq0ss(f1,f2,f3,f4) 
     &			 + fac*9d0*abs(wmats0)**2 *qch2 
 	  endif

200	continue

c Type: dbar-u -> dbar-u
c-----------------------

c*** generic amplitudes
	do i1=-1,1,2
	do i2=-1,1,2
	  msH4f(i1,i2) = 0d0
	  mtH4f(i1,i2) = 0d0
	enddo
	enddo

c dbar-u -W-> sbar-c    (k3->-p1,k4->-p2,k1->p4,k2->p3)   
	if (qw.eq.1) call M_WW_hvv(msH4f,q4,q3,q1,q2,
     &	                  a1hwws,a2hwws,a3hwws)

c dbar-c -Z-> dbar-c    (k3->-p1,k4->p3,k1->p4,k2->-p2)   
	if (qz.eq.1) call M_ZZ_hvv(mtH4f,guu,gdd,qu,qd,q4,q2,q1,q3,
     &	                  a1haat,a2haat,a3haat,a1hazt,a2hazt,a3hazt,
     &	                  a1hzzt,a2hzzt,a3hzzt)

 	do 205 g1=1,2+qb
 	do 205 g2=1,2
 	do 205 g3=1,2+qb
 	do 205 g4=1,2
 	  f1 = 2*g1-1
 	  f2 = 2*g2
 	  f3 = 2*g3-1
 	  f4 = 2*g4
 	  m2aqq0ss(f1,f2,f3,f4) = 0d0
 	  m2aqq0tt(f1,f2,f3,f4) = 0d0
 	  m2aqq0st(f1,f2,f3,f4) = 0d0

	do 205 i1=-1,1,2
	do 205 i2=-1,1,2
	do 205 i3=-1,1,2
	  mats0  = 0d0
	  matt0  = 0d0

	   if (i1.eq.-i2) mats0  = msH4f(i3,i2)
	   if (i1.eq.-i3) matt0  = mtH4f(i2,i3)

 	  wmats0 = v(g2,g1)*cv(g4,g3)*mats0

 	  if ((f1.eq.f3).and.(f2.eq.f4)) then
 	    m2aqq0ss(f1,f2,f3,f4) = m2aqq0ss(f1,f2,f3,f4) 
     &		         + fac*9d0*abs(wmats0)**2 *qch2 
 	    m2aqq0tt(f1,f2,f3,f4) = m2aqq0tt(f1,f2,f3,f4) 
     &		         + fac*9d0*abs(matt0)**2 *qch2 
 	    m2aqq0st(f1,f2,f3,f4) = m2aqq0st(f1,f2,f3,f4) 
     &	          - fac*3d0*2d0*dreal(wmats0*dconjg(matt0)) *qchint
 	  else
 	    m2aqq0ss(f1,f2,f3,f4) = m2aqq0ss(f1,f2,f3,f4) 
     &			 + fac*9d0*abs(wmats0)**2 *qch2 
 	  endif
205	continue


c Type: ubar-u -> ubar-u
c----------------------- 

c*** generic amplitudes
	do i1=-1,1,2
	do i2=-1,1,2
	  msH4f(i1,i2) = 0d0
	  mtH4f(i1,i2) = 0d0
	enddo
	enddo

c ubar-u -Z-> cbar-c  (k1->-p1,k2->-p2,k3->p4,k4->p3)
	if (qz.eq.1) call M_ZZ_hvv(msH4f,guu,guu,qu,qu,q1,q2,q4,q3,
     &	                  a1haas,a2haas,a3haas,a1hazs,a2hazs,a3hazs,
     &	                  a1hzzs,a2hzzs,a3hzzs)

c ubar-c -Z-> ubar-c  (k1->-p1,k2->p3,k3->p4,k4->-p2)
	if (qz.eq.1) call M_ZZ_hvv(mtH4f,guu,guu,qu,qu,q1,q3,q4,q2,
     &	                  a1haat,a2haat,a3haat,a1hazt,a2hazt,a3hazt,
     &	                  a1hzzt,a2hzzt,a3hzzt)

	m2aqq0ss(2,2,2,2) = 0d0
	m2aqq0tt(2,2,2,2) = 0d0
	m2aqq0st(2,2,2,2) = 0d0
	m2aqq0ss(2,2,4,4) = 0d0
	m2aqq0tt(2,2,4,4) = 0d0
	m2aqq0st(2,2,4,4) = 0d0
	m2aqq0ss(2,4,2,4) = 0d0
	m2aqq0tt(2,4,2,4) = 0d0
	m2aqq0st(2,4,2,4) = 0d0

	do 201 i1=-1,1,2
	do 201 i2=-1,1,2
	do 201 i3=-1,1,2
	  mats0  = 0d0
	  matt0  = 0d0

	  if (i1.eq.-i2) mats0  = msH4f(i2,i3)
	  if (i1.eq.-i3) matt0  = mtH4f(i3,i2)

c*** ubar-u -ZZ-> ubar-u    
	m2aqq0ss(2,2,2,2) = m2aqq0ss(2,2,2,2) 
     &			    + fac*9d0*abs(mats0)**2 *qch2
	m2aqq0tt(2,2,2,2) = m2aqq0tt(2,2,2,2) 
     &			    + fac*9d0*abs(matt0)**2 *qch2 
	m2aqq0st(2,2,2,2) = m2aqq0st(2,2,2,2) 
     &	      - fac*3d0*2d0*dreal(mats0*dconjg(matt0)) *qchint
c*** ubar-u -Z-> cbar-c 
	m2aqq0ss(2,2,4,4) = m2aqq0ss(2,2,4,4) 
     &			    + fac*9d0*abs(mats0)**2 *qch2 
c*** ubar-c -Z-> ubar-c 
	m2aqq0tt(2,4,2,4) = m2aqq0tt(2,4,2,4) 
     &			    + fac*9d0*abs(matt0)**2 *qch2 
201	continue
	m2aqq0ss(4,4,4,4) = m2aqq0ss(2,2,2,2) 
	m2aqq0tt(4,4,4,4) = m2aqq0tt(2,2,2,2) 
	m2aqq0st(4,4,4,4) = m2aqq0st(2,2,2,2) 
	m2aqq0ss(4,4,2,2) = m2aqq0ss(2,2,4,4)
	m2aqq0tt(4,4,2,2) = m2aqq0tt(2,2,4,4)
	m2aqq0st(4,4,2,2) = m2aqq0st(2,2,4,4)
	m2aqq0ss(4,2,4,2) = m2aqq0ss(2,4,2,4)
	m2aqq0tt(4,2,4,2) = m2aqq0tt(2,4,2,4)
	m2aqq0st(4,2,4,2) = m2aqq0st(2,4,2,4)


c Type: ubar-u -> dbar-d
c----------------------- 

c*** generic amplitudes
	do i1=-1,1,2
	do i2=-1,1,2
	  msH4f(i1,i2) = 0d0
	  mtH4f(i1,i2) = 0d0
	enddo
	enddo

c ubar-u -Z-> sbar-s  (k1->-p1,k2->-p2,k3->p4,k4->p3)
	if (qz.eq.1) call M_ZZ_hvv(msH4f,guu,gdd,qu,qd,q1,q2,q4,q3,
     &	                  a1haas,a2haas,a3haas,a1hazs,a2hazs,a3hazs,
     &	                  a1hzzs,a2hzzs,a3hzzs)

c ubar-u -W-> dbar-d  (k1->-p1,k2->p3,k3->p4,k4->-p2)
	if (qw.eq.1) call M_WW_hvv(mtH4f,q1,q3,q4,q2,
     &	                  a1hwwt,a2hwwt,a3hwwt)

 	do 202 g1=1,2
 	do 202 g2=1,2
 	do 202 g3=1,2+qb
 	do 202 g4=1,2+qb
 	  f1 = 2*g1
 	  f2 = 2*g2
 	  f3 = 2*g3-1
 	  f4 = 2*g4-1
 	  m2aqq0ss(f1,f2,f3,f4) = 0d0
 	  m2aqq0tt(f1,f2,f3,f4) = 0d0
 	  m2aqq0st(f1,f2,f3,f4) = 0d0

	do 202 i1=-1,1,2
	do 202 i2=-1,1,2
	do 202 i3=-1,1,2
	  mats0  = 0d0
	  matt0  = 0d0

	  if (i1.eq.-i2) mats0  = msH4f(i2,i3)
	  if (i1.eq.-i3) matt0  = mtH4f(i3,i2)

 	  wmatt0 = cv(g1,g3)*v(g2,g4)*matt0

 	  if ((f1.eq.f2).and.(f3.eq.f4)) then
 	    m2aqq0ss(f1,f2,f3,f4) = m2aqq0ss(f1,f2,f3,f4) 
     &		          + fac*9d0*abs(mats0)**2 *qch2 
 	    m2aqq0tt(f1,f2,f3,f4) = m2aqq0tt(f1,f2,f3,f4) 
     &		          + fac*9d0*abs(wmatt0)**2 *qch2
 	    m2aqq0st(f1,f2,f3,f4) = m2aqq0st(f1,f2,f3,f4) 
     &	          - fac*3d0*2d0*dreal(wmatt0*dconjg(mats0)) *qchint
 	  else
 	    m2aqq0tt(f1,f2,f3,f4) = m2aqq0tt(f1,f2,f3,f4) 
     &			 + fac*9d0*abs(wmatt0)**2 *qch2 
 	  endif
202	  continue


c Type: dbar-d -> dbar-d
c-----------------------

c*** generic amplitudes
	do i1=-1,1,2
	do i2=-1,1,2
	  msH4f(i1,i2) = 0d0
	  mtH4f(i1,i2) = 0d0
	enddo
	enddo

c ubar-u -Z-> cbar-c  (k1->-p1,k2->-p2,k3->p4,k4->p3)
	if (qz.eq.1) call M_ZZ_hvv(msH4f,gdd,gdd,qd,qd,q1,q2,q4,q3,
     &	                  a1haas,a2haas,a3haas,a1hazs,a2hazs,a3hazs,
     &	                  a1hzzs,a2hzzs,a3hzzs)

c ubar-c -Z-> ubar-c  (k1->-p1,k2->p3,k3->p4,k4->-p2)
	if (qz.eq.1) call M_ZZ_hvv(mtH4f,gdd,gdd,qd,qd,q1,q3,q4,q2,
     &	                  a1haat,a2haat,a3haat,a1hazt,a2hazt,a3hazt,
     &	                  a1hzzt,a2hzzt,a3hzzt)

	m2aqq0ss(1,1,1,1) = 0d0
	m2aqq0tt(1,1,1,1) = 0d0
	m2aqq0st(1,1,1,1) = 0d0
	m2aqq0ss(1,1,3,3) = 0d0
	m2aqq0tt(1,1,3,3) = 0d0
	m2aqq0st(1,1,3,3) = 0d0
	m2aqq0ss(1,3,1,3) = 0d0
	m2aqq0tt(1,3,1,3) = 0d0
	m2aqq0st(1,3,1,3) = 0d0

	do 203 i1=-1,1,2
	do 203 i2=-1,1,2
	do 203 i3=-1,1,2
	  mats0  = 0d0
	  matt0  = 0d0

	  if (i1.eq.-i2) mats0  = msH4f(i2,i3)
	  if (i1.eq.-i3) matt0  = mtH4f(i3,i2)

c*** dbar-d -ZZ-> dbar-d    
	m2aqq0ss(1,1,1,1) = m2aqq0ss(1,1,1,1) 
     &			    + fac*9d0*abs(mats0)**2 *qch2
	m2aqq0tt(1,1,1,1) = m2aqq0tt(1,1,1,1) 
     &			    + fac*9d0*abs(matt0)**2 *qch2
	m2aqq0st(1,1,1,1) = m2aqq0st(1,1,1,1) 
     &	       - fac*3d0*2d0*dreal(mats0*dconjg(matt0)) *qchint
c*** dbar-d -Z-> sbar-s 
	m2aqq0ss(1,1,3,3) = m2aqq0ss(1,1,3,3) 
     &			    + fac*9d0*abs(mats0)**2 *qch2 
c*** dbar-s -Z-> dbar-s 
	m2aqq0tt(1,3,1,3) = m2aqq0tt(1,3,1,3) 
     &			    + fac*9d0*abs(matt0)**2 *qch2 
203	continue
	m2aqq0ss(3,3,3,3) = m2aqq0ss(1,1,1,1) 
	m2aqq0tt(3,3,3,3) = m2aqq0tt(1,1,1,1) 
	m2aqq0st(3,3,3,3) = m2aqq0st(1,1,1,1) 
	m2aqq0ss(3,3,1,1) = m2aqq0ss(1,1,3,3)
	m2aqq0tt(3,3,1,1) = m2aqq0tt(1,1,3,3)
	m2aqq0st(3,3,1,1) = m2aqq0st(1,1,3,3)
	m2aqq0ss(3,1,3,1) = m2aqq0ss(1,3,1,3)
	m2aqq0tt(3,1,3,1) = m2aqq0tt(1,3,1,3)
	m2aqq0st(3,1,3,1) = m2aqq0st(1,3,1,3)
	if (qb.eq.1) then
        m2aqq0ss(5,5,5,5) = m2aqq0ss(1,1,1,1)
        m2aqq0tt(5,5,5,5) = m2aqq0tt(1,1,1,1)
        m2aqq0st(5,5,5,5) = m2aqq0st(1,1,1,1)
        m2aqq0ss(1,1,5,5) = m2aqq0ss(1,1,3,3)
        m2aqq0tt(1,1,5,5) = m2aqq0tt(1,1,3,3)
        m2aqq0st(1,1,5,5) = m2aqq0st(1,1,3,3)
        m2aqq0ss(3,3,5,5) = m2aqq0ss(1,1,3,3)
        m2aqq0tt(3,3,5,5) = m2aqq0tt(1,1,3,3)
        m2aqq0st(3,3,5,5) = m2aqq0st(1,1,3,3)
        m2aqq0ss(5,5,1,1) = m2aqq0ss(1,1,3,3)
        m2aqq0tt(5,5,1,1) = m2aqq0tt(1,1,3,3)
        m2aqq0st(5,5,1,1) = m2aqq0st(1,1,3,3)
        m2aqq0ss(5,5,3,3) = m2aqq0ss(1,1,3,3)
        m2aqq0tt(5,5,3,3) = m2aqq0tt(1,1,3,3)
        m2aqq0st(5,5,3,3) = m2aqq0st(1,1,3,3)
        m2aqq0ss(1,5,1,5) = m2aqq0ss(1,3,1,3)
        m2aqq0tt(1,5,1,5) = m2aqq0tt(1,3,1,3)
        m2aqq0st(1,5,1,5) = m2aqq0st(1,3,1,3)
        m2aqq0ss(3,5,3,5) = m2aqq0ss(1,3,1,3)
        m2aqq0tt(3,5,3,5) = m2aqq0tt(1,3,1,3)
        m2aqq0st(3,5,3,5) = m2aqq0st(1,3,1,3)
        m2aqq0ss(5,1,5,1) = m2aqq0ss(1,3,1,3)
        m2aqq0tt(5,1,5,1) = m2aqq0tt(1,3,1,3)
        m2aqq0st(5,1,5,1) = m2aqq0st(1,3,1,3)
        m2aqq0ss(5,3,5,3) = m2aqq0ss(1,3,1,3)
        m2aqq0tt(5,3,5,3) = m2aqq0tt(1,3,1,3)
        m2aqq0st(5,3,5,3) = m2aqq0st(1,3,1,3)
	endif

c Type: dbar-d -> ubar-u
c-----------------------
	do i1=-1,1,2
	do i2=-1,1,2
	  msH4f(i1,i2) = 0d0
	  mtH4f(i1,i2) = 0d0
	enddo
	enddo

c*** generic amplitudes

c dbar-d -Z-> cbar-c  (k3->-p1,k4->-p2,k1->p4,k2->p3)
	if (qz.eq.1) call M_ZZ_hvv(msH4f,guu,gdd,qu,qd,q4,q3,q1,q2,
     &	                  a1haas,a2haas,a3haas,a1hazs,a2hazs,a3hazs,
     &	                  a1hzzs,a2hzzs,a3hzzs)

c dbar-s -W-> ubar-c  (k3->-p1,k4->p3,k1->p4,k2->-p2)
	if (qw.eq.1) call M_WW_hvv(mtH4f,q4,q2,q1,q3,
     &	                  a1hwwt,a2hwwt,a3hwwt)

 	do 204 g1=1,2+qb
 	do 204 g2=1,2+qb
 	do 204 g3=1,2
 	do 204 g4=1,2
 	  f1 = 2*g1-1
 	  f2 = 2*g2-1
 	  f3 = 2*g3
 	  f4 = 2*g4
 	  m2aqq0ss(f1,f2,f3,f4) = 0d0
 	  m2aqq0tt(f1,f2,f3,f4) = 0d0
 	  m2aqq0st(f1,f2,f3,f4) = 0d0

	do 204 i1=-1,1,2
	do 204 i2=-1,1,2
	do 204 i3=-1,1,2
	  mats0  = 0d0
	  matt0  = 0d0

	  if (i1.eq.-i2) mats0  = msH4f(i3,i2)
	  if (i1.eq.-i3) matt0  = mtH4f(i2,i3)

 	  wmatt0 = v(g3,g1)*cv(g4,g2)*matt0

 	  if ((f1.eq.f2).and.(f3.eq.f4)) then
 	    m2aqq0ss(f1,f2,f3,f4) = m2aqq0ss(f1,f2,f3,f4) 
     &		         + fac*9d0*abs(mats0)**2 *qch2
 	    m2aqq0tt(f1,f2,f3,f4) = m2aqq0tt(f1,f2,f3,f4) 
     &		         + fac*9d0*abs(wmatt0)**2 *qch2
 	    m2aqq0st(f1,f2,f3,f4) = m2aqq0st(f1,f2,f3,f4) 
     &	          - fac*3d0*2d0*dreal(wmatt0*dconjg(mats0)) *qchint
 	  else
 	    m2aqq0tt(f1,f2,f3,f4) = m2aqq0tt(f1,f2,f3,f4) 
     &			 + fac*9d0*abs(wmatt0)**2 *qch2 
 	  endif
204	  continue

	end

************************************************************************
        subroutine fixahvv(q1,q2,
     &             a1hwwf,a2hwwf,a3hwwf,a1haaf,a2haaf,a3haaf,
     &             a1hazf,a2hazf,a3hazf,a1hzzf,a2hzzf,a3hzzf)
************************************************************************
*       rescale anomalous HVV couplings by form factor
*	q1,q2 = virtualities of vector bosons
*-----------------------------------------------------------------------
*       14.6.12 Stefan Dittmaier
************************************************************************
        implicit real*8 (a-z)
	complex*16 v,cv
	integer qhvv

        common/param/pi,el,alpha,alpha0,alphaz,GF,alphas,
     &         v(3,3),cv(3,3),cw,cw2,sw,sw2,mw,mw2,gw,mz,mz2,gz,mh,mh2,
     &         ml(3),ml2(3),mqp(3),mqp2(3),mqm(3),mqm2(3)
        common/hvv/rsm,d,db,dt,dtb,lambdahvv,
     &             a1hww,a2hww,a3hww,a1haa,a2haa,a3haa,
     &             a1haz,a2haz,a3haz,a1hzz,a2hzz,a3hzz,qhvv

	if (lambdahvv.gt.0d0) then
	  formfac_ww = lambdahvv**4/(lambdahvv**2+abs(q1))
     &	                           /(lambdahvv**2+abs(q2))
	else
	  formfac_ww = 1d0
	endif

	formfac_zz = formfac_ww
	formfac_aa = formfac_ww 
	formfac_az = formfac_ww 

	a1hwwf = a1hww 
	a2hwwf = a2hww * formfac_ww
	a3hwwf = a3hww * formfac_ww
	a1haaf = a1haa 
	a2haaf = a2haa * formfac_aa
	a3haaf = a3haa * formfac_aa
	a1hazf = a1haz 
	a2hazf = a2haz * formfac_az
	a3hazf = a3haz * formfac_az
	a1hzzf = a1hzz 
	a2hzzf = a2hzz * formfac_zz
	a3hzzf = a3hzz * formfac_zz

	end

************************************************************************
        subroutine M_WW_hvv(mH4f,fa,fb,fc,fd,a1,a2,a3)
************************************************************************
*       Helicity amplitudes for
*       H -> WW -> fa + fb-bar + fc + fd-bar 
*       with anomalous HVV couplings
*-----------------------------------------------------------------------
*       14.6.12 Stefan Dittmaier
************************************************************************
        implicit real*8 (a-z)
	complex*16 sp(-4:4,-4:4),csp(-4:4,-4:4)
        complex*16 ccw,ccw2,csw,csw2,cmw,cmw2,cmz,cmz2,cmh,cmh2
        complex*16 pwab,pwcd
        complex*16 Amm,mH4f(-1:1,-1:1)
        complex*16 fac
        complex*16 v,cv
        complex*16 xcw,xcw2,xsw,xsw2,xmw,xmw2,xmz,xmz2,xmh,xmh2,null
        complex*16 xml,xml2,xmqp,xmqp2,xmqm,xmqm2
        real*8 vp(-4:4,-4:4)
        integer fa,fb,fc,fd,ia,ib,ic,id,i1,i2,i3
        integer qborn,qw,qz,qschan,qtchan,qch2,qchint,
     &                   qbini,qbfin,qwidth,qfact,qbos,qferm,qsoft,qhh2,
     &                   qqcddiag,qqcdnondiag,qqcdgsplit,qqcdggfus,qcp
        logical check

        common/prods3/sp,csp,vp
        common/cparam/ccw,ccw2,csw,csw2,cmw,cmw2,cmz,cmz2
        common/xparam/xcw,xcw2,xsw,xsw2,xmw,xmw2,xmz,xmz2,xmh,xmh2,null,
     &       xml(3),xml2(3),xmqp(3),xmqp2(3),xmqm(3),xmqm2(3)
        common/param/pi,el,alpha,alpha0,alphaz,GF,alphas,
     &         v(3,3),cv(3,3),cw,cw2,sw,sw2,mw,mw2,gw,mz,mz2,gz,mh,mh2,
     &         ml(3),ml2(3),mqp(3),mqp2(3),mqm(3),mqm2(3)
        common/rcoptions/qborn,qw,qz,qschan,qtchan,qch2,qchint,
     &                   qbini,qbfin,qwidth,qfact,qbos,qferm,qsoft,qhh2,
     &                   qqcddiag,qqcdnondiag,qqcdgsplit,qqcdggfus,qcp

	Amm(ia,ib,ic,id,a1,a2,a3) = 2d0*a1*sp(ia,ic)*csp(ib,id)
     &	  + dcmplx(a2,+a3)*sp(ic,ia)*csp(ib,ia)*sp(ia,ic)*csp(id,ic)
     &	  + dcmplx(a2,-a3)*sp(ia,ib)*csp(id,ib)*sp(ic,id)*csp(ib,id)

        sab  = vp(fa,fb)
        scd  = vp(fc,fd)

        if ((qwidth.eq.0).or.((qwidth.eq.2).and.(sab.lt.0d0))) then
          pwab = 1d0/(sab-mw2)
        else
          pwab = 1d0/(sab-cmw2)
        endif
        if ((qwidth.eq.0).or.((qwidth.eq.2).and.(scd.lt.0d0))) then
          pwcd = 1d0/(scd-mw2)
        else
          pwcd = 1d0/(scd-cmw2)
        endif

        check = ((qschan.eq.1).and.(sab.gt.0d0).and.(scd.gt.0d0))
     &     .or. ((qtchan.eq.1).and.(sab.lt.0d0).and.(scd.lt.0d0))

        if (.not.check) pwab = 0d0

	mH4f(+1,+1) = 0d0
	mH4f(+1,-1) = 0d0
	mH4f(-1,+1) = 0d0
	mH4f(-1,-1) = pwab*pwcd/2d0/xsw2*Amm(fa,fb,fc,fd,a1,a2,a3)

	end

************************************************************************
        subroutine M_ZZ_hvv(mH4f,gaz,gcz,qa,qc,fa,fb,fc,fd,
     &	             a1aa,a2aa,a3aa,a1az,a2az,a3az,a1zz,a2zz,a3zz)
************************************************************************
*       Helicity amplitudes for
*       H -> AA/AZ/ZZ -> fa + fb-bar + fc + fd-bar 
*       with anomalous HVV couplings
*-----------------------------------------------------------------------
*       14.6.12 Stefan Dittmaier
************************************************************************
        implicit real*8 (a-z)
	complex*16 sp(-4:4,-4:4),csp(-4:4,-4:4)
        complex*16 ccw,ccw2,csw,csw2,cmw,cmw2,cmz,cmz2,cmh,cmh2
        complex*16 pzab,pzcd
        complex*16 App,Apm,Amp,Amm,mH4f(-1:1,-1:1)
        complex*16 fac
        complex*16 gaz(-1:1),gcz(-1:1)
        complex*16 v,cv
        complex*16 xcw,xcw2,xsw,xsw2,xmw,xmw2,xmz,xmz2,xmh,xmh2,null
        complex*16 xml,xml2,xmqp,xmqp2,xmqm,xmqm2
        real*8 vp(-4:4,-4:4)
        integer fa,fb,fc,fd,ia,ib,ic,id,i1,i2,i3
        integer qborn,qw,qz,qschan,qtchan,qch2,qchint,
     &                   qbini,qbfin,qwidth,qfact,qbos,qferm,qsoft,qhh2,
     &                   qqcddiag,qqcdnondiag,qqcdgsplit,qqcdggfus,qcp
        logical check

        common/prods3/sp,csp,vp
        common/cparam/ccw,ccw2,csw,csw2,cmw,cmw2,cmz,cmz2
        common/xparam/xcw,xcw2,xsw,xsw2,xmw,xmw2,xmz,xmz2,xmh,xmh2,null,
     &       xml(3),xml2(3),xmqp(3),xmqp2(3),xmqm(3),xmqm2(3)
        common/param/pi,el,alpha,alpha0,alphaz,GF,alphas,
     &         v(3,3),cv(3,3),cw,cw2,sw,sw2,mw,mw2,gw,mz,mz2,gz,mh,mh2,
     &         ml(3),ml2(3),mqp(3),mqp2(3),mqm(3),mqm2(3)
        common/rcoptions/qborn,qw,qz,qschan,qtchan,qch2,qchint,
     &                   qbini,qbfin,qwidth,qfact,qbos,qferm,qsoft,qhh2,
     &                   qqcddiag,qqcdnondiag,qqcdgsplit,qqcdggfus,qcp

	App(ia,ib,ic,id,a1,a2,a3) = 2d0*a1*csp(ia,ic)*sp(ib,id)
     &	  + dcmplx(a2,-a3)*csp(ic,ia)*sp(ib,ia)*csp(ia,ic)*sp(id,ic)
     &	  + dcmplx(a2,+a3)*csp(ia,ib)*sp(id,ib)*csp(ic,id)*sp(ib,id)

	Apm(ia,ib,ic,id,a1,a2,a3) = 2d0*a1*csp(ia,id)*sp(ib,ic)
     &	  + dcmplx(a2,-a3)*csp(id,ia)*sp(ib,ia)*csp(ia,id)*sp(ic,id)
     &	  + dcmplx(a2,+a3)*csp(ia,ib)*sp(ic,ib)*csp(id,ic)*sp(ib,ic)

	Amm(ia,ib,ic,id,a1,a2,a3) = 2d0*a1*sp(ia,ic)*csp(ib,id)
     &	  + dcmplx(a2,+a3)*sp(ic,ia)*csp(ib,ia)*sp(ia,ic)*csp(id,ic)
     &	  + dcmplx(a2,-a3)*sp(ia,ib)*csp(id,ib)*sp(ic,id)*csp(ib,id)

	Amp(ia,ib,ic,id,a1,a2,a3) = 2d0*a1*sp(ia,id)*csp(ib,ic)
     &	  + dcmplx(a2,+a3)*sp(id,ia)*csp(ib,ia)*sp(ia,id)*csp(ic,id)
     &	  + dcmplx(a2,-a3)*sp(ia,ib)*csp(ic,ib)*sp(id,ic)*csp(ib,ic)

        sab  = vp(fa,fb)
        scd  = vp(fc,fd)

        paab = 1d0/sab
        pacd = 1d0/scd
        if ((qwidth.eq.0).or.((qwidth.eq.2).and.(sab.lt.0d0))) then
          pzab = 1d0/(sab-mz2)
        else
          pzab = 1d0/(sab-cmz2)
        endif
        if ((qwidth.eq.0).or.((qwidth.eq.2).and.(scd.lt.0d0))) then
          pzcd = 1d0/(scd-mz2)
        else
          pzcd = 1d0/(scd-cmz2)
        endif

        check = ((qschan.eq.1).and.(sab.gt.0d0).and.(scd.gt.0d0))
     &     .or. ((qtchan.eq.1).and.(sab.lt.0d0).and.(scd.lt.0d0))

        if (.not.check) paab = 0d0
        if (.not.check) pzab = 0d0

	mH4f(+1,+1) = ((-qa)*(-qc)*paab*pacd)
     &	              *App(fa,fb,fc,fd,a1aa,a2aa,a3aa)
     &	             +((-qa)*gcz(+1)*paab*pzcd+gaz(+1)*(-qc)*pzab*pacd)
     &	              *App(fa,fb,fc,fd,a1az,a2az,a3az)
     &	             +gaz(+1)*gcz(+1)*pzab*pzcd
     &	              *App(fa,fb,fc,fd,a1zz,a2zz,a3zz)	

	mH4f(+1,-1) = ((-qa)*(-qc)*paab*pacd)
     &	              *Apm(fa,fb,fc,fd,a1aa,a2aa,a3aa)
     &	             +((-qa)*gcz(-1)*paab*pzcd+gaz(+1)*(-qc)*pzab*pacd)
     &	              *Apm(fa,fb,fc,fd,a1az,a2az,a3az)
     &	             +gaz(+1)*gcz(-1)*pzab*pzcd
     &	              *Apm(fa,fb,fc,fd,a1zz,a2zz,a3zz)	

	mH4f(-1,+1) = ((-qa)*(-qc)*paab*pacd)
     &	              *Amp(fa,fb,fc,fd,a1aa,a2aa,a3aa)
     &	             +((-qa)*gcz(+1)*paab*pzcd+gaz(-1)*(-qc)*pzab*pacd)
     &	              *Amp(fa,fb,fc,fd,a1az,a2az,a3az)
     &	             +gaz(-1)*gcz(+1)*pzab*pzcd
     &	              *Amp(fa,fb,fc,fd,a1zz,a2zz,a3zz)	

	mH4f(-1,-1) = ((-qa)*(-qc)*paab*pacd)
     &	              *Amm(fa,fb,fc,fd,a1aa,a2aa,a3aa)
     &	             +((-qa)*gcz(-1)*paab*pzcd+gaz(-1)*(-qc)*pzab*pacd)
     &	              *Amm(fa,fb,fc,fd,a1az,a2az,a3az)
     &	             +gaz(-1)*gcz(-1)*pzab*pzcd
     &	              *Amm(fa,fb,fc,fd,a1zz,a2zz,a3zz)	

	end
       subroutine setprods3(p)
************************************************************************                                                                                   
*       Weyl-van der Waerden and Minkowski products                                                                                                        
*-----------------------------------------------------------------------                                                                                   
*       23.8.06 Stefan Dittmaier                                                                                                                           
************************************************************************                                                                                   
        implicit real*8 (a-z)
        complex*16 eiph1,eiph2,sp(-4:4,-4:4),csp(-4:4,-4:4),v,cv
        real*8 p(6,0:3),vp(-4:4,-4:4)
        integer i,j

        common/prods3/sp,csp,vp
        common/param/pi,el,alpha,alpha0,alphaz,GF,alphas,
     &         v(3,3),cv(3,3),cw,cw2,sw,sw2,mw,mw2,gw,mz,mz2,gz,mh,mh2,
     &         ml(3),ml2(3),mqp(3),mqp2(3),mqm(3),mqm2(3)

c energies and angles                                                                                                                                      
        eb1  = p(1,0)
        eb2  = p(2,0)
        e1   = p(3,0)
        e2   = p(4,0)
        cth1 = p(3,3)/p(3,0)
        cth2 = p(4,3)/p(4,0)
        sth1 = sqrt(p(3,1)**2+p(3,2)**2)/p(3,0)
        sth2 = sqrt(p(4,1)**2+p(4,2)**2)/p(4,0)
        eiph1= dcmplx(p(3,1),p(3,2))/p(3,0)/sth1
        eiph2= dcmplx(p(4,1),p(4,2))/p(4,0)/sth2
        if (p(1,3).lt.0d0) then
          cth1 = -cth1
          cth2 = -cth2
          eiph1= dconjg(eiph1)
          eiph2= dconjg(eiph2)
        endif
        r2eb1= sqrt(2d0*eb1)
        r2eb2= sqrt(2d0*eb2)
        r2e1 = sqrt(2d0*e1)
        r2e2 = sqrt(2d0*e2)
        if (cth1.gt.0d0) then
          c2th1 = sqrt((1d0+cth1)/2d0)
          s2th1 = sth1/2d0/c2th1
        else
          s2th1 = sqrt((1d0-cth1)/2d0)
          c2th1 = sth1/2d0/s2th1
        endif
        if (cth2.gt.0d0) then
          c2th2 = sqrt((1d0+cth2)/2d0)
          s2th2 = sth2/2d0/c2th2
        else
          s2th2 = sqrt((1d0-cth2)/2d0)
          c2th2 = sth2/2d0/s2th2
        endif

c spinor products                                                                                                                                          
        sp(1,2) = r2eb1*r2eb2
        sp(1,3) = r2eb1*r2e1*(s2th1)
        sp(1,4) = r2eb1*r2e2*(s2th2)
        sp(2,3) = r2eb2*r2e1*(-1d0/eiph1*c2th1)
        sp(2,4) = r2eb2*r2e2*(-1d0/eiph2*c2th2)
        sp(3,4) = r2e1 *r2e2*(c2th1/eiph1*s2th2-c2th2/eiph2*s2th1)

        do 101 i=1,3
        do 101 j=i+1,4
           sp(j,i) = -sp(i,j)
101     continue

        do 100 i=1,4
        do 100 j=1,4
          csp(i,j)  = dconjg(sp(i,j))
           sp(-i, j) =  sp(i,j)
           sp( i,-j) =  sp(i,j)
           sp(-i,-j) =  sp(i,j)
          csp(-i, j) =-csp(i,j)
          csp( i,-j) =-csp(i,j)
          csp(-i,-j) = csp(i,j)
           vp( i, j) =  sp(i,j)*csp(i,j)
           vp(-i, j) = -vp(i,j)
           vp( i,-j) = -vp(i,j)
           vp(-i,-j) =  vp(i,j)
100     continue

        end

C*********************************************************************************************************
C pdfstuff
      subroutine getpdf(x,mfact,pdf)
C our interface to pdfs
      implicit none
      real*8 x, mfact, mfact2, pdf(-6,6)
      mfact2=mfact*mfact 
      call pdgrv(1,x,mfact2,pdf)
      return
      end
C
      SUBROUTINE PDGRV(ISET,X,Q2,XPDF)
 
C...ISET = 1 - LO,          Lambda_4=0.20 GeV, N_f=5
C...       2 - NLO, MS_bar, Lambda_4=0.20 GeV, N_f=5
C...X          - Bjorken x
C...Q2         - square of the momentum scale  (in GeV**2)
C...XPDF(-6:6) - matrix containing  x*p(x,Q2)
C...     IPDF = -6 ,  -5 ,  -4 ,  -3 ,  -2 ,  -1 ,0 ,1,2,3,4,5,6
C...          t_bar,b_bar,c_bar,s_bar,u_bar,d_bar,gl,d,u,s,c,b,t
C...range of validity:
C...     D-5  < X  < 1
C...      0.3 < Q2 < D8  GeV^2
C...REAL*8 version
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XPDF(-6:6)
 
C...LO PARAMETRIZATIONS :
      IF (ISET.EQ.1) THEN
        AMU2  = 0.25
        ALAM2 = 0.232 * 0.232
        S  = LOG (LOG(Q2/ALAM2)/LOG(AMU2/ALAM2))
        S2 = S * S
        S3 = S2 * S
 
C...X * (UV + DV) :
        DNUD = 0.663 + 0.191 * S - 0.041 * S2 + 0.031 * S3
        AKUD = 0.326
        AGUD = -1.97 +  6.74 * S -  1.96 * S2
        BUD  =  24.4 -  20.7 * S +  4.08 * S2
        DUD  =  2.86 +  0.70 * S -  0.02 * S2
        UDV  = PDFV (X, DNUD, AKUD, AGUD, BUD, DUD)
 
C...X * DV :
        DND  = 0.579 + 0.283 * S + 0.047 * S2
        AKD = 0.523 - 0.015 * S
        AGD =  2.22 -  0.59 * S -  0.27 * S2
        BD  =  5.95 -  6.19 * S +  1.55 * S2
        DD  =  3.57 +  0.94 * S -  0.16 * S2
        DV  = PDFV (X, DND, AKD, AGD, BD, DD)
 
C...X * G :
        ALG =  0.558
        BEG =  1.218
        AKG =   1.00 -  0.17 * S
        BKG =   0.0
        AGG =   0.0  + 4.879 * S - 1.383 * S2
        BGG =  25.92 - 28.97 * S + 5.596 * S2
        CG  = -25.69 + 23.68 * S - 1.975 * S2
        DG  =  2.537 + 1.718 * S + 0.353 * S2
        EG  =  0.595 + 2.138 * S
        ESG =  4.066
        GL =PDFW (X, S, ALG, BEG, AKG, BKG, AGG, BGG, CG, DG, EG, ESG)
 
C...X * UBAR = X * DBAR :
        ALU =  1.396
        BEU =  1.331
        AKU =  0.412 - 0.171 * S
        BKU =  0.566 - 0.496 * S
        AGU =  0.363
        BGU = -1.196
        CU  =  1.029 + 1.785 * S - 0.459 * S2
        DU  =  4.696 + 2.109 * S
        EU  =  3.838 + 1.944 * S
        ESU =  2.845
        UDB=PDFW (X, S, ALU, BEU, AKU, BKU, AGU, BGU, CU, DU, EU, ESU)
 
C...X * SBAR = X * S :
        SS  =   0.0
        ALS =  0.803
        BES =  0.563
        AKS =  2.082 - 0.577 * S
        AGS = -3.055 + 1.024 * S **  0.67
        BS  =   27.4 -  20.0 * S ** 0.154
        DS  =   6.22
        EST =   4.33 + 1.408 * S
        ESS =   8.27 - 0.437 * S
        SB =PDFWS (X, S, SS, ALS, BES, AKS, AGS, BS, DS, EST, ESS)
 
C...X * CBAR = X * C :
        SC  =  0.888
        ALC =   1.01
        BEC =   0.37
        AKC =   0.0
        AGC =   0.0
        BC  =   4.24 - 0.804 * S
        DC  =   3.46 + 1.076 * S
        EC  =   4.61 + 1.490 * S
        ESC =  2.555 + 1.961 * S
        CB =PDFWS (X, S, SC, ALC, BEC, AKC, AGC, BC, DC, EC, ESC)
 
C...X * BBAR = X * B :
        SBO =  1.351
        ALB =   1.00
        BEB =   0.51
        AKB =   0.0
        AGB =   0.0
        BBO =  1.848
        DB  =  2.929 + 1.396 * S
        EB  =   4.71 + 1.514 * S
        ESB =   4.02 + 1.239 * S
        BB =PDFWS (X, S, SBO, ALB, BEB, AKB, AGB, BBO, DB, EB, ESB)
 
C...HO parametrization:
      ELSEIF(ISET.EQ.2) THEN
        AMU2  = 0.3
        ALAM2 = 0.248 * 0.248
        S  = LOG (LOG(Q2/ALAM2)/LOG(AMU2/ALAM2))
        DS = SQRT (S)
        S2 = S * S
        S3 = S2 * S
 
C...X * (UV + DV) :
        DNUD  = 0.330 + 0.151 * S - 0.059 * S2 + 0.027 * S3
        AKUD = 0.285
        AGUD = -2.28 + 15.73 * S -  4.58 * S2
        BUD  =  56.7 -  53.6 * S + 11.21 * S2
        DUD  =  3.17 +  1.17 * S -  0.47 * S2 +  0.09 * S3
        UDV  = PDFV (X, DNUD, AKUD, AGUD, BUD, DUD)
 
C...X * DV :
        DND  = 0.459 + 0.315 * DS + 0.515 * S
        AKD = 0.624              - 0.031 * S
        AGD =  8.13 -  6.77 * DS +  0.46 * S
        BD  =  6.59 - 12.83 * DS +  5.65 * S
        DD  =  3.98              +  1.04 * S  -  0.34 * S2
        DV  = PDFV (X, DND, AKD, AGD, BD, DD)
 
C...X * G :
        ALG =  1.128
        BEG =  1.575
        AKG =  0.323 + 1.653 * S
        BKG =  0.811 + 2.044 * S
        AGG =   0.0  + 1.963 * S - 0.519 * S2
        BGG =  0.078 +  6.24 * S
        CG  =  30.77 - 24.19 * S
        DG  =  3.188 + 0.720 * S
        EG  = -0.881 + 2.687 * S
        ESG =  2.466
        GL =PDFW (X, S, ALG, BEG, AKG, BKG, AGG, BGG, CG, DG, EG, ESG)
 
C...X * UBAR = X * DBAR :
        ALU =  0.594
        BEU =  0.614
        AKU =  0.636 - 0.084 * S
        BKU =   0.0
        AGU =  1.121 - 0.193 * S
        BGU =  0.751 - 0.785 * S
        CU  =   8.57 - 1.763 * S
        DU  =  10.22 + 0.668 * S
        EU  =  3.784 + 1.280 * S
        ESU =  1.808 + 0.980 * S
        UDB=PDFW (X, S, ALU, BEU, AKU, BKU, AGU, BGU, CU, DU, EU, ESU)
 
C...X * SBAR = X * S :
        SS  =   0.0
        ALS =  0.756
        BES =  0.101
        AKS =  2.942 - 1.016 * S
        AGS =  -4.60 + 1.167 * S
        BS  =   9.31 - 1.324 * S
        DS  =  11.49 - 1.198 * S + 0.053 * S2
        EST =  2.630 + 1.729 * S
        ESS =   8.12
        SB =PDFWS (X, S, SS, ALS, BES, AKS, AGS, BS, DS, EST, ESS)
 
C...X * CBAR = X * C :
        SC  =  0.820
        ALC =   0.98
        BEC =   0.0
        AKC = -0.625 - 0.523 * S
        AGC =   0.0
        BC  =  1.896 + 1.616 * S
        DC  =   4.12 + 0.683 * S
        EC  =   4.36 + 1.328 * S
        ESC =  0.677 + 0.679 * S
        CB =PDFWS (X, S, SC, ALC, BEC, AKC, AGC, BC, DC, EC, ESC)
 
C...X * BBAR = X * B :
        SBO =  1.297
        ALB =   0.99
        BEB =   0.0
        AKB =   0.0  - 0.193 * S
        AGB =   0.0
        BBO =   0.0
        DB  =  3.447 + 0.927 * S
        EB  =   4.68 + 1.259 * S
        ESB =  1.892 + 2.199 * S
        BB =PDFWS (X, S, SBO, ALB, BEB, AKB, AGB, BBO, DB, EB, ESB)
      ELSE
       WRITE(*,*) ' error in PDGRV: wrong ISET value'
      ENDIF
 
C...final results
      XPDF(0)=GL
      XPDF(1)=DV+UDB
      XPDF(2)=UDV-DV+UDB
      XPDF(3)=SB
      XPDF(4)=CB
      XPDF(5)=BB
      XPDF(6)=0.
      XPDF(-1)=UDB
      XPDF(-2)=UDB
      XPDF(-3)=SB
      XPDF(-4)=CB
      XPDF(-5)=BB
      XPDF(-6)=0.
 
      RETURN
      END
      FUNCTION PDFV (X, DN, AK, AG, B, D)
 
C...functional forms for ho and lo parametrizations :
      IMPLICIT REAL*8 (A-H,O-Z)
       DX = SQRT (X)
       PDFV = DN * X**AK * (1.+ AG*DX + B*X) * (1.- X)**D
      RETURN
      END
C
      FUNCTION PDFW (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
      IMPLICIT REAL*8 (A-H,O-Z)
       ALX = LOG (1./X)
       PDFW = (X**AK * (AG + X * (BG + X*C)) *ALX**BK + S**AL
     1      * EXP (-E + SQRT (ES * S**BE *ALX))) * (1.- X)**D
      RETURN
      END
C-----------------------------------------------------
      FUNCTION PDFWS (X, S, ST, AL, BE, AK, AG, B, D, E, ES)
      IMPLICIT REAL*8 (A-H,O-Z)
       DX = SQRT (X)
       ALX = LOG (1./X)
       IF (S .LE. ST) THEN
         FWS = 0.0
       ELSE
         FWS = (S-ST)**AL / ALX**AK * (1.+ AG*DX + B*X) * (1.- X)**D
     1          * EXP (-E + SQRT (ES * S**BE *ALX))
       ENDIF
       PDFWS=FWS
      RETURN
      END
      
