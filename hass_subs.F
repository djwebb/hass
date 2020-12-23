!
!=======================================================================
!23456789012345678901234567890123456789012345678901234567890123456789012
!
!  File: hass_subs.F
!
! ***  Version 2.1   (December 2020)
! ***  For unix systems running Fortran 90.
! ***
! ***  Copyright 1978, 2012, 2020 D.J. Webb
! ***  Released under licence GPL-3.0-or-later
! ***
!=======================================================================
!
      function hass_f(kx, ky, pcon, fmax, angle)
!
!23456789012345678901234567890123456789012345678901234567890123456789012
!=======================================================================
!
!  Function to calculate Hassleman Integral for a Pierson-Moskowitz
!  spectrum of sea waves.
!
!  Input:
!    kx, ky   : x and y components of wavenumber of calculate
!    pcon     ; Maximum amplitude (Phillips constant) of spectrum.
!    fmax     : Frequency (Hz) at which maximum occurs.
!    angle    : angle of wave field (angle in radians and
!               anti-clockwise relative to x-axis).
!
!  Output:
!
!    hass_f   : rate of change of action at wavenumber (kx,ky).
!
!  Notes:
!
!    1. The action equals e(kx,ky)/w(kx,ky) where:
!          w is the angular velocity of wave (kx,ky)
!          e is the wave height variance per unit area of wavenumber space.
!    2.  All variables are 64-bit real quantities.
!
!=======================================================================
!
      use iso_fortran_env
      implicit none
!
      real(real64)            :: hass_f
      real(real64),intent(in) :: kx, ky, pcon, fmax, angle
!
      real(real64) :: kx1,ky1,kx2,ky2,kx3,ky3,kx4,ky4,k1,k2,k3,k4,
     &                w1,w2,w3,w4,f1,f2,f3,f4,rel1,rel2,rel3,
     &                abs1,abs2,abs3,stp1,stp2,stp3,scalei,rkmax,
     &                th,al,fm,gam,rhoa,rhob,swx,swy,
     &                sww,swh
      real(real64) :: pi,pi2,radian,degree,g,den
      integer      :: iprint, nextl, icall, i, j
!
      real(real64) :: dummy,x1,y1,epdf,pump,diff,smof,time,tim,tot,tt
      real(real64) :: ans(4)
      logical      :: logi,logp
      real(real64) ::  acc3, ss, cc
!
      common /c_wav/kx1,ky1,kx2,ky2,kx3,ky3,kx4,ky4,k1,k2,k3,k4,
     &       w1,w2,w3,w4,f1,f2,f3,f4
      common /c_int/rel1,rel2,rel3,abs1,abs2,abs3,stp1,stp2,stp3,
     &       scalei,rkmax,th
      common /c_con/ pi,pi2,radian,degree,g,den
      common /jonswp/al,fm,gam,rhoa,rhob,swx,swy,sww,swh
      common /c_test/iprint,icall
!
!=======================================================================
!
      pi  = 4d0*atan(1d0)
      pi2 = pi+pi
      radian = pi/180d0
      degree = 180d0/pi
      g = 9.81
      den = 1d3
      iprint = 0
      icall  = 0
!
      rel1 = 1d-2
      abs1 = 1d-7
      stp1 = 1d-3
      rel2 = 3d-3
      abs2 = 5d-8
      stp2 = 1d-3
      rel3 = 1d-3
      abs3 = 5d-8
      stp3 = 1d-3
!
      scalei = 0.3
      rkmax = 50d0
!
      al  = pcon
      fm  = fmax
!  Johnswap constants
      gam  = 1.0d0
      rhoa = 0.07d0
      rhob = 0.09d0
!  Swell option
      swx = 0d0
      swy = 0d0
      swh = 0d0
      sww = 0.01
!
!  Test input
!
      hass_f = 0d0
      if(kx.eq.0d0.and.ky.eq.0d0)return
      if(pcon.eq.0d0.or.fmax.eq.0d0)return
!
!  Rotate (kx,ky) by angle (this avoids modifying function snpm)
!
      cc = cos(angle)
      ss = sin(angle)
      kx1 = cc*kx + ss*ky
      ky1 = cc*ky - ss*kx
!
!  Call main routine
!
      call hass(ans,kx1,ky1)
!
      hass_f = ans(1)
      return
      end function hass_f

!
!  The remainder of this file consists of the following routines
!
!  hass   -  level 1 - calls quad41 to integrate over angle of k3
!  calc1  -  level 2 - calls quad41 to integrate over k3 radially
!  calc2  -  level 3 - calls quad41 to integrate around loop (k2 or k4)
!  calc3  -  integrand at a single point on the loop
!  quad41 -  integration routine (re-entrant)
!  snpm   -  wave spectrum routine
!

      subroutine hass(ans,x1,y1)
      use iso_fortran_env
      implicit none
!
!=======================================================================
!
!
!=======================================================================
!
      real(real64),intent(in)  ::x1,y1
      real(real64),intent(out) ::ans(4)
!
      real(real64) ::kx1,ky1,kx2,ky2,kx3,ky3,kx4,ky4,k1,k2,k3,k4,
     &               w1,w2,w3,w4,f1,f2,f3,f4
      real(real64) ::rel1,rel2,rel3,abs1,abs2,abs3,stp1,stp2,stp3,
     &               scalei,rkmax,th
      real(real64) ::acc1
      real(real64) ::pi,pi2,radian,degree,g,den
      external calc1
!
      common /c_wav/kx1,ky1,kx2,ky2,kx3,ky3,kx4,ky4,k1,k2,k3,k4,
     &       w1,w2,w3,w4,f1,f2,f3,f4
      common /c_int/rel1,rel2,rel3,abs1,abs2,abs3,stp1,stp2,stp3,
     &       scalei,rkmax,th
      common /c_con/ pi,pi2,radian,degree,g,den
!
      kx1 = x1
      ky1 = y1
      call quad41(ans,0d0,pi2,rel1,abs1,acc1,stp1,calc1)
      return
      end subroutine hass


      subroutine calc1(z,ans)
      use iso_fortran_env
      implicit none
!
!=======================================================================
!
!
!=======================================================================
!
      real(real64),intent(in)  ::z
      real(real64),intent(out) ::ans(4)
!
      real(real64) ::rel1,rel2,rel3,abs1,abs2,abs3,stp1,stp2,stp3,
     &               scalei,rkmax,th
      real(real64) ::pi,pi2,radian,degree,g,den
      integer      ::iprint,icall
!
      real(real64) ::acc2
      logical logp3
      external :: calc2
!
      common /c_int/rel1,rel2,rel3,abs1,abs2,abs3,stp1,stp2,stp3,
     &       scalei,rkmax,th
      common /c_con/ pi,pi2,radian,degree,g,den
      common /c_test/iprint,icall
!
      logp3 = iprint.ge.3
      if(logp3)print 21," Enter calc1 with z = ",z
  21  format(a,2f10.5)
      th=z
      call quad41(ans,0d0,1d0,rel2,abs2,acc2,stp2,calc2)
      if(logp3)write(6,7)z,ans
   7  format(5x,1pe12.4,5x,4e12.4)
      return
      end subroutine calc1


      subroutine calc2(z,ans)
      use iso_fortran_env
      implicit none
!
!=======================================================================
!
!
!=======================================================================
!
      real(real64),intent(in)  ::z
      real(real64),intent(out) ::ans(4)
!
      real(real64) :: r,d,tst4(4)
      real(real64) :: acc3
      logical      :: logi, logp5, logp6
      external :: calc3
!
      real(real64) ::kx1,ky1,kx2,ky2,kx3,ky3,kx4,ky4,k1,k2,k3,k4,
     &               w1,w2,w3,w4,f1,f2,f3,f4
      real(real64) ::rel1,rel2,rel3,abs1,abs2,abs3,stp1,stp2,stp3,
     &               scalei,rkmax,th
      real(real64) ::pi,pi2,radian,degree,g,den
      integer      ::iprint,icall
!
      common /c_wav/kx1,ky1,kx2,ky2,kx3,ky3,kx4,ky4,k1,k2,k3,k4,
     &       w1,w2,w3,w4,f1,f2,f3,f4
      common /c_int/rel1,rel2,rel3,abs1,abs2,abs3,stp1,stp2,stp3,
     &       scalei,rkmax,th
      common /c_con/ pi,pi2,radian,degree,g,den
      common /c_test/iprint,icall
!
      logp5 = iprint.ge.5
      logp6 = iprint.ge.6
!
      if(logp5)print 21,"  Enter calc2 with z = ",z
  21  format(a,6(f12.4,1x))
      if(z.le.0d0.or.z.ge.1d0)goto 200
      r=scalei*z/(1d0-z)
      d=scalei*r/(1d0-z)**2
      kx3 = kx1+r*cos(th)
      ky3 = ky1+r*sin(th)
      k3 = sqrt(kx3**2+ky3**2)
!  Set integral to zero if k3 is very large (>= rkmax)
      if(k3.ge.rkmax)goto 200
!
      if(logp6)print 21,"       calc2.  kx1, ky1 = ",kx1,ky1,r,d
      if(logp6)print 21,"               kx3, ky3 = ",kx3,ky3
!
      call set(logi)
      if(logp6)print *, "       calc2.  logi     = ",logi
      if(logi)then
        call quad41(ans,-1d0,1d0,rel3,abs3,acc3,stp3,calc3)
        ans = ans*d
      else
        ans = 0d0
      endif
!
!  If clean return - carry out integral
!
      if(logp5)print 21,"  Exit calc2. ans = ",ans
      return
!
 200  ans = 0d0
      return
      end subroutine calc2


      recursive subroutine quad41(ANSWER,ALIM,BLIM,RELE,ABSE,ACC,HMI,
     &                            calc1)
      use iso_fortran_env
      implicit none
!
!23456789012345678901234567890123456789012345678901234567890123456789012
!=======================================================================
!
!  MODIFIED FORM OF QA02AD - OCT 1976 - DJW
!  ========================================
!  Updated 2012
!
!=======================================================================
!
      real(real64) :: ANSWER(4),ALIM,BLIM,RELE,ABSE,ACC,HMI
      interface
        subroutine calc1(z,a)
        use iso_fortran_env
        real(real64),intent(in)  :: z
        real(real64),intent(out) :: a(4)
        end subroutine calc1
      end interface
!
      real(real64) :: F(4,9),X(9),G(4,64),Y(64)
      real(real64) :: RANGE0,RERR,H,H1,HH,HMAX,HCRIT,HMIN,HMINN,
     &                ERROR,ASTEP,DEC,EST,AMODV,AMOD1,SUM0,X1,X2,
     &                Z1,Z2,ESTL,P,QACC
      INTEGER  :: I,J,L,IA,IERRC,KERR,K4D,NUM,LTOP,MTOP
      LOGICAL  :: TEST
      DATA    NUM/4/
      TEST(Z1,Z2)=MTOP.GT.0.AND.ABS(Z1-Z2).LT.HMINN
!
!  INITIALIZE LOCAL VARIABLES
!
      RANGE0=BLIM-ALIM
      IF(RANGE0.LE.0D0)GOTO 40
      H=RANGE0/16D0
      SUM0=0D0
      DO I=1,17
        Y(I)=BLIM+(1-I)*H
        CALL CALC1(Y(I),G(1,I))
        SUM0=SUM0+ABS(G(1,I))
      ENDDO
      RERR=(SUM0-0.5D0*(ABS(G(1,I))+ABS(G(1,17))))*H*RELE
      ACC=MAX(ABSE,RERR)
      HMAX=H
      HCRIT=RANGE0*1.0D-5
      HMINN=MAX(HMI,RANGE0*1.0D-15)
      HMIN=HMINN
      ASTEP=8.D0
      ERROR=180D0*ACC/RANGE0
      DEC=ERROR/20D0
      QACC=ACC*0.06D0
      EST=0D0
      IERRC=0
      K4D=0
      AMOD1=0.0D0
      H1=1.0D0
      DO I=1,NUM
        ANSWER(I)=0D0
      ENDDO
      DO I=1,5
        X(I)=Y(18-I)
        DO L=1,NUM
          F(L,I)=G(L,18-I)
        ENDDO
      ENDDO
      MTOP=12
      X1=X(5)
      GOTO 8
!
!  START
!
   2  X1=ALIM
      IA=1
      IF(LTOP.EQ.5)THEN
        DO J=1,4
          IF(MTOP.GE.64)GOTO 502
          MTOP=MTOP+1
          Y(MTOP)=X(6-J)
          DO L=1,NUM
            G(L,MTOP)=F(L,6-J)
          ENDDO
        ENDDO
      ENDIF
 502  DO I=2,5
        X1=X1+H
        X(I)=X1
        IF(TEST(X(I),Y(MTOP)))THEN
          DO L=1,NUM
            F(L,I)=G(L,MTOP)
          ENDDO
          MTOP=MTOP-1
        ELSE
          CALL CALC1(X(I),F(1,I))
        ENDIF
      ENDDO
   8  LTOP=5
      AMODV=ABS(F(1,1)+F(1,5)-4D0*(F(1,2)+F(1,4))+6D0*F(1,3))
      KERR=0
      IF(AMODV.LE.ERROR)GOTO 11
      kERR=1
      IF(H.LE.HCRIT.AND.H*AMODV.LE.5D0*QACC)GOTO 11
      IF(H.GE.HCRIT.OR.AMODV.LE.H1*AMOD1)GOTO 72
      K4D=K4D+1
      IF(K4D.GT.1)GOTO 7
  72  AMOD1=AMODV
      H1=1D0
   6  HH=H/2D0
      IF(HH.LT.HMIN)GOTO 7
      H=HH
      IA=2
      H1=0.3D0*H1
      AMODV=0.0625D0*AMODV
      ASTEP=2D0*ASTEP
      IF(AMODV.LE.ERROR.OR.K4D.GT.0.OR.H*AMOD1.LE.5D0*QACC)GOTO 2
      GOTO 6
!  DECREASE ERROR REQUIREMENTS
   7  ERROR=ERROR*16
      DEC=ERROR/20D0
      IERRC=IERRC+1
      IF(AMODV.GT.ERROR)GOTO 7
      K4D=0
      GOTO (11,2,21,20,17),IA
!
  11  IF(LTOP.NE.9) GOTO 512
      DO J=1,4
        IF(MTOP.GE.64) GOTO 512
        MTOP=MTOP+1
        Y(MTOP)=X(10-J)
        DO L=1,NUM
          G(L,MTOP)=F(L,10-J)
        ENDDO
      ENDDO
 512  DO I=6,9
        X1=X1+H
        X(I)=X1
        IF(TEST(X(I),Y(MTOP)))THEN
          DO L=1,NUM
            F(L,I)=G(L,MTOP)
          ENDDO
          MTOP=MTOP+1
        ELSE
          CALL CALC1(X(I),F(1,I))
        ENDIF
      ENDDO
      LTOP=9
      AMOD1=16D0*AMODV
      GOTO 13
!
!  HALVE STEP
!
  14  IA=3
      HH=H/2D0
      IF(HH.LT.HMIN)GOTO 7
 141  AMOD1=AMODV
      H=HH
      ASTEP=2D0*ASTEP
      IF(LTOP.NE.9)GOTO 522
      DO J=1,4
        IF(MTOP.GE.64)GOTO 522
        MTOP=MTOP+1
        Y(MTOP)=X(10-J)
        DO L=1,NUM
          G(L,MTOP)=F(L,10-J)
        ENDDO
      ENDDO
 522  I=5
      J=8
      X2=X1-H
  16  X(J+1)=X(I)
      X(J)=X2
      DO L=1,NUM
        F(L,J+1)=F(L,I)
      ENDDO
      CALL CALC1(X(J),F(1,J))
      J=J-2
      I=I-1
      X2=X2-H-H
      IF(J.GT.0)GOTO 16
  13  LTOP=9
      AMODV=ABS(F(1,3)+F(1,7)-4D0*(F(1,4)+F(1,6))+6D0*F(1,5))
      IF(AMODV.LE.ERROR)GOTO 17
      IF(IA.LE.2)THEN
        IF(H.LE.HCRIT.AND.H*AMODV.LE.2D0*QACC)GOTO 171
        KERR=0
      ENDIF
  18  IF(H.LT.HCRIT.AND.AMODV.GT.0.3D0*AMOD1)GOTO 181
      IA=5
      HH=H/2D0
      IF(HH.LT.HMIN)GOTO 7
      X1=X1-4D0*H
      GOTO 141
 181  H=2D0*H
      HMIN=H
      ASTEP=ASTEP/2D0
      J=3
      DO I=2,5
        DO L=1,NUM
          F(L,I)=F(L,J)
        ENDDO
        X(I)=X(J)
        J=J+2
      ENDDO
      LTOP=5
      IA=4
      P=H/3D0
      GOTO 7
!
!  ENTER MAIN LOOP
!
 171  KERR=0
  17  P=H/3D0
      DO L=1,NUM
        ANSWER(L)=ANSWER(L)+P*(4D0*(F(L,2)+F(L,4))+2D0*F(L,3)
     &                          +F(L,1)+F(L,5))
      ENDDO
      ASTEP=ASTEP-2D0
      AMODV=MAX(ABS(F(1,1)+F(1,5)-4D0*(F(1,2)+F(1,4))
     &                         +6D0*F(1,3)),AMODV)
      ESTL=2D0*H*AMODV
      IF(IERRC.NE.0.OR.KERR.NE.0)ESTL=ESTL*15D0
      EST=EST+ESTL
      DO I=1,5
        X(I)=X(I+4)
        DO L=1,NUM
          F(L,I)=F(L,I+4)
        ENDDO
      ENDDO
      LTOP=5
!
!  MAIN LOOP
!
  20  AMODV=ABS(F(1,1)+F(1,5)-4D0*(F(1,2)+F(1,4))+6D0*F(1,3))
      KERR=0
      IF(AMODV.GT.ERROR)THEN
        IF(H.GT.HCRIT.OR.H*AMODV.GT.QACC)GOTO 14
        KERR=1
        QACC=QACC-5D-3*H*AMODV
      ENDIF
  21  DO L=1,NUM
        ANSWER(L)=ANSWER(L)+P*(F(L,1)+4D0*F(L,2)+F(L,3))
      ENDDO
      ASTEP=ASTEP-1D0
      ESTL=H*AMODV
      IF(IERRC.NE.0.OR.KERR.NE.0)ESTL=ESTL*15D0
      EST=EST+ESTL
      IF(H+H.GT.HMAX.OR.AMODV.GT.DEC)GOTO 26
      IF(IERRC.GT.0.AND.H.GT.4D0*HMIN)GOTO 33
      IF(MOD(ASTEP,2D0).NE.0D0)GOTO 26
      H=H+H
      ASTEP=ASTEP/2D0
      X1=X1+H
      P=2D0*P
      DO L=1,NUM
        F(L,4)=F(L,5)
      ENDDO
      X(4)=X(5)
      X(5)=X1
      IF(TEST(X(5),Y(MTOP)))THEN
        DO L=1,NUM
          F(L,5)=G(L,MTOP)
        ENDDO
        MTOP=MTOP-1
      ELSE
        CALL CALC1(X(5),F(1,5))
      ENDIF
  26  IF(ASTEP.LT.1.1D0)GOTO 30
      DO L=1,NUM
        F(L,1)=F(L,3)
      ENDDO
      X(1)=X(3)
      DO I=4,5
        DO L=1,NUM
          F(L,I-2)=F(L,I)
        ENDDO
        X(I-2)=X(I)
        X1=X1+H
        X(I)=X1
        IF(TEST(X(I),Y(MTOP)))THEN
          DO L=1,NUM
            F(L,I)=G(L,MTOP)
          ENDDO
          MTOP=MTOP-1
        ELSE
          CALL CALC1(X(I),F(1,I))
        ENDIF
      ENDDO
      GOTO 20
  33  ERROR=ERROR/16D0
      DEC=0.05D0*ERROR
      IERRC=IERRC-1
      IF(IERRC.LE.0)HMIN=HMINN
      GOTO 26
!
!  EXIT
!
  30  DO L=1,NUM
        ANSWER(L)=ANSWER(L)+P*(F(L,3)+4D0*F(L,4)+F(L,5))
      ENDDO
      ACC=(EST+ESTL)/90D0
      RETURN
!
!  ERROR EXIT
!
  40  DO L=1,NUM
        ANSWER(L)=0D0
      ENDDO
      ACC=0D0
      RETURN
      END subroutine quad41


      subroutine calc3(z,soln)
      implicit none
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c=======================================================================
c
!  Subroutine calc3
!
!  Subroutines calc3 and set are designed to calculate the
!  Hassleman Coefficient for a point on the inner integral loop.
!
!  input:
!    k1  :  wavenumber for which dn/dt is calculated
!    k3  :  nearby wavenumber
!    z   ;  angle (0-2pi) around loop giving positions
!           of k2 and k3
!  output
!    soln(4)  : soln(1) contains integrand.
!    solna(8) : coefficients
!
c=======================================================================
c
      real*8, intent(in)  :: z
      real*8, intent(out) :: soln(4)
!
      real*8 phi,cc,c5,c1,c6,c2,s,s1,s2,s11,s22,f,sa,sb,
     &       fd,h,cost,sint,delta,e1,e2,e3,e4,e5,e6,e7,e8,e9,
     &       e10,e11,con1,qq,b1,b2,dph
      real*8 w12,w23,wp12,wm23,k12,k14,k23,k24,k34
      real*8 px,py,p,cosp,sinp,q,qp,wm13,w13,k13,ec3,d1,d2,
     &       u,ur,si,c3,c4,ax,ay,a,cosa,sina,
     &       sqg,clast,scalea,anga, phi0, t1, t2
      real*8 snpm
      integer is,it,in,iloop, iret
      logical logi, logp7, logp8
c
!      save s,si,sqg,clast,in,px,py,p,cosp,sinp,q,qp,wm13,w13,k13,
!     &     u,ur,ec3,c1,c2,c3,c4,c5,c6,cc,b1,b2,scalea,anga
      save s,si,sqg,clast,in,wm13,w13,k13,
     &     ec3,c1,c2,c3,c4,c5,c6,cc,b1,b2,scalea
c
      real*8 kx1,ky1,kx2,ky2,kx3,ky3,kx4,ky4,k1,k2,k3,k4,
     &       w1,w2,w3,w4,f1,f2,f3,f4
      real*8 :: solna(21)
      real*8 rel1,rel2,rel3,abs1,abs2,abs3,stp1,stp2,stp3,scalei,
     &       rkmax,th
      real*8 pi,pi2,radian,degree,g,den
      integer iprint,icall
c
      common /c_wav/kx1,ky1,kx2,ky2,kx3,ky3,kx4,ky4,k1,k2,k3,k4,
     &       w1,w2,w3,w4,f1,f2,f3,f4
      common /c_hass/solna
      common /c_int/rel1,rel2,rel3,abs1,abs2,abs3,stp1,stp2,stp3,
     &       scalei,rkmax,th
      common /c_con/ pi,pi2,radian,degree,g,den
      common /c_test/iprint,icall
      common /c_calc3/ px,py,p,cosp,sinp,q,qp,u,ur,ax,ay,a,
     &                 cosa,sina,anga
      data   in/0/
c
!      print *," Enter calc3 ",icall, iprint
      icall = icall+1
      logp7 = iprint.ge.7
      logp8 = iprint.ge.8
      iret = 1
      if(.false.)then
        print *,icall,z
        soln = 0d0
        soln(1) = z*z
        return
      endif
c
      soln = 0d0
      solna = 0d0
c
      phi0 = pi*z
      phi  = phi0 + anga
      do while(phi.gt.pi)
        phi = phi - pi2
      enddo
      do while(phi.le.-pi)
        phi = phi + pi2
      enddo
      dph  = pi
c
  21  format(a,12f12.4)
      if(logp7)then
        print 21
        print 21," Enter calc3.       z = ",z
      endif
      if(logp8)then
        print 21,"  phi0, phi, anga, dph = ",phi0, phi,anga,dph
        print 21,"         k1 = ",kx1,ky1,k1,w1,f1
        print 21,"         k3 = ",kx3,ky3,k3,w3,f3
      endif
c
c  start with guess that s = u
c
      s = ur
      cc = cos(phi)
      c5 = cc*b1
      c1 = c5+c5
      c6 = cc*u
      c2 = c6+c6
      if(logp8)print *," c1, c3, c2, c4 = ",c1,c3,c2,c4
      if(logp8)print *," Loop           s          s1         s2",
     &        "          s11          s22         f         sa",
     &        "          sb          fd          h"
      is = 0
  10  it = 0
  20  s1 = sqrt((s+c1)*s+c3)
      s2 = sqrt((s+c2)*s+c4)
      s11 = sqrt(s1)
      s22 = sqrt(s2)
      f=s11-s22-qp
      sa = s1*s11
      sb = s2*s22
      if(sa.eq.0d0.or.sb.eq.0d0)goto 70
      fd = ((s+c5)/sa-(s+c6)/sb)*0.5d0
      if(fd.eq.0d0)goto 70
      if(logp8)print *," p,q,f,fd = ",p,q,f,fd
      if(fd.gt.0)goto 70
      h=f/fd
      if(logp8)print 21," Loop ",s,s1,s2,s11,s22,f,sa,sb,fd,h
      it = it+1
      s=s-h
      if(s.le.0d0) goto 70
      if(abs(h).gt.s*1d-10.and.it.le.40)goto 20
c
      cost=(s2**2+u**2-s**2)/(2d0*s2*u)
      if(abs(cost).gt.1d0)cost = sign(1d0,cost)
      sint = sqrt(1d0-cost*cost)
      if(phi.lt.0d0)sint = -sint
      delta = abs(1d0/fd)
      if(q.ge.0d0)then
        kx2 = s2*(cost*cosp-sint*sinp)
        ky2 = s2*(sint*cosp+cost*sinp)
        kx4 = kx2+px
        ky4 = ky2+py
        k2 = s2
        w2 = sqrt(s2)
        w4 = w2+q
        k4 = w4*w4
      else
        kx4 =-s2*(cost*cosp-sint*sinp)
        ky4 =-s2*(sint*cosp+cost*sinp)
        kx2 = kx4-px
        ky2 = ky4-py
        k4 = s2
        w4 = sqrt(k4)
        w2 = w4-q
        k2 = w2*w2
      endif
c
c  My form of integral is zero if |k1-k3| > |k1-k4|
c                   equivalent to       p > |k3-k2| .
c
      w23=sqrt((kx2-kx3)**2+(ky2-ky3)**2)
      if(p.gt.w23.and..true.)then
        iret = 21
        goto 200
      endif
c
c  Integral set to zero if k2 or k4 is very large (>= rkmax)
c      if(k2.ge.rkmax.or.k4.ge.rkmax)goto 200
c
      f2=snpm(k2,w2,kx2,ky2)
      f4=snpm(k4,w4,kx4,ky4)
      if(logp8)then
        print 21,"         k1 = ",kx1,ky1,k1,w1,f1
        print 21,"         k3 = ",kx3,ky3,k3,w3,f3
        print 21,"         k4 = ",kx4,ky4,k4,w4,f4
        print 21,"         k2 = ",kx2,ky2,k2,w2,f2
        if(f2.eq.0d0.or.f4.eq.0d0)iret = 11
      endif
c  Integral set to zero if k2 or k4 is very large (>= rkmax)
      if((k2.ge.rkmax.or.k4.ge.rkmax).and..true.)then
        iret = 22
        goto 200
      endif
c
      wp12 = (w1+w2)**2
      wm23 = (w2-w3)**2
      w12=sqrt((kx1+kx2)**2+(ky1+ky2)**2)
      w23=sqrt((kx2-kx3)**2+(ky2-ky3)**2)
      k12=kx1*kx2+ky1*ky2
      k14=kx1*kx4+ky1*ky4
      k23=kx2*kx3+ky2*ky3
      k24=kx2*kx4+ky2*ky4
      k34=kx3*kx4+ky3*ky4
c
      e1 = w12-wp12
      e2 = w23-wm23
      if(e1.ne.0d0)e1=2d0*wp12*(k1*k2-k12)*(k3*k4-k34)/e1
      if(e2.ne.0d0)e2=2d0*wm23*(k2*k3+k23)*(k1*k4+k14)/e2
      e3=ec3*(k2*k4+k24)
      e4=0.5d0*k14*k23
      e5=0.5d0*k13*k24
      e6=0.5d0*k12*k34
      e7=+0.25d0*(k14+k23)*wm23**2
      e8=+0.25d0*(k13+k24)*wm13**2
      e9=-0.25d0*(k12+k34)*wp12**2
      e10=2.5d0*k1*k2*k3*k4
      e11=wp12*wm13*wm23*(k1+k2+k3+k4)
      con1=(e1+e2+e3+e4+e5+e6+e7+e8+e9+e10+e11)**2
     &      /(9d0*w1*w2*w3*w4)
!      soln(5) = con1
c
      qq=2d0*dph*s*clast*delta*con1
! Revised numbering
!#if 1
      soln(2)=-qq*f1*f3*f2
      soln(3)= qq*f2*f4*f3
      soln(4)=-qq*f2*f4*f1
      soln(1)= qq*f1*f3*f4+soln(2)+soln(3)+soln(4)
!#else
!      soln(2)= qq*f1*f3*f4
!      soln(3)=-qq*f1*f2*f3
!      soln(4)=-qq*f1*f2*f4
!      soln(1)= qq*f2*f3*f4+soln(2)+soln(3)+soln(4)
!#endif
      if(logp7.or..false.)print 64,iret,z,phi*degree,kx4,ky4,soln
  64  format(" soln ",i2," = ",4f12.6,1p4e16.6)
  65  format(" soln ",i2," = ",48x,1p4e16.6)
!
!  Save details in common c_hass
!
      solna = 0d0
      solna(1) = qq*f1*f3*f4
      solna(2) = soln(2)
      solna(3) = soln(3)
      solna(4) = soln(4)
      solna(5) = con1
      solna(6) = delta
      solna(7) = 2d0*dph*s*clast
      solna(9) = qq
      solna(10) = 9d0*w1*w2*w3*w4
      solna(11) = e1
      solna(12) = e2
      solna(13) = e3
      solna(14) = e4
      solna(15) = e5
      solna(16) = e6
      solna(17) = e7
      solna(18) = e8
      solna(19) = e9
      solna(20) = e10
      solna(21) = e10
!
!     print *," Exit calc3 ..."
      return
c
  70  print 21,"  calc3 failed to converge"
      print *,"  k1 = ",kx1,ky1,k1,w1
      print *,"  k3 = ",kx3,ky3,k3,w3
      print *,"  p, q = ",p,q
      print 21,"  phi  = ",phi,phi*degree
      print *, "  s, h  = ",s+h,h
      print *, "  f, fd = ",f,fd
      print *, "  sa, sb = ",sa, sb
c      print 21," program stopping ..."
c      stop 70
      soln = 0d0
      iret = 31
      if(logp7)print 65,iret,z,soln
      return
c
 200  soln = 0d0
      if(logp7)print 64,iret,z,phi*degree,kx4,ky4,soln
      return
c
      entry set(logi)
c
c  entry point to initialise variables for inner integration loop
c  set logi to false if integral is zero
c
      if(iprint.ge.5)print *," Enter set "
!
      logp7 = iprint.ge.7
      logp8 = iprint.ge.8
      if(logp7)print *,"   Enter set.  Icall = ",icall
      if(in.eq.0)then
        sqg=sqrt(g)
        clast=9d0*pi/(4d0*den**2*sqg)
        in=1
      endif
      if(logp8)print 22,"        set. sqg, pi    = ",sqg,pi
      if(logp8)print *,"         set. clast, den = ",clast,den
  22  format(a,6f10.5)
c
      px=kx1-kx3
      py=ky1-ky3
      p = sqrt(px**2+py**2)
c
c  if p is zero, inner integral is zero, set flag
c
      logi = p.ne.0d0
      if(.not.logi)return
c
c  normal case
c
      cosp=px/p
      sinp=py/p
      k1=sqrt(kx1**2+ky1**2)
      k3=sqrt(kx3**2+ky3**2)
      w1=sqrt(k1)
      w3=sqrt(k3)
      q=w1-w3
      logi = q.ne.0d0
      if(.not.logi)return
c
      qp=abs(q)
      wm13=(w1-w3)**2
      w13=sqrt((kx1-kx3)**2+(ky1-ky3)**2)
      k13=kx1*kx3+ky1*ky3
      ec3=2d0*wm13*(k1*k3+k13)/(w13-wm13)
      f1=snpm(k1,w1,kx1,ky1)
      f3=snpm(k3,w3,kx3,ky3)
c
#if 0
      if(p-q**2.lt.20d0*qp)then
        d1=0.25d0*(p/qp-qp)**2
      else
        d1 = 100d0
      endif
#else
      if(0.25d0*(p-qp**2)**2.gt.100d0*qp**2)then
        d1 = 100d0
      else
        d1=0.25d0*(p/qp-qp)**2
      endif
#endif
      d2=0.5d0*(p-qp*sqrt(2d0*p-qp**2))
      u  = 0.5d0*(d1-d2)
      ur = 0.5d0*(d1+d2)
c
      if(logp8)print 21,"        set. px,py,p,q    = ",px,py,p,q
      if(logp8)print 21,"        set. p/q, p/q-q   = ",(p/q),(p/q-q),
     &                                        (p/q-q)**2
      if(logp8)print 21,"        set. d1,d2,u = ",d1,d2,u

      b1=p+u
      b2=u+d2
      c3=b1**2
      c4=u**2
      si=scalei/b2
c
      if(q.ge.0d0)then
        ax=kx3-u*cosp
        ay=ky3-u*sinp
        a=sqrt(ax**2+ay**2)
        cosa=(ax*px+ay*py)/(p*a)
        sina=(ay*px-ax*py)/(p*a)
      else
        ax=kx1+u*cosp
        ay=ky1+u*sinp
        a=sqrt(ax**2+ay**2)
        cosa=-(ax*px+ay*py)/(p*a)
        sina=-(px*ay-py*ax)/(p*a)
      endif
      s = ur
      if(logp8)print 21,"        set. ax,ay,a,cosa,sina = ",
     &                       ax,ay,a,cosa,sina
c
c  calculate scaling for inner integral over angle
c
      anga = atan2(sina,cosa)
c      scalea = scalei/a
c      scalea = 2d0 - 4d0*scalea/pi
c      if(scalea.lt.0d0)scalea = 0d0
c      s = a
c
      if(logp8)print 21,"        set. anga,scalea = ",anga,scalea
c
      if(iprint.ge.5)print *," Exit set "
      return
c
      end


      function snpm(k,w,kx,ky)
!
!=======================================================================
!23456789012345678901234567890123456789012345678901234567890123456789012
!
!  Subroutine to return action density of a Pierson-Moskowitz
!  or Jonswap spectrum with swell assuming a cos^2 angluar distribution
!  around the x-axis  The units used by the routine are SI (or m.k.s.)
!
!  Input:
!    kx  - x component of wavenumber (inverse/metres)
!    ky  - y component of wavenumber (ditto)
!    k   - equals sqrt(kx*kx + ky*ky)
!    w   - equals sqrt(k)
!
!    Note: w is not angular velocity.  Angular velocity is w times
!          sqrt(g) where g = 9.81 metres/sec^2.
!
!    The wave field is specified by the constants in common jonswp.
!    These should be set before the routine is first entered.  Any
!    changes made later are ignored.
!
!=======================================================================
!
      use iso_fortran_env
      implicit none
!
      real(real64) :: snpm
!
      real(real64),intent(in) ::k,w,kx,ky
!
      real(real64) ::r2,r6,pi,sqg,c1,c2,gamlog,c3,c4a,c4b,c5,rmax,rmin
      real(real64) :: g,den
      logical lin,lsw,lg
!
      save   pi,sqg,c1,c2,gamlog,c3,c4a,c4b,c5,rmax,rmin,
     &       g, den, lin, lsw, lg
!
      real(real64) ::al,fm,gam,rhoa,rhob,swx,swy,sww,swh
      common /jonswp/al,fm,gam,rhoa,rhob,swx,swy,sww,swh
!
      data lin,g,den /.true.,9.81d0,1.0d3/
!
      if(lin)goto 6
!
!  Pierson-Moskowitz Spectrum Calculation
!
   1  if(kx.le.0d0.or.k.gt.rmax.or.k.lt.rmin)goto 4
      r2 = k*k
      r6 = r2*r2*r2
      if(lg)goto 2
      snpm = exp(c2/r2)*c1*kx*kx/(r6*w)
      if(lsw)goto 5
      return
!
!  Johnswap Spectrum Calculation
!
   2  if(w.gt.c3)goto 3
      snpm=exp(c2/r2 + gamlog*exp((w-c3)**2*c4a))*c1*kx*kx/(r6*w)
      if(lsw)goto 5
      return
!
   3  snpm=exp(c2/r2 + gamlog*exp((w-c3)**2*c4b))*c1*kx*kx/(r6*w)
      if(lsw)goto 5
      return
!
!  Johnswap and pierson-Moskowitz are zero for negative angles
!
   4  snpm = 0d0
      if(lsw.and.k.gt.RMAX)goto 5
      return
!
!  Add swell component
!
   5  snpm = snpm+swh*exp(-((kx-swx)**2+(ky-swy)**2)*c5)
      return
!
!  On first entry calculate constants used to improve efficiecy
!
   6  pi = 4d0*datan(1d0)
      sqg=sqrt(g)
      c1 = sqg*den*al/pi
      c2 = -1.25d0*((2d0*pi*fm)**2/g)**2
      gamlog=log(gam)
      c3=2d0*pi*fm/sqg
      c4a=-0.5d0/(c3*rhoa)**2
      c4b=-0.5d0/(c3*rhob)**2
      c5 = 0d0
      if(sww.gt.0d0)c5=1d0/(sww**2)
      lin=.false.
      lsw = swh.gt.0d0
      lg  = gam.ne.1d0
!      write(6,7)
!   7  format(" Jonswop spectrum with swell")
      rmax = 1000d0*sqrt(-c2)
      rmin = 0.1d0*sqrt(-c2)
      goto 1
      end function snpm

