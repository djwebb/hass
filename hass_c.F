      program hass_c
!
!23456789012345678901234567890123456789012345678901234567890123456789012
!=======================================================================
!
!  File: hass_c.F
!
!  Program to calculate the rate of change of action density due to
!  wave-wave interactions (Hassleman Integral) at a value of kx and ky
!  given on the command line.
!
! ***  Version 2.1   (December 2020)
! ***  For unix systems running Fortran 90.
! ***
! ***  Copyright 1978, 2012, 2020 D.J. Webb
! ***  Released under licence GPL-3.0-or-later
! ***
!
!  The command has the form
!
!    hass_c  0.3  0.0      : where kx= 0.3 and ky=0.0
!
!  It is also possible to use five parameters i.e.
!
!    hass_c 0.3 0.0 0.01 0.3 0.0
!
!  where the additional parameters are the Phillips Constant (giving
!  the amplitude of the wave spectrum), the frequency (Hz) at the
!  peak of the spectrum and the angle (degrees anticlockwise from
!  east) of the spectral peak.
!
!  The results are printed one per line in the form:
!
!     kx   ky  dn/dt
!
!   Check values
!
!   pcon = 0.01, fmax = 0.3, angl = 0.0, kx = 0.4, ky = 0.0
!
!   Result : 0.400     0.000       1.7000E-02
!
!  Compile program in a directory containing both "hass_c.F" and
!  "hass_subs.F".  With the gnu fortran compiler the command is:
!
!     gfortran -O3 -o hass_c hass_c.F
!
!=======================================================================
!
!
      use iso_fortran_env
      implicit none
!
      integer       :: nn
      real(real64)  :: kx, ky, dndt, pcon, fmax, angl, pi, radian
      character(len=80) :: string80
!
      interface
        function hass_f(kx,ky,pcon,fmax,angl)
        use iso_fortran_env
        real(real64) :: hass_f
        real(real64), intent(in) :: kx, ky, pcon, fmax, angl
        end function hass_f
      end interface
!
!  Initialise
!
      pcon = 0.01d0
      fmax = 0.3d0
      angl = 0.0d0
      pi = 4d0*atan(1d0)
      radian = pi/180d0
!
!  Read command arguments
!
      nn = command_argument_count()
      if(nn.ne.2.and.nn.ne.5)then
        print *," Use: hass_c kx ky"
        print *," or:  hass_c kx ky pcon fmax angl(degrees)"
        stop
      endif
      call get_command_argument(1,string80)
      read(string80,*)kx
      call get_command_argument(2,string80)
      read(string80,*)ky
      if(nn.eq.5)then
        call get_command_argument(3,string80)
        read(string80,*)pcon
        call get_command_argument(4,string80)
        read(string80,*)fmax
        call get_command_argument(5,string80)
        read(string80,*)angl
      endif
!
      dndt = hass_f(kx,ky,pcon,fmax,angl*radian)
      print 21,kx,ky,dndt
  21  format(2f10.3,5x,es12.4)
!
      stop
      end

#include "hass_subs.F"
