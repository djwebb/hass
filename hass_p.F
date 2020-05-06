      program hass_p
!
!23456789012345678901234567890123456789012345678901234567890123456789012
!=======================================================================
!
!  File: hass_p.F
!
!  Program to calculate the rate of change of action density due to
!  wave-wave interactions (Hassleman Integral) over a range of
!  wavenumbers for a Pierson Moskowitz spectrum with:
!    amplitude (Phillips constant) = 0.01
!    spectrum peak                 = 0.3 Hz
!
!  If the 2-D wavenumber (kx,ky) is given by (i*dx,j*dy) then
!  the range of values output are specified by the parameters:
!
!  i1, i2 :  first and last value of integer 'i'
!  j1, j2 :  first and alst value of integer 'j'
!  dx, dy ;  size of step in x and y directions
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
!  Compile program in a directory containing both "hass_p.F" and
!  "hass_subs.F". With the gnu fortran compiler the command is:
!
!     gfortran -O3 -o hass_p hass_p.F
!
!=======================================================================
!
!
      use iso_fortran_env
      implicit none
!
      integer, parameter      :: i1 = 0        ! Starting i index
      integer, parameter      :: i2 = 10       ! Ending i index
      real(real64), parameter :: dx = 0.1d0    ! Interval in kx direction
      integer, parameter      :: j1 = 0        ! Starting j index
      integer, parameter      :: j2 = 10       ! Ending j index
      real(real64), parameter :: dy = 0.1d0    ! Interval in ky direction
!
      real(real64), parameter :: pcon = 0.01d0   ! Phillips Constant
      real(real64), parameter :: fmax = 0.3d0    ! Spectral Peak
      real(real64), parameter :: angl = 0.0d0    ! Wind sea angle (radians)
!
      integer       :: i,j
      real(real64)  :: kx, ky, dndt
!
      interface
        function hass_f(kx,ky,pcon,fmax,angl)
        use iso_fortran_env
        real(real64) :: hass_f
        real(real64), intent(in) :: kx, ky, pcon, fmax, angl
        end function hass_f
      end interface
!
      do j=j1,j2
        ky = j*dy
        do i=i1,i2
          kx = i*dx
          dndt = hass_f(kx,ky,pcon,fmax,angl)
          print 21,kx,ky,dndt
        enddo
      enddo
  21  format(2f10.3,5x,es12.4)
!
      stop
      end

#include "hass_subs.F"
