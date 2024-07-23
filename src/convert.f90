! This file is part of irtools.
!
! MIT License
!   
! Copyright (c) 2024 Philipp Pracht
! 
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
! 
!==========================================================================!
!> Partially adapted from the CREST and xtb codes under GNU LGPL3
!==========================================================================!

module irtools_convert
  use iso_fortran_env,only:wp => real64
  implicit none
  public

  real(wp),parameter,public :: bohr = 0.52917726_wp
  real(wp),parameter,public :: angstrom = 1.0_wp/bohr
  real(wp),parameter,public :: autoaa = bohr
  real(wp),parameter,public :: aatoau = angstrom

  real(wp),parameter :: au_to_rcm = 219474.63_wp
  real(wp),parameter :: amutokg = 1.660539040e-27_wp
  real(wp),parameter :: metokg  = 9.10938356e-31_wp
  real(wp),parameter :: kgtome = 1.0_wp / metokg
  real(wp),parameter :: amutoau = amutokg*kgtome
  real(wp),parameter :: amutoau2 = amutoau**2

  real(wp),parameter :: kB = 3.166808578545117e-06_wp

  !> lightspeed
  real(wp),public,parameter :: lightspeed = 137.0359990740_wp
  !> femtosectons to atomic time units
  real(wp),public,parameter :: fstoau = 41.3413733365614_wp
  !> Coulomb to atomic charge units (electrons)
  real(wp),public,parameter :: autoc = 1.6021766208e-19_wp
  !> Debye to atomic units
  real(wp),public,parameter :: autod = autoc*lightspeed*autoaa**2*fstoau*1.0e+16_wp
  real(wp),public,parameter :: dtoau = 1.0_wp/autod

!> ----- DIPOLE DERIVATIVE UNITS -----
!  Dipole derivatives along mass-weighted normal mode coordinates (a.u.) to km/mol (IR int.)
!  autokmmol is the factor to convert dipole derivatives along
!  (mass-weighted) normal mode coordinates from a.u. into
!  the Naperian absorption coefficient
!
!        (pi)*(Navogadro)*(dipole derivative along normal mode)**2
!   A =  ---------------------------------------------------------  ,
!                       3*(4*(pi)*(epsilon0))*c**2
!
!  measured in km/mol. Cf. International Union of Pure and Applied
!  Chemistry (IUPAC), Quantities, Units and Symbols in Physical
!  Chemistry, Recommendations 1993, Reprinted with Corrections 1995.
!  Blackwell Scientific Publications.
  real(wp),public,parameter :: au_to_kmmol = 1.7770969e+6_wp

end module irtools_convert
