! This file is part of vibtools.
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

module vibtools_statmech
  use iso_fortran_env,only:wp => real64,stdout=>output_unit
  use vibtools_io_mod
  use vibtools_convert
  use vibtools_atmasses
  use vibtools_maths
  use vibtools_thermo
  use vibtools_axis
  implicit none
  private

  public :: compute_thermodynamics

  real(wp),parameter :: autorcm = 219474.63067_wp
  real(wp),parameter :: rcmtoau = 1.0_wp/autorcm
  real(wp),parameter :: autocal = 627.50947428_wp*1000.0_wp
  real(wp),parameter :: MHztoinvcm = 1.0_wp/2.99792458d+4

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine compute_thermodynamics(nat,at,xyz,nfreq,freq,T,rotnum, &
  &                                 ithr,sthr, zpve,et,ht,ts,g)
!******************************************************************************
!* Wrapper for thermo calculations.
!* Requires input structure, frequencies, a symmetry factor, and a temperature.
!* Returns a lot of stuff...
!*
!******************************************************************************
    implicit none
    !> INPUT
    integer,intent(in)  :: nat
    integer,intent(in)  :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in)  :: nfreq
    real(wp),intent(in) :: freq(nfreq)
    real(wp),intent(in) :: T
    integer,intent(in)  :: rotnum
    real(wp),intent(in) :: ithr
    real(wp),intent(in) :: sthr

    real(wp),intent(out) ::  zpve
    real(wp),intent(out) :: et,ht,ts,g
    !> LOCAL
    real(wp) :: A,B,C,rabc(3),avmom
    real(wp) :: symnum,molmass
    logical :: atom,linear,pr
    integer :: n3,nvib,nimag
    real(wp),allocatable :: vibs(:)
    integer :: i,j,k,l

    real(wp),parameter :: vibthr = 1.0_wp

    pr = .false.
    atom = .false.
    linear = .false.

    molmass = molweight(nat,at)
    symnum = real(rotnum)

!>--- rotational constants calculation
    call axis(nat,at,xyz,rabc(1:3),avmom)
    rabc = rabc*MHztoinvcm !> MHz to cm⁻¹
    A = rabc(3)
    B = rabc(2)
    C = rabc(1)
    rabc(1) = A
    rabc(3) = C

!>--- vib. preparation
    n3 = nfreq
    allocate (vibs(n3),source=0.0_wp)

!>--- structure check
    if (c .lt. 1.d-10) linear = .true.
    if (a+b+c .lt. 1.d-6) then
      atom = .true.
      nvib = 0
    end if

    nvib = 0
    do i = 1,n3
      if (abs(freq(i)) .gt. vibthr) then
        nvib = nvib+1
        vibs(nvib) = freq(i)
      end if
    end do

!>--- invert imaginary modes
    nimag = 0
    do i = 1,nvib
      if (vibs(i) .lt. 0.0_wp) then
        nimag = nimag+1
        if (vibs(i) .gt. ithr) then
          vibs(i) = -vibs(i)
          !if (pr) write (stdout,'(a,i5," :",f10.2)') 'Inverting frequency',i,vibs(i)
        end if
      end if
    end do
 
    vibs = vibs*rcmtoau   !> convert vibs from cm⁻¹ to Eh

!>--- calculate the (harmonic) ZPVE
    zpve = 0.5_wp*sum(vibs(1:nvib))

!>--- evaluate partition functions and return enthalpy, entropy, heat capacity
    call thermodyn(stdout,A,B,C,avmom,linear,atom,symnum,molmass, &
      &              vibs,nvib,T,sthr,et,ht,g,ts,zpve,pr)


    deallocate(vibs)
  end subroutine compute_thermodynamics

!========================================================================================!
!========================================================================================!
end module vibtools_statmech
