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

module irtools_c_interface
  use iso_c_binding
  use iso_fortran_env,only:wp => real64
  use irtools_io_mod
  use irtools_convert
  use irtools_atmasses
  use irtools_maths
  use irtools_core
  implicit none
  private

  public :: c_computespec_core

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine c_computespec_core(c_nat,c_at,c_xyz,c_hess,c_dipd, &
    &                           c_ams, c_fscal, c_freq, c_ints) &
    &                        bind(C,name="c_computespec_core")
    implicit none
    !> Input arguments from C
    integer(c_int),value,intent(in) :: c_nat
    integer(c_int),target,intent(in) :: c_at(*)
    real(c_double),target,intent(in) :: c_xyz(3,*) !> NOTE Fortran/C matrix orders
    real(c_double),target,intent(inout) :: c_hess(3*c_nat,3*c_nat)
    real(c_double),target,intent(in) :: c_dipd(3,*)  !> dimension(3,nat*3)
    real(c_double),target,intent(in) :: c_ams(118)  !> atomic masses
    real(c_double),value,intent(in) :: c_fscal
    !> Output arguments to C
    real(c_double),target,intent(out) :: c_freq(3*c_nat)
    real(c_double),target,intent(out) :: c_ints(3*c_nat) 

    !> Local variables
    integer :: nat
    integer,pointer :: at(:)
    real(wp),pointer :: xyz(:,:)
    real(wp),pointer :: hess(:,:)
    real(wp),pointer :: dipd(:,:)
    real(wp),pointer :: ams(:)
    real(wp) :: fscal
    real(wp),pointer :: freq(:)
    real(wp),pointer :: ints(:)

    nat = c_nat
    fscal = c_fscal

    call c_f_pointer(c_loc(c_at),at, [c_nat])
    call c_f_pointer(c_loc(c_xyz),xyz, [3,c_nat]) !> Assumes xyz[nat][3] in C     
    call c_f_pointer(c_loc(c_hess),hess, [3*c_nat,3*c_nat])
    call c_f_pointer(c_loc(c_dipd),dipd, [3,3*c_nat])
    call c_f_pointer(c_loc(c_ams),ams, [118])
    call c_f_pointer(c_loc(c_freq),freq, [3*c_nat])  
    call c_f_pointer(c_loc(c_ints),ints, [3*c_nat])

    call computespec_core(nat,at,xyz,hess,dipd,ams,fscal,freq,ints)

    c_hess(1:3*nat,1:3*nat) = hess(1:3*nat,1:3*nat)
    c_freq(1:3*nat) = freq(1:3*nat)
    c_ints(1:3*nat) = ints(1:3*nat)

  end subroutine c_computespec_core

!========================================================================================!
!========================================================================================!
end module irtools_c_interface
