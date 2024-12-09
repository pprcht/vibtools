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

module vibtools_c_interface
  use iso_c_binding
  use iso_fortran_env,only:wp => real64
  use vibtools_io_mod
  use vibtools_convert
  use vibtools_atmasses
  use vibtools_maths
  use vibtools_core
  use vibtools_plot
  use vibtools_statmech
  implicit none
  private

  public :: c_computespec_core
  public :: c_print_vib_spectrum_stdout
  public :: c_lorentzian_broadening
  public :: c_compute_thermodynamics

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
    real(wp),allocatable :: freqtmp(:),intstmp(:)

    nat = c_nat
    fscal = c_fscal

    call c_f_pointer(c_loc(c_at),at, [c_nat])
    call c_f_pointer(c_loc(c_xyz),xyz, [3,c_nat]) !> Assumes xyz[nat][3] in C     
    call c_f_pointer(c_loc(c_hess),hess, [3*c_nat,3*c_nat])
    call c_f_pointer(c_loc(c_dipd),dipd, [3,3*c_nat])
    call c_f_pointer(c_loc(c_ams),ams, [118])
    call c_f_pointer(c_loc(c_freq),freq, [3*c_nat])  
    call c_f_pointer(c_loc(c_ints),ints, [3*c_nat])

    call computespec_core(nat,at,xyz,hess,dipd,ams,fscal,freqtmp,intstmp)
    
    freq(:) = freqtmp(:)
    ints(:) = intstmp(:)

    c_hess(1:3*nat,1:3*nat) = hess(1:3*nat,1:3*nat)
    c_freq(1:3*nat) = freq(1:3*nat)
    c_ints(1:3*nat) = ints(1:3*nat)

    if(allocated(intstmp)) deallocate(intstmp)
    if(allocated(freqtmp)) deallocate(freqtmp)
  end subroutine c_computespec_core

!========================================================================================!

subroutine c_print_vib_spectrum_stdout(c_nat3, c_freq, c_intens) &
    &            bind(C, name="c_print_vib_spectrum_stdout")
    implicit none
    integer(c_int),value,intent(in) :: c_nat3
    real(c_double),target, intent(in) :: c_freq(*), c_intens(*)
    real(wp),pointer :: freq(:)
    real(wp),pointer :: intens(:)
    integer :: nat3
    nat3= c_nat3
    call c_f_pointer(c_loc(c_freq),freq, [c_nat3])
    call c_f_pointer(c_loc(c_intens),intens, [c_nat3])
    !> Call the original Fortran subroutine
    call print_vib_spectrum_stdout(nat3, freq, intens)
end subroutine c_print_vib_spectrum_stdout

!========================================================================================!

  subroutine c_lorentzian_broadening(c_nmodes,c_freq,c_intens,c_npoints, &
  &                c_plt,c_xmin,c_xmax,c_dx,c_fwhm) bind(C, name='c_lorentzian_broadening')
    implicit none
    !> C variables
    integer(c_int),value,intent(in) :: c_nmodes
    real(c_double),target,intent(in) :: c_freq(*),c_intens(*)
    integer(c_int),value,intent(in) :: c_npoints
    real(c_double),target,intent(inout) :: c_plt(*)
    real(c_double),value,intent(in) :: c_xmin,c_xmax,c_dx,c_fwhm
    !> fortran variables
    integer :: nmodes,npoints
    real(wp) :: xmin,xmax,dx,fwhm
    real(wp),pointer :: freq(:),intens(:),plt(:)
    nmodes = c_nmodes
    xmin = c_xmin
    xmax = c_xmax
    dx = c_dx
    fwhm = c_fwhm
    npoints = c_npoints
    call c_f_pointer(c_loc(c_freq),freq, [c_nmodes])
    call c_f_pointer(c_loc(c_intens),intens, [c_nmodes])
    call c_f_pointer(c_loc(c_plt),plt, [c_npoints])
    call lorentzian_broadening(nmodes,freq,intens,npoints,plt,xmin,xmax,dx,fwhm)
  end subroutine c_lorentzian_broadening

!========================================================================================!


  subroutine c_compute_thermodynamics(c_nat,c_at,c_xyz,c_nfreq,c_freq, &
    &                                 c_T,c_sthr,c_ithr, c_rotnum, &
    &                                 c_zpve, c_et, c_ht, c_ts, c_cp, c_g) &
    &                        bind(C,name="c_compute_thermodynamics")
    implicit none
    !> Input arguments from C
    integer(c_int),value,intent(in)  :: c_nat
    integer(c_int),target,intent(in) :: c_at(*)
    real(c_double),target,intent(in) :: c_xyz(3,*) !> NOTE Fortran/C matrix orders
    integer(c_int),value,intent(in)  :: c_nfreq
    real(c_double),target,intent(in) :: c_freq(c_nfreq)
    real(c_double),value,intent(in)  :: c_T
    real(c_double),value,intent(in)  :: c_sthr
    real(c_double),value,intent(in)  :: c_ithr
    integer(c_int),value,intent(in)  :: c_rotnum
    !> Output arguments to C
    real(c_double),target,intent(out) :: c_zpve
    real(c_double),target,intent(out) :: c_et     
    real(c_double),target,intent(out) :: c_ht
    real(c_double),target,intent(out) :: c_ts
    real(c_double),target,intent(out) :: c_cp
    real(c_double),target,intent(out) :: c_g


    !> Local variables
    integer :: nat,nfreq,rotnum
    real(wp) :: T,sthr,ithr
    integer,pointer :: at(:)
    real(wp),pointer :: xyz(:,:)
    real(wp),pointer :: freq(:)
    real(wp) :: zpve,et,ht,ts,cp,g

    nat = c_nat
    nfreq = c_nfreq
    T = c_T
    sthr = c_sthr
    ithr = c_ithr
    rotnum = c_rotnum

    call c_f_pointer(c_loc(c_at),at, [c_nat])
    call c_f_pointer(c_loc(c_xyz),xyz, [3,c_nat]) !> Assumes xyz[nat][3] in C     
    call c_f_pointer(c_loc(c_freq),freq, [c_nfreq])
  
    call compute_thermodynamics(nat,at,xyz,nfreq,freq,T,rotnum, &
    &                           ithr,sthr, zpve,et,ht,ts,cp,g)

    c_zpve = zpve
    c_et = et
    c_ht = ht
    c_ts = ts
    c_cp = cp
    c_g = g


  end subroutine c_compute_thermodynamics

!========================================================================================!
!========================================================================================!
end module vibtools_c_interface
