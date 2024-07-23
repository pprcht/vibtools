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

module irtools_plot
  use iso_fortran_env,only:wp => real64
  use irtools_io_mod
  use irtools_convert
  use irtools_atmasses
  use irtools_maths
  implicit none
  private

  public :: plotvibspec
  interface plotvibspec
    module procedure :: plotvibspec_wrapper
  end interface plotvibspec

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine plotvibspec_wrapper(vibname,xmin,xmax,dx,fwhm)
!******************************************************************************
!* Wrapper for plotting a vibspectrum file
!* Takes all the filenames and does the I/O, Lorentzian broadening
!* and finally writes a plain-text file with frequency/intensity pairs
!*
!******************************************************************************
    implicit none
    !> INPUT
    character(len=*),intent(in) :: vibname
    real(wp),intent(in) :: xmin,xmax,dx,fwhm
    !> LOCAL
    integer :: nat,nat3
    real(wp),allocatable :: freq(:)
    real(wp),allocatable :: intens(:)
    real(wp) :: zpve,temp,Svib,Hvib
    character(len=:),allocatable :: line
    integer :: npoints,plotfile
    real(wp),allocatable :: plt(:)
    real(wp) :: nfac,wref,wcurrent,addval,wtmp
    integer :: i,j

    line = firstline(vibname)
    if (line .ne. '$vibrational spectrum') then
      write (*,*) 'ERROR: '//line
      stop
    end if
    write (*,*) repeat('-',50)
    write (*,'(1x,a,a)') 'Input vibspectrum : ',trim(vibname)
    call read_vib_spectrum(vibname,nat3,freq,ints=intens)
    write (*,'(1x,a,i10)') 'Number of modes   :',nat3
    write (*,*) repeat('-',50)

    npoints = nint(abs(xmin-xmax)/dx)
    write (*,'(1x,a,f10.1)') 'Xmin (cm^-1)          :',xmin
    write (*,'(1x,a,f10.1)') 'Xmax (cm^-1)          :',xmax
    write (*,'(1x,a,f10.1)') 'dx (cm^-1)            :',dx
    write (*,'(1x,a,i10)') 'Number of plot points :',npoints
    allocate (plt(npoints),source=0.0_wp)
    write (*,'(1x,a,f10.1)') 'FWHM (cm^-1)          :',fwhm

    call lorentzian_broadening(nat3,freq,intens,npoints,plt,xmin,xmax,dx,fwhm)

    call squarenorm(npoints,plt,nfac)

    open (newunit=plotfile,file='spectrum.txt')
    wcurrent = xmin
    do i = 1,npoints
      write (plotfile,'(3F25.15)') wcurrent,plt(i),0.0_wp
!>--- Add the calculated frequencies as "sticks" to the thrid coloumn. Same normalization factor
      do j = 1,nat3
        if (freq(j) .ge. wcurrent.and.freq(j) .lt. wcurrent+dx) then
          wtmp = freq(j)
          addval = sumphi(wtmp,nat3,freq,intens,fwhm) * nfac
          write (plotfile,'(3F25.15)') wtmp,addval,intens(j)*nfac
        end if
      end do
      wcurrent = wcurrent+dx
    end do
    close (plotfile)

    write(*,'(a)') 'Data written to spectrum.txt'

  end subroutine plotvibspec_wrapper

!========================================================================================!

  subroutine lorentzian_broadening(nmodes,freq,intens,npoints,plt,xmin,xmax,dx,fwhm)
    implicit none
    integer,intent(in) :: nmodes
    real(wp),intent(in) :: freq(nmodes),intens(nmodes)
    integer,intent(in) :: npoints
    real(wp),intent(inout) :: plt(npoints)
    real(wp),intent(in) :: xmin,xmax,dx,fwhm

    real(wp) :: wcurrent,val
    integer :: i,j

    wcurrent = xmin
    do i = 1,npoints
      val = 0.0_wp
      !do j = 1,nmodes
      !  val = val+philorentz(wcurrent,freq(j),intens(j),fwhm)
      !end do
      val = sumphi(wcurrent,nmodes,freq,intens,fwhm)
      plt(i) = val
      wcurrent = wcurrent+dx
    end do

  end subroutine lorentzian_broadening

  function sumphi(nu,nmodes,freq,intens,fwhm) result(sphi)
    implicit none
    real(wp),intent(in) :: nu
    integer,intent(in)  :: nmodes
    real(wp),intent(in) :: freq(nmodes),intens(nmodes)
    real(wp),intent(in) :: fwhm
    real(wp) :: sphi
    integer :: j
    sphi = 0.0_wp
    do j = 1,nmodes
      sphi = sphi+philorentz(nu,freq(j),intens(j),fwhm)
    end do
  end function sumphi

  function philorentz(nu,nu0,i0,w) result(phi)
!*******************************************************************************
!* Function to calculate the value of φ(ν) = I_0/(1 + ((ν_0 - ν)/0.5ω)²)  at ν
!* I.e., φ(ν) is a Lorentzian (Cauchy) line shape function.
!*
!* on Input:  ν         - input frequency
!*            ν_0       - peak position
!*            I_0       - peak intensity
!*            ω         - full width at half maximum (FWHM)
!*
!* on Output: φ(ν)      - function value
!******************************************************************************
    implicit none
    real(wp) :: phi
    real(wp),intent(in)  :: nu
    real(wp),intent(in)  :: nu0
    real(wp),intent(in)  :: i0
    real(wp),intent(in)  :: w
    real(wp) :: x
    phi = 0.0_wp
    x = subfunc(nu,nu0,w)
    phi = i0*(1.0_wp/(1.0_wp+x**2))
    return
  end function philorentz

  function subfunc(nu,nu0,w) result(sf)
!************************************************************
!* Subsidiary function s(ν) =  (ν_0 - ν)/0.5ω)
!*
!* on Input:  ν         - input frequency
!*            ν_0       - peak position
!*            ω         - full width at half maximum (FWHM)
!*
!* on Output: s(ν)      - function value
!***********************************************************
    implicit none
    real(wp) :: sf
    real(wp),intent(in)  :: nu
    real(wp),intent(in)  :: nu0
    real(wp),intent(in)  :: w
    sf = 0.0_wp
    sf = (nu0-nu)/(w/2.0_wp)
    return
  end function subfunc

!========================================================================================!

  subroutine squarenorm(npoints,plt,nfac)
    implicit none
    integer,intent(in) :: npoints
    real(wp),intent(inout) :: plt(npoints)
    real(wp),intent(out) :: nfac
    real(wp) :: tmpa
    integer :: i,j
    nfac = 1.0_wp
    tmpa = 0.0_wp
    do i = 1,npoints
      tmpa = tmpa+plt(i)
    end do
    tmpa = sqrt(tmpa)
    nfac = 1.0_wp/tmpa
    plt = plt*nfac
  end subroutine squarenorm

!========================================================================================!
!========================================================================================!
end module irtools_plot
