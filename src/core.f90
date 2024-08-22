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

module irtools_core
  use iso_fortran_env,only:wp => real64
  use irtools_io_mod
  use irtools_convert
  use irtools_atmasses
  use irtools_maths
  implicit none
  private

  public :: computespec, computespec_core
  interface computespec
    module procedure :: computespec_wrapper
    module procedure :: computespec_core
  end interface computespec

  public :: computethermo
  interface computethermo
    module procedure :: compute_thermo_wrapper
  end interface computethermo

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine computespec_wrapper(xyzname,hessname,dipname,massname,fscal,outfile)
!********************************************
!* Wrapper for the IR spectra calculation
!* Takes all the filenames and does the I/O
!* before calling the actual calculator
!********************************************
    implicit none
    !> INPUT
    character(len=*),intent(in) :: xyzname
    character(len=*),intent(in) :: hessname
    character(len=*),intent(in) :: dipname
    character(len=*),intent(in) :: massname
    character(len=*),intent(in) :: outfile
    real(wp),intent(in) :: fscal
    !> LOCAL
    integer :: nat,nat3
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    real(wp),allocatable :: hess(:,:)
    real(wp),allocatable :: dipd(:,:)
    real(wp),allocatable :: amass(:)
    real(wp),allocatable :: freq(:)
    real(wp),allocatable :: intens(:)
    real(wp) :: zpve,temp,Svib,Hvib

!>--- Info printout
    write (*,*) '**************************************************************'
    write (*,*) 'PLEASE MAKE SURE HESSIAN AND DIPOLE FILES ARE IN ATOMIC UNITS!'
    write (*,*) '**************************************************************'
    write (*,*) 
    write (*,*) repeat('-',50)
    write (*,'(1x,a,a)') 'Input structure : ',trim(xyzname)
    write (*,'(1x,a,a)') 'Input Hessian   : ',trim(hessname)
    write (*,'(1x,a,a)') 'Input dμ/dR     : ',trim(dipname)
    if (exists(massname)) write (*,'(1x,a,a)') 'Input dμ/dR     : ',trim(massname)
    write (*,'(1x,a,f10.5)') 'Scaling factor  :',fscal
    write (*,'(1x,a,a)') 'Output file     : ',trim(outfile)
    write (*,*) repeat('-',50)

!>--- read files (the subroutines should do the allocation)
    !> Cartesian coordinates
    call readxyz(xyzname,nat,at,xyz)
    xyz = xyz*aatoau
    nat3 = nat*3

    !> Hessian
    call readhess(hessname,nat,hess)

    !> molecular dipole moment derivatives
    call readdipgrad(dipname,nat,dipd)

    !> Atomic masses
    allocate (amass(118))
    amass(:) = ams(:)
    call readnewmass(massname,amass)

    write (*,*) 'Allocated data:'
    write (*,'(3x,a,i0)') 'nat ',nat
    write (*,'(3x,a,i0,a)') 'at(',size(at,1),')'
    write (*,'(3x,a,i0,a,i0,a)') 'xyz(',size(xyz,1),',',size(xyz,2),')'
    write (*,'(3x,a,i0,a,i0,a)') 'hess(',size(hess,1),',',size(hess,2),')'
    write (*,'(3x,a,i0,a,i0,a)') 'dipd(',size(dipd,1),',',size(dipd,2),')'
    write (*,*) repeat('-',50)

!>--- pass on to the calculator
    call computespec_core(nat,at,xyz,hess,dipd,amass,fscal,freq,intens)

!>--- ZPVE (in Hartree)
    zpve = 0.0_wp
    zpve = 0.5_wp*sum(freq(:)/au_to_rcm)
    write (*,'(1x,a,f25.15)') 'ZPVE / Eh        :',zpve
!>--- vib. entropy and enthalpy
    temp = 298.15_wp
    call vibthermo(freq(7:),temp,Svib,Hvib)
    write (*,'(1x,a,f8.2,a)') 'thermodynamics @',temp,'K :'
    write (*,'(1x,a,f25.15)') 'H_vib / kcal/mol :',Hvib
    write (*,'(1x,a,f25.15)') 'S_vib / cal/molK :',Svib
    write (*,*) repeat('-',50)

!>--- create the output file
    call print_vib_spectrum(nat,at,nat3,xyz,freq,intens,outfile)
    write (*,'(1x,a,a)') trim(outfile),' written.'

  end subroutine computespec_wrapper

!========================================================================================!

  subroutine computespec_core(nat,at,xyz,hess,dipd,amass,fscal,freq,intens)
!*****************************************************
!* Core routine that processes the structure, Hessian
!* and dipole moment derivatives to calculate
!* frequencies and IR intensities according to
!* the double harmonic approximation
!*****************************************************
    implicit none
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    real(wp),intent(inout) :: hess(3*nat,3*nat)
    real(wp),intent(in) :: dipd(3,nat*3)
    real(wp),intent(in) :: amass(118)
    real(wp),intent(in) :: fscal
    !> OUTPUT
    real(wp),allocatable,intent(out) :: freq(:)
    real(wp),allocatable,intent(out) :: intens(:)
    !> LOCAL
    integer :: nat3,i,j,k

    nat3 = 3*nat
    allocate (freq(nat3),intens(nat3),source=0.0_wp)

    !> Project translation and rotation, apply mass weighting
    call prj_mw_hess(nat,at,nat3,xyz,amass,hess)

    !> Diagonalize
    call diagH(nat3,hess,freq)
    !> convert to cm^-1
    do i = 1,nat3
      if (freq(i) < 0.0_wp) then
        freq(i) = -sqrt(abs(freq(i)))*au_to_rcm
      else
        freq(i) = sqrt(abs(freq(i)))*au_to_rcm
      end if
    end do
    !> Apply linear frequency scaling factor
    freq = freq*fscal

    !> calculat intensities from dipole derivatives
    call IR_intensities(nat,at,amass,hess,dipd,intens)

  end subroutine computespec_core

!========================================================================================!
  
  subroutine compute_thermo_wrapper(vibname,fscal)
!******************************************************************************
!* Wrapper for thermo calculations if the input is a vibspectrum
!* Takes all the filenames and does the I/O
!*
!* BE CAREFUL ABOUT THE fscal PARAMETER, VIBSPECTRUM MIGHT ALREADY INCLUDE IT!
!*
!******************************************************************************
    implicit none
    !> INPUT
    character(len=*),intent(in) :: vibname
    real(wp),intent(in) :: fscal
    !> LOCAL
    integer :: nat,nat3
    real(wp),allocatable :: freq(:)
    real(wp),allocatable :: intens(:)
    real(wp) :: zpve,temp,Svib,Hvib
    character(len=:),allocatable :: line
    
    line=firstline(vibname)
    if(line.ne.'$vibrational spectrum')then
      write(*,*) 'ERROR: '//line
      stop 
    endif
    write (*,*) repeat('-',50)
    write (*,'(1x,a,a)') 'Input vibspectrum : ',trim(vibname)
    write (*,'(1x,a,f10.5)') 'Scaling factor    :',fscal
    call read_vib_spectrum(vibname,nat3,freq,ints=intens)
    write (*,'(1x,a,i10)') 'Number of modes   :',nat3  
    write (*,*) repeat('-',50)

    freq = freq*fscal

!>--- ZPVE (in Hartree)
    zpve = 0.0_wp
    zpve = 0.5_wp*sum(freq(:)/au_to_rcm)
    write (*,'(1x,a,f25.15)') 'ZPVE / Eh         :',zpve
!>--- vib. entropy and enthalpy
    temp = 298.15_wp
    call vibthermo(freq(7:),temp,Svib,Hvib)
    write (*,'(1x,a,f8.2,a)') 'thermodynamics @',temp,'K :'
    write (*,'(1x,a,f25.15)') 'H_vib / kcal/mol  :',Hvib
    write (*,'(1x,a,f25.15)') 'S_vib / cal/molK  :',Svib
    write (*,*) repeat('-',50)
    write (*,*) 'done.'

  end subroutine compute_thermo_wrapper

!========================================================================================!
!> HESSIAN HANDLING
!========================================================================================!

  subroutine prj_mw_hess(nat,at,nat3,xyz,amass,hess)
!***************************************************************
!* Processing of the Hessian matrix
!* Projection of the translational and rotational DOF out of
!* the numerical Hessian plus the mass-weighting of the Hessian
!***************************************************************
    implicit none
    integer,intent(in) :: nat,nat3
    integer,intent(in) :: at(nat)
    real(wp),intent(in) ::  xyz(3,nat),amass(118)
    real(wp),intent(inout) :: hess(nat3,nat3)
    !real(wp) ::  hess_ut(nat3*(nat3+1)/2),pmode(nat3,1)
    real(wp),allocatable ::  hess_ut(:),pmode(:,:)
    real(wp),allocatable :: invmass(:)
    integer :: ii,jj,k

    allocate (hess_ut(nat3*(nat3+1)/2),source=0.0_wp)
    allocate (pmode(nat3,1),source=0.0_wp)

    !> Transforms matrix of the upper triangle vector
    call dsqtoh(nat3,hess,hess_ut)

    !> Projection (in module irtools_maths)
    call trproj(nat,nat3,xyz,hess_ut,.false.,0,pmode,1)

    !> Transforms vector of the upper triangle into matrix
    call dhtosq(nat3,hess,hess_ut)

    !> Mass weighting
    call get_invsqmass(nat,at,amass,invmass)
    do ii = 1,nat3
      do jj = 1,nat3
        hess(jj,ii) = hess(jj,ii)*invmass(ii)*invmass(jj)
      end do
    end do

    deallocate (invmass)
    deallocate (pmode,hess_ut)
  end subroutine prj_mw_hess

!========================================================================================!

  subroutine dsqtoh(n,a,b)
!**********************************************************
!* converts upper triangle of a matrix into a vector
!* automatically symmetrizes the Hessian, for good measure
!**********************************************************
    implicit none
    integer,intent(in)  :: n
    real(wp),intent(in) :: a(n,n)
    real(wp),intent(out)  :: b(n*(n+1)/2)
    integer :: i,j,k
    k = 0
    do i = 1,n
      do j = 1,i
        k = k+1
        b(k) = (a(i,j)+a(j,i))*0.5_wp
      end do
    end do
  end subroutine dsqtoh

  subroutine dhtosq(n,a,b)
!*********************************************************
!* converts upper triangle vector into a symmetric matrix
!*********************************************************
    implicit none
    integer,intent(in)  :: n
    real(wp),intent(out) :: a(n,n)
    real(wp),intent(in)  :: b(n*(n+1)/2)
    integer :: i,j,k
    k = 0
    do i = 1,n
      do j = 1,i
        k = k+1
        a(j,i) = b(k)
        a(i,j) = b(k)
      end do
    end do
  end subroutine dhtosq

!========================================================================================!

  subroutine get_invsqmass(nat,at,amass,invmass)
!*************************************************
!* Creates an inverse squared atomic mass tensor
!* for mass weighting Hessian etc.
!*************************************************
    implicit none
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: amass(118)
    !> OUTPUT
    real(wp),allocatable :: invmass(:)
    integer :: i,j,k,ii
    allocate (invmass(nat*3),source=0.0_wp)
    do i = 1,nat
      k = at(i)
      do j = 1,3
        ii = (i-1)*3+j
        invmass(ii) = 1.0_wp/sqrt(amass(k)*amutoau)
      end do
    end do
  end subroutine get_invsqmass

!========================================================================================!

  subroutine diagH(nat3,hess,freq)
!********************************************************
!* Diagonalization of the Hessian to obtain
!* frequencies and normal modes (which overwrite hess !)
!* Everything is returned in a.u.
!********************************************************
    implicit none
    !> INPUT
    integer,intent(in) :: nat3
    !> OUTPUT
    real(wp),intent(inout) :: hess(nat3,nat3)
    real(wp),intent(inout) :: freq(nat3)
    !> LOCAL
    integer,allocatable :: iwork(:)
    real(wp),allocatable :: work(:)
    integer :: lwork,liwork,info,i
    !> LAPACK
    external :: dsyevd

    !> LAPACK workspace parameters
    lwork = 1+6*nat3+2*nat3**2
    liwork = 3+5*nat3
    allocate (work(lwork),iwork(liwork))

    !> Diagonalization
    call dsyevd('V','U',nat3,hess,nat3,freq,work,lwork,iwork,liwork,info)

    deallocate (work,iwork)
  end subroutine diagH

!========================================================================================!

  subroutine IR_intensities(nat,at,amass,hess,dipd,intens)
!************************************************************************
!* The actual part referring to the double harmonic approximation
!* project the dipole derivatives via the modes and get the IR intensity
!* Intensities are returned in IR units km/mol
!************************************************************************
    implicit none
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: amass(118)
    real(wp),intent(in) :: hess(nat*3,nat*3)
    real(wp),intent(in) :: dipd(3,nat*3)
    !> OUTPUT
    real(wp),allocatable,intent(out) :: intens(:)
    !> LOCAL
    real(wp),allocatable :: invmass(:)
    integer :: i,j,k,nat3
    real(wp) :: sum2,trdip(3)
    call get_invsqmass(nat,at,amass,invmass)
    nat3 = nat*3
    allocate (intens(nat3),source=0.0_wp)
    do i = 1,nat3
      do k = 1,3
        sum2 = 0.0_wp
        do j = 1,nat3
          sum2 = sum2+dipd(k,j)*(hess(j,i)*invmass(j))
        end do
        trdip(k) = sum2
      end do
      intens(i) = au_to_kmmol*(trdip(1)**2+trdip(2)**2+trdip(3)**2)
    end do

    deallocate (invmass)
  end subroutine IR_intensities

!========================================================================================!

  subroutine Raman_intensities(nat,at,amass,hess,dalphadr,intens)
!************************************************************************
!* The actual part referring to the double harmonic approximation
!* project the polarizability derivatives via the modes and
!* get the Raman intensity
!************************************************************************
    implicit none
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: amass(118)
    real(wp),intent(in) :: hess(nat*3,nat*3)
    real(wp),intent(in) :: dalphadr(6,nat*3)
    !> OUTPUT
    real(wp),allocatable,intent(out) :: intens(:)
    !> LOCAL
    real(wp),allocatable :: invmass(:)
    real(wp),allocatable :: dalphadq(:,:)
    integer :: i,j,k,nat3
    real(wp) :: sum2,asq,gamsq
    call get_invsqmass(nat,at,amass,invmass)
    nat3 = nat*3
    allocate (intens(nat3),source=0.0_wp)
    allocate (dalphadq(6,nat3),source=0.0_wp)
    do i = 1,nat3
      do k = 1,6
        sum2 = 0.0_wp
        do j = 1,nat3
          sum2 = sum2+(hess(j,i)*invmass(j))*dalphadr(k,j)
        end do
        dalphadq(k,i) = sum2
      end do
      asq = (dalphadq(1,i)+dalphadq(3,i)+dalphadq(6,i))**2/9.0_wp
      gamsq = ((dalphadq(1,i)-dalphadq(3,i))**2+(dalphadq(3,i)-dalphadq(6,i))**2+ &
         & (dalphadq(6,i)-dalphadq(1,i))**2+6.0_wp*(dalphadq(2,i)**2+ &
         & dalphadq(5,i)**2+dalphadq(4,i)**2))*0.5_wp
      intens(i) = (45.0_wp*asq+7.0_wp*gamsq)
    end do
    !intens = intens * autoaa4byamu()

    deallocate (dalphadq)
    deallocate (invmass)
  end subroutine Raman_intensities

!========================================================================================!

  subroutine vibthermo(vibs,temp,Svib,Hvib)
!******************************************************
!* Get the harmonic oscillator entropy and enthalpy
!* for the molecular vibrations at a given temperature
!******************************************************
    implicit none
    real(wp),intent(in) :: vibs(:)  !> vib. frequencies in cm^-1
    real(wp),intent(in) :: temp     !> temperature in K
    real(wp),intent(out) :: Svib    !> vibrational entropy in cal/molK
    real(wp),intent(out) :: Hvib    !> vibrational enthalpy in kcal/mol
    !> LOCAL
    integer :: i,j,k
    integer :: nvibs
    real(wp) :: omega,beta,ewj,mu,avmom
    real(wp) :: q_vib,h_vib,cp_ho,cpvib,s_vib,sv_ho
    real(wp),parameter :: R = 1.98726D0    ! GAS CONSTANT IN CALORIES/MOLE

    nvibs = size(vibs,1)

    beta=1.0_wp/kB/Temp ! beta in 1/Eh
    q_vib = 1.0_wp
    h_vib = 0.0_wp
    cp_ho = 0.0_wp
    s_vib = 0.0_wp
    sv_ho  = 0.0_wp
    cpvib = 0.0_wp
 
    do i = 1,nvibs
      omega = vibs(i)/au_to_rcm
      ! omega in Eh, beta in 1/Eh
      ewj = exp(-omega*beta)
      q_vib = q_vib/(1.0_wp-ewj)
      ! h_vib in Eh
      h_vib = h_vib+omega*ewj/(1.0_wp-ewj)
      ! cp_ho in Eh²
      cp_ho = omega**2*ewj/(1.0_wp-ewj)/(1.0_wp-ewj)
      if (omega .gt. 0.0_wp) then
        ! sv is S/R which is dimensionless
        ! harm. osc. entropy
        sv_ho = (vibs(i)/au_to_rcm)*beta*ewj/(1.0_wp-ewj)-log(1.0_wp-ewj)
      else
        sv_ho = 0.0_wp
      end if
      ! heat capacity, all in Eh²
      cpvib = cpvib+cp_ho
      ! entropy s_vib is converted to cal/mol/K... by multiplying with R
      !write(*,*)sv_ho
      s_vib = s_vib+R*sv_ho
    end do
    ! h_vib in Eh, beta is in 1/Eh, T is in K, R is in cal/mol/K,
    Hvib=h_vib*R*beta*Temp/1000.0_wp
    ! same here
    ! cpvib is in Eh², beta in 1/Eh, R in cal/mol/K
    cpvib=cpvib*R*beta**2
    ! Entropy already in correct units
    Svib = s_vib 
    
  end subroutine vibthermo

!========================================================================================!
!========================================================================================!
end module irtools_core
