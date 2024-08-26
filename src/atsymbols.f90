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

module vibtools_atsymbols
  use iso_fortran_env,only:wp => real64
  implicit none
  private

!&<
  !> Element symbols
  character(len=2),parameter,public :: PSE(118) = [ &
   & 'H ',                                                                                'He', &
   & 'Li','Be',                                                  'B ','C ','N ','O ','F ','Ne', &
   & 'Na','Mg',                                                  'Al','Si','P ','S ','Cl','Ar', &
   & 'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
   & 'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe', &
   & 'Cs','Ba','La',                                                                            &
   &                'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',      &
   &                'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn', &
   & 'Fr','Ra','Ac',                                                                            &
   &                'Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',      &
   &                'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og' ]
!&>

  public :: i2e          !> function to convert atomic number to element symbol
  public :: asym         !> "
  interface asym         !> "
    module procedure i2e !> "
  end interface asym
  public :: e2i          !> function to convert element symbol into atomic number

!========================================================================================!
!========================================================================================!
contains ! MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  integer function e2i(cin)
!*********************************************************
!* e2i is used to map the element (as a string) to integer
!*********************************************************
    implicit none
    character(len=*),intent(in) :: cin
    character(len=:),allocatable :: c
    integer :: iout
    integer :: i,j,k,ich,io,Z
    logical :: ex
    c = trim(convertlable(cin))
    read (cin,*,iostat=io) j
    if (io == 0) Z = j
    if (any(PSE(:) .eq. c)) then
      do i = 1,118
        if (trim(PSE(i)) .eq. c) then
          iout = i
          exit
        end if
      end do
    else if (io == 0.and.Z <= 118) then
      iout = Z
    else !> special cases
      select case (trim(c))
      case ('D'); iout = 1
      case ('T'); iout = 1
      case default; iout = 0
      end select
    end if
    e2i = iout
  end function e2i

!========================================================================================!

  character(len=2) function i2e(iin,oformat)
!************************************************************
!* i2e is used to map the element (as a integer) to a string
!************************************************************
    implicit none
    integer,intent(in) :: iin
    character(len=:),allocatable :: c
    character(len=*),optional :: oformat
    if (iin <= 118) then
      c = uppercase(PSE(iin))
    else
      c = 'XX'
    end if
    i2e = trim(c)
    if (present(oformat)) then
      select case (oformat)
      case ('lc','lowercase')
        i2e = lowerCase(trim(c))
      case ('nc','nicecase')
        if (len_trim(c) .gt. 1) then
          c(2:2) = lowerCase(c(2:2))
          i2e = trim(c)
        end if
      case default
        continue
      end select
    end if
  end function i2e

!========================================================================================!

  function upperCase(s)
    implicit none
    character(len=*),intent(in) :: s
    character(len=:),allocatable :: sout
    character(len=:),allocatable :: upperCase
    integer :: ic,i
    character(26),Parameter :: high = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26),Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
    sout = s
    do i = 1,LEN_TRIM(s)
      ic = INDEX(low,s(i:i))
      if (ic > 0) sout(i:i) = high(ic:ic)
    end do
    call move_alloc(sout,upperCase)
  end function upperCase

  function lowerCase(s)
    implicit none
    character(len=*),intent(in) :: s
    character(len=:),allocatable :: sout
    character(len=:),allocatable :: lowerCase
    integer :: ic,i
    character(26),Parameter :: high = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26),Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
    sout = s
    do i = 1,LEN_TRIM(s)
      ic = INDEX(high,s(i:i))
      if (ic > 0) sout(i:i) = low(ic:ic)
    end do
    call move_alloc(sout,lowerCase)
  end function lowerCase

  function convertlable(s)
    implicit none
    character(len=*),intent(in) :: s
    character(len=:),allocatable :: sout
    character(len=:),allocatable :: convertlable
    integer :: ic,i
    character(14),parameter :: lab = '0123456789*_+-'
    character(26),parameter :: high = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26),parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
    sout = s
    do i = 1,len_trim(s)
      ic = index(lab,s(i:i))
      if (ic > 0) sout(i:i) = ' '
      ic = index(low,s(i:i))
      if (ic > 0) sout(i:i) = high(ic:ic)
    end do
    sout = trim(adjustl(sout))
    if (len_trim(sout) .gt. 1) then
      sout(2:2) = lowerCase(sout(2:2))
    else
      sout = sout//' '
    end if
    call move_alloc(sout,convertlable)
  end function convertlable

!========================================================================================!
!========================================================================================!
end module vibtools_atsymbols
