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

module irtools_io_mod
  use iso_fortran_env,only:wp => real64
  use irtools_atsymbols
  implicit none
  private

  public :: exists
  public :: firstline
  public :: readxyz
  public :: readdipgrad
  public :: readnewmass
  public :: readhess

  public :: print_vib_spectrum
  public :: read_vib_spectrum

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES BEGIN HERE
!========================================================================================!
!========================================================================================!

  function exists(filename) result(yesno)
!***************************************************
!* one-line wrapper to check for a files' existance
!***************************************************
    character(len=*),intent(in) :: filename
    logical :: yesno
    inquire (file=trim(filename),exist=yesno)
  end function exists

!========================================================================================!

  function firstline(filename) result(first)
!********************************************
!* returns the TRIMMED first line of a file
!********************************************
    implicit none
    character(len=*),intent(in) :: filename
    character(len=:),allocatable :: first
    integer :: ich,io
    character(len=500) :: atmp
    if(.not.(exists(filename)))then
      first='file does not exists' 
    else
      open(newunit=ich,file=filename)
      read(ich,'(a)',iostat=io) atmp
      close(ich)
      if(io/=0)then
        first='file read error'
      else
        first=trim(adjustl(atmp))
      endif
    endif
  end function firstline

!========================================================================================!

  subroutine readxyz(filename,nat,at,xyz)
!*************************************************
!* A simple reader for xyz files
!* Takes the file name as input and outputs
!* the number of atoms (nat), the atom types (at)
!* and coordinates in Angström
!*************************************************
    implicit none
    character(len=*),intent(in) :: filename
    integer,intent(out) :: nat
    integer,allocatable,intent(out) :: at(:)
    real(wp),allocatable,intent(out) :: xyz(:,:)
    integer :: i,j,k,ich,io
    character(len=80) :: atmp
    character(len=2) :: symb
!>--- reset
    nat = 0
    if (.not.exists(filename)) return
!>--- open and read file
    open (newunit=ich,file=filename)
    read (ich,*,iostat=io) j
    if (io == 0) then
      nat = j
    else
      return
    end if
    allocate (at(nat),source=0)
    allocate (xyz(3,nat),source=0.0_wp)
    read (ich,'(a)') atmp
    do i = 1,nat
      read (ich,*,iostat=io) symb,xyz(1:3,i)
      if (io == 0) then
        at(i) = e2i(symb)
      else
        error stop 'error in readxyz()'
      end if
    end do
    close (ich)
  end subroutine readxyz

!========================================================================================!

  subroutine readdipgrad(filename,nat,dipgrad)
!****************************************************************
!* A minimal reader for dipole gradients from a plain text file
!* Takes the file name  and number of atoms (nat) as input and
!* returns 3*3N entries for dipgrad
!****************************************************************
    implicit none
    character(len=*),intent(in) :: filename
    integer,intent(in) :: nat
    real(wp),allocatable,intent(out) :: dipgrad(:,:)
    integer :: i,j,k,l,ich,io

    allocate (dipgrad(3,nat*3),source=0.0_wp)
    if (.not.exists(filename)) return
!>--- open and read file
    open (newunit=ich,file=filename)
    do i = 1,nat*3
      read (ich,*,iostat=io) dipgrad(1:3,i)
      if (io .ne. 0) then
        error stop 'error in readdipgrad()'
      end if
    end do
    close (ich)
  end subroutine readdipgrad

!========================================================================================!

  subroutine readnewmass(filename,masses)
!****************************************************************
!* A reader that reads a plain text two-column file:
!* the first column specifies an atomic number, the second
!* column specifies the associated atomic mass
!****************************************************************
    implicit none
    character(len=*),intent(in) :: filename
    real(wp),intent(inout) :: masses(118)
    integer :: i,j,k,l,ich,io
    real(wp) :: dum
    if (.not.exists(filename)) return

!>--- open and read file
    write (*,'(a,a)') 'Reading atomic masses from file ',trim(filename)
    open (newunit=ich,file=filename)
    do
      read (ich,*,iostat=io) k,dum
      if (io .ne. 0) exit !> EOF
      masses(k) = dum
    end do
    close (ich)
  end subroutine readnewmass

!========================================================================================!

  subroutine readhess(filename,nat,hess)
!*******************************************************************
!* A minimal reader for the plain (i.e. non-mass-weighted) Hessian
!* Takes the file name  and number of atoms (nat) as input and
!* reads 3N lines with 3N entries each.
!* EXPECTED IN ATOMIC UNITS !! (Hartree and Bohr)
!*******************************************************************
    implicit none
    character(len=*),intent(in) :: filename
    integer,intent(in) :: nat
    real(wp),allocatable,intent(out) :: hess(:,:)
    integer :: i,j,k,l,ich,io,n3
    character(len=80) :: atmp
    if (.not.exists(filename)) return

    n3 = 3*nat

    allocate (hess(n3,n3),source=0.0_wp)
!>--- open and read file
    open (newunit=ich,file=filename)
    read (ich,'(a)') atmp
    if (index(trim(atmp),'$hessian') .ne. 0) then
      close (ich)
      call rdhess(n3,hess,filename)
    else
      rewind (ich)
      !write(*,'(a,a)') 'Reading Hessian from plain text file ',trim(filename)
      do i = 1,n3
        read (ich,*,iostat=io) hess(1:n3,i)
        if (io .ne. 0) then
          error stop 'error in readhess()'
        end if
      end do
      close (ich)
    end if
  end subroutine readhess

  subroutine rdhess(nat3,h,filename)
!******************************************
!* Hessian reader in the Turbomole format
!******************************************
    integer,intent(in)  :: nat3
    real(wp),intent(out) :: h(nat3,nat3)
    character(len=*),intent(in) :: filename
    integer  :: iunit,i,j,mincol,maxcol
    character(len=5)  :: adum
    character(len=80) :: a80

    !write(*,'(a,a)') 'Reading Hessian from Turbomole format file ',trim(filename)
    open (newunit=iunit,file=filename)
50  read (iunit,'(a)') a80
    if (index(a80,'$hessian') .ne. 0) then
      do i = 1,nat3
        maxcol = 0
200     mincol = maxcol+1
        maxcol = min(maxcol+5,nat3)
        read (iunit,*) (h(j,i),j=mincol,maxcol)
        if (maxcol .lt. nat3) goto 200
      end do
      close (iunit)
      goto 300
    end if
    goto 50

300 return
  end subroutine rdhess

!========================================================================================!

  subroutine print_vib_spectrum(nat,at,nat3,xyz,freq,intens,outname)
!*****************************************************
!* creates a file in Turbomoles "vibspectrum" format
!*****************************************************
    integer,intent(in) :: nat,nat3
    integer,intent(in) :: at(nat)
    real(wp),intent(in) ::  xyz(3,nat)
    real(wp),intent(in) ::  freq(nat3),intens(nat3)
    character(len=*),intent(in) :: outname
    integer :: i,ich
    real(wp),parameter :: thr = 0.01_wp
    open (newunit=ich,file=outname)
    write (ich,'("$vibrational spectrum")')
    write (ich,'("#  mode    symmetry    wave number    IR intensity    selection rules")')
    write (ich,'("#                       1/cm              km/mol         IR    RAMAN")')
    do i = 1,nat3
      if (abs(freq(i)) .lt. thr) then
        write (ich,'(i6,9x,    f18.2,f16.5,7x," - ",5x," - ")') &
          i,freq(i),0.0_wp
      else
        write (ich,'(i6,8x,"a",f18.2,f16.5,7x,"YES",5x,"YES")') &
          i,freq(i),intens(i)
      end if
    end do
    write (ich,'("$end")')
    close (ich)
  end subroutine print_vib_spectrum

!========================================================================================!

  subroutine read_vib_spectrum(fname,nmodes,freq,ints)
!********************************************************************
!* read vibspectrum file in TM format
!* determines the number of modes and allocates arrays automatically
!********************************************************************
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(out)   :: nmodes
    real(wp),allocatable,intent(out) :: freq(:)          !> frequencies
    real(wp),allocatable,intent(out),optional :: ints(:) !> intensities (optional)
    integer :: k,ich,io,n
    character(len=256) :: atmp
    real(wp) :: floats(10)
    logical :: ex
    nmodes=0
    if (.not.exists(fname)) return
!>--- first, count modes and allocate
    k = 1 ! mode counter
    open (file=fname,newunit=ich)
    rdfile: do
      read (ich,'(a)',iostat=io) atmp
      if (io < 0) exit
      if (index(atmp,'$vibrational spectrum') .ne. 0) then
        rdblock: do
          read (ich,'(a)',iostat=io) atmp
          if (io < 0) exit rdfile
          if (index(atmp,'$end') .ne. 0) exit rdfile
          if (index(atmp,'#') .ne. 0) cycle rdblock !skip comment lines
          k = k+1
        end do rdblock
      end if
    end do rdfile
    nmodes=k-1
    allocate(freq(nmodes),source=0.0_wp)
    if(present(ints))then
      allocate(ints(nmodes), source=0.0_wp)
    endif

    rewind(ich)
!>--- then, READ frequencies (and intensities, if applicable)
    k=1
    rdfile2: do
      read (ich,'(a)',iostat=io) atmp
      if (io < 0) exit
      if (index(atmp,'$vibrational spectrum') .ne. 0) then
        rdblock2: do
          read (ich,'(a)',iostat=io) atmp
          if (io < 0) exit rdfile2
          if (index(atmp,'$end') .ne. 0) exit rdfile2
          if (index(atmp,'#') .ne. 0) cycle rdblock2 !skip comment lines
          call sreadl(atmp,floats,n)
          freq(k) = floats(2)
          if(present(ints))then
             ints(k) = floats(3)
          endif
          k = k+1
        end do rdblock2
      end if
    end do rdfile2

    close (ich)
    return
  end subroutine read_vib_spectrum

!=================================================================================!

  subroutine sreadl(a,xx,nn)
!*****************************************************************
!* sreadl is a helper function that splits a string a into floats
!* This is an replacement of an older F77 routine.
!*****************************************************************
    implicit none
    character(len=*) :: a
    real(wp) :: xx(*)
    integer,intent(out) :: nn
    integer :: al,bl
    integer :: i,j,k,l,io
    character(len=:),allocatable :: b
    character(len=:),allocatable :: dum
    character(len=1) :: c
    real(wp) :: dval
    dum = ''
    nn = 0
    b = trim(adjustl(a))
    al = len_trim(b)+1
    if (al < 1) return
    do i = 1,al
      bl = len_trim(b)
      if (bl >= 1) then
        c = b(1:1)
      end if
      if (c == ' '.or.i == al.or.bl < 1) then
        dum = trim(dum)
        read (dum,*,iostat=io) dval
        if (io == 0) then !if we've just read a float
          nn = nn+1
          xx(nn) = dval
        end if
        dum = '' !reset
        b = trim(adjustl(b)) !advance to the next
        if (bl < 1) exit
      else
        dum = dum//c !add to dum
        bl = len(b)
        b = b(2:bl) !truncate b
      end if
    end do

    deallocate (dum,b)
    return
  end subroutine sreadl

!========================================================================================!
!========================================================================================!
end module irtools_io_mod
