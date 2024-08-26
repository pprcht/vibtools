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

program vibtools
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use vibtools_io_mod
  use vibtools_core
  use vibtools_plot
  implicit none
  integer :: i,j,k,l
  integer :: nargs
  intrinsic :: iargc,getarg
  character(len=128) :: tmp
  character(len=:),allocatable :: tmp2
  real(wp) :: rdum
  integer :: io
!> PARAMETER
  character(len=:),allocatable :: xyzname
  character(len=:),allocatable :: outfile
  character(len=:),allocatable :: hessname
  character(len=:),allocatable :: dipname
  character(len=:),allocatable :: massname
  logical :: printvers
  real(wp) :: fscal
  integer :: RUNTYPE
  real(wp) :: xmin,xmax,dx,fwhm 
!> DEFAULTS
  RUNTYPE = 0
  printvers = .false.
  outfile = 'vibspectrum'
  massname = 'none given'
  fscal = 1.0_wp
  xmin=100.0_wp
  xmax=5000.0_wp
  dx=1.0_wp
  fwhm=30.0_wp

!> ARGPARSER BEGIN +++++++++++++++++++++++
!>--- parse arguments
  nargs = iargc()
  do i = 1,nargs
    call getarg(i,tmp)
    tmp2 = trim(tmp)

    !>-- first arg *can* be input file
    if (i == 1) then
      if (exists(tmp2)) xyzname = tmp2
    end if

    select case (tmp2)
    case ('-h','--help')
      call help()
    case ('-v','--v','--version')
      printvers = .true.
    case ('-i','--file')
      if (i+1 > nargs) cycle
      call getarg(i+1,tmp)
      xyzname = trim(tmp)
    case ('-o','--output')
      if (i+1 > nargs) cycle
      call getarg(i+1,tmp)
      outfile = trim(tmp)
    case ('-hess','--hessian')
      if (i+1 > nargs) cycle
      call getarg(i+1,tmp)
      hessname = trim(tmp)
      RUNTYPE=1
    case ('-dip','--dipole')
      if (i+1 > nargs) cycle
      call getarg(i+1,tmp)
      dipname = trim(tmp)
      RUNTYPE=1
    case ('-m','--mass')
      if (i+1 > nargs) cycle
      call getarg(i+1,tmp)
      massname = trim(tmp)
    case ('-s','--scal')
      if (i+1 > nargs) cycle
      call getarg(i+1,tmp)
      read (tmp,*,iostat=io) rdum
      if (io == 0) fscal = rdum
    case ('--plot')
      RUNTYPE = 3
    end select
  end do
!> ARGPARSER END +++++++++++++++++++++++
  call printout(printvers)

!>--- Check for all required arguments
  if (.not.allocated(xyzname)) then
    error stop 'No input specified but this is required! See --help'
  end if
  if (RUNTYPE == 0) then
    if (firstline(xyzname) .eq. '$vibrational spectrum') then
      RUNTYPE = 2
    end if
  end if

  select case (RUNTYPE)
  case (1)
!>--- vibspectrum generator from Hessian and dipole derivatives
    if (.not.allocated(hessname)) then
      error stop 'No Hessian file specified but this is required! See --help'
    end if
    if (.not.allocated(dipname)) then
      error stop 'No dipole derivative file specified but this is required! See --help'
    end if
!>--- pass on to actual calculator
    call computespec(xyzname,hessname,dipname,massname,fscal,outfile)

  case (2)
!>--- calculate ZPVE for a given vibspectrum
    call computethermo(xyzname,fscal)

  case(3)
!>--- Create a plot-able file from a vibspectrum one. Applies Lorentzian broadening
    call plotvibspec(xyzname,xmin,xmax,dx,fwhm)

  case default
    write (*,*) 'No (valid) runtype selected! See --help!'
  end select

end program vibtools

!========================================================================================!

subroutine printout(printvers)
  implicit none
  logical :: printvers
  include 'vibtools_metadata.fh'
  if (printvers) then
    write (*,'(1x,a,1x,"v",a)') 'vibtools',trim(version)
    stop
  else
    write (*,'(1x,a,1x,"v",a,1x,a )') 'vibtools',trim(version),trim(date)
    write (*,'(1x,"commit (",a,") compiled by ",a)') commit,author
    call pr_disclaim()
    write (*,*)
  end if
end subroutine printout

subroutine pr_disclaim
  write (*,*)
  write (*,*) 'This program is distributed in the hope that it will be useful,'
  write (*,*) 'but WITHOUT ANY WARRANTY; without even the implied warranty of'
  write (*,*) 'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.'
end subroutine pr_disclaim

!========================================================================================!

subroutine help()
  write (*,*) "-------"
  write (*,*) "vibtools"
  write (*,*) "-------"
  call printout(.false.)

  write (*,*) "A tool for creating an IR spectrum via the double harmonic approximation"
  write (*,*) "Two runtypes are available:" 
  write (*,*) "In the first, the tool requires a molecular structure in xyz format as <input>," 
  write (*,*) "a plain text file with the plain Hessian (in atomic units, 3N*3N entries)"
  write (*,*) "and a file containing molecular dipole direvatives (3N entries)."
  write (*,*) "On output the IR spectrum is created in Turbomoles vibspectrum format."
  write (*,*) "In the second mode, <input> can be a vibspectrum file itself and the"
  write (*,*) "tool will only calculate the ZPVE."
  write (*,*)
  write (*,*) "Usage:"
  write (*,*) "  vibtools <input> [arguments]"
  write (*,*) ""
  write (*,*) "Available options are:"
  write (*,*) ""
  write (*,*) "  -h,--help          Display this message."
  write (*,*) "  -i,--file <str>    Specify input file (alternatively, the input file "
  write (*,*) "                     can be provided as the first argument)"
  write (*,*) "  -hess <str>        Specify Hessian input file (REQUIRED ARGUMENT)"
  write (*,*) "  -dip <str>         Specify dipole derivative input file (REQUIRED ARGUMENT)"
  write (*,*) "  -m,--mass <str>    Specify input file containing user-defined atomic masses,"
  write (*,*) "                     which can be used for atomic mass scaling as discussed"
  write (*,*) "                     in https://dx.doi.org/10.1021/acs.jctc.0c00877"
  write (*,*) "                     Needs two columns: the atomic number and associated mass."
  write (*,*) "  -s,--scal <float>  Specify a linear frequency scaling factor"
  write (*,*) "  -o,--output <str>  Specify the output file (if not provided, the output"
  write (*,*) "                     is simply written to 'vibspectrum')"
  write (*,*) "  --plot             Apply Lorentzian line shapes to a vibspectrum file that"
  write (*,*) "                     was provided as <input> and write plain-text spectrum.txt"
  write (*,*) "                     which can be plotted e.g. with gnuplot."
  write (*,*)

  write (*,*)
  write (*,*) '--help exit.'
  stop
end subroutine help
