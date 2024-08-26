## vibtools
![CI workflow](https://github.com/pprcht/vibtools/actions/workflows/build.yml/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-coral.svg)](./LICENSE)
<h3 align="center">Molecular thermodynamics evaluation and IR spectra generation<br>via the doubly-harmonic approximation</h3>


This repository contains tools to evaluate molecular Hessian matrices and generate IR spectra from the doubly-harmonic approximation.

### Installation

While some portion the program backend is written in Fortran, the Python and C++ bindings enable to install the package via `pip`, which will automatically compile the code. Use the `-v` option to keep track of the compilation process:
```bash
git clone git@github.com:pprcht/vibtools.git
cd vibtools
pip install . -v
```

---

### Running `pyvibtools` in a Python script

See the [`example/`](example/) subdirectory for an use-case example.
The directory contrains the following required files:
- [`pyvibtools.ipynb`](example/pyvibtools.ipynb), a notebook demonstrating the use of `pyvibtools`, also using the following files
- [`struc.xyz`](example/struc.xyz), the input geometry at which dipole derivatives and the Hessian were calculated
- [`dipgrad`](example/dipgrad), a plain-text file with the 9N<sub>at</sub> entries (3x3N<sub>at</sub>) that are the Cartesian dipole derivatives in atomic units
- [`numhess`](example/numhess), the seminumerical Hessian (non-massweighted, in Hartree and Bohr) here in the Turbomole output format. A plain-text format with 3N<sub>at</sub> lines á 3N<sub>at</sub> entries is also valid input.

Additionally the [`vibspectrum.ref`](example/vibspectrum.ref) file provides a higher-level reference spectrum.

---

### Running `pyvibtools` as a command line tool

The `pip install` will generate a cli application of the same name as the package.
The help menu (`pyvibtools -h`) will output:

<details>
<summary>cli output (dropdown tab, click here)</summary>

```
usage: pyvibtools [-h] -i FILE [-hess FILE] [-dip FILE] [-s <float>] [-o FILE] [--plot] [-msc FILE2]

A program to process vibrational spectra.

options:
  -h, --help            show this help message and exit
  -i FILE, --file FILE  Specify input file. Depending on the application, this file can be variety
                        of file types, e.g. a molecular structure or a vibspectrum file.

General Options:
  General options for input and output operations

  -hess FILE, --hessian FILE
                        Specify non-massweighted Hessian input file. The Hessian must be in atomic
                        units (Hartree, Bohr)
  -dip FILE, --dipole-gradient FILE
                        Specify Cartesian dipole derivative input file. Dipole gradient must be in
                        atomic units (charge, Bohr)
  -s <float>, --scal <float>
                        Specify a linear frequency scaling factor
  -o FILE, --output FILE
                        Specify the output file (written in TM vibspectrum format)

Processing/Plotting Options:
  These options will be applied to the spectrum and are useful for additional processing, like
  plotting or comparing spectra via match scores

  --plot                Apply Lorentzian line shapes to a the computed/read spectrumand plot via
                        matplotlib.
  -msc FILE2, --matchscore FILE2
                        Specify vibrational spectrum file for matchscore calculation comparing the
                        computed/read-in (-i) spectrum with FILE2
```
</details>


Typical use cases are (taking files from the [`example/`](example/) subdirectory for demo purposes):
1. The generation of `vibspectrum` files in the Turbomole format
   ```
   pyvibtools -i struc.xyz -hess numhess -dip dipgrad -o vibspectrum.new
   ```
   which will generate a new file called `vibspectrum.new`

2. Plotting a vibrational spectrum with `matplotlib`:
   ```
   pyvibtools -i struc.xyz -hess numhess -dip dipgrad --plot
   ```
   or (with `vibspectrum.new` generated in the previous step)
   ```
   pyvibtools -i vibspectrum.new --plot
   ```

3. Computing matchscores between two spectra:
   ```
   pyvibtools -i vibspectrum.new -msc vibspectrum.ref
   ```


---


### Background information

The program uses the unmodified Hessian and structure to first calculate the mass-weighted Hessian and project out translation and rotation of the molecule.
By diagonalization of the mass-weighted Hessian $\mathbf{F}^{(m)}\mathbf{Q}=\mathbf{\epsilon{Q}}$, vibrational frequencies $\tilde{\nu_p}=\frac{1}{2\pi}\sqrt{\epsilon_p}$ and the corresponding normal modes $Q_p^{(m)}$ are determined. 
Finally, the Cartesian dipole derivatives $\frac{\partial \mu}{\partial R}$ are projected along the modes to determine IR intensities $A_{\tilde{\nu_p}}$ according to the doubly-harmonic approximation:
```math
A_{\tilde{\nu_p}} \propto \left(\frac{\partial{\mu}}{\partial{Q^{(m)}_p}}\right)^2= \sum_{\alpha} \left( \sum^{3N_\mathrm{at}}_{j} \frac{\partial{\mu_{\alpha}}}{\partial{R_j}} \frac{\partial{{R_j}}}{\partial{Q}^{(m)}_{p,j}}  \right)^2
```

The for visualization, the spectrum is processed further by applying Lorentzian line shape functions 
```math
\phi_p(\nu) = I_{_p} \left( 1 + \frac{{\nu_p} - \nu }{0.5w}  \right)^{-1}
```
to each frequency/intensity pair (here $\nu_p = \tilde{\nu_p}$ and $I_p = A_{\tilde{\nu_p}}$ with a FWHM of $w$ = 30 cm<sup>-1</sup>) which is then plotted via `matplotlib`:

![IR spectrum for the example](assets/spectrum.png)



---

<details>
<summary><h4>Building and using the standalone Fortran program</h4> (dropdown tab)</summary>

The following setup is optional and ***not required*** if using the `pip` install.
In fact, the Python command line tool `pyvibtools` has some more functionalities than the Fortran version, which is kept around only for legacy purposes.
To build the `vibtools` Fortran binary with CMake use the following chain of commands (in this example with `gfortran/gcc` compilers)
```bash
FC=gfortran CC=gcc cmake -B _build -Dbuild_exe=true
```
and then followed by
```bash
make -C _build
```
The program can then be found within the `_build` directory and only needs to be added to the program path.

With the compiled program run

```bash
vibtools struc.xyz -dip dipgrad -hess numhess -s 0.9606 
```

which should produce the following output:

```
 vibtools v0.1 Wed, 15 May 14:43:21, 05/15/2024
 commit (5896188) compiled by 'philipp@xps15'
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 
 **************************************************************
 PLEASE MAKE SURE HESSIAN AND DIPOLE FILES ARE IN ATOMIC UNITS!
 **************************************************************
 
 --------------------------------------------------
 Input structure : struc.xyz
 Input Hessian   : numhess
 Input dμ/dR     : dipgrad
 Scaling factor  :   0.96060
 Output file     : vibspectrum
 --------------------------------------------------
 Allocated data:
   nat 14
   at(14)
   xyz(3,14)
   hess(42,42)
   dipd(3,42)
 --------------------------------------------------
 ZPVE / Eh        :        0.088380841694293
 thermodynamics @  298.15K :
 H_vib / kcal/mol :        3.731781748516586
 S_vib / cal/molK :       22.992769821068045
 --------------------------------------------------
 vibspectrum written.
```

The produced `vibspectrum` file can be processed further by the command

```bash
vibtools vibspectrum --plot
```
to write a plain-text `spectrum.txt`

</details>


## TODO

[x] Interface Fortran/C to Python
[x] Create `ipynb` example
[ ] Implement thermodynamics evaluation
[ ] Implement Hessian/dipole gradient readers for other programs
[ ] Implement other types of doubly-harmonic spectra (Raman?) 
