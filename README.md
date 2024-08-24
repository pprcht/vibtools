## IRtools
![CI workflow](https://github.com/pprcht/IRtools/actions/workflows/build.yml/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-coral.svg)](./LICENSE)
<h3 align="center">Molecular thermodynamics evaluation and IR spectra generation<br>via the doubly-harmonic approximation</h3>


This repository contains tools to evaluate molecular Hessian matrices and generate IR spectra from the doubly-harmonic approximation.

### Installation

While some portion the program backend is written in Fortran, the Python and C++ bindings enable to install the package via `pip`, which will automatically compile the code. Use the `-v` option to keep track of the compilation process:
```bash
git clone git@github.com:pprcht/IRtools.git
cd IRtools
pip install . -v
```

---

### Running the program

See the [`example/`](example/) subdirectory for an use-case example.
The directory contrains the following required files:
- [`pyirtools.ipynb`](example/pyirtools.ipynb), a notebook demonstrating the use of `pyirtools`, also using the following files
- [`struc.xyz`](example/struc.xyz), the input geometry at which dipole derivatives and the Hessian were calculated
- [`dipgrad`](example/dipgrad), a plain-text file with the 9N<sub>at</sub> entries (3x3N<sub>at</sub>) that are the Cartesian dipole derivatives in atomic units
- [`numhess`](example/numhess), the seminumerical Hessian (non-massweighted, in Hartree and Bohr) here in the Turbomole output format. A plain-text format with 3N<sub>at</sub> lines á 3N<sub>at</sub> entries is also valid input.

Additionally the [`vibspectrum.ref`](example/vibspectrum.ref) file provides a higher-level reference spectrum.


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
To build the `irtools` Fortran binary with CMake use the following chain of commands (in this example with `gfortran/gcc` compilers)
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
irtools struc.xyz -dip dipgrad -hess numhess -s 0.9606 
```

which should produce the following output:

```
 irtools v0.1 Wed, 15 May 14:43:21, 05/15/2024
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
irtools vibspectrum --plot
```
to write a plain-text `spectrum.txt`

</details>
