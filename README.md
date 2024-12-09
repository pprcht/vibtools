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
- [`thermo.ipynb`](example/thermo.ipynb) An example for the use of the `thermo` routine for thermostatistical evaluation of the molecular thermodynamics (see below).


Additionally the [`vibspectrum.ref`](example/vibspectrum.ref) file provides a higher-level reference spectrum.

---

### Running `pyvibtools` as a command line tool

The `pip install` will generate a cli application of the same name as the package.
The help menu (`pyvibtools -h`) will output:

<details>
<summary>cli output (dropdown tab, click here)</summary>

```
usage: pyvibtools [-h] -i FILE [-hess FILE] [-dip FILE] [-xyz FILE] [-o FILE] [--plot]
                  [--multiplot FILE [FILE ...]] [-ocsv] [-msc FILE2] [--thermo] [-s <float>]
                  [-xmin <float>] [-xmax <float>] [--autox] [-dx <int>]

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
  -xyz FILE             Specify the used structure explicitly. For use in combination with -i to
                        read in a vibspectrum.
  -o FILE, --output FILE
                        Specify the output file (written in TM vibspectrum format)

Processing/Plotting Options:
  These options will be applied to the spectrum and are useful for additional processing, like
  plotting or comparing spectra via match scores.

  --plot                Apply Lorentzian line shapes to a the computed/read spectrumand plot via
                        matplotlib.
  --multiplot FILE [FILE ...], -mp FILE [FILE ...]
                        One or more files to plot together with -i
  -ocsv                 Write the spectrum (additionally) to vibspectrum.csv
  -msc FILE2, --matchscore FILE2
                        Specify vibrational spectrum file for matchscore calculation comparing the
                        computed/read-in (-i) spectrum with FILE2
  --thermo              Evaluate molecular partition functions and calculate ZPVE, enthalpy,
                        entropy and free energy.

Spectral Data Processing Parameters:
  These options are used to define parameters for the spectra processing via the CLI. These
  options will be primarily applied to the spectrum defined via the -i argument.

  -s <float>, --scal <float>
                        Specify a linear frequency scaling factor
  -xmin <float>         Specify minimum frequency value considered for spectra processing
  -xmax <float>         Specify maximum frequency value considered for spectra processing
  --autox               Allow automatic determination of frequency range if at least one spectrum
                        is experimental
  -dx <int>             Specify frequency interval considered for spectra processing. The minimum
                        value is 1 cm⁻¹
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


### Thermodynamics evaluation

`pyvibtools` includes functions for thermostatistical evaluation of the molecular modes via the **`thermo`** routine. The implementation includes Grimme's rigid-rotor/harmonic-oscillator interpolation of the vibrational entropy (see [S. Grimme, *Chem. Eur. J.*, **2012**, *18*, 9955-9964](https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/chem.201200497) and [P. Pracht, S. Grimme, *Chem. Sci.*, **2021**, *12*,
6551-6568](https://pubs.rsc.org/en/content/articlelanding/2021/sc/d1sc00621e)).

A use-case example can be found in the dedicated Jupyter notebook [`example/thermo.ipynb`](example/thermo.ipynb).

The `thermo`-evaluation provides the following output:
```
   ...................................................
   :                  THERMO SETUP                   :
   :.................................................:
   :  # frequencies                          36      :
   :  # imaginary                            0       :
   :  temperature                        298.15 K    :
   :  symmetry                               c1      :
   :  rotational number                       1      :
   :  scaling factor                  1.0000000      :
   :  rotor cutoff                   25.0000000 cm⁻¹ :
   :  imag. cutoff                  -50.0000000 cm⁻¹ :
   :.................................................:

THERMO Results
==============
 Heat capacity [Cp(T)]      5.23251990e-05 Eh, 3.28345595e+01 cal/mol/K

 ZPVE                       9.20042114e-02 Eh  
 Enthalpy [H(0)-H(T)+PV]    9.47053317e-03 Eh
 -------------------------------------------
                            1.01474745e-01 Eh

 Entropy [S]                1.49488492e-04 Eh, 9.38054486e+01 cal/mol/K
 -T*S                      -4.45699938e-02 Eh

 Free energy [H-T*S]        5.69047508e-02 Eh, 3.57082717e+01 kcal/mol
```

The free energy contribution calculated here is **additive** to the total (electronic) energy.


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

 * [x] Interface Fortran/C to Python
 * [x] Create `ipynb` example
 * [x] Implement thermodynamics evaluation
 * [ ] Implement Hessian/dipole gradient readers for other programs
 * [ ] Implement other types of doubly-harmonic spectra (Raman?) 
