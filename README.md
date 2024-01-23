# SuAVE 2.0.0
## DESCRIPTION: 

 The `Surface Assessment Via grid Evaluation (SuAVE)` software was developed to account 
 for the effect of curvature in the calculations of structural properties of chemical 
 interfaces regardless of the chemical composition, asymmetry, and level of atom coarseness. 
 It employs differential geometry techniques, enabling the representation of chemical 
 surfaces as fully differentiable. 


 `Denys E. S. Santos` wrote the code and conceived the algorithm and mathematical formalism. 
 `Thereza A. Soares` conceived the idea and supervised code development, and `Kaline Coutinho`
 revised the code and provided expertise with the equations and discussion.  
 
 ## CITATION: 

 Surface Assessment via Grid Evaluation (SuAVE) for Every Surface Curvature and Cavity 
 Shape. (2022) Denys E. S. Santos, Kaline Coutinho, Thereza A. Soares. J. Chem. Inf. Model.,
 v. 62, p. 4690–4701. [![DOI](https://zenodo.org/badge/71169051.svg)](https://doi.org/10.1021/acs.jcim.2c00673) 


 SuAVE: A Tool for Analyzing Curvature-Dependent Properties in Chemical Interfaces
 (2020) Denys E. S. Santos, Frederico J. S. Pontes, Roberto D. Lins, Kaline Coutinho, 
 Thereza A. Soares. J. Chem. Inf. Model., v. 60(2), p. 473-484. [![DOI](https://zenodo.org/badge/71169051.svg)](https://doi.org/10.1021/acs.jcim.9b00569)

 ## INSTALLATION:

 SuAVE can be compiled in any operational system, is distributed free of charge and
 it is maintained by a single developer. If you experience code issues, we kindly 
 request you to contact us at `e-mail: suave.biomat@gmail.com`


 Prerequisites for installation:
 
 -> `cmake installed`
   
 -> `gfortran installed`

 -> `library libquadmath installed` 


 In order to install the programs and enjoy your analysis time, please follow the steps 
 below:


 **1- Download the code:**

-> From [GitHub web Page](https://github.com/SuAVE-Software/source) 
   clone the repository through the following commmand:


```bash
   git clone https://github.com/SuAVE-Software/source.git
```


-> Dowload the compiled version from [SuAVE Web Page](https://www.biomatsite.net/suave-software)


 **2- Compile the source code:**

 If you downloaded the compiled version you just need to insert the files on a 
 convenient path in order to be run. If you downloaded the source code, please compile
 it by the use of any `FORTRAN` compiler and the `MakeFile`. 

 
 ->> Enter the directory where you have downloaded the code


```bash
   cd ~/PATH_TO_SRC/
```


 ->> Edit makefile in order to proceed with the installation process. Edit `INSTALL_PATH`
     content to update the `PATH` where you want to place the compiled source, and also
     the `FCFLAGS` to adapt it to your needs.


```bash
   FCFLAGS = -O2 (DEFAULT)
   INSTALL_PATH = /usr/local/suave (DEFAULT)
```


   (The use of flags `-O2` or `-O3` is well accepted by this compiler and the code. It will be 
   helpful for extracting the best performance of SuAVE)
   (SuAVE version 2.0.0 is also able to use `OpenMP parallelization`. To do so, the user has to 
   insert the flag `-fopenmp on FCFLAGS directive`)

 
 ->> Once with the makefile updated, run make!


```bash
   make
```

```bash
   sudo make install
```


 ->> Insert the following directives in the .bash_profile file:


```bash
   export SUAVE=/usr/local/suave
   export PATH=$SUAVE:$PATH
```


 ->> Update the bash


```bash
   source .bash_profile
```


   or


```bash
   source .bashrc
```


 ->> Have fun with your analysis !


 ## USAGE:

 A classical use for most of SuAVE tools is exemplified below for s_densph:

 Usage: s_densph -in file.pdb -ind1 file1.ndx -ind2 file2.ndx -dens dens.ndx
 
 file.pdb ---- atomic coordinates in PDB format

 file.ndx ---- index file containing user-selected atoms used to fit the grid points to 
 the chemical surface.

 dens.ndx ---- index file containing user-selected atoms used to calculate density profile.
 

 Options:
 
 `-bin`             defines the number of rectangular partition bins along the phi- and 
                  psi-angles
 
 `-grid`            generates a PDB file containing the grid points used in the fitting 
                  for the last frame in the trajectory file.
 
 `-rmsd`            calculates the RMSD between the fitted grid and the selected atoms in the 
                  index files. This estimates how precisely is the grid surface fitted to the
                  chemical surface throughout the trajectory file.
 
 `-coarse`          generates a coarse grid over the surface index atoms from which a finer grid 
 		  will be generated. This is recommended for surfaces defined by atoms which 
		  greatly fluctuate throughout the trajectory. 
 
 `-begin`           first frame to use in the calculations
 
 `-end`             last frame to use in the calculations
 
 `-skip`            number of trajectory frames to be skipped during the analysis 
 
 `-slices`          defines the number of slices along the axis normal to the system used 
 		  to calculate the density profile
 
 `-rough`           percentage of the original surface roughness the user wants to keep 
 		  (1.0 by default, meaning 100%)
 
 `-help`            prints HELP information and quits


 For more information as to the usage of each code, please use -help flag, or see the Tutorial on [SuAVE
 Web page](https://www.biomatsite.net/suave-software) 
