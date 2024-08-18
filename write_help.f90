subroutine options_cart()
  
  write(*, *) ""
  write(*, *) "-bin             defines the number of rectangular partition "
  write(*, *) "bins along the x- and y-axes"
  write(*, *) ""
  write(*, *) "-rmsd            calculates the RMSD between the fitted" 
  write(*, *) "grid and the selected atoms in the index files. This" 
  write(*, *) "estimates how precisely is the grid surface fitted to the" 
  write(*, *) "chemical surface throughout the trajectory file."
  write(*, *) ""
  write(*, *) "-coarse          generates a coarse grid over the surface" 
  write(*, *) "index atoms from which a finer grid will be generated. This" 
  write(*, *) "is recommended for surfaces defined by atoms which greatly" 
  write(*, *) "fluctuate throughout the trajectory. "
  write(*, *) ""
  write(*, *) "-begin           first frame to use in the calculations"
  write(*, *) ""
  write(*, *) "-end             last frame to use in the calculations"
  write(*, *) ""
  write(*, *) "-skip            number of trajectory frames to be skipped" 
  write(*, *) "during the analysis "
  write(*, *) ""
  write(*, *) "-rough           percentage of the original surface roughness"
  write(*, *) "the user wants to keep (1.0 by default, meaning 100%)"
  write(*, *) ""
  write(*, *) "-help            prints HELP information and quits"
     
end subroutine options_cart

subroutine options_sphe()

  write(*, *) ""
  write(*, *) "-bin             defines the number of rectangular partition "
  write(*, *) "bins along the phi- and psi-angles"
  write(*, *) ""
  write(*, *) "-rmsd            calculates the RMSD between the fitted"
  write(*, *) "grid and the selected atoms in the index files. This"
  write(*, *) "estimates how precisely is the grid surface fitted to the"
  write(*, *) "chemical surface throughout the trajectory file."
  write(*, *) ""
  write(*, *) "-coarse          generates a coarse grid over the surface"
  write(*, *) "index atoms from which a finer grid will be generated. This"
  write(*, *) "is recommended for surfaces defined by atoms which greatly"
  write(*, *) "fluctuate throughout the trajectory. "
  write(*, *) ""
  write(*, *) "-begin           first frame to use in the calculations"
  write(*, *) ""
  write(*, *) "-end             last frame to use in the calculations"
  write(*, *) ""
  write(*, *) "-skip            number of trajectory frames to be skipped"
  write(*, *) "during the analysis "
  write(*, *) ""
  write(*, *) "-rough           percentage of the original surface roughness"
  write(*, *) "the user wants to keep (1.0 by default, meaning 100%)"
  write(*, *) ""
  write(*, *) "-help            prints HELP information and quits"
  
end subroutine options_sphe

subroutine write_help(name)

  character(len=8) :: name

  select case (trim(name))

     case ('area')
     
        write(*, *) ""
        write(*, *) ""
        write(*, *) "s_area calculates the time-dependent area per molecule of the" 
        write(*, *) "interface."
        write(*, *) ""
        write(*, *) "Usage: s_area -in file.pdb -ind1 file1.ndx -ind2 file2.ndx"
        write(*, *) "              -lipid N_lipids"
        write(*, *) ""
        write(*, *) "file.pdb ---- atomic coordinates in PDB format"
        write(*, *) "file.ndx ---- index file containing user-selected atoms used" 
        write(*, *) "to fit the grid points to the chemical surface."
        write(*, *) "N_lipids ---- number of lipids composing the leaflet"
        write(*, *) ""
        write(*, *) "Options:"
        write(*, *) ""
        write(*, *) "-grid            generates a PDB file containing the grid"
        write(*, *) "points used in the fitting for the last frame in the "
        write(*, *) "trajectory file."
        call options_cart()
        stop

     case('dens') 

        write(*, *) ""
        write(*, *) ""
        write(*, *) "s_dens calculates the density profile for the system taking" 
        write(*, *) "into account the curvature of the surface/interface."
        write(*, *) ""
        write(*, *) "Usage: s_dens -in file.pdb -ind1 file1.ndx -ind2 file2.ndx"
        write(*, *) "-dens dens.ndx"
        write(*, *) ""
        write(*, *) "file.pdb ---- atomic coordinates in PDB format"
        write(*, *) "file.ndx ---- index file containing user-selected atoms used" 
        write(*, *) "to fit the grid points to the chemical surface."
        Write(*, *) "dens.ndx ---- index file containing user-selected atoms used "
        Write(*, *) "to calculate density profile."  
        write(*, *) ""
        write(*, *) "Options:"
        write(*, *) ""
        write(*, *) "-grid            generates a PDB file containing the grid"
        write(*, *) "points used in the fitting for the last frame in the "
        write(*, *) "trajectory file."
        write(*, *) ""
        write(*, *) "-inside           used to calculate number of atoms of a" 
        write(*, *) "specified molecule inside the structure"
        write(*, *) ""
        write(*, *) "-map             used to generate a 2D density map for a "
        write(*, *) "specified group of atoms" 
        write(*, *) ""
        write(*, *) "-slices          defines the number of slices along the axis"
        write(*, *) "normal to the system used to calculate the density profile"
        call options_cart()
        stop

     case ('order')

        write(*, *) ""
        write(*, *) ""
        Write(*, *) "s_order calculates the distribution of angles between the "
        Write(*, *) "vectors normal to the rectangular partition surfaces and the z" 
        Write(*, *) "axis of the system. s_order also calculates the curvature" 
        Write(*, *) "order parameter P(theta) or SC."
        write(*, *) ""
        write(*, *) "Usage: s_order -in file.pdb -ind file.ndx"
        write(*, *) ""
        write(*, *) "file.pdb ---- atomic coordinates in PDB format"
        write(*, *) "file.ndx ---- index file containing user-selected atoms used" 
        write(*, *) "to fit the grid points to the chemical surface."
        write(*, *) ""
        write(*, *) "Options:"
        write(*, *) ""
        write(*, *) "-grid            generates a PDB file containing the grid"
        write(*, *) "points used in the fitting for the last frame in the "
        write(*, *) "trajectory file."
        write(*, *) ""
        write(*, *) "-range           defines the range specified in the XPM file"
        call options_cart()
        stop

     case('thick')

        write(*, *) ""
        write(*, *) ""
        write(*, *) "s_thick calculates the average bilayer thickness and volume"
        write(*, *) "per molecule along the simulation time"
        write(*, *) ""
        write(*, *) "Usage: s_thick -in file.pdb -ind1 file1.ndx -ind2 file2.ndx"
        write(*, *) "               -lipid N_lipids"
        write(*, *) ""
        write(*, *) "file.pdb ---- atomic coordinates in PDB format"
        write(*, *) "file.ndx ---- index file containing user-selected atoms used"
        write(*, *) "to fit the grid points to the chemical surface."
        write(*, *) "N_lipids ---- number of lipids composing the leaflet"
        write(*, *) ""
        write(*, *) "Options: "
        write(*, *) ""
        write(*, *) "-grid            generates a PDB file containing the grid"
        write(*, *) "points used in the fitting for the last frame in the "
        write(*, *) "trajectory file."
        write(*, *) ""
        write(*, *) "-range           defines the range specified in the XPM file"
        call options_cart()
        stop
        
     case('topog')

        write(*, *) ""
        write(*, *) ""
        write(*, *) "s_topog calculates the level surface of the chemical interface"
        write(*, *) ""
        write(*, *) "Usage: s_topog -in file.pdb -ind1 file1.ndx -ind2 file2.ndx"
        write(*, *) ""
        write(*, *) "file.pdb ---- atomic coordinates in PDB format"
        write(*, *) "file.ndx ---- index file containing user-selected atoms used"
        write(*, *) "to fit the grid points to the chemical surface."
        write(*, *) ""
	write(*, *) "Options:"
        write(*, *) ""
        write(*, *) "-grid            generates a PDB file containing the grid"
        write(*, *) "points used in the fitting for the last frame in the "
        write(*, *) "trajectory file."
        write(*, *) ""
        write(*, *) "-range           defines the range specified in the XPM file"
        call options_cart()
        stop
        
     case('grid')

        write(*, *) ""
        write(*, *) ""
        write(*, *) "s_grid builds a grid per frame throughout a trajectory file."
        write(*, *) "Its output is useful to verify how accurate is the fitting of"
        write(*, *) "the calculated grid on the chemical surface.  "
        write(*, *) ""
        write(*, *) "Usage: s_grid -in file.pdb -ind file.ndx"
        write(*, *) ""
        write(*, *) "file.pdb ---- atomic coordinates in PDB format"
        write(*, *) "file.ndx ---- index file containing user-selected atoms used"
        write(*, *) "to fit the grid points to the chemical surface."
        write(*, *) ""
	write(*, *) "Options:"
        write(*, *) ""
        call options_cart()
        stop
        
     case('index')

        write(*, *) ""
        write(*, *) ""
        write(*, *) "s_index creates the input file containing indexes assigned to"
        write(*, *) "the atoms, which will be used to define the surface/interface."
        write(*, *) "The user should select atoms distributed along the "
        write(*, *) "full surface, and preferentially, with low atomic "
        write(*, *) "fluctuation along time."
        write(*, *) " "
        write(*, *) "Usage: s_index -in file.pdb"
        write(*, *) ""
        write(*, *) "file.pdb ---- atomic coordinates in PDB format"
        write(*, *) ""
        write(*, *) "Options: "
        write(*, *) ""
        write(*, *) "-residue         selects a full residue to fit the grid or"
        write(*, *) "to be used for the calculation of the density profile."
        write(*, *) "Otherwise the s_index will ask for specific atoms."
        write(*, *) ""
        write(*, *) "-sphere          defines compact systems (micelles,"
        write(*, *) "vesicles, etc.)"
        write(*, *) ""
        write(*, *) "-join            joins all indexes into a single file."
        write(*, *) "Creates the index files for the calculation of the"
        write(*, *) "density profile. This option must be used to calculate"
        write(*, *) "monolayers as opposed to bilayers."
        write(*, *) ""
        write(*, *) "-gromacs         converts density index files from GROMACS"
        write(*, *) "format to SuAVE format"
        write(*, *) ""
        write(*, *) "-bound           defines the position of placement of a"
        write(*, *) "divisor plane between the two leaflets of a bilayer. If the"
        write(*, *) "flag is not used, s_index will assume the boundary at the"
        write(*, *) "average distance between the two surfaces defined by the"
        write(*, *) "atoms listed in the index file."
        write(*, *) ""
        write(*, *) "-help            prints HELP information and quits."
        stop

     case('spher')

        write(*, *) ""
        write(*, *) ""
        write(*, *) "s_spher calculates the time-dependent area per molecule, the"
        write(*, *) "average radius, the radius of gyration and sphericity of "
        write(*, *) "compact interfaces."
        write(*, *) ""
        write(*, *) "Usage: s_spher -in file.pdb -ind1 file1.ndx -ind2 file2.ndx"
        write(*, *) "               -lipid N_lipids"
        write(*, *) ""
        write(*, *) "file.pdb ---- atomic coordinates in PDB format"
        write(*, *) "file.ndx ---- index file containing user-selected atoms used"
        write(*, *) "to fit the grid points to the chemical surface."
        write(*, *) "N_lipids ---- number of lipids composing the leaflet"
        write(*, *) ""
        write(*, *) "Options:"
        write(*, *) ""
        write(*, *) "-grid            generates a PDB file containing the grid"
        write(*, *) "points used in the fitting for the last frame in the "
        write(*, *) "trajectory file."
        call options_sphe()
        stop

     case('shell')

        write(*, *) ""
        write(*, *) ""
        write(*, *) "s_shell calculates the average thickness and volume "
        write(*, *) "for spherical shells along the simulation time."
        write(*, *) ""
        write(*, *) "Usage: s_shell -in file.pdb -ind1 file1.ndx -ind2 file2.ndx"
        write(*, *) ""
        write(*, *) "file.pdb ---- atomic coordinates in PDB format"
        write(*, *) "file.ndx ---- index file containing user-selected atoms used"
        write(*, *) "to fit the grid points to the chemical surface."
        write(*, *) ""
        write(*, *) "Options:"
        write(*, *) ""
        write(*, *) "-grid            generates a PDB file containing the grid"
        write(*, *) "points used in the fitting for the last frame in the "
        write(*, *) "trajectory file."
	write(*, *) ""
        write(*, *) "-range           defines the range specified in the XPM file"
        call options_sphe()
        stop

     case('bend')

        write(*, *) ""
        write(*, *) ""
        Write(*, *) "s_bend calculates the distribution of angles between the "
        Write(*, *) "vectors normal to the partition surfaces and the gradient "
        Write(*, *) "vector of the sphere circumscribing the compact surface. "
        write(*, *) "s_bend also calculates the angular deflection among these "
        write(*, *) "two vectors as an order parameter P(theta)"
        Write(*, *) "or SC."
        write(*, *) ""
        write(*, *) "Usage: s_bend -in file.pdb -ind file.ndx "
        write(*, *) ""
        write(*, *) "file.pdb ---- atomic coordinates in PDB format"
        write(*, *) "file.ndx ---- index file containing user-selected atoms used"
        write(*, *) "to fit the grid points to the chemical surface."
        write(*, *) ""
        write(*, *) "Options:"
        write(*, *) ""
        write(*, *) "-grid            generates a PDB file containing the grid"
        write(*, *) "points used in the fitting for the last frame in the "
        write(*, *) "trajectory file."
        write(*, *) ""
        write(*, *) "-range           defines the range specified in the XPM file"
        call options_sphe()
        stop

     case('densph')

        write(*, *) ""
        write(*, *) ""
        write(*, *) "s_densph calculates the density profile for compact surfaces"
        write(*, *) "taking into account their curvature"
        write(*, *) ""
        write(*, *) "Usage: s_densph -in file.pdb -ind1 file1.ndx -ind2 file2.ndx"
        write(*, *) "-dens dens.ndx"
        write(*, *) ""
        write(*, *) "file.pdb ---- atomic coordinates in PDB format"
        write(*, *) "file.ndx ---- index file containing user-selected atoms used" 
        write(*, *) "to fit the grid points to the chemical surface."
        Write(*, *) "dens.ndx ---- index file containing user-selected atoms used "
        Write(*, *) "to calculate density profile."
        write(*, *) ""
        write(*, *) "Options:"
        write(*, *) ""
        write(*, *) "-grid            generates a PDB file containing the grid"
        write(*, *) "points used in the fitting for the last frame in the "
        write(*, *) "trajectory file."
        write(*, *) ""
        write(*, *) "-slices          defines the number of slices along the axis"
        write(*, *) "normal to the system used to calculate the density profile"
        call options_sphe()
        stop

     case('gridsph')

        write(*, *) ""
        write(*, *) ""
        write(*, *) "s_gridsph builds a grid per frame throughout a trajectory file"
        write(*, *) "of compact surfaces."
        write(*, *) "Its output is useful to verify how accurate is the fitting of"
        write(*, *) "the calculated grid on the chemical surface.  "
        write(*, *) ""
        write(*, *) "Usage: s_gridsph -in file.pdb -ind file.ndx "
        write(*, *) ""
        write(*, *) "file.pdb ---- atomic coordinates in PDB format"
        write(*, *) "file.ndx ---- index file containing user-selected atoms used"
        write(*, *) "to fit the grid points to the chemical surface."
        write(*, *) ""
        write(*, *) "Options:"
        call options_sphe()
        stop
        
     case('stat')

        write(*, *) ""
        write(*, *) ""
        write(*, *) "s_stat performs central tendency and dispersion measurements,"
        write(*, *) "as well as assimetry and kurtosis upon the Probability "
        write(*, *) "Density Function describing the input data."
        write(*, *) ""
        write(*, *) "Usage: s_stat -in file.xvg"
        write(*, *) ""
        write(*, *) "file.xvg ---- contains a SuAVE output"
        write(*, *) ""
        write(*, *) "Options: "
        write(*, *) ""
        write(*, *) "-help            prints HELP information and quits"
        stop

     case('access')

        write(*, *) ""
        write(*, *) ""
        write(*, *) "s_access calculates accessible volumes within structures of "
        write(*, *) "interest using Monte Carlo integration. "
        write(*, *) ""
        write(*, *) "Usage: s_access -in file.pdb  -ind1 sys.ndx "
        write(*, *) "                -param par.dat" 
        write(*, *) ""
        write(*, *) "file.pdb ---- atomic coordinates in PDB format"
        Write(*, *) "sys.ndx  ---- index file containing the atoms composing the"
        Write(*, *) "system of interest"
        write(*, *) "par.dat  ---- file containing van der Waals sigma parameter"
        write(*, *) "[nm] for each atom being considered in sys.ndx file"
        write(*, *) ""
        write(*, *) "Options:"
        write(*, *) ""
        write(*, *) "-bin             defines the number of rectangular partition "
        write(*, *) "bins along the phi- and psi-angles"
        write(*, *) ""
        write(*, *) "-grid            generates a PDB file containing the grid"
        write(*, *) "points used in the fitting for the last frame in the "
        write(*, *) "trajectory file."
        write(*, *) ""
        write(*, *) "-coarse          generates a coarse grid over the surface"
        write(*, *) "index atoms from which a finer grid will be generated. This"
        write(*, *) "is recommended for surfaces defined by atoms which greatly"
        write(*, *) "fluctuate throughout the trajectory. "
        write(*, *) ""
        write(*, *) "-begin           first frame to use in the calculations"
        write(*, *) ""
        write(*, *) "-end             last frame to use in the calculations"
        write(*, *) ""
        write(*, *) "-skip            number of trajectory frames to be skipped"
        write(*, *) "during the analysis "
        write(*, *) ""
        write(*, *) "-atempt          defines number of atempts to calculate volume"
        write(*, *) ""
        write(*, *) "-const           defines the constant multiplying the radius"
        write(*, *) ""
        write(*, *) "-rough           percentage of the original surface roughness"
        write(*, *) "the user wants to keep (1.0 by default, meaning 100%)"
        write(*, *) ""
        write(*, *) "-help            prints HELP information and quits"
        stop

     case('filter')

        write(*, *) ""
        write(*, *) ""
        write(*, *) "s_filter performs a Discrete Fourier Transform (DFT)"
        write(*, *) "on the signal described in the input file and performs"
        write(*, *) "filtering treatments"
        write(*, *) ""
        write(*, *) "Usage: s_filter -in file.xvg"
        write(*, *) ""
        write(*, *) "file.xvg ---- contains a SuAVE output"
        write(*, *) ""
        write(*, *) "Options: "
        write(*, *) ""
        write(*, *) "-hfreq           defines the highest angular frequency to be"
        write(*, *) "used in the Inverse Fourier Transform (default: no filter)"
        write(*, *) ""
        write(*, *) "-lfreq           defines the lowest angular frequency to be"
        write(*, *) "used in the Inverse Fourier Transform (default: no filter)"
        write(*, *) ""
        write(*, *) "-lowint          defines the lowest intensity to be used in the"
        write(*, *) "Inverse Fourier Transform (default: -1)"
        write(*, *) ""
        write(*, *) "-param           defines the parameter used on the unidimensional" 
        write(*, *) "SuAVE fitting process. In this case the DFT is not applied"
        write(*, *) ""
        write(*, *) "-help            prints HELP information and quits"
        stop

     case('inertia')
        
        write(*, *) ""
        write(*, *) ""
        Write(*, *) "s_inertia calculates the Principal Moments of Inertia of the"
        write(*, *) "closed surface"
        write(*, *) ""
        write(*, *) "Usage: s_inertia -in file.pdb -ind file.ndx "
        write(*, *) ""
        write(*, *) "file.pdb ---- atomic coordinates in PDB format"
        write(*, *) "file.ndx ---- index file containing user-selected atoms used"
        write(*, *) "to fit the grid points to the chemical surface."
        write(*, *) ""
        write(*, *) "Options:"
        call options_sphe()
        stop
        
     case('count')
        
        write(*, *) ""
        write(*, *) ""
        write(*, *) "s_count calculates the number of molecules inside a specified"
        write(*, *) "structure surrounded by SuAVE grid."
        write(*, *) ""
        write(*, *) "Usage: s_count -in file.pdb -ind file1.ndx  -dens dens.ndx "
        write(*, *) ""
        write(*, *) "file.pdb ---- atomic coordinates in PDB format"
        write(*, *) "file.ndx ---- index file containing the atoms composing the"
        write(*, *) "the structure "
        Write(*, *) "dens.ndx ---- index file indicating the molecules that have"
        Write(*, *) "to be counted"
        write(*, *) ""
        write(*, *) "Options:"
        write(*, *) ""
        write(*, *) "-const           defines the constant multiplying the radius"
        write(*, *) ""
        write(*, *) "-molec           generates a PDB file containing the molecules"
        write(*, *) "observed inside the structure in each frame of the trajectory"
        call options_sphe()
        stop
        
     case('gauss')

        write(*, *) ""
        write(*, *) ""
        Write(*, *) "s_gauss calculates the Gaussian and Mean Curvatures for"
        write(*, *) "the interface" 
        write(*, *) ""
        write(*, *) "Usage: s_gauss -in file.pdb -ind file.ndx"
        write(*, *) ""
        write(*, *) "file.pdb ---- atomic coordinates in PDB format"
        write(*, *) "file.ndx ---- index file containing user-selected atoms used" 
        write(*, *) "to fit the grid points to the chemical surface."
        write(*, *) ""
        write(*, *) "Options:"
        write(*, *) ""
        write(*, *) "-grid            generates a PDB file containing the grid" 
        write(*, *) "points used in the fitting for the last frame in the "
        write(*, *) "trajectory file."
        call options_cart()
        stop
        
     end select

end subroutine write_help

