README

Instructions on how to calculate total stress, including gravitational and boundary stresses, from topography, crustal structure, and geodetic velocities using a sequence of Matlab programs.

Codes written by Hamish Hirschberg.
Citation: Hirschberg, H. P., S. Lamb, M. K. Savage (2018). Strength of an obliquely convergent plate boundary: lithospheric stress magnitudes and viscosity in New Zealand. Geophysical Journal International.
Based on the method of Flesch et al. (2001).

INTRODUCTION TO CODES:
moho2gpe.m calculates the vertically averaged gravitational potential energy (GPE) from topography and crustal structure.
grav_stress.m calculates the gravitational stress from the GPE.
basis_fns.m calculates the stress basis functions from the region geometry.
tot_stress.m calculates the total deviatoric stress.
stress_calc.m runs the above code using a script intended to be edited by the user.

HOW TO USE CODES:
The codes are written as functions in Matlab and are designed to be used as functions in the script 'stress_calc.m' or in the Matlab command window. The inputs and outputs are designed to work with either GMT or Matlab.
The example inputs provided are:
topo.xyz - topographic elevation in m
moho.xyz - crustal thickness/depth to Moho in km
vel.xyz  - velocities (east and north) in mm/yr

NOTES ON CONVENTIONS:
Extension is taken as positive.
The magnitude for stress and strain rate is taken as the square-root of the second invariant, i.e. sqrt(sum(s_ij*s_ij)/2) summed over i and j.
Effective viscosity is calculated as total deviatoric stress magnitude divided by strain rate magnitude, i.e. T/E.

************************************************************************
FULL WORKFLOW SCRIPT - stress_calc.m

This script provides a workflow for the above programs in Matlab, calling them in order and passing the required variables between them.

The script allows the user to change the name of the input and output files and input constants. It also allows the user to set a model name that is suffixed to the output file names to keep them distinguishable, along with setting a directory for the output different from the working directory (where the programs are).

How to use:
Open 'stress_calc.m' in the Matlab editor.
Change input file names, constants, working directory, and the model name to suit.
Run 'stress_calc.m'.
If the plotting script 'map.txt' does not work on your system, replace it with the plotting method of your choice.

INPUT/OUTPUT CONSIDERATIONS
The following should be taken into account in regards to the input and output:
- Input files are assumed to have a constant spacing (e.g. 0.5 deg long, 0.4 deg lat) on a regular, rectangular grid.
- All input files must use the same grid. Output files use the same grid.
- No smoothing (pre or post) is done by any of the functions. If smoothing is desired, it can be done separately.
- Note the conventions used in the outputs (see "Notes on Conventions" above).

EXAMPLE MAP - map.txt
The script 'map.txt' is provided as an example of how the output can be plotted using GMT. It is written in bash. It is set up to plot the total deviatoric stress calculated from the sample data but can be modified to plot other output results and other regions. Its output is 'map_model.ps' for a model named '_model'.

************************************************************************
INDIVIDUAL PROGRAMS

************************************************************************
moho2gpe.m - GPE from topography and crustal structure

Input Variables: (use as many in order as desired; use [] to skip)
1. topography in m as xyz file name or Matlab variable (required)
2. output name of output xyz file (requires >=1 input xyz file)
3. depth in km GPE averaged down to as a scalar
4. crustal density in g/cm3 as scalar, xyz file name, or Matlab variable
5. Moho depth in km as xyz file name or Matlab variable
6. mantle density in g/cm3 as scalar
7. points at sea (1) vs points on land (0) as xyz file or Matlab variable (use if region contains land below sea level)

Output Variable (optional):
1. GPE in MPa as Matlab variable
To output as a text file, use input variable 2 which specifies output file name.

Default values:
Depth averaged to: 25 km
Crustal density: 2.76 g/cm3
Moho depth: default is to consider only crust
Mantle density: 3.2 g/cm3
Offshore points: default is all points with negative elevation (below sea level)

Examples of how to run:
moho2gpe('topo.xyz','gpe.xyz') - calculate GPE based on just topography from the file 'topo.xyz' with a constant density crust of 2.76 g/cm3 averaged over top 25 km and outputted to the file 'gpe.xyz'.
moho2gpe('topo.xyz','gpe.xyz',30,2.8) - calculate GPE as above but averaged over top 30 km with a constant density crust of 2.8 g/cm3.
moho2gpe('topo.xyz','gpe.xyz',100,[],'moho.xyz',3.3) - calculate GPE averaged over top 100 km using the default crustal density, the Moho depth given by the file 'moho.xyz', a mantle density of 3.3 g/cm3, and outputted to the file 'gpe.xyz'.
gpe=moho2gpe('topo.xyz',[],100,[],'moho.xyz',3.3,'sea.txt') - calculate GPE as above but outputting to the Matlab variable 'gpe', no output text file, and defining offshore points using the file 'sea.txt' (use if region contains land below sea level).

************************************************************************
grav_stress.m - gravitational stress from GPE

Input Variables (choose a below combination)
1. xyz file or Matlab variable containing long, lat, and GPE
OR
1. xyz file or Matlab variable containing long, lat, and GPE AND
2. name of output file in xyz-based format
OR
1. Matlab vector or m x n matrix containing long AND
2. Matlab vector or m x n matrix containing lat AND
3. Matlab vector or m x n matrix containing GPE, each with same dimensions
OR
1. - 3. as above AND
4. name of output file in xyz-based format

The input longitude, latitude, and GPE need to be for grid points in a rectangular grid with constant longitudinal and latitudinal spacing at the desired calculation spacing. Constructing a gridfile (.grd) in GMT and converting it to xyz is useful for getting this.

Output Variables (optional)
1. m x n x 3 Matlab variable containing xx, yy, and xy components of deviatoric gravitational stress

To output as a text file, provide the name of the file in the input after the other variable(s).
The output file columns follow the follow format (designed for use with GMT's psvelo -Sx):
Column 1: Longitude
Column 2: Latitude
Column 3: Most extensional stress (SHmin)
Column 4: Most compressional stress (SHmax)
Column 5: Azimuth of SHmax
Column 6: Stress magnitude - sqrt(0.5*tau_ij*tau_ij)

Examples of how to run:
grav_stress('gpe.xyz','grav.xyz') - calculate gravitational stress from input GPE 'gpe.xyz' and outputting to text file 'grav.xyz'.
grav=grav_stress(long,lat,gpe) - calculate gravitational stress using grid points at 'long' and 'lat' with GPE of 'gpe', outputting the result as the Matlab variable 'grav'.

************************************************************************
basis_fns.m - basis functions from region geometry

Input: EITHER
1. region for basis functions as xy(z) file or Matlab variable
OR
1. Matlab vector or m x n matrix containing long AND
2. Matlab vector or m x n matrix containing lat of region

The input longitudes and latitudes need to be for grid points in a rectangular grid with constant longitudinal and latitudinal spacing. For later steps, these grid points need to match with the grid points for which the gravitational stress was calculated and match that used for the velocity data. An easy way to ensure this is to use the output file of grav_stress.m as input.

Output:
1. m x n x k x 3 Matlab variable containing k=6*(m+n)-9 stress basis functions

The output of this function is designed to be used by tot_stress.m as a Matlab variable.

Examples of how to run:
bases=basis_fns('gpe.xyz') - calculate basis functions at the longitudinal and latitudinal spacing defined by the file 'gpe.xyz', storing them as the Matlab variable 'bases'.
bases=basis_fns(long,lat) - calculate basis functions at the longitudinal and latitudinal spacing of 'long' and 'lat', storing them as the Matlab variable 'bases'.

************************************************************************
tot_stress.m - calculate total deviatoric stress

This program requires grav_stress.m and basis_fn.m to have been run previously, with the gravitational stress and basis functions stored as Matlab variables (see example workflow).

Requires gravitational stress, basis functions, and velocities to be gridded according to the same rectangular grid.

Input:
1. m x n x 3 Matlab variable of gravitationally-induced stresses
2. m x n x k x 3 Matlab variable of stress basis functions
AND EITHER
3. xyz-based file with velocities in mm/yr
OR
3. m x n x 3 Matlab variable of strain rates in units of 10^-12/s
WITH OPTIONAL (use [] to skip if using later options; requires xyz-based velocity file)
4. name of output xyz-based file of total deviatoric stress
5. name of output xyz file of misfit
6. name of output xyz file of effective viscosity
7. name of output xyz-based file of boundary stress
8. name of output xyz-based file of strain rates
AND OPTIONAL (use [] to skip if using later options)
9. minimum strain rate second invariant in units of 10^-12/s which is considered reliable (default: 0, i.e. no minimum)
10. maximum misfit considered reliable (default: 1, i.e. no maximum)

Output Variables (Optional):
1. m x n x 3 Matlab variable of total deviatoric stress
2. m x n matrix of misfit
3. m x n matrix of effective viscosity
4. m x n x 3 Matlab variable of boundary stress
5. m x n x 3 Matlab variable of strain rates

Output text file format for total stress, boundary stress, and strain rates:
Column 1: Longitude
Column 2: Latitude
Column 3: Most extensional stress/strain rate (SHmin)
Column 4: Most compressive stress/strain rate (SHmax)
Column 5: Azimuth of SHmax
Column 6: Magnitude of stress/strain rate - sqrt(sum(s_ij*s_ij)/2) summed over i and j

Examples of how to run:
tot_stress(grav,bases,'vel.txt','tot.txt') - calculate the total deviatoric stress with the gravitational stresses 'grav', the stress basis functions 'bases', and velocity file 'vel.txt', writing the result to 'tot.txt'.
tot_stress(grav,bases,'vel.txt','tot.txt','mis.txt','visc.txt') - as above, but also writing the stress-strain rate orientations misfit to 'mis.txt' and the effective viscosity to 'visc.txt'.
[tot,mis,visc]=tot_stress(grav,bases,strain) - calculate the total deviatoric stress from 'grav' and 'bases' as above and using strain rates 'strain', storing as Matlab variables the total stress 'tot', the misfit 'mis', and the viscosity 'visc'.

************************************************************************
SUPPORTING FUNCTIONS

Other, supporting functions are used. These are generally not called directly by the user (although they can be); instead they are called by the main functions. Information on them is supplied below.

***
totlincom.m - function to calculate the functional that is minimised in tot_stress.m
This functional is equivalent to the difference between the total deviatoric stress and strain rate orientations integrated over the region of interest. By minimising this, tot_stress.m aims to find the closest solution to an isotropic effective viscosity.
L=totlincom(grav,bases,strain,coeffs) - calculate the value of the functional 'L' based on the gravitational stress 'grav', the basis functions 'bases' with coefficients 'coeffs', and the strain rates 'strain'.

***
princaxes.m - function to calculate the principal stresses from components
The magnitudes and orientations of the principal axes of stress are calculated for plotting using GMT's psvelo. They are converted from either an m x n x 3 Matlab variable or 3 m x n Matlab variables with the xx, yy, and xy horizontal components of stress for m x n points. They can be outputted as magnitudes and orientation (for plotting in GMT) or x- and y-components (for plotting in Matlab).
Works for any symmetric 2 x 2 tensor, such as horizontal stress or strain.
axes=princaxes(stress) - calculate the horizontal principal axes of stress and orientation 'axes' for the stresses 'stress'
[axes,comps]=princaxes(xx,yy,xy) - calculate the principal axes, outputted as magnitude and orientation 'axes' and as x- and y-components 'comps'.

















