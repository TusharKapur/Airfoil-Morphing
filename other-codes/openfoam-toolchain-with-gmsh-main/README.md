For Ubuntu 18.04 bionicbeaver necessary programs:

    // install openscad with //
        sudo add-apt-repository ppa:openscad/releases
        sudo apt-get update
        sudo apt-get install openscad

    // install matlab //

    // install openfoam with snappyHexmesh and paraView //


GOOD TO KNOW:

    - add an alias into .bashrc to automatically start the toolchain in the 
      right directory or just always click Run_Toolchain.m in the folder:
      alias mat_tool='cd /home/intern/OpenFoamToolchain/'
    
    - Paths are without '/' in the End, and usually in the complete form. 
      paths can be used straight away, but when concenating paths add '/'
      Try to always specify the whole path.

    - opened in the working directory 'Matlab_Toolchain',
      the Header of 'Run_Toolchain' recursively adds all necessary subfolders to searchpath 

    - no file with 'airfoilXYZX.txt' in the directories pls

    - When choosing mesher and solvers -> first meshertype and then SetupXY 
      (Setup Template has to be stored in Matlab_Toolchain/simulation_templates)

    - Warnings with "... not found in path" can be ignored

    - Dont rotate/translate with OpenSCAD

PREPARATION:

    - Set up your simulation Scenarios/Directories in the 
      /Matlab_Toolchain/simulation_templates in the corresponding folder

    - Use the Simulation_inputs/airfoil_coordinates folder to store relevant
      airfoil geometries in x/y coordinates

    - Use the Simulation_inputs/simulation_parameters folder to change the 
      simulation Parameters for the mesher and solver

    - If you want to run multiple Simulations set general.RunMultipleSimulations = 1; 
      Otherwise it must be 0.

    - For MULTIPLE SIMULATIONS specify the Simulation parameters in:
      Simulation_Inputs/simulation_parameters/MultipleSimulations.xlsx

HOW TO RUN:

    - Run in Terminal: mat_tool 
      (or your specified alias or run Matlab and the startup script in the OpenFoamToolchain folder)

    - Specify everything of the Toolchain in the first Chapter
        - Parameters to Translate or Rotate with OpenSCAD are stored in 
          Simulation_inputs/airfoil_coordinates/OpenSCAD_flapParam.mat

    - For MULTIPLE SIMULATIONS: use Run_multipleSimulations.m

INFO FOR OPTION SELECTION:

    - In Simulation_inputs/simulation_parameters/SHM_param_case* the mesh is defined. 
      Every sheet is a different mesh.  
      The mesh can be selected either in the Excel-sheet MultipleSimulation "Sheet" field, 
      or in the RunToolchain.mat "InputSheet" variable. 

    - Similarly, also different solver options are defined in the Excel-Sheet SimpleFoam_param_case*, in each sheet. 

    - Important: definition of airfoil name must be consistent across Matlab, SnappyHexMesh and SimpleFoam. 

    - Until now the following options are investigated:

	  - SnappyHexMesh-ExcelSheet
		  * computational domain 
		  * maximum refinement level
		  * layers (number, thickness) 

	  - SimpleFoam-ExcelSheet
		  * boundary types and conditions 
		  * aoA and U_freestream

OPEN TO DOS:

    - Implement the Matlab triangulation(translate,rotate) and stl export

    - Post processing:
        * Residuals (k, w, Ux, Uy, Uz, p, nut)
        * Pressure coefficient
        * Residuals (k, w, Ux, Uy, Uz, p, nut)
        * Pressure coefficient

ADDITIONAL INSTRUCTIONS:

    - ADDITIONAL PREREQUISITES
        * Install Gmsh (version 4.4.1)
        * Openfoam version - v2012 from https://www.openfoam.com/
    - FOR SWITCHING BETWEEN 2D AND 3D
        * In RunToolchain.m - directions.solv.Type = 'simpleFoam-2D' or 'simpleFoam-3D' 
        * In RunToolchain.m - general.meshDimension = 2D or 3D
        * In gmsh_v1.xls - change extrude_thickness and no_layers. In case of 2D, no_layers should be set = 1.
    - Notes
        * In case of Gmsh, the mesh parameters can be changed in the gmsh_v1.xls. It contains values for variables that are replaced in "airfoil.geo" (the configuration file for Gmsh), which is created by the function "dattogeo".
        * There are two different directories for the solver: "simpleFoam-2D" and "simpleFoam-3D", for the 2D and 3D cases respectively. This is because of the difference in front and back boundary conditions.
        * The code has currently not been implemented to be run on MAC OS. 
    




