%% initialize
% bugfix for library problems
setenv('LD_LIBRARY_PATH', ['/usr/lib/x86_64-linux-gnu:',getenv('LD_LIBRARY_PATH')]);

% DO YOU WANT TO RUN MULTIPLE SIMULATIONS?
    general.RunMultipleSimulations = 0;                     % Input needed
% IS THIS A MAC (MAC OS)? For MACOS 1=yes 0=no
    general.switchMAC = 0;                                  % Input needed 
    
%% Check for multiple runs
if ~general.RunMultipleSimulations
%% Specify Input
clear ans directions flapParam mainParam meshParam solvParam

% SIMULATION NAME %
    directions.simName = 'Optimization_testParallel';                 % Input needed
    
% MESHER %
    % currently only 'snappyHexMesh' available
    directions.mesh.Type = 'snappyHexMesh';                 % Input needed [SHM, cfMesh, gMesh]
    directions.mesh.Szenario = 'Sebastiano_1';              % Input needed [refinement box, distance, simplegrading]
    directions.mesh.InputFile = 'SHM_param_case1.xlsx';     % Input needed
    % Sheet, Start and End of Mesher ExcelImport
    directions.mesh.InputSheet = 'Mesher';                  % Input needed [coarse, medium, fine, etc..]
    directions.mesh.InputStart = 'A2';                      % Input needed
    directions.mesh.InputEnd = 'F55';                       % Input needed
    
% SOLVER %
    % currently only 'simpleFoam' available
    directions.solv.Type = 'simpleFoam';                    % Input needed
    directions.solv.Szenario = 'kwSST_modified';              % Input needed [kwSST, SpalartAllmaras]
    directions.solv.InputFile = 'SimpleFoam_param_case1.xlsx';     % Input needed
    % Sheet, Start and End of Solver ExcelImport
    directions.solv.InputSheet = 'solver';                  % Input needed
    directions.solv.InputStart = 'A2';                      % Input needed
    directions.solv.InputEnd = 'F55';                       % Input needed   

% OPTIMIZATION %
    directions.solv.InputFileOpt = 'Optim_param_case1.xlsx';
    directions.solv.InputSheetOpt = 'kwSST';
    directions.solv.InputStartOpt = 'A2';
    directions.solv.InputEndOpt = 'F30';
    
% AIRFOIL FILES - dont add data ending / no .txt etc
    general.twoElemAirfoil = 1;   % if not just use mainAirfoil  % Input needed
    MainAirfoilName = 'NLR7301';                            % Input needed
    FlapAirfoilName = 'NLR7301FLAP';                        % Input needed
    % Do you want to have the two airfoil files separated? (1 = yes, 0 = no)
    general.flapSeparated = 1;                              % Input needed    
% MESHING OPTIONS
    % Option1-4: snappyHexMesh
    % Option1: blockMesh + surfaceFetExtr + snappyHex (overwr) + checkMesh
    % Option2: blockMesh + snappyHex
    % Option3: Remesh
    % Option4: Mesh in Parallel
    general.meshingOption = '4';                            % Input needed 
    general.meshParallel  = 0;                              % Input needed
%     general.meshCores     = 8;                              % Input needed
    general.meshDimension = '2D';                           % Input needed [2D, 3D] (3D not implemented)
% SOLVER OPTIONS
    % Option1: Solve
    general.solverOption  = '1';                            % Input needed 
    % Do YOU WANT TO OPEN THE PARAVIEW GUI?
    general.openParaView = 0;                               % Input needed
end
%% directories

% MAC                                  
    if general.switchMAC
        directionsMAC.path.homeDirectory='/Users/sgb';          % Input needed
        directionsMAC.OpenFoam.workingDirectory='workingDir';   % Input needed
        generalMAC.command.docker = '/usr/local/bin/docker exec -it -u ofuser of_v2012 sh -c ';    % Input needed for usr/local/bin
        generalMAC.command.cd  = append(generalMAC.command.docker, '"cd '); % Pay attention to closing parenthesis (") at end of terminal command
    end

% set directory paths
    directions.path.openfoamtool = pwd;
    directions.path.simulation = append(pwd,'/Simulation_Outputs/',directions.simName);
    directions.path.mesher = append(directions.path.simulation,'/',directions.mesh.Type);
    directions.path.geometry = append(directions.path.mesher,'/constant/triSurface');
    directions.path.solver = append(directions.path.simulation,'/',directions.solv.Type);
    directions.path.inputsGeo = append(pwd,'/Simulation_Inputs/airfoil_coordinates');
    directions.path.inputsSim = append(pwd,'/Simulation_Inputs/simulation_parameters');
    directions.path.inputsExp = append(pwd,'/Simulation_Inputs/experimental_data');
    
        % MAC
    if general.switchMAC
        directionsMAC.path.mesher = erase(directions.path.mesher, directionsMAC.path.homeDirectory);  
        directionsMAC.path.mesher = [directionsMAC.OpenFoam.workingDirectory, directionsMAC.path.mesher];  

        directionsMAC.path.solver = erase(directions.path.solver, directionsMAC.path.homeDirectory);
        directionsMAC.path.solver = [directionsMAC.OpenFoam.workingDirectory, directionsMAC.path.solver];
    end
    
% include work and input folders
    addpath (genpath(append(pwd,'/Mat_OF_library')));

% Set up new working directories - Create new simulation directory
    system(['mkdir ',directions.path.simulation]);
% Copy mesher and solver files into folder - rename them
if ~general.switchMAC
    
    system(['cp -ar ',directions.path.openfoamtool,'/Templates/mesher/',directions.mesh.Type,'/',directions.mesh.Szenario,' ',directions.path.simulation]);
    system(['mv ',directions.path.simulation,'/',directions.mesh.Szenario,' ',directions.path.mesher]);
    system(['cp -ar ',directions.path.openfoamtool,'/Templates/solver/',directions.solv.Type,'/',directions.solv.Szenario,' ',directions.path.simulation]);
    system(['mv ',directions.path.simulation,'/',directions.solv.Szenario,' ',directions.path.solver]);
    
else
    
    system(['cp -r ',directions.path.openfoamtool,'/Templates/mesher/',directions.mesh.Type,'/',directions.mesh.Szenario,' ',directions.path.simulation]);
    system(['mv ',directions.path.simulation,'/',directions.mesh.Szenario,' ',directions.path.mesher]);
    system(['cp -r ',directions.path.openfoamtool,'/Templates/solver/',directions.solv.Type,'/',directions.solv.Szenario,' ',directions.path.simulation]);
    system(['mv ',directions.path.simulation,'/',directions.solv.Szenario,' ',directions.path.solver]);

end
% exclude old simulations, include new simulation
    rmpath (genpath(append(directions.path.openfoamtool,'/Simulation_Outputs')));
    rmpath (genpath(append(directions.path.openfoamtool,'/Simulation_Inputs')));
    addpath (genpath(directions.path.simulation)); 

%% save command window output
diary (append(directions.path.simulation,'/CommandHistory.txt'));

%% Conversion from .dat to .stl

% Ask for Extrusion Parameters
    mainPath = append(directions.path.inputsGeo,'/OpenSCAD_mainParam.mat');
    load(mainPath);
    mainParam.airfoilName = MainAirfoilName;

if general.twoElemAirfoil
    flapPath = append(directions.path.inputsGeo,'/OpenSCAD_flapParam.mat');
    load(flapPath);
    flapParam.airfoilName = FlapAirfoilName;
end

    clear mainPath flapPath MainAirfoilName FlapAirfoilName;

% run conversion into stl
    % main airfoil conversion
    dattostl(mainParam, directions.path.inputsGeo, directions.path.geometry, general.switchMAC);

    % flap conversion
    dattostl(flapParam, directions.path.inputsGeo, directions.path.geometry, general.switchMAC);
    
    % stitch airfoils together, if there are two airfoils and wanted
    if ~general.twoElemAirfoil            % check for two Element airfoil
       general.flapSeparated = 0;
    end
    
    if general.flapSeparated == 0
        % read flap file and save to flap vaiable
        % flapID = fopen(append(flapParam.airfoilName,'.stl'),'r');
        flap = fileread(append(flapParam.airfoilName,'.stl'));

        % open mainAirfoil
        mainID = fopen(append(directions.path.mesher,'/constant/triSurface/',mainParam.airfoilName,'.stl'),'a');

        % append flap stl to end of main Airfoil
        if mainID > 0 
            fprintf(mainID,flap);
        end
    
        % and then delete the old file
        system(['rm -f ',directions.path.mesher,'/constant/triSurface/',flapParam.airfoilName,'.stl']);
   
        % clear local variables
    end

    clear main mainID flap flapID flapParam mainParam

    
%% Set Mesher input
    
    % import Input from Excel file
    meshParam = importSearchNReplaceSheet(directions.path.inputsSim,directions.mesh.InputFile,directions.mesh.InputSheet,directions.mesh.InputStart,directions.mesh.InputEnd);

    % write mesher input into files
    varNr = height(meshParam);                  % number of variables to be replaces

for i = 1:varNr                                 % run through different file names
    
    if meshParam{i,4} ~= ""
    % DataType does not need to be set, since system data does not have file ending
        searchnreplace (meshParam {i,2},directions.path.mesher,'',meshParam{i,3},meshParam{i,4});
    end
end

 clear i varNr

 %% Set solver inputs

    % import Input from Excel file
    solvParam = importSearchNReplaceSheet(directions.path.inputsSim,directions.solv.InputFile,directions.solv.InputSheet,directions.solv.InputStart,directions.solv.InputEnd);

    % write solver input into files
    varNr = height(solvParam);               % number of variables to be replaces

for i = 1:varNr                              % run through different file names
     
    if solvParam{i,4}~= ""
    % DataType does not need to be set, since system data does not have file ending
        searchnreplace (solvParam {i,2},directions.path.solver,'',solvParam{i,3},solvParam{i,4});
    end    
end

 clear i varNr 

%% run mesher

switch general.meshingOption
    
    case '1'
        % MESHING OPTION 1 - Sebastiano Option
        if general.switchMAC         
            !/usr/local/bin/docker start of_v2012;
            system([generalMAC.command.cd, directionsMAC.path.mesher, ' && openfoam blockMesh"']);
            system([generalMAC.command.cd, directionsMAC.path.mesher, ' && openfoam surfaceFeatureExtract"']);
            system([generalMAC.command.cd, directionsMAC.path.mesher, ' && openfoam snappyHexMesh -overwrite"']);
            system([generalMAC.command.cd, directionsMAC.path.mesher, ' && openfoam checkMesh -latestTime | tee checkMesh3D.log"']);
        else
            system(['cd ',directions.path.mesher,' && blockMesh']);
            system(['cd ',directions.path.mesher,' && surfaceFeatureExtract']);
            system(['cd ',directions.path.mesher,' && snappyHexMesh -overwrite']);
            system(['cd ',directions.path.mesher,' && checkMesh -allGeometry -latestTime']);
        end

    case '2'
        % MESHING OPTION 2
        system(['blockMesh && snappyHexMesh – overwrite && checkMesh -allGeometry -allTopology -writeAllFields']);
        
    case '3'
        % MESHING OPTION 3
        system(['rm -r constant/polyMesh/ && blockMesh && snappyHexMesh -overwrite && checkMesh -allGeometry -allTopology -writeAllFields']);

    case '4'    
        % MESHING OPTION 4 - Calum Douglas, Youtube Tutorial
        % decompose, mesh parallel, add together, delete processors folders
        system(['cd ',directions.path.mesher,' && blockMesh']);
        system(['cd ',directions.path.mesher,' && surfaceFeatureExtract']);
        system(['cd ',directions.path.mesher,' && decomposePar']);
        system(['cd ',directions.path.mesher,' && mpirun -np 4 snappyHexMesh -overwrite -parallel']);
        system(['cd ',directions.path.mesher,' && reconstructParMesh -constant']);
        system(['cd ',directions.path.mesher,' && rm -r processor*']);

end

%% Switch 2D - 3D
switch general.meshDimension
    case '2D'
        if general.switchMAC   
            system([generalMAC.command.cd, directionsMAC.path.solver, ' && openfoam extrudeMesh"']);
            system([generalMAC.command.cd, directionsMAC.path.solver, ' && openfoam changeDictionary"']);
            system([generalMAC.command.cd, directionsMAC.path.solver, ' && openfoam checkMesh -latestTime | tee checkMesh2D.log"']);    
        else
            system(['cd ',directions.path.solver,' && extrudeMesh']);
            system(['cd ',directions.path.solver,' && changeDictionary']);
            system(['cd ',directions.path.solver,' && checkMesh -allGeometry -latestTime']);  
        end
    case '3D'
        % to be defined 
end
 
%% run solver

switch general.solverOption
    
    case '1'
        % Solver OPTION 1 - Sebastiano1
        if general.switchMAC
            !/usr/local/bin/docker start of_v2012;
            system([generalMAC.command.cd, directionsMAC.path.solver, ' && openfoam simpleFoam"']);
        else
            run cmaes.m
        end

end

%% run post

% import 
% addpath(genpath(append(directions.path.solver,'/postProcessing')));
% run importForceCoeffs.m
% run cp_distribution.m
% run importYPlus.m

% TimeArr = table2array(forceCoeffs(101:end,1));
% Cm = table2array(forceCoeffs(101:end,2));
% Cd = table2array(forceCoeffs(101:end,3));
% Cl = table2array(forceCoeffs(101:end,4));
% Clf = table2array(forceCoeffs(101:end,5));
% Clr = table2array(forceCoeffs(101:end,6));
% 
% figure
% 
% hold on
% 
% plot (TimeArr,Cm);
% plot (TimeArr,Cd);
% plot (TimeArr,Cl);
% plot (TimeArr,Clf);
% plot (TimeArr,Clr);
% 
% legend('Cm','Cd','Cl','Clf','Clr');
% title(directions.simName);
% saveas(figure(2),[directions.path.simulation,'/',directions.simName, '.fig']);
% hold off

if general.openParaView
    system(['cd ',directions.path.solver, ' && paraFoam']);
end

clear Cd Cl Clf Clr Cm TimeArr