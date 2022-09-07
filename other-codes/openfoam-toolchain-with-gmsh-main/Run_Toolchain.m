%% initialize
% bugfix for library problems
setenv('LD_LIBRARY_PATH', ['/usr/lib/x86_64-linux-gnu:',getenv('LD_LIBRARY_PATH')]);

% DO YOU WANT TO RUN MULTIPLE SIMULATIONS?
    general.RunMultipleSimulations = 0;                     % Input needed
% IS THIS A MAC (MAC OS)? For MACOS 1=yes 0=no
% Check if docker deamon is running before running the script
    general.switchMAC = 0;                                  % Input needed 
    
%% Check for multiple runs
if ~general.RunMultipleSimulations
    %% Specify Input
    clc
    clear ans directions flapParam mainParam meshParam solvParam

    % SIMULATION NAME %
        directions.simName = 'test_2';                 % Input needed

    % MESHER %
        % 'snappyHexMesh' and 'gmsh' available
        directions.mesh.Type = 'gmsh';                 % Input needed [SHM, cfMesh, gMesh]
        directions.mesh.Szenario = 'v6'; % Input needed in case of snappyHexMesh [refinement box, distance, simplegrading]
        directions.mesh.InputFile = 'gmsh_v1.xlsx';     % Input needed 
        % Sheet, Start and End of Mesher ExcelImport
        directions.mesh.InputSheet = 'Sheet1';                  % Input needed [coarse, medium, fine, etc..]
        directions.mesh.InputStart = 'A2';                      % Input needed
        directions.mesh.InputEnd = 'F31';                       % Input needed
        
    % SOLVER %
        % currently only 'simpleFoam' available
        directions.solv.Type = 'simpleFoam-2D';                    % Input needed
        directions.solv.Szenario = '4_freestream_kwSST_cp';              % Input needed [kwSST, SpalartAllmaras, freestream, etc..]
        directions.solv.InputFile = 'SimpleFoam_param_case2_nlf26.xlsx';     % Input needed
        % Sheet, Start and End of Solver ExcelImport
        directions.solv.InputSheet = 'solver6';                  % Input needed
        directions.solv.InputStart = 'A2';                      % Input needed
        directions.solv.InputEnd = 'F30';                       % Input needed 
        
    % AIRFOIL FILES - dont add data ending / no .txt etc
        general.twoElemAirfoil = 1;   % if not just use mainAirfoil  % Input needed
        MainAirfoilName = 'bladenlf26';                            % Input needed 
        FlapAirfoilName = 'bladenlf26FLAP';                        % Input needed 
        % Do you want to have the two airfoil files separated? (1 = yes, 0 = no)
        general.flapSeparated = 0;                              % Input needed  
        
    % MESHING OPTIONS
        % Option1-4: snappyHexMesh
        % Option1: blockMesh + surfaceFetExtr + snappyHex (overwr) + checkMesh
        % Option2: blockMesh + snappyHex
        % Option3: Remesh
        % Option4: Mesh in Parallel
        % Option5: gmsh
        general.meshingOption = '5';                            % Input needed 
        general.meshParallel  = 0;                              % Input needed
%         general.meshCores     = 8;                              % Input needed
        general.meshDimension = '2D';                           % Input needed [2D, 3D] (3D not implemented)
        
    % SOLVER OPTIONS
        % Option1: Solve
        % Option2: Run solver in Parallel
        general.solverOption  = '2';                            % Input needed 
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
    
% MAC
if general.switchMAC
    directionsMAC.path.mesher = erase(directions.path.mesher, directionsMAC.path.homeDirectory);  
    directionsMAC.path.mesher = [directionsMAC.OpenFoam.workingDirectory, directionsMAC.path.mesher];  

    directionsMAC.path.solver = erase(directions.path.solver, directionsMAC.path.homeDirectory);
    directionsMAC.path.solver = [directionsMAC.OpenFoam.workingDirectory, directionsMAC.path.solver];
end
    
% include work and input folers
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
if general.meshingOption ~= '5'
    % main airfoil conversion
    dattostl(mainParam, directions.path.inputsGeo, directions.path.geometry, general.switchMAC);

    % flap conversion
    dattostl(flapParam, directions.path.inputsGeo, directions.path.geometry, general.switchMAC);
elseif general.meshingOption == '5'
    dattogeo(mainParam, flapParam, directions.path.inputsGeo, directions.path.mesher, general.switchMAC);
    % system(['rm ',directions.path.mesher,' airfoil.geo airfoil.msh']);
    % system(['touch ',directions.path.mesher,' airfoil.geo']);
    % dattogeo
    % system(['rm dattogeo.m']);
    % system(['cp dattogeo-old.m dattogeo.m']);
end
    % stitch airfoils together, if there are two airfoils and wanted
    if ~general.twoElemAirfoil            % check for two Element airfoil
       general.flapSeparated = 0;
    end


    if general.flapSeparated == 0 & general.meshingOption ~= '5'
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

    % write mesher input into files
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
            system([generalMAC.command.cd, directionsMAC.path.mesher, ' && openfoam checkMesh -constant -allGeometry"']);
        else
            system(['cd ',directions.path.mesher,' && blockMesh']);
            system(['cd ',directions.path.mesher,' && surfaceFeatureExtract']);
            system(['cd ',directions.path.mesher,' && snappyHexMesh -overwrite']);
            system(['cd ',directions.path.mesher,' && checkMesh -constant -allGeometry']);
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
    
    case '5'    
        % MESHING OPTION 5 - Using Gmsh
        % system(['mv ',directions.path.mesher,'/airfoil.geo ',directions.path.simulation]);
        system(['rm ',directions.path.mesher,'/airfoil.msh']);
        % system(['cp dattogeo-old.m dattogeo.m']);
        system(['cd ',directions.path.mesher,' && gmsh airfoil.geo -3 airfoil.msh -format msh2']); % -format msh2
        system(['cd ',directions.path.mesher,' && gmshToFoam airfoil.msh']);
        system(['cd ',directions.path.mesher,' && checkMesh']);
        % system(['cd ',directions.path.mesher,' && changeDictionary']);
        system(['cd ',directions.path.solver,'/constant/',' && mkdir polyMesh']);
        system(['cd ',directions.path.mesher,' && cp constant/polyMesh/* ',directions.path.solver,'/constant/polyMesh/.']);


end

%% Switch 2D - 3D
    switch general.meshDimension
        case '2D'
            if general.switchMAC   
                if  general.meshingOption ~= '5' 
                    system([generalMAC.command.cd, directionsMAC.path.solver, ' && openfoam extrudeMesh"']);
                end
                system([generalMAC.command.cd, directionsMAC.path.solver, ' && openfoam changeDictionary"']);
                system([generalMAC.command.cd, directionsMAC.path.solver, ' && openfoam checkMesh -latestTime | tee checkMesh2D.log"']);    
            else
                if  general.meshingOption ~= '5' 
                    system(['cd ',directions.path.solver,' && extrudeMesh']);
                end
                system(['cd ',directions.path.solver,' && changeDictionary']);
                system(['cd ',directions.path.solver,' && checkMesh -allGeometry -latestTime']);  
            end
        case '3D'
            if general.switchMAC   
                system([generalMAC.command.cd, directionsMAC.path.solver, ' && openfoam changeDictionary"']);
                system([generalMAC.command.cd, directionsMAC.path.solver, ' && openfoam checkMesh -latestTime | tee checkMesh2D.log"']);    
            else
                system(['cd ',directions.path.solver,' && changeDictionary']);
                system(['cd ',directions.path.solver,' && checkMesh -allGeometry -latestTime']);  
            end
    end

%% run solver

switch general.solverOption
    
    case '1'
        % Solver OPTION 1 - Sebastiano1
        if general.switchMAC
            !/usr/local/bin/docker start of_v2012;
            system([generalMAC.command.cd, directionsMAC.path.solver, ' && openfoam simpleFoam"']);
        else
            system(['cd ',directions.path.solver,' && simpleFoam']);
            system(['cd ',directions.path.solver,' && postProcess -func sample -latestTime']);
            system(['cd ',directions.path.solver,'/postProcessing/sample && ls >> list.txt']);
            addpath(genpath(append(directions.path.solver,'/postProcessing/sample')));
            sampleFolder = fopen('list.txt');
            listSampleFolder = textscan(sampleFolder,'%f');
            nameSampleFolder = listSampleFolder{1};
            if nameSampleFolder ~= 0 
                system(['mv ',directions.path.solver,'/postProcessing/sample/[1-9]* ',directions.path.solver,'/postProcessing/sample/0']);
            end
        end
        
    case '2'
        % Solver OPTION 2 - Run the solver simpleFoam in parallel
        system(['cd ',directions.path.solver,' && decomposePar']);
        system(['cd ',directions.path.solver,' && mpirun -np 4 simpleFoam -parallel']);
        system(['cd ',directions.path.solver,' && reconstructPar']);
        system(['cd ',directions.path.solver,' && rm -r processor*']);
        system(['cd ',directions.path.solver,' && postProcess -func sample -latestTime']);
        system(['cd ',directions.path.solver,'/postProcessing/sample && ls >> list.txt']);
        addpath(genpath(append(directions.path.solver,'/postProcessing/sample')));
        sampleFolder = fopen('list.txt');
        listSampleFolder = textscan(sampleFolder,'%f');
        nameSampleFolder = listSampleFolder{1};
        if nameSampleFolder ~= 0 
            system(['mv ',directions.path.solver,'/postProcessing/sample/[1-9]* ',directions.path.solver,'/postProcessing/sample/0']);
        end
        
end

%% run post

% import 
addpath(genpath(append(directions.path.solver,'/postProcessing')));
run importForceCoeffs.m

if  ~general.switchMAC
    run cp_distribution.m   % run to calculate and plot the cp distribution
    run plotResiduals.m    
end

% run importYPlus.m

TimeArr = table2array(forceCoeffs(101:end,1));
Cm = table2array(forceCoeffs(101:end,2));
Cd = table2array(forceCoeffs(101:end,3));
Cl = table2array(forceCoeffs(101:end,4));
Clf = table2array(forceCoeffs(101:end,5));
Clr = table2array(forceCoeffs(101:end,6));

h = figure(2);

hold on

plot (TimeArr,Cm);
plot (TimeArr,Cd);
plot (TimeArr,Cl);
plot (TimeArr,Clf);
plot (TimeArr,Clr);

legend('Cm','Cd','Cl','Clf','Clr');
title('Force Coefficients');
saveas(h,[directions.path.simulation,'/',directions.simName, '.fig']);
saveas(h,[directions.path.simulation,'/',directions.simName, '.eps'],'epsc');
hold off

if general.openParaView
    system(['cd ',directions.path.solver, ' && paraFoam']);
end

clear Cd Cl Clf Clr Cm TimeArr

if general.RunMultipleSimulations
    pause(1)
    close all; 
end

%% Save Output 

% Save workspace of current simulation in directory 
save([directions.path.simulation,'/',directions.simName, '.mat']);

% Save Simulation Inputs and Outputs to cell variable as archive
clear cfd_archive
load('Simulation_Outputs/cfd_archive.mat')
[r, ~] = size(cfd_archive);
cfd_archive{r+1,1} = directions.simName;
cfd_archive{r+1,2} = table2array(forceCoeffs(end,:));   % saves last line of coefficients
cfd_archive{r+1,3} = meshParam;
cfd_archive{r+1,4} = solvParam;
cfd_archive{r+1,5} = forceCoeffs;
save('Simulation_Outputs/cfd_archive.mat','cfd_archive')
clear ans r 