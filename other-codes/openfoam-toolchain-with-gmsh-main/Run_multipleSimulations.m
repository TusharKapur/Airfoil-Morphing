%% Script for running multiple Simulations in a row
% Input through
% Simulation_Inputs/simulation_parameters/MultipleSimulations.xlsx
% Script reads out everything from there
clear;
clc;
%% Initialize
% NUMBER OF SIMULATIONS (From top of Excel list)
    general.NumberOfSimulations = 2;            % Input needed
    general.switchMAC = 1;                      % Input needed (DockerDeamon must run)
% Import Multiple Simulation Options
	addpath(genpath(pwd));
    importMultipleSimulationOptions; 
    
for k=1:general.NumberOfSimulations
%% Specify Inputs

% SIMULATION NAME %
    directions.simName = char(MultipleSimulations.SimulationNames(k));
    
% MESHER %
    % currently only 'snappyHexMesh' available
    directions.mesh.Type       = char(MultipleSimulations.MesherType(k));
    directions.mesh.Szenario   = char(MultipleSimulations.MesherSzenario(k));
    directions.mesh.InputFile  = char(MultipleSimulations.MesherInputFile(k));
    % Sheet, Start and End of Mesher ExcelImport
    directions.mesh.InputSheet = char(MultipleSimulations.MesherSheet(k));
    directions.mesh.InputStart = char(MultipleSimulations.MesherStart(k));
    directions.mesh.InputEnd   = char(MultipleSimulations.MesherEnd(k));
% SOLVER %
    % currently only 'simpleFoam' available
    directions.solv.Type       = char(MultipleSimulations.SolverType(k));
    directions.solv.Szenario   = char(MultipleSimulations.SolverSzenario(k));
    directions.solv.InputFile  = char(MultipleSimulations.SolverInputFile(k));
    % Sheet, Start and End of Solver ExcelImport
    directions.solv.InputSheet = char(MultipleSimulations.SolverSheet(k)); 
    directions.solv.InputStart = char(MultipleSimulations.SolverStart(k));
    directions.solv.InputEnd   = char(MultipleSimulations.SolverEnd(k));
% AIRFOIL FILES - dont add data ending / no .txt etc
    general.twoElemAirfoil     = char(MultipleSimulations.TwoElementAirfoil(k));   % if not just use mainAirfoil
    MainAirfoilName            = char(MultipleSimulations.MainAirfoil(k));
    FlapAirfoilName            = char(MultipleSimulations.FlapAirfoil(k));
    % Do you want to have the two airfoil files separated? (1 = yes, 0 = no)
    general.flapSeparated      = char(MultipleSimulations.FlapSeparated(k));  
% MESHING OPTIONS
    % Option1-4: snappyHexMesh
    % Option1: blockMesh + surfaceFetExtr + snappyHex (overwr) + checkMesh
    % Option2: blockMesh + snappyHex
    % Option3: Remesh
    % Option4: Mesh in Parallel
    general.meshingOption = char(MultipleSimulations.MeshingOption(k));
    general.meshParallel  = char(MultipleSimulations.MeshParallel(k));
    general.meshCores     = char(MultipleSimulations.MeshCores(k));
    general.meshDimension = char(MultipleSimulations.MeshDimension(k));    % [2D, 3D] (3D not implemented)
% SOLVER OPTIONS
    % Option1: Solve
    general.solverOption  = char(MultipleSimulations.SolverOption(k)); 
    % Do YOU WANT TO OPEN THE PARAVIEW GUI?
    general.openParaView = 0; 
   
    
    
   run Run_Toolchain.m
    
end
 