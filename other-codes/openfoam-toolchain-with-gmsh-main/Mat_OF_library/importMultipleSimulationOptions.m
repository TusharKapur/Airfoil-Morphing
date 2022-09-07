% Import Simulation Options for running several Simulations
% Excel Table in OpenFoamToolchain/Simulation_Inputs/simulation_parameters directions.inputPath with 
 
%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 16);

% Specify sheet and range
opts.Sheet = "Inputs";
importEnd = int2str(general.NumberOfSimulations + 1);
opts.DataRange = append('A2:P', importEnd);

% Specify column names and types
opts.VariableNames = ["SimulationNames", "MesherType", "MesherSzenario", "MesherInputFile", "MesherSheet", "MesherStart", "MesherEnd", "SolverType", "SolverSzenario", "SolverInputFile", "SolverSheet", "SolverStart", "SolverEnd", "TwoElementAirfoil", "MainAirfoil", "FlapAirfoil", "FlapSeparated", "MeshingOption", "MeshParallel", "MeshCores", "MeshDimension", "SolverOption"];
opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "int8", "string", "string", "int8", "string", "int8", "int8", "string", "string"];

% Specify variable properties
opts = setvaropts(opts, ["SimulationNames", "SolverType", "SolverInputFile"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["SimulationNames", "MesherType", "MesherSzenario", "MesherInputFile", "MesherSheet", "MesherStart", "MesherEnd", "SolverType", "SolverSzenario", "SolverInputFile", "SolverSheet", "SolverStart", "SolverEnd", "MainAirfoil", "FlapAirfoil", "MeshDimension"], "EmptyFieldRule", "auto");

% Import the data
if ~general.switchMAC
    MultipleSimulations = readtable("/home/intern/OpenFoamToolchain/Simulation_Inputs/simulation_parameters/MultipleSimulations.xlsx", opts, "UseExcel", false);
else
    MultipleSimulations = readtable("/sgb/OpenFoamToolchain/Simulation_Inputs/simulation_parameters/MultipleSimulations.xlsx", opts, "UseExcel", false);
end

%% Clear temporary variables
clear opts importEnd
