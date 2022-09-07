function output = importSearchNReplaceSheet (importPath, importFileName, sheet, rangeStart, rangeEnd)
% Import data from spreadsheet inside directions.inputPath with 
% importExcelSheet (importFileName, sheet, rangeStart, rangeEnd, output). 
%  
% E.g.: importExcelsheet ('SHM_param_case1.xlsx','Mesher','A2','F14','meshParam')

    %% Set up the Import Options and import the data
    opts = spreadsheetImportOptions("NumVariables", 6);

    % Specify sheet and range
    opts.Sheet = sheet;
    opts.DataRange = append(rangeStart,':', rangeEnd);
    
    % Specify column names and types
    opts.VariableNames = ["name", "file", "alias", "value", "type", "description"];
    opts.VariableTypes = ["char", "string", "string", "string", "string", "string"];

    % Specify variable properties
    opts = setvaropts(opts, ["name", "alias", "value", "type", "description"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["name", "file", "alias", "value", "type", "description"], "EmptyFieldRule", "auto");

    % Import the data
    output = readtable(append(importPath,'/',importFileName), opts, "UseExcel", false);
    %% Clear temporary variables
    clear opts
end