
%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [7, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Time", "p", "Ux", "Uy", "k", "nut", "omega"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";


residuals = readtable(append(directions.path.solver,'/postProcessing/residuals/0/residuals.dat'), opts);

%% Clear temporary variables
clear opts

%% plot residuals

TimeArr = table2array(residuals(101:end,1));
p = table2array(residuals(101:end,2));
Ux = table2array(residuals(101:end,3));
Uy = table2array(residuals(101:end,4));
k = table2array(residuals(101:end,5));
nut = table2array(residuals(101:end,6));
omega = table2array(residuals(101:end,7));



figure

hold on

plot (TimeArr,p);
plot (TimeArr,Ux);
plot (TimeArr,Uy);
plot (TimeArr,k);
plot (TimeArr,nut);
plot (TimeArr,omega);

legend('p','Ux','Uy','k','nut','omega');
title('residuals');
saveas(figure(2),[directions.path.simulation,'/','residualsPlot.fig']);
hold off