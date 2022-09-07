function [NRMSE,oldValue] = cp_distribution(oldValue,newValue,solvParam,directions,optmParam,nDesVar)
% clc
% clear all

% addpath(genpath(append(directions.path.solver,'/postProcessing/sample/0')));
% addpath(genpath(append(directions.path.inputsExp)));

%% run solver
newValue = arrayfun(@num2str,newValue,'UniformOutput',false);
for i = 1:nDesVar 
    searchnreplace(optmParam{i,2},directions.path.solver,'',oldValue(i),newValue(i));
    oldValue(i) = newValue(i);
end

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

clear i

addpath(genpath(append(directions.path.solver,'/postProcessing/sample/0')));

%% Import the data
% open the text files where the data is stored
OFmain_ID = fopen('p_main.raw'); % simulation data
OFflap_ID = fopen('p_flap.raw');
cpmain_ID = fopen('cp_main.txt'); % experimental data
cpflap_ID = fopen('cp_flap.txt');

% specify the format of the data 
OFmain_data = textscan(OFmain_ID,'%f%f%f%f','headerLines',2);
OFflap_data = textscan(OFflap_ID,'%f%f%f%f','headerLines',2);
cpmain_data = textscan(cpmain_ID,'%f%f');
cpflap_data = textscan(cpflap_ID,'%f%f');

clear OFmain_ID OFflap_ID cpmain_ID cpflap_ID

%% extract the data from the cell arrays

xcp_main = cell2mat(cpmain_data); 
xcp_flap = cell2mat(cpflap_data);
OF_xcp_main = cell2mat(OFmain_data([1 4]));
OF_xcp_flap = cell2mat(OFflap_data([1 4]));
Uinf = str2num(table2array(solvParam(16,4)));
% Uinf = 63;
xcp_flap(:,1) = xcp_flap(:,1) + 1;
OF_xcp_main(:,2) = OF_xcp_main(:,2)/(0.5*Uinf^2)*(-1);
OF_xcp_flap(:,2) = OF_xcp_flap(:,2)/(0.5*Uinf^2)*(-1);

%% Main element: divide the data into approximately lower and upper surface

% search the point with the highest "x-value" and get the corresponding cp
[OFxmax_main,i_xmax_main] = max(OF_xcp_main(:,1)); % OpenFOAM
cp_xmax_OFmain = OF_xcp_main(i_xmax_main,2);
cp_xmax_main = xcp_main(1,2); % Experimental data
xmax_main = xcp_main(1,1);

% first get points with higher cp than T.E.
i_cpup_OFmain = OF_xcp_main(:,2) >= cp_xmax_OFmain; %OpenFOAM 
OFmain_xcp_cpupper = OF_xcp_main(i_cpup_OFmain,:);
i_cpup_main = xcp_main(:,2) >= cp_xmax_main; % Experimental data
main_xcp_cpupper = xcp_main(i_cpup_main,:);

% from that points, separate just the ones from the TE to the point with
% max cp
% OpenFOAM
[OFmax_cpmain,OFindexcp] = max(OFmain_xcp_cpupper(:,2));   
OFx_cpmax = OFmain_xcp_cpupper(OFindexcp,1); 
OFindex1_main = OFmain_xcp_cpupper(:,1) >= OFx_cpmax & OFmain_xcp_cpupper(:,1) <= OFxmax_main;
OFmain_xcp_upper1 = OFmain_xcp_cpupper(OFindex1_main,:); % from TE to cpmax
OFmain_xcp_upper_1 = sortrows(OFmain_xcp_upper1,-1); % from TE to cpmax ordered
OFmain_xcp_cpupper2 = OFmain_xcp_cpupper(OFindex1_main == 0,:); % points from cpmax to LE
% Experimental data
[max_cpmain,indexcp] = max(main_xcp_cpupper(:,2));   
x_cpmax = main_xcp_cpupper(indexcp,1); 
index1_main = main_xcp_cpupper(:,1) >= x_cpmax & main_xcp_cpupper(:,1) <= xmax_main;
main_xcp_upper1 = main_xcp_cpupper(index1_main,:); % from TE to cpmax
main_xcp_upper_1 = sortrows(main_xcp_upper1,-1); % from TE to cpmax ordered
main_xcp_cpupper2 = main_xcp_cpupper(index1_main == 0,:); % points from cpmax to LE

% get points from xmin to cpmax
% OpenFOAM
[OFxmin_main,OFi_xmin_main] = min(OFmain_xcp_cpupper2(:,1));
OFcp_xmin_main = OFmain_xcp_cpupper2(OFi_xmin_main,2);
OFindex2_main = OFmain_xcp_cpupper2(:,2) >= OFcp_xmin_main & OFmain_xcp_cpupper2(:,2) <= OFmax_cpmain;
OFmain_xcp_upper2 = OFmain_xcp_cpupper2(OFindex2_main,:); % points from LE to cpmax
OFmain_xcp_upper_2 = sortrows(OFmain_xcp_upper2,-1); % from cpmax to LE ordered
OFmain_upper = [OFmain_xcp_upper_1;OFmain_xcp_upper_2];
% Experimental data
[xmin_main,i_xmin_main] = min(main_xcp_cpupper2(:,1));
cp_xmin_main = main_xcp_cpupper2(i_xmin_main,2);
index2_main = main_xcp_cpupper2(:,2) >= cp_xmin_main & main_xcp_cpupper2(:,2) <= max_cpmain;
main_xcp_upper2 = main_xcp_cpupper2(index2_main,:); % points from LE to cpmax
main_xcp_upper_2 = sortrows(main_xcp_upper2,-1); % from cpmax to LE
main_upper = [main_xcp_upper_1;main_xcp_upper_2];

% get points from lower surface
% OpenFOAM
OFmain_xcp_lower1 = OF_xcp_main(i_cpup_OFmain == 0,:); % from TE 
OFmain_xcp_lower_1 = sortrows(OFmain_xcp_lower1,-1);
OFmain_xcp_lower2 = OFmain_xcp_cpupper2(OFindex2_main == 0,:); % from LE
OFmain_xcp_lower_2 = sortrows(OFmain_xcp_lower2,-1);
OFmain_lower = [OFmain_xcp_lower_1;OFmain_xcp_lower_2];
% Experimental data
main_xcp_lower1 = xcp_main(i_cpup_main == 0,:); % from TE 
main_xcp_lower_1 = sortrows(main_xcp_lower1,-1);
main_xcp_lower2 = main_xcp_cpupper2(index2_main == 0,:); % from LE
main_xcp_lower_2 = sortrows(main_xcp_lower2,-1);
main_lower = [main_xcp_lower_1;main_xcp_lower_2];

%% Flap: divide the data into approximately lower and upper surface
% OpenFOAM
[OFxmax_flap,i_xmax_flap] = max(OF_xcp_flap(:,1)); % point at the TE
cp_xmax_OFflap = OF_xcp_flap(i_xmax_flap,2); % corresponding cp
OFxmin_flap = min(OF_xcp_flap(:,1)); % point at the LE

index_OFflap = OF_xcp_flap(:,2) >= cp_xmax_OFflap; 
OFflap_xcp_upper = OF_xcp_flap(index_OFflap,:); % points with higher cp than TE
OFflap_xcp_lower = OF_xcp_flap(index_OFflap == 0,:); % points with lower cp than TE

OFflap_upper = sortrows(OFflap_xcp_upper,-1); 
OFflap_lower = sortrows(OFflap_xcp_lower,-1);

% Experimental data
cp_xmax_flap = xcp_flap(1,2); % cp of the point at the TE
index_expFlap = xcp_flap(:,2) >= cp_xmax_flap;
flap_upper = xcp_flap(index_expFlap,:);
flap_lower = xcp_flap(index_expFlap ==0,:);

%% Move flap data
OFxmin_lowerF = min(OFflap_lower(:,1));
xmin_lowerF = min(flap_lower(:,1));

dif = xmin_lowerF - OFxmin_lowerF;
OFflap_upper(:,1) = OFflap_upper(:,1) + dif;
OFflap_lower(:,1) = OFflap_lower(:,1) + dif;

OFxmin_flap2 = min(OFflap_lower(:,1));
OFxmax_flap2 = max(OFflap_upper(:,1));

%% interpolation
x_main = linspace(OFxmin_main,OFxmax_main,3000);
cp_main_up = interp1(OFmain_upper(:,1),OFmain_upper(:,2),x_main,'pchip');
cp_main_low = interp1(OFmain_lower(:,1),OFmain_lower(:,2),x_main,'pchip');
OFmain_up = [x_main(:),cp_main_up(:)];
OFmain_low = [x_main(:),cp_main_low(:)];

x_flap = linspace(OFxmin_flap2,OFxmax_flap2,1500);
cp_flap_up = interp1(OFflap_upper(:,1),OFflap_upper(:,2),x_flap,'pchip');
cp_flap_low = interp1(OFflap_lower(:,1),OFflap_lower(:,2),x_flap,'pchip');
OFflap_up = [x_flap(:),cp_flap_up(:)];
OFflap_low = [x_flap(:),cp_flap_low(:)];

% figure(1)
% plot(x_main,cp_main_up,'r-',x_main,cp_main_low,'ro','DisplayName','Simulation data')
% hold on 
% plot(x_flap,cp_flap_up,'r-',x_flap,cp_flap_low,'ro','DisplayName','Simulation data')
% hold on
% plot(main_upper(:,1),main_upper(:,2),'ko',main_lower(:,1),main_lower(:,2),'ko','DisplayName','Experimental data')
% hold on 
% plot(flap_upper(:,1),flap_upper(:,2),'ko',flap_lower(:,1),flap_lower(:,2),'ko','DisplayName','Experimental data')
% hold on
% 
% title('cp distribution');
% ylabel('-cp');
% xlabel('x/chord');
% grid on;

%% determine the error between experimental and simulation results
% Find the closest simulation data to the measurement points
counter_upmain = 1;
numPoints_upmain = length(main_upper);
indexMin_upmain = zeros(length(main_upper),1); 
cpUpMain_distance = length(main_upper);
for counter_upmain = 1 : numPoints_upmain
    
    x_distances_upmain = abs(main_upper(counter_upmain,1) - OFmain_up(:,1));
    [min_xDist_upmain,indexMin_upmain(counter_upmain)] = min(x_distances_upmain);
    cpUpMain_distance(counter_upmain) = OFmain_up(indexMin_upmain(counter_upmain),2) - main_upper(counter_upmain,2);
%     plot(OFmain_up(indexMin_upmain(counter_upmain),1),OFmain_up(indexMin_upmain(counter_upmain),2),'*r','DisplayName','closest data')
    counter_upmain = counter_upmain + 1;

end

counter_upflap = 1;
numPoints_upflap = length(flap_upper);
indexMin_upflap = zeros(length(flap_upper),1); 
cpUpFlap_distance = length(flap_upper);
for counter_upflap = 1 : numPoints_upflap
    
    x_distances_upflap = abs(flap_upper(counter_upflap,1) - OFflap_up(:,1));
    [min_xDist_upflap,indexMin_upflap(counter_upflap)] = min(x_distances_upflap);
    cpUpFlap_distance(counter_upflap) = OFflap_up(indexMin_upflap(counter_upflap),2) - flap_upper(counter_upflap,2);
%     plot(OFflap_up(indexMin_upflap(counter_upflap),1),OFflap_up(indexMin_upflap(counter_upflap),2),'*r','DisplayName','closest data')
    counter_upflap = counter_upflap + 1;

end

counter_lowmain = 1;
numPoints_lowmain = length(main_lower);
indexMin_lowmain = zeros(length(main_lower),1);
cpLowMain_distance = length(main_lower);
for counter_lowmain = 1 : numPoints_lowmain
    
    x_distances_lowmain = abs(main_lower(counter_lowmain,1) - OFmain_low(:,1));
    [min_xDist_lowmain,indexMin_lowmain(counter_lowmain)] = min(x_distances_lowmain);
    cpLowMain_distance(counter_lowmain) = OFmain_low(indexMin_lowmain(counter_lowmain),2) - main_lower(counter_lowmain,2);
%     plot(OFmain_low(indexMin_lowmain(counter_lowmain),1),OFmain_low(indexMin_lowmain(counter_lowmain),2),'*r','DisplayName','closest data')
    counter_lowmain = counter_lowmain + 1;

end

counter_lowflap = 1;
numPoints_lowflap = length(flap_lower);
indexMin_lowflap = zeros(length(flap_lower),1);
cpLowFlap_distance = length(flap_lower);
for counter_lowflap = 1 : numPoints_lowflap
    
    x_distances_lowflap = abs(flap_lower(counter_lowflap,1) - OFflap_low(:,1));
    [min_xDist_lowflap,indexMin_lowflap(counter_lowflap)] = min(x_distances_lowflap);
    cpLowFlap_distance(counter_lowflap) = OFflap_low(indexMin_lowflap(counter_lowflap),2) - flap_lower(counter_lowflap,2);
%     plot(OFflap_low(indexMin_lowflap(counter_lowflap),1),OFflap_low(indexMin_lowflap(counter_lowflap),2),'*r','DisplayName','closest data')
    counter_lowflap = counter_lowflap + 1;

end

% legend(legendUnq());
% saveas(figure(1),[directions.path.simulation,'/','cp_distribution.fig']);

numPoints = numPoints_upmain + numPoints_lowmain + numPoints_upflap + numPoints_lowflap;
RMSE = sqrt(1/numPoints*(sum(cpUpMain_distance.^2) + sum(cpLowMain_distance.^2) + sum(cpUpFlap_distance.^2) + sum(cpLowFlap_distance.^2)));
NRMSE = RMSE/(max(xcp_main(:,2)) - min(xcp_main(:,2)));

clearvars -except NRMSE RMSE oldValue OFmain_upper main_upper OFmain_lower main_lower OFflap_upper OFflap_lower flap_upper flap_lower OFmain_up OFmain_low OFflap_up OFflap_low

%% Plot cp distribution

% figure (2)
% 
% plot(OFmain_upper(:,1),OFmain_upper(:,2),'ok','DisplayName','simulation data')
% hold on
% plot(OFmain_lower(:,1),OFmain_lower(:,2),'or','DisplayName','simulation data')
% hold on
% plot(OFflap_upper(:,1),OFflap_upper(:,2),'og','DisplayName','experimental data')
% hold on
% plot(OFflap_lower(:,1),OFflap_lower(:,2),'ob','DisplayName','experimental data')
% hold on
% plot(main_upper(:,1),main_upper(:,2),'ok','DisplayName','simulation data')
% hold on
% plot(main_lower(:,1),main_lower(:,2),'or','DisplayName','simulation data')
% hold on
% plot(flap_upper(:,1),flap_upper(:,2),'ok','DisplayName','experimental data')
% hold on
% plot(flap_lower(:,1),flap_lower(:,2),'or','DisplayName','experimental data')
% hold on
% 
% title('cp distribution');
% ylabel('-cp');
% xlabel('x/chord');
% grid on;
% legend(legendUnq());
end