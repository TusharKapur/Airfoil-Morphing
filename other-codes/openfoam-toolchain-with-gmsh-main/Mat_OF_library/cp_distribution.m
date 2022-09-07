% clc
% clear all

addpath(genpath(append(directions.path.solver,'/postProcessing/sample/0')));

%% Import the data
% open the text files where the data is stored
OFmain_ID = fopen('p_main.raw'); % simulation data
OFflap_ID = fopen('p_flap.raw');

% specify the format of the data 
OFmain_data = textscan(OFmain_ID,'%f%f%f%f','headerLines',2);
OFflap_data = textscan(OFflap_ID,'%f%f%f%f','headerLines',2);

%% extract the data from the cell arrays and calculate cp

OF_xcp_main = cell2mat(OFmain_data([1 4]));
OF_xcp_flap = cell2mat(OFflap_data([1 4]));
Uinf = str2num(table2array(solvParam(23,4)));
OF_xcp_main(:,2) = OF_xcp_main(:,2)/(0.5*Uinf^2)*(-1);
OF_xcp_flap(:,2) = OF_xcp_flap(:,2)/(0.5*Uinf^2)*(-1);

%% Main element: divide the data into approximately lower and upper surface

% search the point with the highest "x-value" and get the corresponding cp
[OFxmax_main,i_xmax_main] = max(OF_xcp_main(:,1)); 
cp_xmax_OFmain = OF_xcp_main(i_xmax_main,2);

% first get points with higher cp than T.E.
i_cpup_OFmain = OF_xcp_main(:,2) >= cp_xmax_OFmain;  
OFmain_xcp_cpupper = OF_xcp_main(i_cpup_OFmain,:);

% from that points, separate just the ones from the TE to the point with
% max cp
[OFmax_cpmain,OFindexcp] = max(OFmain_xcp_cpupper(:,2));   
OFx_cpmax = OFmain_xcp_cpupper(OFindexcp,1); 
OFindex1_main = OFmain_xcp_cpupper(:,1) >= OFx_cpmax & OFmain_xcp_cpupper(:,1) <= OFxmax_main;
OFmain_xcp_upper1 = OFmain_xcp_cpupper(OFindex1_main,:); % from TE to cpmax
OFmain_xcp_upper_1 = sortrows(OFmain_xcp_upper1,-1); % from TE to cpmax ordered
OFmain_xcp_cpupper2 = OFmain_xcp_cpupper(OFindex1_main == 0,:); % points from cpmax to LE

% get points from xmin to cpmax
[OFxmin_main,OFi_xmin_main] = min(OFmain_xcp_cpupper2(:,1));
OFcp_xmin_main = OFmain_xcp_cpupper2(OFi_xmin_main,2);
OFindex2_main = OFmain_xcp_cpupper2(:,2) >= OFcp_xmin_main & OFmain_xcp_cpupper2(:,2) <= OFmax_cpmain;
OFmain_xcp_upper2 = OFmain_xcp_cpupper2(OFindex2_main,:); % points from LE to cpmax
OFmain_xcp_upper_2 = sortrows(OFmain_xcp_upper2,-1); % from cpmax to LE ordered
OFmain_upper = [OFmain_xcp_upper_1;OFmain_xcp_upper_2];

% get points from lower surface
OFmain_xcp_lower1 = OF_xcp_main(i_cpup_OFmain == 0,:); % from TE 
OFmain_xcp_lower_1 = sortrows(OFmain_xcp_lower1,-1);
OFmain_xcp_lower2 = OFmain_xcp_cpupper2(OFindex2_main == 0,:); % from LE
OFmain_xcp_lower_2 = sortrows(OFmain_xcp_lower2,-1);
OFmain_lower = [OFmain_xcp_lower_1;OFmain_xcp_lower_2];

%% Flap: divide the data into approximately lower and upper surface

[OFxmax_flap,i_xmax_flap] = max(OF_xcp_flap(:,1)); % point at the TE
cp_xmax_OFflap = OF_xcp_flap(i_xmax_flap,2); % corresponding cp
OFxmin_flap = min(OF_xcp_flap(:,1)); % point at the LE

index_OFflap = OF_xcp_flap(:,2) >= cp_xmax_OFflap; 
OFflap_xcp_upper = OF_xcp_flap(index_OFflap,:); % points with higher cp than TE
OFflap_xcp_lower = OF_xcp_flap(index_OFflap == 0,:); % points with lower cp than TE

OFflap_upper = sortrows(OFflap_xcp_upper,-1); 
OFflap_lower = sortrows(OFflap_xcp_lower,-1);

%% Plot cp distribution

figure (1)

plot(OFmain_upper(:,1),OFmain_upper(:,2),'ok')
hold on
plot(OFmain_lower(:,1),OFmain_lower(:,2),'ok')
hold on
plot(OFflap_upper(:,1),OFflap_upper(:,2),'or')
hold on
plot(OFflap_lower(:,1),OFflap_lower(:,2),'or')
hold on

title('cp distribution');
ylabel('-cp');
xlabel('x/chord');
grid on;
% legend(legendUnq());