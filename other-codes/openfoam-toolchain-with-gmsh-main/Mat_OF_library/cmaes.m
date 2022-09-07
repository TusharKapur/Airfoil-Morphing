%
% Copyright (c) 2015, Mostapha Kalami Heris & Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "LICENSE" file for license terms.
%
% Project Code: YPEA108
% Project Title: Covariance Matrix Adaptation Evolution Strategy (CMA-ES)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Cite as:
% Mostapha Kalami Heris, CMA-ES in MATLAB (URL: https://yarpiz.com/235/ypea108-cma-es), Yarpiz, 2015.
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

% clc;
% clear;
% close all;

%% set directories

addpath (genpath(append(pwd,'/Mat_OF_library')));
addpath (genpath(directions.path.simulation)); 
addpath(genpath(append(directions.path.solver,'/postProcessing')));
addpath(genpath(append(directions.path.inputsExp)));

 %% Set initial optimization inputs

    % import Input from Excel file
    optmParam = importSearchNReplaceSheet(directions.path.inputsSim,directions.solv.InputFileOpt,directions.solv.InputSheetOpt,directions.solv.InputStartOpt,directions.solv.InputEndOpt);

    % determine number of design variables
    loc = cellfun('isempty', optmParam{:,'name'} );
    nDesVar = height(optmParam) - sum(loc == 1);    % number of variables to be replaced
    
    % define lower and upper boundaries of the design variables
    lb = optmParam.type;
    lb(cellfun('isempty',lb)) = [];
    ub = optmParam.description;
    ub(cellfun('isempty',ub)) = [];
    
    oldValue = optmParam.alias;
    oldValue(cellfun('isempty',oldValue)) = [];

 clear i 
 
%% Problem Settings

% CostFunction = cp_distribution(oldValue,solvParam,directions,optmParam);   % Cost Function

nVar = 5;                % Number of Unknown (Decision) Variables

VarSize = [1 nVar];       % Decision Variables Matrix Size

VarMin = str2double(lb);
VarMax = str2double(ub);

% VarMin = zeros(nDesVar,1);
% VarMax = zeros(nDesVar,1);
% VarMin(1,1) = 0.08;             % Lower Bound of Decision Variables for first design variable
% VarMin(2,1) = 0.8;             % Lower Bound of Decision Variables for second design variable
% VarMax(1,1) = 1;             % Upper Bound of Decision Variables for first design variable
% VarMax(2,1) = 0.9;             % Upper Bound of Decision Variables for second design variable

%% CMA-ES Settings

% Maximum Number of Iterations
MaxIt = 5;

% Population Size (and Number of Offsprings)
lambda = (4+round(3*log(nVar)))*10;

% Number of Parents
mu = round(lambda/2);

% Parent Weights
w = log(mu+0.5)-log(1:mu);
w = w/sum(w);

% Number of Effective Solutions
mu_eff = 1/sum(w.^2);

% Step Size Control Parameters (c_sigma and d_sigma);
sigma0 = zeros(nDesVar,1);
for desVar = 1:nDesVar
    sigma0(desVar,1) = 0.3*(VarMax(desVar,1) - VarMin(desVar,1));
end    
cs = (mu_eff+2)/(nVar+mu_eff+5);
ds = 1+cs+2*max(sqrt((mu_eff-1)/(nVar+1))-1, 0);
ENN = sqrt(nVar)*(1-1/(4*nVar)+1/(21*nVar^2));

% Covariance Update Parameters
cc = (4+mu_eff/nVar)/(4+nVar+2*mu_eff/nVar);
c1 = 2/((nVar+1.3)^2+mu_eff);
alpha_mu = 2;
cmu = min(1-c1, alpha_mu*(mu_eff-2+1/mu_eff)/((nVar+2)^2+alpha_mu*mu_eff/2));
hth = (1.4+2/(nVar+1))*ENN;

%% Initialization

ps = cell(MaxIt, 1);
pc = cell(MaxIt, 1);
C = cell(MaxIt, 1);
sigma = cell(MaxIt, 1);

ps{1} = zeros(VarSize);
pc{1} = zeros(VarSize);
C{1} = eye(nVar);
sigma{1} = sigma0;

empty_individual.Position = [];
empty_individual.Mean = [];
empty_individual.Step = [];
empty_individual.Cost = [];

M = repmat(empty_individual, MaxIt, 1);
M(1).Position = zeros(nDesVar,nVar);
for k = 1:nDesVar
    M(1).Position(k,:) = unifrnd(VarMin(k,1), VarMax(k,1), VarSize);
end
M(1).Mean = zeros(nDesVar,1);
for k = 1:nDesVar
    M(1).Mean(k,1) = mean(M(1).Position(k,:),2);
end
M(1).Step = zeros(VarSize);
[M(1).Cost,oldValue] = cp_RMSE(oldValue,M(1).Mean,solvParam,directions,optmParam,nDesVar);

BestSol = M(1);

BestCost = zeros(MaxIt, 1);

%% Move directories from simulation 

system(['mkdir ',directions.path.solver,'/optimization']);
system(['mv ',directions.path.solver,'/postProcessing/ ',directions.path.solver,'/[1-9]* ',directions.path.solver,'/optimization']);

%% CMA-ES Main Loop

for g = 1:MaxIt
    
    % Generate Samples
    pop = repmat(empty_individual, lambda, 1);
    for i = 1:lambda

        % Generating Sample
        pop(i).Step = mvnrnd(zeros(VarSize), C{g});
        for k = 1:nDesVar
            pop(i).Position(k,:) = M(g).Position(k,:) + sigma{g}(k,1)*pop(i).Step;
        end
        
        % Applying Bounds
        for k = 1:nDesVar 
            pop(i).Position(k,:) = max(pop(i).Position(k,:), VarMin(k,1));
            pop(i).Position(k,:) = min(pop(i).Position(k,:), VarMax(k,1));
        end
        
        % Obtain mean values
        for k = 1:nDesVar
            pop(i).Mean(k,1) = mean(pop(i).Position(k,:),2);
        end
        
        % Evaluation
        
        [pop(i).Cost,oldValue] = cp_RMSE(oldValue,pop(i).Mean,solvParam,directions,optmParam,nDesVar);
        
        % Update Best Solution Ever Found
        if pop(i).Cost<BestSol.Cost
            BestSol = pop(i);
            system(['rm -rf ',directions.path.solver,'/optimization/*']);
            system(['mv ',directions.path.solver,'/postProcessing/ ',directions.path.solver,'/[1-9]* ',directions.path.solver,'/optimization']);
        else
            system(['rm -rf ',directions.path.solver,'/postProcessing ',directions.path.solver,'/[1-9]*']);
        end
    end
    
    % Sort Population
    Costs = [pop.Cost];
    [Costs, SortOrder] = sort(Costs);
    pop = pop(SortOrder);
  
    % Save Results
    BestCost(g) = BestSol.Cost;
    
    % Display Results
    disp(['Iteration ' num2str(g) ': Best Cost = ' num2str(BestCost(g))]);
    
    % Exit At Last Iteration
    if g == MaxIt
        break;
    end
    
    % Update Mean
    M(g+1).Step = 0;
    for j = 1:mu
        M(g+1).Step = M(g+1).Step+w(j)*pop(j).Step;
    end
    for k = 1:nDesVar
        M(g+1).Position(k,:) = M(g).Position(k,:) + sigma{g}(k,1)*M(g+1).Step;
    end
        
    % Applying Bounds
    for k = 1:nDesVar
        M(g+1).Position(k,:) = max(M(g+1).Position(k,:), VarMin(k,1));
        M(g+1).Position(k,:) = min(M(g+1).Position(k,:), VarMax(k,1));
    end
    
    % Obtain mean values (CELIA)
    for k = 1:nDesVar
        M(g+1).Mean(k,1) = mean(M(g+1).Position(k,:),2);
    end   
        
    % Evaluation
    [M(g+1).Cost,oldValue] = cp_RMSE(oldValue,M(g+1).Mean,solvParam,directions,optmParam,nDesVar);
    
    % Update Best Solution Ever Found
    if M(g+1).Cost < BestSol.Cost
        BestSol = M(g+1);
        system(['rm -rf ',directions.path.solver,'/optimization/*']);
        system(['mv ',directions.path.solver,'/postProcessing/ ',directions.path.solver,'/[1-9]* ',directions.path.solver,'/optimization']);
    else
        system(['rm -rf ',directions.path.solver,'/postProcessing ',directions.path.solver,'/[1-9]*']);
    end
    
    % Update Step Size
    ps{g+1} = (1-cs)*ps{g}+sqrt(cs*(2-cs)*mu_eff)*M(g+1).Step/chol(C{g})';
    for k = 1:nDesVar
        sigma{g+1}(k,1) = sigma{g}(k,1)*exp(cs/ds*(norm(ps{g+1})/ENN-1))^0.3;
    end
    
    % Update Covariance Matrix
    if norm(ps{g+1})/sqrt(1-(1-cs)^(2*(g+1)))<hth
        hs = 1;
    else
        hs = 0;
    end
    delta = (1-hs)*cc*(2-cc);
    pc{g+1} = (1-cc)*pc{g}+hs*sqrt(cc*(2-cc)*mu_eff)*M(g+1).Step;
    C{g+1} = (1-c1-cmu)*C{g}+c1*(pc{g+1}'*pc{g+1}+delta*C{g});
    for j = 1:mu
        C{g+1} = C{g+1}+cmu*w(j)*pop(j).Step'*pop(j).Step;
    end
    
    % If Covariance Matrix is not Positive Defenite or Near Singular
    [V, E] = eig(C{g+1});
    if any(diag(E)<0)
        E = max(E, 0);
        C{g+1} = V*E/V;
    end
    
end

%% Display Results

figure;
saveas(figure,[directions.path.simulation,'/','cmaes.fig']);

% plot(BestCost, 'LineWidth', 2);
semilogy(BestCost, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;

