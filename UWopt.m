%% Script Summary:
% This script is designed to optimize the mods value (modifiers) using a cost function from runSim.m.
% It integrates GA, Patternsearch, and fminsearch. Users can choose which one they want to use.

% Created by Feng Gu
% Last modified: 03/20/2025

% The current strategy combines hybrid GA and pattern search,  
% followed by running fminsearch until convergence.  
% The goal is to make the parameters as unique as possible  
% to create a true digital twin.  

load('UWcohort.mat')

tic
MRI_flag = 0;
[~, ~, mods] = targetVals_UW(UWpatients,PATIENT_NO,MRI_flag);
m = ones(1,length(mods));
idx = find(ismember(mods, {'Amref_RV','Vw_RV','R_tPA','R_tSA','R_t_c','R_p_c','R_m_c','R_a_c'}));

%% Set up options for GA, Patternsearch, and fminsearch

[ub, lb] = m_bounds(mods); % use a simple function to set up boundary conditions for all modifiers
maxStallGen = 1;
maxGen = 3;
popSize = 100*length(idx);
numVars = length(lb);
halfPop = round(popSize / 2);

% 50% uniform distribution  
popUniform = lb + rand(halfPop, numVars) .* (ub - lb);  

% 50% normal distribution (mean at the center of the search space,  
% standard deviation set to 1/6 of the search range)  
mu = 1;  
sigma = (1.5 - 0.5) / 6;  
popNormal = normrnd(mu, sigma, halfPop, numVars);  

% Combine both parts of the population  
InitialPopulationMatrix = [popUniform; popNormal];  

% Enforce bounds to ensure all values stay within lb and ub  
InitialPopulationMatrix = max(min(InitialPopulationMatrix, ub), lb);  

ga_options = optimoptions("ga", ...
    'Display', 'iter', ...
    'MaxStallGenerations', maxStallGen, ...
    'UseParallel', true, ...
    'MaxGenerations', maxGen, ...
    'PopulationSize', popSize, ...
    'MutationFcn', {@mutationadaptfeasible, 0.2}, ... % Increase mutation rate 
    'EliteCount', ceil(0.05 * popSize), ... % Preserve the best solutions  
    'InitialPopulationMatrix', InitialPopulationMatrix); % Randomly initialize the population  


options_PS_first = optimoptions('patternsearch',...
    'SearchFcn', {@searchga, idx,ga_options},...
    'Display', 'iter',...
    'MeshTolerance', 1e-6, ...
    'StepTolerance', 1e-6, ...
    'FunctionTolerance', 1e-6, ...
    'MaxIterations', 5,...
    'UseCompletePoll', true,...
    'UseParallel', true);

options_PS_following = optimoptions('patternsearch',...
    ...'SearchFcn', {@searchga, idx,ga_options},...
    'Display', 'iter',...
    'MeshTolerance', 1e-6, ...
    'StepTolerance', 1e-6, ...
    'FunctionTolerance', 1e-6, ...
    'MaxIterations', 98,...
    'UseCompletePoll', true,...
    'UseParallel', true);

options_Fmin = optimset('Display','iter','TolFun', 1e-4, 'TolX', 1e-3, 'MaxIter',98); % reduce maxiter if you think it's getting stuck

%% Optimize parameters

[m, ~, ~, ~] = patternsearch(@(m)evaluateModelUW(m,UWpatients,PATIENT_NO), m, [], [], [], [],lb, ub, [], options_PS_first);
m = fminsearch(@(m)evaluateModelUW(m,UWpatients,PATIENT_NO), m, options_Fmin);
CurrentBestCost  = evaluateModelUW(m,UWpatients,PATIENT_NO);

maxIter = 9;
iter = 0;
while iter < maxIter
    iter = iter + 1;
    [mPS, ~, ~, ~] = patternsearch(@(m)evaluateModelUW(m,UWpatients,PATIENT_NO), m, [], [], [], [],lb, ub, [], options_PS_following);
    mFmin = fminsearch(@(m)evaluateModelUW(m,UWpatients,PATIENT_NO), mPS, options_Fmin);
    costnew = evaluateModelUW(mFmin,UWpatients,PATIENT_NO,Geo_Opt);
    if abs(CurrentBestCost - costnew) <= 10 || isinf(costnew)
        break;
    else
        CurrentBestCost = costnew;
        m = mFmin;
    end
end

% Save the output structure
[targets, inputs, mods] = targetVals_UW(UWpatients,PATIENT_NO,MRI_flag);
[INIparams, INIinit] = estiminiParams(targets,inputs);
[params, init] = optParams(INIparams, INIinit, mods,m);
try
    runSim
    output.mods = mods;
    output.m = m;
    output.params = params;
    output.init = init;
    output.targets = targets;
    output.inputs = inputs;
    save(sprintf('SimsUWwithoutCMR0/P_NO%d',PATIENT_NO), "output");
catch
    error('The optimization has already reached convergence')
end
toc


%%%% Delete the parallel pool
delete(gcp("nocreate"));

