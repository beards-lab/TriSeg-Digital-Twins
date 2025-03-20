%% Script Summary:
% This script is specifically designed to optimize the mods value (modifiers)  
% for PKU patients using only RHC and 2D echo and a cost function from runSim.m.  
% It integrates GA, Pattern Search, and fminsearch, allowing users to select  
% which optimization method to use.  

% Created by Feng Gu  
% Last modified: 03/19/2024  
m = 1*ones(1,length(mods)); % if the predefined modifiers do not exist
[~, ~, mods] = targetVals_PKU(MRI_flag);
CurrentBestCost  = evaluateModelPKU(m,MRI_flag);
% Set up options for patternsearch and fminsearch
idx = find(ismember(mods, {'Amref_RV', 'Vw_RV','R_tPA','R_tSA','R_t_c','R_p_c','R_m_c','R_a_c'}));

%% Optimization
parpool(6)

% Set boundary conditions for all modifiers using a simple function
[ub, lb] = m_bounds(mods); 

% Set optimization parameters
maxStallGen = 1;  % Maximum stall generations
maxGen = 3;       % Maximum generations
popSize = 100 * length(idx);  % Population size

numVars = length(lb);  % Number of variables
halfPop = round(popSize / 2);  % Half of the population size

% 50% Uniform distribution
popUniform = lb + rand(halfPop, numVars) .* (ub - lb);

% 50% Normal distribution (mean at the midpoint of the search space, standard deviation set to 1/6 of the search range)
mu = 1;
sigma = (1.5 - 0.5) / 6;
popNormal = normrnd(mu, sigma, halfPop, numVars);

% Combine both populations
InitialPopulationMatrix = [popUniform; popNormal];

% Ensure all values are within the boundary range (lb and ub)
InitialPopulationMatrix = max(min(InitialPopulationMatrix, ub), lb);

% Set GA optimization options
ga_options = optimoptions("ga", ...
    'Display', 'iter', ...
    'MaxStallGenerations', maxStallGen, ...
    'UseParallel', true, ...
    'MaxGenerations', maxGen, ...
    'PopulationSize', popSize, ...
    'MutationFcn', {@mutationadaptfeasible, 0.2}, ...  % Increase mutation rate
    'EliteCount', ceil(0.05 * popSize), ...  % Retain the best solutions
    'InitialPopulationMatrix', InitialPopulationMatrix);  % Initialize population randomly

% Set Pattern Search options for the first optimization phase
options_PS_first = optimoptions('patternsearch', ...
    'SearchFcn', {@searchga, idx, ga_options}, ...
    'Display', 'iter', ...
    'MeshTolerance', 1e-6, ...
    'StepTolerance', 1e-6, ...
    'FunctionTolerance', 1e-6, ...
    'MaxIterations', 98, ...
    'UseCompletePoll', true, ...
    'UseParallel', true);

% Set fminsearch options
options_Fmin = optimset('Display', 'iter', 'TolFun', 1e-4, 'TolX', 1e-3, 'MaxIter', 98);  % Reduce max iterations if optimization is stuck

% Perform the first phase of optimization using patternsearch and fminsearch
[m, ~, ~, ~] = patternsearch(@(m)evaluateModelPKU(m, MRI_flag), m, [], [], [], [], lb, ub, [], options_PS_first);
m = fminsearch(@(m)evaluateModelPKU(m, MRI_flag), m, options_Fmin);

% Evaluate the cost of the optimized solution
costnew = evaluateModelPKU(m, MRI_flag);

% Set Pattern Search options for the subsequent optimization phases
options_PS_following = optimoptions('patternsearch', ...
    'Display', 'iter', ...
    'MeshTolerance', 1e-6, ...
    'StepTolerance', 1e-6, ...
    'FunctionTolerance', 1e-6, ...
    'MaxIterations', 98, ...
    'UseCompletePoll', true, ...
    'UseParallel', true);

% Maximum number of iterations for the following optimization loop
maxIter = 9;
iter = 0;

% Perform additional optimization until convergence or the maximum number of iterations is reached
while iter < maxIter
    iter = iter + 1;
    
    % Optimize using patternsearch and fminsearch
    [m, ~, ~, ~] = patternsearch(@(m)evaluateModelPKU(m, MRI_flag), m, [], [], [], [], lb, ub, [], options_PS_following); 
    m = fminsearch(@(m)evaluateModelPKU(m, MRI_flag), m, options_Fmin);
    
    % Evaluate the new cost
    costnew = evaluateModelPKU(m, MRI_flag);
    
    % Check for convergence
    if abs(CurrentBestCost - costnew) <= 10 
        break;
    else
        CurrentBestCost = costnew;
    end
end

% Save the final optimized results
output.mods = mods;
output.modifiers = m;
save PKUmFengRV.mat output

% Delete the parallel pool
delete(gcp("nocreate"))


