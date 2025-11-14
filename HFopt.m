%% Script Summary:
% This script is designed to optimize the mods value (modifiers) using a cost function from runSim.m.
% It integrates GA, Patternsearch, and fminsearch. Users can choose which one they want to use.
% Note that GA will take a lot of time.

% Created by Andrew Meyer and Feng Gu
% Last modified: 11/14/2025

idx = 1:length(mods);
m = zeros(1, length(mods));

%% Set up boundary conditions and GA basic parameters
[ub, lb] = m_bounds(mods);                 % Parameter bounds for all modifiers
maxStallGen = 5;                           % Early stopping for GA
maxGen = inf;                              % Max GA generations
popSize = 150;                             % Total population size for GA

numVars = length(lb);
halfPop = round(popSize / 2);

%% Generate initial population (50% uniform + 50% normal)

% 50% uniform distribution within bounds
popUniform = lb + rand(halfPop, numVars) .* (ub - lb);

% 50% normal distribution (mean = midpoint of scaled space, std = range/6)
mu = 0;
sigma = (1.5 - 0.5) / 6;
popNormal = normrnd(mu, sigma, halfPop, numVars);

% Combine the two distributions
InitialPopulationMatrix = [popUniform; popNormal];

% Clamp to ensure all samples fall within [lb, ub]
InitialPopulationMatrix = max(min(InitialPopulationMatrix, ub), lb);

%% GA options
ga_options = optimoptions("ga", ...
    'Display', 'iter', ...
    'MaxStallGenerations', maxStallGen, ...
    'UseParallel', true, ...
    'MaxGenerations', maxGen, ...
    'PopulationSize', popSize, ...
    'MutationFcn', {@mutationadaptfeasible, 0.2}, ...   % Increased mutation rate
    'EliteCount', ceil(0.05 * popSize), ...             % Preserve top 5% individuals
    'InitialPopulationMatrix', InitialPopulationMatrix); % Custom hybrid population

%% Patternsearch (first stage): uses GA-based search function
options_PS_first = optimoptions('patternsearch', ...
    'SearchFcn', {@searchga, idx, ga_options}, ...      % Hybrid: GA inside patternsearch
    'Display', 'iter', ...
    'MeshTolerance', 1e-6, ...
    'StepTolerance', 1e-6, ...
    'FunctionTolerance', 1e-6, ...
    'MaxIterations', 3, ...                             % Limited iterations for first stage
    'UseCompletePoll', true, ...
    'UseParallel', true);

%% fminsearch options (local refinement)
options_Fmin = optimset( ...
    'Display', 'iter', ...
    'TolFun', 1e-4, ...
    'TolX', 1e-3, ...
    'MaxIter', 98);                                     % Reduce if convergence stalls

%% First round: patternsearch → fminsearch
[m, ~, ~, ~] = patternsearch(@(m)evaluateModelUmich(m,patients,PATIENT_NO,ModelWin,MRI_flag), ...
                             m, [], [], [], [], lb, ub, [], options_PS_first);

m = fminsearch(@(m)evaluateModelUmich(m,patients,PATIENT_NO,ModelWin,MRI_flag), ...
               m, options_Fmin);

CurrentBestCost = evaluateModelUmich(m,patients,PATIENT_NO,ModelWin,MRI_flag);

%% Iterative refinement loop
options_PS_following = optimoptions('patternsearch', ...
    'Display', 'iter', ...
    'MeshTolerance', 1e-6, ...
    'StepTolerance', 1e-6, ...
    'FunctionTolerance', 1e-6, ...
    'MaxIterations', 98, ...                            % Full patternsearch iterations
    'UseCompletePoll', true, ...
    'UseParallel', true);

maxIter = 3;
iter = 0;

while iter < maxIter
    iter = iter + 1;

    % Patternsearch refinement
    [m, ~, ~, ~] = patternsearch(@(m)evaluateModelUmich(m,patients,PATIENT_NO,ModelWin,MRI_flag), ...
                                 m, [], [], [], [], lb, ub, [], options_PS_following);

    % Local refinement
    m = fminsearch(@(m)evaluateModelUmich(m,patients,PATIENT_NO,ModelWin,MRI_flag), ...
                   m, options_Fmin);

    % Evaluate improvement
    costnew = evaluateModelUmich(m,patients,PATIENT_NO,ModelWin,MRI_flag);

    % Stop if improvement is small or solution diverges
    if abs(CurrentBestCost - costnew) <= 10 || isinf(costnew)
        break;
    else
        CurrentBestCost = costnew;
    end
end

%% Construct and save output structure
[~, targets, inputs, mods] = targetVals_HF(patients, PATIENT_NO, ModelWin, MRI_flag);
[INIparams, INIinit] = estiminiParams(targets, inputs);
[params, init] = optParams(INIparams, INIinit, mods, m);

try
    runSimOnGL;

    output.mods    = mods;
    output.m       = m;
    output.params  = params;
    output.init    = init;
    output.targets = targets;
    output.inputs  = inputs;
    output.o       = o;
    output.y       = y;
    output.t       = t;

    save(sprintf('Sims0707UM/P_NO%d', PATIENT_NO), "output");

catch
    % If runSimOnGL fails, save minimal outputs for debugging
    output.mods = mods;
    output.m    = m;

    save(sprintf('Sims0707UM/P_NO%d_failedGL', PATIENT_NO), "output");
    error('Simulation failed in runSimOnGL.');
end