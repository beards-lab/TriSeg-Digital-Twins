%% Script Summary:
% This script is designed to optimize the mods value (modifiers) only for canonical subjects using a
% cost function from runSim.m. It integrates GA and fminsearch. Users can choose which one they want
% to use.

% Created by Andrew Meyer and Feng Gu
% Opt in log space
% Last modified: 11/14/2025

m = zeros(1,length(mods)); % First time run
idx = (1:length(mods));
% if GENDER ==1
%     m = readmatrix('modifiers_male.csv');
% else
%     m = readmatrix('modifiers_female.csv');
% end

% Set up options for patternsearch and fminsearch

[ub, lb] = m_bounds(mods); % use a simple function to set up boundary conditions for all modifiers
maxStallGen = 5;
maxGen = inf;
popSize = 150;


numVars = length(lb);
halfPop = round(popSize / 2);

% Generate 50% of the population from a uniform distribution
popUniform = lb + rand(halfPop, numVars) .* (ub - lb);

% Generate 50% of the population from a normal distribution
% Mean = midpoint of the search range; Std = (range)/6 so 99.7% lies within bounds
mu = 0;
sigma = (1.5 - 0.5) / 6;
popNormal = normrnd(mu, sigma, halfPop, numVars);

% Combine uniform and normal samples to form the initial population
InitialPopulationMatrix = [popUniform; popNormal];

% Clamp values to ensure all individuals remain within [lb, ub]
InitialPopulationMatrix = max(min(InitialPopulationMatrix, ub), lb);

% GA options
ga_options = optimoptions("ga", ...
    'Display', 'iter', ...                          % Show iteration details
    'MaxStallGenerations', maxStallGen, ...         % Early stopping
    'UseParallel', true, ...                        % Enable parallel evaluation
    'MaxGenerations', maxGen, ...                   % Max GA generations
    'PopulationSize', popSize, ...                  % GA population size
    'MutationFcn', {@mutationadaptfeasible, 0.2}, ... % Higher mutation rate
    'EliteCount', ceil(0.05 * popSize), ...         % Preserve top 5% individuals
    'InitialPopulationMatrix', InitialPopulationMatrix); % Custom initial population

options_PS_first = optimoptions('patternsearch',...
    'SearchFcn', {@searchga, idx,ga_options},...
    'Display', 'iter',...
    'MeshTolerance', 1e-6, ...
    'StepTolerance', 1e-6, ...
    'FunctionTolerance', 1e-6, ...
    'MaxIterations', 10,...
    'UseCompletePoll', true,...
    'UseParallel', true);

options_Fmin = optimset('Display','iter','TolFun', 1e-4, 'TolX', 1e-3, 'MaxIter',98); % reduce maxiter if you think it's getting stuck

[m, ~, ~, ~] = patternsearch(@(m)SrdevaluateModel(m,GENDER), m, [], [], [], [],lb, ub, [], options_PS_first);
m = fminsearch(@(m)SrdevaluateModel(m,GENDER), m, options_Fmin);
CurrentBestCost  = SrdevaluateModel(m,GENDER);

%%
options_PS_following = optimoptions('patternsearch',...
    ...'SearchFcn', {@searchga, idx,ga_options},...
    'Display', 'iter',...
    'MeshTolerance', 1e-6, ...
    'StepTolerance', 1e-6, ...
    'FunctionTolerance', 1e-6, ...
    'MaxIterations', 98,...
    'UseCompletePoll', true,...
    'UseParallel', true);


maxIter = 20;
iter = 0;
while iter < maxIter
    iter = iter + 1;
    [m, ~, ~, ~] = patternsearch(@(m)SrdevaluateModel(m,GENDER), m, [], [], [], [],lb, ub, [], options_PS_following);
    m = fminsearch(@(m)SrdevaluateModel(m,GENDER), m, options_Fmin);
    costnew = SrdevaluateModel(m,GENDER);
    if abs(CurrentBestCost - costnew) <= 10
        break;
    else
        CurrentBestCost = costnew;
    end
end


m = num2cell(m);
M = [mods;m];
if GENDER == 1
    writecell(M,'modifiers_male.csv');
else
    writecell(M,'modifiers_female.csv');
end