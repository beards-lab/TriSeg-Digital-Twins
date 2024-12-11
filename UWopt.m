%% Script Summary:
% This script is designed to optimize the mods value (modifiers) using a cost function from runSim.m.
% It integrates GA, Patternsearch, and fminsearch. Users can choose which one they want to use.
% Note that GA will take a lot of time.

% Created by Feng Gu
% Last modified: 11/19/2024
% 
load(sprintf('UWSims/P_NO%d.mat',PatID)); % if it does exist
m = output.modifiers;
% m = 1*ones(1,length(mods)); % if the predefined modifiers do not exist
cost = evaluateModelUW(m,UWpatients,PatID); % call cost function in runSim.m

%% GA
% [ub, lb] = m_bounds(mods); % use a simple function to set up boundary conditions for all modifiers
% ub_0 = 4 .* ones(1,length(ub));
% lb_0 = 0.25 .* ones(1,length(lb));
% maxStallGen = 8;
% maxGen = Inf;
% popSize = 100;
% options = optimoptions("ga",'Display','iter', 'MaxStallGenerations', maxStallGen, 'UseParallel',true,'MaxGenerations',maxGen,'PopulationSize',popSize,'InitialPopulationRange',[lb_0; ub_0]);
% [m,fcost,~,ga_out,fpop,fscores] = ga(@(m)evaluateModelUW(m,UWpatients,PatID), length(m),[],[],[],[],lb, ub,[], options);

%% patternsearch
% Set up options for patternsearch
% options = optimoptions('patternsearch',...
%     'Display', 'iter',...
%     'PlotFcn', @psplotbestf,...
%     'FunctionTolerance', 1e-4,...
%     'StepTolerance', 1e-3,...
%     'MaxIterations', 18,...
%     'UseCompletePoll', true,...
%     'UseParallel', true);
% [m, fval, exitflag, output] = patternsearch(@(m)evaluateModelUW(m,UWpatients,PatID), m, [], [], [], [], [], [], [], options);

%% Fminsearch
while true
    options = optimset('Display','iter','PlotFcns',@optimplotfval, 'TolFun', 1e-4, 'TolX', 1e-3, 'MaxIter',100); % reduce maxiter if you think it's getting stuck
    m = fminsearch(@(m)evaluateModelUW(m,UWpatients,PatID), m, options);
    costnew = evaluateModelUW(m,UWpatients,PatID);
    if  cost - costnew <=1
        break;
    else
        cost = costnew;
    end
end

%% Save the output structure
output.mods = mods;
output.modifiers = m;
save(sprintf('UWSims/P_NO%d.mat',PatID), "output");
