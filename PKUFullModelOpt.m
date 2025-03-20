%% Script Summary:
% This script is specifically designed to optimize the mods value (modifiers)  
% for PKU patients using 3D echo data and a cost function from runSim.m.  
% It integrates GA, Pattern Search, and fminsearch, allowing users to select  
% which optimization method to use.  

% Created by Feng Gu  
% Last modified: 03/19/2024  


m = 1*ones(1,length(mods)); % if the predefined modifiers do not exist
cost = evaluateModelPKU(m,MRI_flag); % call cost function in runSim.m

%% GA
[ub, lb] = m_bounds(mods); % use a simple function to set up boundary conditions for all modifiers
ub_0 = 4 .* ones(1,length(ub));
lb_0 = 0.25 .* ones(1,length(lb));
maxStallGen = 3;
maxGen = Inf;
popSize = 100;
options = optimoptions("ga",'Display','iter', 'MaxStallGenerations', maxStallGen, 'UseParallel',true,'MaxGenerations',maxGen,'PopulationSize',popSize,'InitialPopulationRange',[lb_0; ub_0]);
[m,fcost,~,ga_out,fpop,fscores] = ga(@(m)evaluateModelPKU(m,MRI_flag), length(m),[],[],[],[],lb, ub,[], options);

%% patternsearch
% Set up options for patternsearch
options = optimoptions('patternsearch',...
    'Display', 'iter',...
    'PlotFcn', @psplotbestf,...
    'FunctionTolerance', 1e-4,...
    'StepTolerance', 1e-3,...
    'MaxIterations', 648,...
    'UseCompletePoll', true,...
    'UseParallel', true);
[m, fval, exitflag, output] = patternsearch(@(m)evaluateModelPKU(m,MRI_flag), m, [], [], [], [], [], [], [], options);

%% Fminsearch
while true
    options = optimset('Display','iter','PlotFcns',@optimplotfval, 'TolFun', 1e-4, 'TolX', 1e-3, 'MaxIter',100); % reduce maxiter if you think it's getting stuck
    m = fminsearch(@(m)evaluateModelPKU(m,MRI_flag), m, options);
    costnew = evaluateModel(m,patients,PatID,ModelWin);
    if  cost - costnew <=1
        break;
    else
        cost = costnew;
    end
end

%% Save the output structure
output.mods = mods;
output.modifiers = m;
save PKUm.mat output
