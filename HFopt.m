%% Script Summary:
% This script is designed to optimize the mods value (modifiers) using a cost function from runSim.m.
% It integrates GA, Patternsearch, and fminsearch. Users can choose which one they want to use.
% Note that GA will take a lot of time.

% Created by Andrew Meyer and Feng Gu
% Last modified: 10/29/2024

% load(sprintf('Sims/P_NO%dWindow%d.mat',PATIENT_NO,ModelWin)); % if it does exist
% m = output.modifiers;
% m = 1*ones(1,length(mods)); % if the predefined modifiers do not exist
load P_NO1.mat
m = output.modifiers;
% m = [modifiers(1:4) 0.25 modifiers(5:end)];
m(2) = 0.9;
m(5) = 0.25;
cost = evaluateModelUmich(m,patients,PATIENT_NO,ModelWin,MRI_flag); % call cost function in runSim.m

% %% GA
% [ub, lb] = m_bounds(mods); % use a simple function to set up boundary conditions for all modifiers
% ub_0 = 4 .* ones(1,length(ub));
% lb_0 = 0.25 .* ones(1,length(lb));
% maxStallGen = 8;
% maxGen = Inf;
% popSize = 100;
% options = optimoptions("ga",'Display','iter', 'MaxStallGenerations', maxStallGen, 'UseParallel',true,'MaxGenerations',maxGen,'PopulationSize',popSize,'InitialPopulationRange',[lb_0; ub_0]);
% [m,fcost,~,ga_out,fpop,fscores] = ga(@(m)evaluateModel(m,patients,PATIENT_NO,ModelWin), length(m),[],[],[],[],lb, ub,[], options);
% 
% %% patternsearch
% % Set up options for patternsearch
% options = optimoptions('patternsearch',...
%     'Display', 'iter',...
%     'PlotFcn', @psplotbestf,...
%     'FunctionTolerance', 1e-4,...
%     'StepTolerance', 1e-3,...
%     'MaxIterations', 648,...
%     'UseCompletePoll', true,...
%     'UseParallel', true);
% [m, fval, exitflag, output] = patternsearch(@(m)evaluateModel(m,patients,PATIENT_NO,ModelWin), m, [], [], [], [], [], [], [], options);

%% Fminsearch
while true
    options = optimset('Display','iter','PlotFcns',@optimplotfval, 'TolFun', 1e-4, 'TolX', 1e-3, 'MaxIter',100); % reduce maxiter if you think it's getting stuck
    m = fminsearch(@(m)evaluateModelUmich(m,patients,PATIENT_NO,ModelWin,MRI_flag), m, options);
    costnew = evaluateModelUmich(m,patients,PATIENT_NO,ModelWin,MRI_flag);
    if  cost - costnew <=1
        break;
    else
        cost = costnew;
    end
end

%% Save the output structure
output.mods = mods;
output.modifiers = m;
save P_NO1.mat output
% save(sprintf('Sims/P_NO%dWindow%d.mat',PATIENT_NO,ModelWin), "output");
