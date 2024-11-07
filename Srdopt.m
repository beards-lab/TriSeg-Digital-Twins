%% Script Summary:
% This script is designed to optimize the mods value (modifiers) only for canonical subjects using a
% cost function from runSim.m. It integrates GA and fminsearch. Users can choose which one they want
% to use.

% Created by Andrew Meyer and Feng Gu
% Last modified: 10/29/2024

% m = ones(1,length(mods)); % First time run
if GENDER ==1
    m = readmatrix('modifiers_male.csv');
else
    m = readmatrix('modifiers_female.csv');
end
cost = SrdevaluateModel(m,GENDER);% Requiring different inputs compared to the real patients

%% GA
[ub, lb] = m_bounds(mods);
ub_0 = 3 .* ones(1,length(ub));
lb_0 = 0.33 .* ones(1,length(lb));
maxStallGen = 8;
maxGen = Inf;
popSize = 100;
options = optimoptions("ga",'Display','iter', 'MaxStallGenerations', maxStallGen, 'UseParallel',true,'MaxGenerations',maxGen,'PopulationSize',popSize,'InitialPopulationRange',[lb_0; ub_0]);
[m,fcost,~,ga_out,fpop,fscores] = ga(@(m)evaluateModel(m,GENDER),length(m),[],[],[],[],lb, ub,[], options);

%% Fminsearch
options = optimset('Display','iter','PlotFcns',@optimplotfval, 'TolFun', 1e-4, 'TolX', 1e-3, 'MaxIter', 138); % reduce maxiter if you think it's getting stuck
m = fminsearch(@(m)SrdevaluateModel(m,GENDER), m, options);
m = num2cell(m);
M = [mods;m];
if GENDER == 1
    writecell(M,'modifiers_male.csv');
else
    writecell(M,'modifiers_female.csv');
end