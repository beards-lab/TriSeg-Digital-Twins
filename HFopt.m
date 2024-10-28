load(sprintf('Sims/P_NO%dWindow%d.mat',PatID,ModelWin));
m = output.modifiers;
% m = 1*ones(1,length(mods)); % the first time running
cost = evaluateModel(m,patients,PatID,ModelWin);
%% GA
% [ub, lb] = m_bounds(mods);
% ub_0 = 4 .* ones(1,length(ub));
% lb_0 = 0.25 .* ones(1,length(lb));
% tic
% maxStallGen = 8;
% maxGen = Inf;
% % maxGen = 2;
% popSize = 100;
% % popSize = 10;
% % timelimit = getenv('SLURM_WALL_CLOCK_LIMIT')
% % jobtime = getenv('SLURM_JOB_TIME')
% options = optimoptions("ga",'Display','iter', 'MaxStallGenerations', maxStallGen, 'UseParallel',true,'MaxGenerations',maxGen,'PopulationSize',popSize,'InitialPopulationRange',[lb_0; ub_0]);
% [m,fcost,~,ga_out,fpop,fscores] = ga(@(m)evaluateModel(m,patients,PatID,ModelWin), length(m),[],[],[],[],lb, ub,[], options);


%% Fminsearch
tic
while true
    options = optimset('Display','iter','PlotFcns',@optimplotfval, 'TolFun', 1e-4, 'TolX', 1e-3, 'MaxIter',18); % reduce maxiter if you think it's getting stuck
    m = fminsearch(@(m)evaluateModel(m,patients,PatID,ModelWin), m, options);
    costnew = evaluateModel(m,patients,PatID,ModelWin);
    if  cost - costnew <=1
        break;
    else
        cost = costnew;
    end
end
output.mods = mods;
output.modifiers = m;
save(sprintf('Sims/P_NO%dWindow%d.mat',PatID,ModelWin),"output");
toc
%% patternsearch
tic
% Set up options for patternsearch
options = optimoptions('patternsearch',...
    'Display', 'iter',...
    'PlotFcn', @psplotbestf,...
    'FunctionTolerance', 1e-4,...
    'StepTolerance', 1e-3,...
    'MaxIterations', 648,...
    'UseCompletePoll', true,...
    'UseParallel', true);
Run patternsearch instead of fminsearch
[m, fval, exitflag, output] = patternsearch(@(m)evaluateModel(m,patients,PatID,ModelWin), m, [], [], [], [], [], [], [], options);

Save the output structure
output.mods = mods;
output.modifiers = m;
save(sprintf('Sims/P_NO%dWindow%d.mat',PatID,ModelWin), "output");
toc
