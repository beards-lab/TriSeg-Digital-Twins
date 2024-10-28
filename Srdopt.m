% m = ones(1,length(mods));
if GENDER ==1
    m = readmatrix('modifiers_male.csv');
else
    m = readmatrix('modifiers_female.csv');
end
cost = evaluateModel(m,GENDER);
%% GA
% [ub, lb] = m_bounds(mods);
% ub_0 = 3 .* ones(1,length(ub));
% lb_0 = 0.33 .* ones(1,length(lb));
% tic
% maxStallGen = 8;
% maxGen = Inf;
% % maxGen = 2;
% popSize = 100;
% % popSize = 10;
% %timelimit = getenv('SLURM_WALL_CLOCK_LIMIT')
% %jobtime = getenv('SLURM_JOB_TIME')
% options = optimoptions("ga",'Display','iter', 'MaxStallGenerations', maxStallGen, 'UseParallel',true,'MaxGenerations',maxGen,'PopulationSize',popSize,'InitialPopulationRange',[lb_0; ub_0]);
% [m,fcost,~,ga_out,fpop,fscores] = ga(@(m)evaluateModel(m,GENDER),length(m),[],[],[],[],lb, ub,[], options);
%% Optimize with fmincon
% opt    = optimset('Display','Iter'); 
% [m,J_fmincon] = fmincon(@(m)evaluateModel(m,GENDER),m,[],[],[],[],lb,ub,[],opt); 

%% Fminsearch
tic
options = optimset('Display','iter','PlotFcns',@optimplotfval, 'TolFun', 1e-4, 'TolX', 1e-3, 'MaxIter', 138); % reduce maxiter if you think it's getting stuck
m = fminsearch(@(m)SrdevaluateModel(m,GENDER), m, options);
m = num2cell(m);
M = [mods;m];
if GENDER == 1
    writecell(M,'modifiers_male.csv');
else
    writecell(M,'modifiers_female.csv');
end
toc