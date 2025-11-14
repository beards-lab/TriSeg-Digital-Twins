%% Script Summary:
% This script is designed to optimize the mods value (modifiers) only for canonical subjects using a
% cost function from runSim.m. It integrates GA and fminsearch. Users can choose which one they want
% to use.

% Created by Andrew Meyer and Feng Gu
% Last modified: 10/29/2024

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
% maxGen = 2;
% popSize = 100*length(idx);
popSize = 150;
% timelimit = getenv('SLURM_WALL_CLOCK_LIMIT')
% jobtime = getenv('SLURM_JOB_TIME')

numVars = length(lb);
halfPop = round(popSize / 2);

% 50% 均匀分布
popUniform = lb + rand(halfPop, numVars) .* (ub - lb);

% 50% 正态分布（均值为搜索空间中点，标准差设为搜索区间的 1/6）
mu = 0;
sigma = (1.5 - 0.5) / 6;
popNormal = normrnd(mu, sigma, halfPop, numVars);

% 合并两部分种群
InitialPopulationMatrix = [popUniform; popNormal];

% 限制值域，确保所有值在 lb 和 ub 之间
InitialPopulationMatrix = max(min(InitialPopulationMatrix, ub), lb);


ga_options = optimoptions("ga", ...
    'Display', 'iter', ...
    'MaxStallGenerations', maxStallGen, ...
    'UseParallel', true, ...
    'MaxGenerations', maxGen, ...
    'PopulationSize', popSize, ...
    'MutationFcn', {@mutationadaptfeasible, 0.2}, ... % 提高变异率
    'EliteCount', ceil(0.05 * popSize), ... % 保留最优解
    'InitialPopulationMatrix', InitialPopulationMatrix); % 随机初始化种群


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

% 
% %% GA
% [ub, lb] = m_bounds(mods);
% ub_0 = 3 .* ones(1,length(ub));
% lb_0 = 0.33 .* ones(1,length(lb));
% maxStallGen = 8;
% maxGen = Inf;
% popSize = 200;
% options = optimoptions("ga",'Display','iter', 'MaxStallGenerations', maxStallGen, 'UseParallel',true,'MaxGenerations',maxGen,'PopulationSize',popSize,'InitialPopulationRange',[lb_0; ub_0]);
% [m,fcost,~,ga_out,fpop,fscores] = ga(@(m)SrdevaluateModel(m,GENDER),length(m),[],[],[],[],lb, ub,[], options);
% 
% %% Fminsearch
% CurrentBestCost = SrdevaluateModel(m,GENDER);
% while true    
%     options = optimset('Display','iter','PlotFcns',@optimplotfval, 'TolFun', 1e-4, 'TolX', 1e-3, 'MaxIter', 98); % reduce maxiter if you think it's getting stuck
%     m = fminsearch(@(m)SrdevaluateModel(m,GENDER), m, options);
%     costnew = SrdevaluateModel(m,GENDER);
%     if  CurrentBestCost - costnew <= 1e-3
%         ValveFlag = true;
%         % try to only use Fminsearch
%         if GENDER ==1
%             [targets, inputs, mods] = targetVals_male();
%         else
%             [targets, inputs, mods] = targetVals_female();
%         end
%         [INIparams, INIinit] = estiminiParams(targets,inputs);
%         [params, init] = optParams(INIparams, INIinit, mods,m);
%         runSim;
%         fields_to_find = {'MVr', 'PVr', 'AVr', 'TVr'};
%         all_fields = fieldnames(targets);
%         field_positions = containers.Map();
%         for i = 1:length(fields_to_find)
%             field_name = fields_to_find{i};
%             if isfield(targets, field_name)
%                 position = find(strcmp(all_fields, field_name));
%                 field_positions(field_name) = position;
%             end
%         end
% 
%         if ~ field_positions.Count == 0
%             if isKey(field_positions,'AVr') && cost(field_positions('AVr')) >=30
%                 ValveFlag = false;
%                 [~,locs] = ismember('R_a_c',mods);
%                 if o_vals.AVr < targets.AVr
%                     m(locs) = m(locs)-0.1;
%                 else
%                     m(locs) = m(locs)+0.1;
%                 end
%             end
% 
%             if isKey(field_positions,'MVr') && cost(field_positions('MVr')) >=30
%                 ValveFlag = false;
%                 [~,locs] = ismember('R_m_c',mods);
%                 if o_vals.MVr < targets.MVr
%                     m(locs) = m(locs)-0.1;
%                 else
%                     m(locs) = m(locs)+0.1;
%                 end
%             end
% 
%             if isKey(field_positions,'PVr') && cost(field_positions('PVr'))>=30
%                 ValveFlag = false;
%                 [~,locs] = ismember('R_p_c',mods);
%                 if o_vals.PVr < targets.PVr
%                     m(locs) = m(locs)-0.1;
%                 else
%                     m(locs) = m(locs)+0.1;
%                 end
%             end
% 
%             if isKey(field_positions,'TVr') && cost(field_positions('TVr'))>=30
%                 ValveFlag = false;
%                 [~,locs] = ismember('R_t_c',mods);
%                 if o_vals.TVr < targets.TVr
%                     m(locs) = m(locs)-0.1;
%                 else
%                     m(locs) = m(locs)+0.1;
%                 end
%             end
%         end
%         costnew = SrdevaluateModel(m,GENDER);
%     end
%     if CurrentBestCost - costnew <= 1e-3 && ValveFlag
%         break;
%     else
%         CurrentBestCost = costnew;
%     end
% end
m = num2cell(m);
M = [mods;m];
if GENDER == 1
    writecell(M,'modifiers_male.csv');
else
    writecell(M,'modifiers_female.csv');
end