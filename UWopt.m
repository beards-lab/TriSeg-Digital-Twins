% %% Script Summary:
% % This script is designed to optimize the mods value (modifiers) using a cost function from runSim.m.
% % It integrates GA, Patternsearch, and fminsearch. Users can choose which one they want to use.
% % Note that GA will take a lot of time.
% 
% % Created by Feng Gu
% % Last modified: 11/19/2024
% %
% % load(sprintf('Sims_GA/P_NO%d.mat',PATIENT_NO)); % if it does exist
% % if MRI_flag == 1
% %     load(sprintf('SimsUWwithCMR/P_NO%d.mat',PATIENT_NO)); % load results from great lake clusters
% % else
% %     load(sprintf('SimsUWwithoutCMR0/P_NO%d.mat',PATIENT_NO)); % load results from great lake clusters
% % end
% % m = output.m;
% % m = [m m(end)];
% % load(sprintf('SimsUWwithoutCMR0/P_NO%d.mat',PATIENT_NO))
% % m = output.modifiers;
% % m = 1*ones(1,length(mods)); % if the predefined modifiers do not exist
% % m(9) = 10;
% % CurrentBestCost = evaluateModelUW(m,UWpatients,PATIENT_NO); % call cost function in runSim.m
% 
% %% GA
% [ub, lb] = m_bounds(mods); % use a simple function to set up boundary conditions for all modifiers
% maxStallGen = 8;
% maxGen = Inf;
% % maxGen = 2;
% popSize = 666;
% % popSize = 10;
% % timelimit = getenv('SLURM_WALL_CLOCK_LIMIT')
% % jobtime = getenv('SLURM_JOB_TIME')
% 
% numVars = length(lb);
% halfPop = round(popSize / 2);
% 
% % 50% 均匀分布
% popUniform = lb + rand(halfPop, numVars) .* (ub - lb);
% 
% % 50% 正态分布（均值为搜索空间中点，标准差设为搜索区间的 1/6）
% mu = 1;
% sigma = (1.5 - 0.5) / 6;
% popNormal = normrnd(mu, sigma, halfPop, numVars);
% 
% % 合并两部分种群
% InitialPopulationMatrix = [popUniform; popNormal];
% 
% % 限制值域，确保所有值在 lb 和 ub 之间
% InitialPopulationMatrix = max(min(InitialPopulationMatrix, ub), lb);
% 
% 
% ga_options = optimoptions("ga", ...
%     'Display', 'iter', ...
%     'MaxStallGenerations', maxStallGen, ...
%     'UseParallel', true, ...
%     'MaxGenerations', maxGen, ...
%     'PopulationSize', popSize, ...
%     'MutationFcn', {@mutationadaptfeasible, 0.2}, ... % 提高变异率
%     'EliteCount', ceil(0.05 * popSize), ... % 保留最优解
%     'InitialPopulationMatrix', InitialPopulationMatrix); % 随机初始化种群
% % [m,fcost,~,ga_out,fpop,fscores] = ga(@(m)evaluateModelUW(m,UWpatients,PATIENT_NO), length(m),[],[],[],[],lb, ub,[], options);
% 
% %% patternsearch
% % Set up options for patternsearch and fminsearch
% idx = find(ismember(mods, {'Amref_RV', 'Vw_RV'}));
% options_PS = optimoptions('patternsearch',...
%     ...'SearchFcn', {@searchga, idx,ga_options},...
%     'Display', 'iter',...
%     'PlotFcn', @psplotbestf,...
%     'MeshTolerance', 1e-6, ...
%     'StepTolerance', 1e-6, ...
%     'FunctionTolerance', 1e-6, ...
%     'MaxIterations', 8,...
%     'UseCompletePoll', true,...
%     'UseParallel', true);
% options_Fmin = optimset('Display','iter','PlotFcns',@optimplotfval, 'TolFun', 1e-4, 'TolX', 1e-3, 'MaxIter',8); % reduce maxiter if you think it's getting stuck
% while true
% [m, fval, exitflag, output] = patternsearch(@(m)evaluateModelUW(m,UWpatients,PATIENT_NO), m, [], [], [], [],lb, ub, [], options_PS);
% offsets = [-0.3 -0.2, -0.1, 0, 0.1, 0.2,0.3];
% PotentialFmincosts = zeros(length(offsets), 1);
% optimized_m = zeros(2, length(offsets));
% parfor i = 1:length(offsets)
%     m_fmin_initial = m(idx);
%     m_fmin_initial(2) = m_fmin_initial(2)+offsets(i)
%     m_fmin = fminsearch(@(m_fmin)evaluateModelSubset(m_fmin,m,idx,UWpatients,PATIENT_NO), m_fmin_initial, options_Fmin);
%     optimized_m(:, i) = m_fmin;
%     PotentialFmincosts(i) = evaluateModelSubset(m_fmin,m,idx,UWpatients,PATIENT_NO);
% end
% [minCost, bestIdx] = min(costs);
% best_m = optimized_m(:, bestIdx);
% m(idx) = best_m;
% costnew = evaluateModelUW(m,UWpatients,PATIENT_NO);
%     if  CurrentBestCost - costnew <= 3
%         ValveFlag = true;
%         [targets, inputs, mods] = targetVals_UW(UWpatients,PATIENT_NO,MRI_flag);
%         [INIparams, INIinit] = estiminiParams(targets,inputs);
%         [params, init] = optParams(INIparams, INIinit, mods,m, targets);
%         runSim;
%         fields_to_find = {'MVr', 'PVr', 'AVr', 'TVr','DNP'};
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
% 
%             if isKey(field_positions,'DNP') && cost(field_positions('DNP'))>=30
%                 ValveFlag = false;
%                 [~,locs] = ismember('R_tPA',mods);
%                 if o_vals.DNP < targets.DNP
%                     m(locs) = m(locs)-0.1;
%                 else
%                     m(locs) = m(locs)*3.33;
%                 end
%             end
%         end
%         costnew = evaluateModelUW(m,UWpatients,PATIENT_NO);
%     end
%     if CurrentBestCost - costnew <= 3 && ValveFlag
%         break;
%     else
%         CurrentBestCost = costnew;
%     end
% end
% 
% %% Fminsearch
% % while true
% %    options = optimset('Display','iter','PlotFcns',@optimplotfval, 'TolFun', 1e-4, 'TolX', 1e-3, 'MaxIter',98); % reduce maxiter if you think it's getting stuck
% %   m = fminsearch(@(m)evaluateModelUW(m,UWpatients,PATIENT_NO), m, options);
% %     costnew = evaluateModelUW(m,UWpatients,PATIENT_NO);
% %     if  CurrentBestCost - costnew <= 10
% %         ValveFlag = true;
% %         % try to only use Fminsearch
% %         [targets, inputs, mods] = targetVals_UW(UWpatients,PATIENT_NO,MRI_flag);
% %         [INIparams, INIinit] = estiminiParams(targets,inputs);
% %         [params, init] = optParams(INIparams, INIinit, mods,m, targets);
% %         runSim;
% %         fields_to_find = {'MVr', 'PVr', 'AVr', 'TVr','DNP'};
% %         all_fields = fieldnames(targets);
% %         field_positions = containers.Map();
% %         for i = 1:length(fields_to_find)
% %             field_name = fields_to_find{i};
% %             if isfield(targets, field_name)
% %                 position = find(strcmp(all_fields, field_name));
% %                 field_positions(field_name) = position;
% %             end
% %         end
% % 
% %         if ~ field_positions.Count == 0
% %             if isKey(field_positions,'AVr') && cost(field_positions('AVr')) >=30
% %                 ValveFlag = false;
% %                 [~,locs] = ismember('R_a_c',mods);
% %                 if o_vals.AVr < targets.AVr
% %                     m(locs) = m(locs)-0.1;
% %                 else
% %                     m(locs) = m(locs)+0.1;
% %                 end
% %             end
% % 
% %             if isKey(field_positions,'MVr') && cost(field_positions('MVr')) >=30
% %                 ValveFlag = false;
% %                 [~,locs] = ismember('R_m_c',mods);
% %                 if o_vals.MVr < targets.MVr
% %                     m(locs) = m(locs)-0.1;
% %                 else
% %                     m(locs) = m(locs)+0.1;
% %                 end
% %             end
% % 
% %             if isKey(field_positions,'PVr') && cost(field_positions('PVr'))>=30
% %                 ValveFlag = false;
% %                 [~,locs] = ismember('R_p_c',mods);
% %                 if o_vals.PVr < targets.PVr
% %                     m(locs) = m(locs)-0.1;
% %                 else
% %                     m(locs) = m(locs)+0.1;
% %                 end
% %             end
% % 
% %             if isKey(field_positions,'TVr') && cost(field_positions('TVr'))>=30
% %                 ValveFlag = false;
% %                 [~,locs] = ismember('R_t_c',mods);
% %                 if o_vals.TVr < targets.TVr
% %                     m(locs) = m(locs)-0.1;
% %                 else
% %                     m(locs) = m(locs)+0.1;
% %                 end
% %             end
% % 
% %             if isKey(field_positions,'DNP') && cost(field_positions('DNP'))>=30
% %                 ValveFlag = false;
% %                 [~,locs] = ismember('R_tPA',mods);
% %                 if o_vals.DNP < targets.DNP
% %                     m(locs) = m(locs)-0.1;
% %                 else
% %                     m(locs) = m(locs)*3.33;
% %                 end
% %             end
% %         end
% %         costnew = evaluateModelUW(m,UWpatients,PATIENT_NO);
% %     end
% %     if CurrentBestCost - costnew <= 10 && ValveFlag
% %         break;
% %     else
% %         CurrentBestCost = costnew;
% %     end
% % end
% 
% 
% 
% %% Save the output structure
% output.mods = mods;
% output.modifiers = m;
% save(sprintf('SimsUWwithoutCMR0/P_NO%d.mat',PATIENT_NO),"output");
% 
load('UWcohort.mat')

UWpatients_const = parallel.pool.Constant(UWpatients);

tic
MRI_flag = 0;
[~, ~, mods] = targetVals_UW(UWpatients,PATIENT_NO,MRI_flag);
m = ones(1,length(mods));
Geo_Opt = 1;
CurrentBestCost  = evaluateModelUW(m,UWpatients,PATIENT_NO,Geo_Opt);

[ub, lb] = m_bounds(mods); % use a simple function to set up boundary conditions for all modifiers

% Set up options for patternsearch and fminsearch
idx = find(ismember(mods, {'Amref_RV', 'Vw_RV'}));
options_PS = optimoptions('patternsearch',...
    ...'SearchFcn', {@searchga, idx,ga_options},...
    'Display', 'iter',...
    'MeshTolerance', 1e-6, ...
    'StepTolerance', 1e-6, ...
    'FunctionTolerance', 1e-6, ...
    'MaxIterations', 5,...
    'UseCompletePoll', true,...
    'UseParallel', true);
options_Fmin = optimset('Display','iter','TolFun', 1e-4, 'TolX', 1e-3, 'MaxIter',98); % reduce maxiter if you think it's getting stuck
while true
[m, ~, ~, ~] = patternsearch(@(m)evaluateModelUW(m,UWpatients,PATIENT_NO,Geo_Opt), m, [], [], [], [],lb, ub, [], options_PS);
offsets = [-0.3 -0.2, -0.1, 0, 0.1, 0.2,0.3];
PotentialFmincosts = zeros(length(offsets), 1);
optimized_m = zeros(2, length(offsets));
parfor i = 1:length(offsets)
    UWpatients_local = UWpatients_const.Value;
    m_fmin_initial = m(idx);
    m_fmin_initial(2) = m_fmin_initial(2)+offsets(i)
    m_fmin = fminsearch(@(m_fmin)evaluateModelSubset(m_fmin,m,idx,UWpatients_local,PATIENT_NO,Geo_Opt), m_fmin_initial, options_Fmin);
    optimized_m(:, i) = m_fmin;
    PotentialFmincosts(i) = evaluateModelSubset(m_fmin,m,idx,UWpatients_local,PATIENT_NO,Geo_Opt);
end
[~, bestIdx] = min(PotentialFmincosts);
best_m = optimized_m(:, bestIdx);
m(idx) = best_m;
costnew = evaluateModelUW(m,UWpatients,PATIENT_NO);
    if  CurrentBestCost - costnew <= 3
        ValveFlag = true;
        [targets, inputs, mods] = targetVals_UW(UWpatients,PATIENT_NO,MRI_flag);
        [INIparams, INIinit] = estiminiParams(targets,inputs);
        [params, init] = optParams(INIparams, INIinit, mods,m, targets);
        runSim;
        fields_to_find = {'MVr', 'PVr', 'AVr', 'TVr','DNP'};
        all_fields = fieldnames(targets);
        field_positions = containers.Map();
        for i = 1:length(fields_to_find)
            field_name = fields_to_find{i};
            if isfield(targets, field_name)
                position = find(strcmp(all_fields, field_name));
                field_positions(field_name) = position;
            end
        end

        if ~ field_positions.Count == 0
            if isKey(field_positions,'AVr') && cost(field_positions('AVr')) >=30
                ValveFlag = false;
                [~,locs] = ismember('R_a_c',mods);
                if o_vals.AVr < targets.AVr
                    m(locs) = m(locs)-0.1;
                else
                    m(locs) = m(locs)+0.1;
                end
            end

            if isKey(field_positions,'MVr') && cost(field_positions('MVr')) >=30
                ValveFlag = false;
                [~,locs] = ismember('R_m_c',mods);
                if o_vals.MVr < targets.MVr
                    m(locs) = m(locs)-0.1;
                else
                    m(locs) = m(locs)+0.1;
                end
            end

            if isKey(field_positions,'PVr') && cost(field_positions('PVr'))>=30
                ValveFlag = false;
                [~,locs] = ismember('R_p_c',mods);
                if o_vals.PVr < targets.PVr
                    m(locs) = m(locs)-0.1;
                else
                    m(locs) = m(locs)+0.1;
                end
            end

            if isKey(field_positions,'TVr') && cost(field_positions('TVr'))>=30
                ValveFlag = false;
                [~,locs] = ismember('R_t_c',mods);
                if o_vals.TVr < targets.TVr
                    m(locs) = m(locs)-0.1;
                else
                    m(locs) = m(locs)+0.1;
                end
            end

            if isKey(field_positions,'DNP') && cost(field_positions('DNP'))>=30
                ValveFlag = false;
                [~,locs] = ismember('R_tPA',mods);
                if o_vals.DNP < targets.DNP
                    m(locs) = m(locs)-0.1;
                else
                    m(locs) = m(locs)*3.33;
                end
            end
        end
        costnew = evaluateModelUW(m,UWpatients,PATIENT_NO,Geo_Opt);
    end
    if CurrentBestCost - costnew <= 3 && ValveFlag
        break;
    else
        CurrentBestCost = costnew;
    end
end

Geo_Opt = 0;
CurrentBestCost  = evaluateModelUW(m,UWpatients,PATIENT_NO,Geo_Opt);
while true
[m, fval, exitflag, output] = patternsearch(@(m)evaluateModelUW(m,UWpatients,PATIENT_NO,Geo_Opt), m, [], [], [], [],lb, ub, [], options_PS);
offsets = [-0.3 -0.2, -0.1, 0, 0.1, 0.2,0.3];
PotentialFmincosts = zeros(length(offsets), 1);
optimized_m = zeros(length(mods), length(offsets));
parfor i = 1:length(offsets)
    UWpatients_local = UWpatients_const.Value;
    m_fmin_initial = m;
    m_fmin_initial(idx(2)) = m_fmin_initial(idx(2))+offsets(i)
    m_fmin = fminsearch(@(m_fmin)evaluateModelUW(m_fmin,UWpatients_local,PATIENT_NO,Geo_Opt), m_fmin_initial, options_Fmin);
    optimized_m(:, i) = m_fmin;
    PotentialFmincosts(i) = evaluateModelUW(m_fmin,UWpatients_local,PATIENT_NO,Geo_Opt);
end
[minCost, bestIdx] = min(PotentialFmincosts);
best_m = optimized_m(:, bestIdx);
m = best_m;
costnew = evaluateModelUW(m,UWpatients,PATIENT_NO);
    if  CurrentBestCost - costnew <= 3
        ValveFlag = true;
        [targets, inputs, mods] = targetVals_UW(UWpatients,PATIENT_NO,MRI_flag);
        [INIparams, INIinit] = estiminiParams(targets,inputs);
        [params, init] = optParams(INIparams, INIinit, mods,m, targets);
        runSim;
        fields_to_find = {'MVr', 'PVr', 'AVr', 'TVr','DNP'};
        all_fields = fieldnames(targets);
        field_positions = containers.Map();
        for i = 1:length(fields_to_find)
            field_name = fields_to_find{i};
            if isfield(targets, field_name)
                position = find(strcmp(all_fields, field_name));
                field_positions(field_name) = position;
            end
        end

        if ~ field_positions.Count == 0
            if isKey(field_positions,'AVr') && cost(field_positions('AVr')) >=30
                ValveFlag = false;
                [~,locs] = ismember('R_a_c',mods);
                if o_vals.AVr < targets.AVr
                    m(locs) = m(locs)-0.1;
                else
                    m(locs) = m(locs)+0.1;
                end
            end

            if isKey(field_positions,'MVr') && cost(field_positions('MVr')) >=30
                ValveFlag = false;
                [~,locs] = ismember('R_m_c',mods);
                if o_vals.MVr < targets.MVr
                    m(locs) = m(locs)-0.1;
                else
                    m(locs) = m(locs)+0.1;
                end
            end

            if isKey(field_positions,'PVr') && cost(field_positions('PVr'))>=30
                ValveFlag = false;
                [~,locs] = ismember('R_p_c',mods);
                if o_vals.PVr < targets.PVr
                    m(locs) = m(locs)-0.1;
                else
                    m(locs) = m(locs)+0.1;
                end
            end

            if isKey(field_positions,'TVr') && cost(field_positions('TVr'))>=30
                ValveFlag = false;
                [~,locs] = ismember('R_t_c',mods);
                if o_vals.TVr < targets.TVr
                    m(locs) = m(locs)-0.1;
                else
                    m(locs) = m(locs)+0.1;
                end
            end

            if isKey(field_positions,'DNP') && cost(field_positions('DNP'))>=30
                ValveFlag = false;
                [~,locs] = ismember('R_tPA',mods);
                if o_vals.DNP < targets.DNP
                    m(locs) = m(locs)-0.1;
                else
                    m(locs) = m(locs)*3.33;
                end
            end
        end
        costnew = evaluateModelUW(m,UWpatients,PATIENT_NO,Geo_Opt);
    end
    if CurrentBestCost - costnew <= 3 && ValveFlag
        break;
    else
        CurrentBestCost = costnew;
    end
end
% Save the output structure
[targets, inputs, mods] = targetVals_UW(UWpatients,PATIENT_NO,MRI_flag);
[INIparams, INIinit] = estiminiParams(targets,inputs);
[params, init] = optParams(INIparams, INIinit, mods,m,targets);
try
runSim
output.mods = mods;
output.m = m;
output.params = params;
output.init = init;
output.targets = targets;
output.inputs = inputs;
save(sprintf('Sims/P_NO%d',PATIENT_NO), "output");
catch
  error('The optimization has already reached convergence')
end
toc


%%%% Delete the parallel pool
delete(gcp("nocreate"));

