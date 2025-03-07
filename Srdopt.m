%% Script Summary:
% This script is designed to optimize the mods value (modifiers) only for canonical subjects using a
% cost function from runSim.m. It integrates GA and fminsearch. Users can choose which one they want
% to use.

% Created by Andrew Meyer and Feng Gu
% Last modified: 10/29/2024

m = ones(1,length(mods)); % First time run
% if GENDER ==1
%     m = readmatrix('modifiers_male.csv');
% else
%     m = readmatrix('modifiers_female.csv');
% end

%% GA
[ub, lb] = m_bounds(mods);
ub_0 = 3 .* ones(1,length(ub));
lb_0 = 0.33 .* ones(1,length(lb));
maxStallGen = 8;
maxGen = Inf;
popSize = 200;
options = optimoptions("ga",'Display','iter', 'MaxStallGenerations', maxStallGen, 'UseParallel',true,'MaxGenerations',maxGen,'PopulationSize',popSize,'InitialPopulationRange',[lb_0; ub_0]);
[m,fcost,~,ga_out,fpop,fscores] = ga(@(m)SrdevaluateModel(m,GENDER),length(m),[],[],[],[],lb, ub,[], options);

%% Fminsearch
CurrentBestCost = SrdevaluateModel(m,GENDER);
while true    
    options = optimset('Display','iter','PlotFcns',@optimplotfval, 'TolFun', 1e-4, 'TolX', 1e-3, 'MaxIter', 98); % reduce maxiter if you think it's getting stuck
    m = fminsearch(@(m)SrdevaluateModel(m,GENDER), m, options);
    costnew = SrdevaluateModel(m,GENDER);
    if  CurrentBestCost - costnew <= 1e-3
        ValveFlag = true;
        % try to only use Fminsearch
        if GENDER ==1
            [targets, inputs, mods] = targetVals_male();
        else
            [targets, inputs, mods] = targetVals_female();
        end
        [INIparams, INIinit] = estiminiParams(targets,inputs);
        [params, init] = optParams(INIparams, INIinit, mods,m);
        runSim;
        fields_to_find = {'MVr', 'PVr', 'AVr', 'TVr'};
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
        end
        costnew = SrdevaluateModel(m,GENDER);
    end
    if CurrentBestCost - costnew <= 1e-3 && ValveFlag
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