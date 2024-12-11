%% Script Summary:
% This script calls the functions: targetVal.m and estimParams.m, as well as the scripts: XXopt.m
% and runSim.m.
% The first section generates simulations for canonical male and female subjects, and the
% second section generates simulations for heart failure patients.

% Created by Feng Gu
% 11/18 Starting to work on dismiss right side infomation.
% Last modified: 10/29/2024

%% Simulation for Canonical Subjects
clear
for GENDER  = 2 % 1 for male and other for female
    if GENDER == 1
        [targets, inputs, mods] = targetVals_male();
        modifiers = readmatrix('modifiers_male.csv','range','2:2');
        % modifiers = ones(1, length(mods));
        [~,params, init] = estimParams(targets,inputs,mods,modifiers);
    else
        [targets, inputs, mods] = targetVals_female();
        modifiers = readmatrix('modifiers_female.csv','range','2:2');
        % modifiers = ones(1, length(mods));
        [~,params, init] = estimParams(targets,inputs,mods,modifiers);
    end
    Sim = 1; % 1 for simulation and other for optimazation
    if Sim ==1
        runSim;
        NplotSrd; % 4-panel figure just for canonical subjects
        GetMovie; % TriSeg model: displacement and stress as functions of time
        See_TriSeg; % slices of GetMoive.m
    else
        Srdopt; % optimize modifiers of canonical subjects
    end
end

%% Simulation for heart failure patients
clear
load AllPatients.mat
secondslot = [17 79 206 256 288 325 352 355 360 361]; % patients with 2 3-month windows
% Shuntlist = [34 41 54 61 83 116 183 231 268 278 312]; % patients with shunt
for PatID = 192 % any number between 1 and 370, example patient in paper is 192
    for ModelWin =  1
        [Windowdate,targets, inputs, mods] = targetVals_HF(patients,PatID,ModelWin);
        RUNOPT = 0; % 0 for simulation and 1 for optimazation
        if RUNOPT ==1
            HFopt; % optimize modifiers of HF patients
        end
        if RUNOPT == 0
            load(sprintf('Sims/P_NO%dWindow%d.mat',PatID,ModelWin)); % load results from great lake clusters
            modifiers = output.modifiers;
            [Error, params, init] = estimParams(targets,inputs,mods, modifiers);
            runSim;
            Print_cost = 1;% 1 for performance, other for patient's pre-condition (PHI sensitive)
            NplotFit; % 6-panel figure for HF patients
            GetMovie;
            See_TriSeg;
        end
    end
end

%% RV info missing try
clear
load('UWcohort.mat');
% TRW = NaN(64,1);
PredictedLVEDV = NaN(64,1);
PredictedLVESV = NaN(64,1);
for PatID = 1:64
    [targets, inputs, mods] = targetVals_UW(UWpatients,PatID);
    RUNOPT = 0; % 0 for simulation and 1 for optimazation
    try
    if RUNOPT ==1
        UWopt; % optimize modifiers of HF patients
    end
    if RUNOPT == 0
        % load(sprintf('UWSims/P_NO%d.mat',PatID)); % load results from great lake clusters
        load(sprintf('Sims_GA/P_NO%d.mat',PatID)); % load results from great lake clusters
        modifiers = output.modifiers;
        % % modifiers(9) = 2;
        % modifiers(16) = 1.1;
        % modifiers(11) = 2;
        % modifiers = ones(1,length(mods));
        [Error, params, init] = estimParams(targets,inputs,mods, modifiers);
        runSim;
        % Print_cost = 1;% 1 for performance, other for patient's pre-condition (PHI sensitive)
        % NplotFit; % 6-panel figure for HF patients
        % GetMovie;
        % See_TriSeg;
        pause(3);
        % TRW(PatID) = o_vals.Hed_RW;
        PredictedLVEDV(PatID) = o_vals.LVEDV;
        PredictedLVESV(PatID) = o_vals.LVESV;
    end
    catch
    end
end
%% extract modifiers from cluster
% clear
for PatID = 1:64
    if PatID < 10
        matFilesDir = sprintf('Runs/p00%d/12-03/',PatID);
    elseif PatID < 100
        matFilesDir = sprintf('Runs/p0%d/12-03/',PatID);
    else
        matFilesDir = sprintf('Runs/p%d/12-03/',PatID);
    end
    matFiles = dir(fullfile(matFilesDir, '*.mat'));
    if ~isempty(matFiles)
        matFileName = matFiles.name;
        matFilePath = fullfile(matFilesDir, matFileName);
        GAresult = load(matFilePath);
        output.modifiers = GAresult.output.m;
        output.mods = GAresult.output.mods;
        save(sprintf('Sims_GA/P_NO%d.mat',PatID),"output");
    end
end
%%
% 去除 NaN 对于 RVEDV vs PredictedRVEDV
validIdx_RVEDV =  ~isnan(PredictedRVEDV);
[r_RVEDV, p_RVEDV] = corr(UWpatients.RVEDV(validIdx_RVEDV), PredictedRVEDV(validIdx_RVEDV), 'Type', 'Pearson');

% 去除 NaN 对于 RVESV vs PredictedRVESV
validIdx_RVESV = ~isnan(PredictedRVESV);
[r_RVESV, p_RVESV] = corr(UWpatients.RVESV(validIdx_RVESV), PredictedRVESV(validIdx_RVESV), 'Type', 'Pearson');

% 去除 NaN 对于 Hed_RW vs TRW
validIdx_RW = ~isnan(TRW);
[r_RW, p_RW] = corr(UWpatients.Hed_RW(validIdx_RW), TRW(validIdx_RW), 'Type', 'Pearson');

% 显示结果
disp(['RVEDV vs PredictedRVEDV: r = ', num2str(r_RVEDV), ', p = ', num2str(p_RVEDV)]);
disp(['RVESV vs PredictedRVESV: r = ', num2str(r_RVESV), ', p = ', num2str(p_RVESV)]);
disp(['Hed_RW vs TRW: r = ', num2str(r_RW), ', p = ', num2str(p_RW)]);

% 可选：在散点图中显示 r 和 p 值
figure(66); scatter(UWpatients.RVEDV(validIdx_RVEDV), PredictedRVEDV(validIdx_RVEDV));
title(['RVEDV: r = ', num2str(r_RVEDV), ', p = ', num2str(p_RVEDV)]);
xlabel('RVEDV'); ylabel('PredictedRVEDV');

figure(77); scatter(UWpatients.RVESV(validIdx_RVESV), PredictedRVESV(validIdx_RVESV));
title(['RVESV: r = ', num2str(r_RVESV), ', p = ', num2str(p_RVESV)]);
xlabel('RVESV'); ylabel('PredictedRVESV');

figure(88); scatter(UWpatients.Hed_RW(validIdx_RW), TRW(validIdx_RW));
title(['Hed_RW: r = ', num2str(r_RW), ', p = ', num2str(p_RW)]);
xlabel('Hed_RW'); ylabel('TRW');
 t1 = [PredictedRVEDV UWpatients.RVEDV];
DI1 = PredictedRVEDV-UWpatients.RVEDV;
figure;
histogram(DI1, 'Normalization', 'probability'); % 归一化为频率
xlabel('DI1 Values');
ylabel('Frequency');
title('Frequency Distribution of DI1');

%%
% 去除 NaN 对于 LVEDV vs PredictedLVEDV
CO = 1000 * UWpatients.CO./ 60;
SV = 60 * CO ./ UWpatients.HR;
PredictedLVEDV =   SV./ (UWpatients.EF.*0.01);
PredictedLVESV =   PredictedLVEDV-SV;
% PredictedLVEDV = UWpatients.LVIDd.^3*0.7851+97.32;
% PredictedLVESV = UWpatients.LVIDs.^3*0.9185+63.92;
validIdx_LVEDV = ~isnan(PredictedLVEDV);
[r_LVEDV, p_LVEDV] = corr(UWpatients.LVEDV(validIdx_LVEDV), PredictedLVEDV(validIdx_LVEDV), 'Type', 'Pearson');

% 去除 NaN 对于 LVESV vs PredictedLVESV
validIdx_LVESV = ~isnan(PredictedLVESV);
[r_LVESV, p_LVESV] = corr(UWpatients.LVESV(validIdx_LVESV), PredictedLVESV(validIdx_LVESV), 'Type', 'Pearson');

% 显示结果
disp(['LVEDV vs PredictedLVEDV: r = ', num2str(r_LVEDV), ', p = ', num2str(p_LVEDV)]);
disp(['LVESV vs PredictedLVESV: r = ', num2str(r_LVESV), ', p = ', num2str(p_LVESV)]);

% 可选：在散点图中显示 r 和 p 值
figure(66); scatter(UWpatients.LVEDV(validIdx_LVEDV), PredictedLVEDV(validIdx_LVEDV));
title(['LVEDV: r = ', num2str(r_LVEDV), ', p = ', num2str(p_LVEDV)]);
xlabel('LVEDV'); ylabel('PredictedLVEDV');

figure(77); scatter(UWpatients.LVESV(validIdx_LVESV), PredictedLVESV(validIdx_LVESV));
title(['LVESV: r = ', num2str(r_LVESV), ', p = ', num2str(p_LVESV)]);
xlabel('LVESV'); ylabel('PredictedLVESV');


% 计算和显示 PredictedLVEDV 与 LVEDV 的差值分布
DI1 = PredictedLVEDV - UWpatients.LVEDV;
figure;
histogram(DI1, 'Normalization', 'probability'); % 归一化为频率
xlabel('DI1 Values');
ylabel('Frequency');
title('Frequency Distribution of DI1');
