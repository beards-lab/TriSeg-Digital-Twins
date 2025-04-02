%% Script Summary:
% This script calls the functions: targetVal.m and estimParams.m, as well as the scripts: XXopt.m
% and runSim.m.
% The first section generates simulations for canonical male and female subjects, and the
% second section generates simulations for heart failure patients.

% Created by Feng Gu
% 03/20/2025 Done with modeling patients without CMR info.
% Till 03/20 UW cohort is purely fit without CMR. Umich cohort is still
% based on full data inputs(including CMR) and PKU cohort only have 1
% example.
% Last modified: 03/20/2025

%% Simulation for heart failure patients from Umich cohort without CMR info
clear
load AllPatients.mat
PredictedRVEDV = NaN(370,1);
PredictedRVEF = NaN(370,1);
secondslot = [17 79 206 256 288 325 352 355 360 361]; % patients with 2 3-month windows
for PATIENT_NO = 1 % any number between 1 and 370, example patient in paper is 192
    try
        for ModelWin =  1
            MRI_flag = 0;
            [Windowdate,targets, inputs, mods] = targetVals_HF(patients,PATIENT_NO,ModelWin,MRI_flag);
            [INIparams, INIinit] = estiminiParams(targets,inputs);
            RUNOPT = 0; % 0 for simulation and 1 for optimazation
            if RUNOPT ==1
                HFopt; % optimize modifiers of HF patients, Current HFopt is still used for full model
            end
            if RUNOPT == 0
                % load(sprintf('UmichSimsWithoutCMR/P_NO%dWindow%d.mat',PATIENT_NO,ModelWin)); % load results from great lake clusters
                load 'P_NO1 (2).mat'
                % load(sprintf('Sims0320/P_NO%d.mat',PATIENT_NO)); % load results from great lake clusters
                modifiers = output.m;
                % modifiers = ones(1,length(mods));
                % modifiers(15) = 0.5;
                % modifiers(16) = 0.5;
                % modifiers = [modifiers(1:4) 1 modifiers(5:end)];
                % modifiers(17) = 100;
                % modifiers(3) = 1;
                [params, init] = optParams(INIparams, INIinit, mods,modifiers);
                % params.K1 = params.K1/3;
                % params.R_PV = testparams.R_PV;
                % params.R_SV = testparams.R_SV;
                % params.k_pas_LV = testparams.k_pas_LV;
                % params.k_pas_RV = testparams.k_pas_RV;
                % params.k_act_LV = testparams.k_act_LV;
                % params.k_act_RV = testparams.k_act_RV;
                % params.Amref_LV = testparams.Amref_LV;
                % params.Amref_RV = testparams.Amref_RV;
                % params.Amref_SEP = testparams.Amref_SEP;
                % params.Vw_LV = testparams.Vw_LV;
                % params.Vw_RV = testparams.Vw_RV;
                % params.Vw_SEP = testparams.Vw_SEP;
                % params.C_SV = testparams.C_SV;
                % params.C_PV = testparams.C_PV;
                % params.C_SA = testparams.C_SA;
                % params.C_PA = testparams.C_PA;
                % params.R_SA = testparams.R_SA;
                % params.R_PA = testparams.R_PA;
                % params.R_tSA = testparams.R_tSA;
                % params.R_tPA = testparams.R_tPA;
                % params.LAV0c = testparams.LAV0c;
                % params.RAV0c = testparams.RAV0c;
                % params.R_PV = testparams.R_PV;
                % params.R_SV = testparams.R_SV;
                % params = testparams;
                % params.K_P = 1.87;
                % params.B_P = 0.12;
                % params.Vh0 = 885.5763;
                % init = testinit;
                % params.K1 = 3.15;
                % params.expPeri = 0.12;
                runSim;
                PredictedRVEDV(PATIENT_NO) = o_vals.RVEDV;
                PredictedRVEF(PATIENT_NO) = EF_RV;
                Print_cost = 1;% 1 for performance, other for patient's pre-condition (PHI sensitive)
                NplotFit; % 6-panel figure for HF patients
                % GetMovie;
                % See_TriSeg;
            end
        end
    catch ME
        disp(['Error: ', ME.message, 'PatID',num2str(PATIENT_NO)])
    end
end


%% Simulation for heart failure patients from UW cohort without CMR info
clear
load('UWcohort.mat');
tic
Cost = NaN(64,1);
RVcost = NaN(64,1);
for PATIENT_NO = 2
    try
        MRI_flag = 0;
        [targets, inputs, mods] = targetVals_UW(UWpatients, PATIENT_NO, MRI_flag);
        [INIparams, INIinit] = estiminiParams(targets,inputs);
        RUNOPT = 0; % 0 for simulation and 1 for optimazation
        if RUNOPT ==1
            UWopt; % optimize modifiers of HF patients without CMR info
        end
        if RUNOPT == 0
            % if MRI_flag == 1
            %     load(sprintf('SimsUWwithCMRLVFromTTE/P_NO%d.mat',PATIENT_NO)); % load results from great lake clusters
            % else
            %     load(sprintf('SimsUWwithoutCMR0/P_NO%d.mat',PATIENT_NO)); % load results from great lake clusters
            % end
            % modifiers = output.m;
            modifiers = ones(1,length(mods));
            [params, init] = optParams(INIparams, INIinit, mods,modifiers);
            evaluateModelUW(modifiers,UWpatients,PATIENT_NO)
            % if isfield(targets, 'EAr')
            %     targets = rmfield(targets, 'EAr'); % 移除字段 'b'
            % end
            runSim;
            Print_cost = 1;% 1 for performance, other for patient's pre-condition (PHI sensitive)
            NplotFit; % 6-panel figure for HF patientsP
            Cost(PATIENT_NO) = total_cost;
            RVcost (PATIENT_NO)  = ((o_vals.RVEDV - UWpatients.RVEDV(PATIENT_NO)) / std(UWpatients.RVEDV)).^2 + ...
                ((o_vals.RVESV - UWpatients.RVESV(PATIENT_NO)) / std(UWpatients.RVESV)).^2 + ...
                ((o_vals.Hed_RW - UWpatients.Hed_RW(PATIENT_NO)) / std(UWpatients.Hed_RW)).^2;
            % GetMovie;
            % See_TriSeg;
            % pause(3);
        end
    catch ME
        disp(['Error: ', ME.message, 'PatID',num2str(PATIENT_NO)])
    end
end
toc

%% Bland-Altman plot
% X = predictedTRW;
% Y = UWpatients.Hed_RW;
X = X(~isnan(X));
Y = Y(~isnan(X));
% Calculate mean and difference
meanXY = (X + Y) / 2;
diffXY = X - Y;

% Calculate limits of agreement (LoA)
meanDiff = mean(diffXY);  % Mean of differences
stdDiff = std(diffXY);    % Standard deviation of differences
LoA_upper = meanDiff + 1.96 * stdDiff;  % Upper limit
LoA_lower = meanDiff - 1.96 * stdDiff;  % Lower limit

% Plot Bland-Altman plot
figure;
scatter(meanXY, diffXY, 'bo');  % Scatter plot
hold on;
yline(meanDiff, 'r-', 'Mean Diff');  % Mean difference line
yline(LoA_upper, 'g--', 'Upper LoA');  % Upper limit line
yline(LoA_lower, 'g--', 'Lower LoA');  % Lower limit line

% Add labels and title
xlabel('Mean of X and Y');
ylabel('Difference (X - Y)');
title('Bland-Altman Plot');
grid on;
hold off;

%% Extract modifiers from cluster
for PATIENT_NO = 1:64
    if PATIENT_NO < 10
        matFilesDir = sprintf('Runs/p00%d/02-11/',PATIENT_NO);
    elseif PATIENT_NO < 100
        matFilesDir = sprintf('Runs/p0%d/02-11/',PATIENT_NO);
    else
        matFilesDir = sprintf('Runs/p%d/02-11/',PATIENT_NO);
    end
    matFiles = dir(fullfile(matFilesDir, '*.mat'));
    if ~isempty(matFiles)
        matFileName = matFiles.name;
        matFilePath = fullfile(matFilesDir, matFileName);
        GAresult = load(matFilePath);
        output.modifiers = GAresult.output.m;
        output.mods = GAresult.output.mods;
        save(sprintf('Sims_GA/P_NO%d.mat',PATIENT_NO),"output");
    end
end
%% 03/19 Simulation for heart failure patients from PKU cohort without CMR info (future...)
tic
for PATIENT_NO = 1
    try
        MRI_flag = 1;
        [targets, inputs, mods] = targetVals_PKU(MRI_flag);
        [INIparams, INIinit] = estiminiParams(targets,inputs);
        RUNOPT = 0; % 0 for simulation and 1 for optimazation
        if RUNOPT ==1
            if MRI_flag == 1
                PKUFullModelOpt; % optimize modifiers of HF patients
            else
                PKUFengRVValidationOpt; % optimize modifiers of HF patients without CMR info
            end
        end
        if RUNOPT == 0
            if MRI_flag == 1
                load PKUm.mat; % load results from great lake clusters
            else
                load PKUmFengRV.mat; % load results from great lake clusters
            end
            modifiers = output.modifiers;
            Geo_Opt = 1;
            [params, init] = optParams(INIparams, INIinit, mods,modifiers);
            runSim;
            Print_cost = 1;% 1 for performance, other for patient's pre-condition (PHI sensitive)
            NplotFit; % 6-panel figure for HF patientsP
            % GetMovie;
            See_TriSeg;
            pause(3);
        end
    catch ME
        disp(['Error: ', ME.message, 'PatID',num2str(PATIENT_NO)])
    end
end
toc

%% 03/31 Simulation for Mia and Pallak's LVAD patients
clear
load('AllPatients_preVAD_data.mat')
tic
for PATIENT_NO = 396:473
    try
        ModelWin = length(patients(PATIENT_NO).snapshots);
        MRI_flag = 0;
        [Windowdate,targets, inputs, mods] = targetVals_VAD(patients,PATIENT_NO);
        [INIparams, INIinit] = estiminiParams(targets,inputs);
        modifiers = ones(1,length(mods));
        [params, init] = optParams(INIparams, INIinit, mods,modifiers);
        runSim;
        % Print_cost = 1;% 1 for performance, other for patient's pre-condition (PHI sensitive)
        % NplotFit; % 6-panel figure for HF patientsP
        % GetMovie;
        % See_TriSeg;
        % pause(3);
    catch ME
        disp(['Error: ', ME.message, 'PatID',num2str(PATIENT_NO)])
        break
    end
end
toc
%%
T = readtable("RVoutcomes.xlsx");

figure(101);clf;
for i = 1:length(patients)
    [~,idx] = ismember(patients(i).snapshots(1).patKey,T.patkey);
    plot([min([patients(i).snapshots(1).RHCDate patients(i).snapshots(1).TTEDate])...
        T.VADOpDate(idx)],i*ones(1,2),'Color',[0.7 0.7 0.7,0.5]);
    hold on;
    scatter(T.VADOpDate(idx),i,'rx');
    numSnapshots = length(patients(i).snapshots);
    colors = jet(numSnapshots); 

    for j = 1:numSnapshots
        plot([min([patients(i).snapshots(j).RHCDate, patients(i).snapshots(j).TTEDate]) ...
            max([patients(i).snapshots(j).RHCDate, patients(i).snapshots(j).TTEDate])], ...
            i * ones(1, 2), 'Color', colors(j, :),'LineWidth',1.5);
        hold on;
    end
end