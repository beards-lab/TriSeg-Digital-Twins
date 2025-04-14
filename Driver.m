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

%% Simulation for Canonical Subjects
clear
MRI_flag = 1;
for GENDER  = 1:2 % 1 for male and other for female
    if GENDER == 1
        [targets, inputs, mods] = targetVals_male();
        modifiers = readmatrix('modifiers_male.csv','range','2:2');
        % modifiers = ones(1,length(mods));
        [INIparams, INIinit] = estiminiParams(targets,inputs);
        % params.HR = 180;
        % params.T = 60/params.HR;
        [params, init] = optParams(INIparams, INIinit, mods,modifiers);
        BSA = 1.9;
        ParamM.C_SA = params.C_SA/BSA;
        ParamM.C_SV = params.C_SV/BSA;
        ParamM.C_PA = params.C_PA/BSA;
        ParamM.C_PV = params.C_PV/BSA;
        ParamM.R_SA = BSA*params.R_SA;
        ParamM.R_PA = BSA*params.R_PA;
        ParamM.R_tSA = BSA*params.R_tSA;
        ParamM.R_tPA = BSA*params.R_tPA;
        ParamM.k_act_LV = params.k_act_LV;
        ParamM.k_act_RV = params.k_act_RV;
        ParamM.k_pas_LV = params.k_pas_LV;
        ParamM.k_pas_RV = params.k_pas_RV;
        ParamM.Vw_LV = params.Vw_LV/BSA;
        ParamM.Vw_RV = params.Vw_RV/BSA;
        ParamM.Vw_SEP = params.Vw_SEP/BSA;
        ParamM.Amref_LV = params.Amref_LV/BSA;
        ParamM.Amref_RV = params.Amref_RV/BSA;
        ParamM.Amref_SEP = params.Amref_SEP/BSA;
        ParamM.R_SV = BSA*params.R_SV;
        ParamM.R_PV = BSA*params.R_PV;
        ParamM.Vh0 = params.Vh0/BSA;ParamM.K1 = params.K1;ParamM.expPeri = params.expPeri;
        ParamM.R_t_o = BSA*params.R_t_o;
        ParamM.R_t_c = BSA*params.R_t_c;
        ParamM.R_p_o = BSA*params.R_p_o;
        ParamM.R_p_c = BSA*params.R_p_c;
        ParamM.R_m_o = BSA*params.R_m_o;
        ParamM.R_m_c = BSA*params.R_m_c;
        ParamM.R_a_o = BSA*params.R_a_o;
        ParamM.R_a_c = BSA*params.R_a_c;
    else
        [targets, inputs, mods] = targetVals_female();
        modifiers = readmatrix('modifiers_female.csv','range','2:2');
        [INIparams, INIinit] = estiminiParams(targets,inputs);
        % params.HR = 180;
        % params.T = 60/params.HR;
        [params, init] = optParams(INIparams, INIinit, mods,modifiers);
        BSA = 1.6;
        ParamF.C_SA = params.C_SA/BSA;
        ParamF.C_SV = params.C_SV/BSA;
        ParamF.C_PA = params.C_PA/BSA;
        ParamF.C_PV = params.C_PV/BSA;
        ParamF.R_SA = BSA*params.R_SA;
        ParamF.R_PA = BSA*params.R_PA;
        ParamF.R_tSA = BSA*params.R_tSA;
        ParamF.R_tPA = BSA*params.R_tPA;
        ParamF.k_act_LV = params.k_act_LV;
        ParamF.k_act_RV = params.k_act_RV;
        ParamF.k_pas_LV = params.k_pas_LV;
        ParamF.k_pas_RV = params.k_pas_RV;
        ParamF.Vw_LV = params.Vw_LV/BSA;
        ParamF.Vw_RV = params.Vw_RV/BSA;
        ParamF.Vw_SEP = params.Vw_SEP/BSA;
        ParamF.Amref_LV = params.Amref_LV/BSA;
        ParamF.Amref_RV = params.Amref_RV/BSA;
        ParamF.Amref_SEP = params.Amref_SEP/BSA;

        ParamF.R_SV = BSA*params.R_SV;
        ParamF.R_PV = BSA*params.R_PV;
        ParamF.Vh0 = params.Vh0/BSA;ParamF.K1 = params.K1;ParamF.expPeri = params.expPeri;
        ParamF.R_t_o = BSA*params.R_t_o;
        ParamF.R_t_c = BSA*params.R_t_c;
        ParamF.R_p_o = BSA*params.R_p_o;
        ParamF.R_p_c = BSA*params.R_p_c;
        ParamF.R_m_o = BSA*params.R_m_o;
        ParamF.R_m_c = BSA*params.R_m_c;
        ParamF.R_a_o = BSA*params.R_a_o;
        ParamF.R_a_c = BSA*params.R_a_c;


        % modifiers = ones(1,length(mods));
    end
    Sim = 1; % 1 for simulation and other for optimazation
    if Sim ==1

        runSim;
        NplotSrd; % 4-panel figure just for canonical subjects
        % GetMovie; % TriSeg model: displacement and stress as functions of time
        % See_TriSeg; % slices of GetMoive.m
    else
        Srdopt; % optimize modifiers of canonical subjects
    end
end
GENDER = [1;2];
standardparamstable = table(GENDER);
for No = 1:length(GENDER)
    paramsname = [fieldnames(ParamM)];
    for j = 1:length(paramsname)
        if No == 1
            if isfield(ParamM,paramsname{j})
                p = ParamM.(paramsname{j});
                % elseif isfield(initM,paramsname{j})
                %     p = initM.(paramsname{j});
                % elseif isfield(PredictionM,paramsname{j})
                %     p = PredictionM.(paramsname{j});
            end
        elseif No ==2
            if isfield(ParamF,paramsname{j})
                p = ParamF.(paramsname{j});
                % elseif isfield(initF,paramsname{j})
                %     p = initF.(paramsname{j});
                % elseif isfield(PredictionF,paramsname{j})
                %     p = PredictionF.(paramsname{j});
            end
        end
        standardparamstable.(paramsname{j})(No) = p;
    end
end
save('standardnormalizationNew',"standardparamstable");


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
                % load 'P_NO1 (1).mat'
                load(sprintf('Sims0401/P_NO%d.mat',PATIENT_NO)); % load results from great lake clusters
                modifiers = output.m;
                % modifiers = ones(1,length(mods));
                % modifiers(15) = 0.5;
                % modifiers(16) = 0.5;
                % modifiers = [modifiers(1:4) 1 modifiers(5:end)];
                % modifiers(17) = 100;
                % modifiers(3) = 1;
                [params, init] = optParamsWrongVersion(INIparams, INIinit, mods,modifiers);
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
                figure(2);clf;plot(V_LV(1:end/2),P_LV(1:end/2));hold on;
                scatter(V_LV(1),P_LV(1));
                scatter(V_LV(Qa_pos_end),P_LV(Qa_pos_end));
                [~,idxESStress] = max(sigma_act_LV(1:end/2));
                scatter(V_LV(idxESStress),P_LV(idxESStress));
                PVSearchingPeriod = [V_LV(idxESStress:Qa_pos_end) P_LV(idxESStress:Qa_pos_end)];
                PVSearchingPeriod = sortrows(PVSearchingPeriod, 1);
                dPdV = gradient(PVSearchingPeriod(:,2),PVSearchingPeriod(:,1));
                [~,idxESLVEmaxindPdV] = min(abs(median(dPdV(dPdV>0))-dPdV));
                idxESLVEmax = Qa_pos_end - idxESLVEmaxindPdV;

                scatter(V_LV(idxESLVEmax),P_LV(idxESLVEmax));
                Vtmax = V_LV(idxESLVEmax);
                tESavc = t(Qa_pos_end);

                tN35 = tESavc*0.35;
                [~, idx_tN35] = min(abs(t - tN35));
                PnTn35 = P_LV(idx_tN35)./P_LV(Qa_pos_end);

                Vtn = V_LV(idx_tN35);
                ENTN35 = 0.4739;
                V0SB35 = (PnTn35*Vtmax-Vtn*ENTN35)/(PnTn35-ENTN35);

                tN30 = tESavc*0.30;
                [~, idx_tN30] = min(abs(t - tN30));
                PnTn30 = P_LV(idx_tN30)./P_LV(Qa_pos_end);

                Vtn = V_LV(idx_tN30);
                ENTN30 = 0.4293;
                V0SB30 = (PnTn30*Vtmax-Vtn*ENTN30)/(PnTn30-ENTN30);

                tN25 = tESavc*0.25;
                [~, idx_tN25] = min(abs(t - tN25));
                PnTn25 = P_LV(idx_tN25)./P_LV(Qa_pos_end);

                Vtn = V_LV(idx_tN25);
                ENTN25 = 0.3802;
                V0SB25 = (PnTn25*Vtmax-Vtn*ENTN25)/(PnTn25-ENTN25);

                V0SB = [V0SB35 V0SB30 V0SB25];
                V0SB = mean(V0SB(V0SB>0));
                if isnan(V0SB)
                    V0SB = 0;
                end
                Linex = [ V_LV(idxESLVEmax) V0SB];
                Liney = [ P_LV(idxESLVEmax) 0];
                plot(Linex,Liney);

                % V0 = V_LV(1) * (0.6 - 0.006 * P_LV(1));
                V0 = V0SB;

                % 计算 V30
                V30 = V0 + ((V_LV(1) - V0) / ((P_LV(1) / 27.78)^(1/2.76)));

                % 计算 V15
                V15 = 0.8 * (V30 - V0) + V0;

                % 计算 α 和 β
                if P_LV(1) <= 22
                    beta = log(P_LV(1) / 30) / log(V_LV(1) / V30);
                    alpha = 30 / (V30^beta);
                else
                    beta = log(P_LV(1) / 15) / log(V_LV(1) / V15);
                    alpha = P_LV(1) / (V_LV(1)^beta);
                end

                % 生成 EDPVR 曲线
                PEDV = linspace(V0, V_LV(1), 100); % 使用 linspace 生成平滑点
                PEDP = alpha * PEDV.^beta;

                plot(PEDV, PEDP);
                xlabel('Volume (mL)');
                ylabel('Pressure (mmHg)');
                title('Pressure-Volume Relationship (PV-Loop)');
                grid on;

                x = [V_LV(1:idxESLVEmax); Linex'; PEDV']; % 所有x坐标
                y = [P_LV(1:idxESLVEmax); Liney'; PEDP']; % 所有y坐标
                Area = polyarea(x, y); % 计算多边形面积
                LVSWandPMW(PatID) = (2.08*Area+0.549)*params.HR ;
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
PredictedRVEDV = NaN(64,1);
PredictedRVEF = NaN(64,1);
PredictedTRW = NaN(64,1);
for PATIENT_NO = 1:64
    try
        MRI_flag = 0;
        [targets, inputs, mods] = targetVals_UW(UWpatients, PATIENT_NO, MRI_flag);
        [INIparams, INIinit] = estiminiParams(targets,inputs);
        RUNOPT = 0; % 0 for simulation and 1 for optimazation
        if RUNOPT ==1
            UWopt; % optimize modifiers of HF patients without CMR info
        end
        if RUNOPT == 0
            if MRI_flag == 1
                load(sprintf('SimsUWwithCMRLVFromTTE/P_NO%d.mat',PATIENT_NO)); % load results from great lake clusters
            else
                load(sprintf('Sims0401UW/P_NO%d.mat',PATIENT_NO)); % load results from great lake clusters
            end
            modifiers = output.m;
            % modifiers = ones(1,length(mods));
            [params, init] = optParamsWrongVersion(INIparams, INIinit, mods,modifiers);
            % if isfield(targets, 'EAr')
            %     targets = rmfield(targets, 'EAr'); % 移除字段 'b'
            % end

            % params = output.params;
            % params.K1 = 1;
            % params.expPeri = 0.4;
            % init = output.init;
            runSim;
            Print_cost = 1;% 1 for performance, other for patient's pre-condition (PHI sensitive)
            NplotFit; % 6-panel figure for HF patientsP
            Cost(PATIENT_NO) = total_cost;
            RVcost (PATIENT_NO)  = ((o_vals.RVEDV - UWpatients.RVEDV(PATIENT_NO)) / std(UWpatients.RVEDV)).^2 + ...
                ((o_vals.RVESV - UWpatients.RVESV(PATIENT_NO)) / std(UWpatients.RVESV)).^2 + ...
                ((o_vals.Hed_RW - UWpatients.Hed_RW(PATIENT_NO)) / std(UWpatients.Hed_RW)).^2;
            PredictedRVEDV(PATIENT_NO) = o_vals.RVEDV;
            PredictedRVEF(PATIENT_NO) = EF_RV;
            PredictedTRW(PATIENT_NO) = o_vals.Hed_RW;
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
X = PredictedRVEDV;
Y = UWpatients.RVEDV;
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
        MRI_flag = 0;
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
DTTable = readtable("LastWindowDTinfo.csv",VariableNamingRule="preserve");

%%
tic
for PATIENT_NO = 25
    try
        if ismember(PATIENT_NO,[40 162 241 349 424])
            ModelWin = length(patients(PATIENT_NO).snapshots)-2;
        elseif ismember(PATIENT_NO,[67 87 111 115 139 161 166 183 189 191 201 211 234 261 278 284 291 294 307 331 346 377 382 389 437 458])
            ModelWin = length(patients(PATIENT_NO).snapshots)-1;
        else
            ModelWin = length(patients(PATIENT_NO).snapshots);
        end
        MRI_flag = 0;
        % windowslot = mean([patients(PATIENT_NO).snapshots.RHCDate;patients(PATIENT_NO).snapshots.TTEDate],1);
        % targetstable = readtable("share_addmisionDate_VADdate.xlsx","VariableNamingRule","preserve");
        % [~,idx] = ismember(patients(PATIENT_NO).snapshots(1).patKey,targetstable.patkey);
        % targetsDate = targetstable.date_admit(idx);
        % [~, ModelWin] = min(abs(windowslot -  targetsDate));
        try
            [Windowdate,targets, inputs, mods] = targetVals_VAD(patients,PATIENT_NO,ModelWin);
            [INIparams, INIinit] = estiminiParams(targets,inputs);
        catch
            ModelWin = ModelWin+1;
            [Windowdate,targets, inputs, mods] = targetVals_VAD(patients,PATIENT_NO,ModelWin);
            [INIparams, INIinit] = estiminiParams(targets,inputs);
        end
        load(sprintf('Sims0402VAD/P_NO%d.mat',PATIENT_NO)); % load results from great lake clusters
        modifiers = output.m;
        % modifiers = ones(1,length(mods));
        [params, init] = optParamsWrongVersion(INIparams, INIinit, mods,modifiers);
        runSim;
        % Prediction.LVpowerIndex_S = LVpower/BSA;
        % Prediction.RVpowerIndex_S = RVpower/BSA;
        % Prename = fieldnames(Prediction);
        % for j = 1:length(Prename)
        %     p = Prediction.(Prename{j});
        %     DTTable.(Prename{j})(PATIENT_NO ) = p;
        % end
        % outputname = fieldnames(o_vals);
        % for j = 1:length(outputname)
        %     o = o_vals.(outputname{j});
        %     DTTable.(outputname{j}+"_S")(PATIENT_NO ) = o;
        % end
        % DTTable.minRVP_S(PATIENT_NO) = o_vals.P_RV_min;
        % DTTable.minLVP_S(PATIENT_NO) = o_vals.P_LV_min;
        % DTTable.LVm_S(PATIENT_NO) = o_vals.LV_m;
        % DTTable.RVm_S(PATIENT_NO) = o_vals.RV_m;
        % DTTable.HedLW_S(PATIENT_NO) = o_vals.Hed_LW;
        % DTTable.HedSW_S(PATIENT_NO) = o_vals.Hed_SW;
        Print_cost = 1;% 1 for performance, other for patient's pre-condition (PHI sensitive)
        NplotFit; % 6-panel figure for HF patients
        % GetMovie;
        See_TriSeg;
        pause(3);
    catch ME
        disp(['Error: ', ME.message, 'PatID',num2str(PATIENT_NO)])
        % break
    end
end
toc
%% Figure for Mia and Pallak
clear
load('AllPatients_preVAD_data.mat')
Table = readtable("RVoutcomes.xlsx");
figure(101);clf;
for i = 2
    [~,idx] = ismember(patients(i).snapshots(1).patKey,Table.patkey);
    plot([min([patients(i).snapshots(1).RHCDate patients(i).snapshots(1).TTEDate])...
        Table.VADOpDate(idx)],i*ones(1,2),'Color',[0.7 0.7 0.7,0.5]);
    hold on;
    scatter(Table.VADOpDate(idx),i,'rx');
    numSnapshots = length(patients(i).snapshots);
    colors = jet(numSnapshots);

    for j = 1:numSnapshots
        plot([min([patients(i).snapshots(j).RHCDate, patients(i).snapshots(j).TTEDate]) ...
            max([patients(i).snapshots(j).RHCDate, patients(i).snapshots(j).TTEDate])], ...
            i * ones(1, 2), 'Color', colors(j, :),'LineWidth',1.5);
        hold on;
    end
end


%%
DTTable = readtable("LastWindowDTinfo.csv");
listorder = (1:473);
PatientNumber =listorder(~isnan(DTTable.C_SA));
TableList = cell(473,1);
for i = 1:473
    PATIENT_NO = PatientNumber(i);
    % if ismember(PATIENT_NO,[40 162 241 349 424])
    %     ModelWin = length(patients(PATIENT_NO).snapshots)-2;
    % elseif ismember(PATIENT_NO,[67 87 111 115 139 161 166 183 189 191 201 211 234 261 278 284 291 294 307 331 346 377 382 389 437 458])
    %     ModelWin = length(patients(PATIENT_NO).snapshots)-1;
    % else
    %     ModelWin = length(patients(PATIENT_NO).snapshots);
    % end
    windowslot = mean([patients(PATIENT_NO).snapshots.RHCDate;patients(PATIENT_NO).snapshots.TTEDate],1);
    targetstable = readtable("share_addmisionDate_VADdate.xlsx","VariableNamingRule","preserve");
    [~,idx] = ismember(patients(PATIENT_NO).snapshots(1).patKey,targetstable.patkey);
    targetsDate = targetstable.date_admit(idx);
    [~, ModelWin] = min(abs(windowslot -  targetsDate));
    DataToConvert = patients(PATIENT_NO).snapshots(ModelWin);
    OutputTable = struct2table(DataToConvert);

    charCols = varfun(@(x) ischar(x), OutputTable, 'OutputFormat', 'uniform');

    for colIdx = find(charCols)
        OutputTable.(colIdx) = string(OutputTable.(colIdx));
    end
    if PATIENT_NO == 1
        DataTable = OutputTable;
    else
        DataTable(PATIENT_NO,:) = OutputTable;
    end
end
DataTable.Age = round(years(mean([DataTable.RHCDate DataTable.TTEDate],2)-DataTable.Birthday));
%%
writetable(DataTable,'Close2AdmitData.csv');
%%
T1 = readtable("Close2AdmitData.csv");
T2 = readtable("LastWindowData.csv");
samewindow = [];
for  i = 1:height(T1)
    if ~isnat(T1.RHCDate(i))
        if T1.rhcId(i) == T2.rhcId(i) &&...
                T1.tteId(i) == T2.tteId(i)
        else
            samewindow = [samewindow i];
        end
    end
end

%% Run Haolin's Virtual paients
clear
Table = readtable("generated_fake_data.csv");
MRI_flag = 1:64;
for PATIENT_NO = 1:64
    try
        FakeParamsT = Table(PATIENT_NO,:);
        FakeParamsT.RAV0u = FakeParamsT.LAV0u*1.1;
        FakeParamsT.expPeri = 0.4;
        [params, init] = ProcessParamsFromVAE(FakeParamsT);
        runSimonFakeDT;
        PlotFakeDTSim;
    catch ME
        disp({ME.message num2str(PATIENT_NO)})
    end
end
