%% Script Summary:
% This script calls the functions: targetVal.m and estimParams.m, as well as the scripts: XXopt.m
% and runSim.m.
% The first section generates simulations for canonical male and female subjects, and the
% second section generates simulations for heart failure patients.

% Created by Feng Gu
% Last modified: 10/29/2024

%% Simulation for Canonical Subjects
clear
<<<<<<< Updated upstream
=======
MRI_flag = 1;
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
for GENDER  = 1 % 1 for male and other for female
    if GENDER == 1
        [targets, inputs, mods] = targetVals_male();
        modifiers = readmatrix('modifiers_male.csv','range','2:2');
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
        [~,params, init] = estimParams(targets,inputs,mods,modifiers);
    else
        [targets, inputs, mods] = targetVals_female();
        modifiers = readmatrix('modifiers_female.csv','range','2:2');
        [~,params, init] = estimParams(targets,inputs,mods,modifiers);
    end
    Sim = 1; % 1 for simulation and other for optimazation
    if Sim ==1
        runSim; 
        NplotSrd; % 4-panel figure just for canonical subjects
        GetMovie; % TriSeg model: displacement and stress as functions of time
        See_TriSeg; % slices of GetMoive.m
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
        % modifiers = zeros(1,length(mods));
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
        % modifiers = zeros(1,length(mods));
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
        % theta = 1;
        params.HR = 200;
        params.T = 60/params.HR;
        params.R_SA = params.R_SA*0.25;
        params.C_SV = params.C_SV*0.25;
        params.R_PA = params.R_PA*0.5;
        params.R_PV = params.R_PV*0.5;
        params.k_act_LV = params.k_act_LV*3;
        params.k_act_RV = params.k_act_RV*3;
        runSimOnGL;
        NplotSrd; % 4-panel figure just for canonical subjects
        % GetMovie; % TriSeg model: displacement and stress as functions of time
        % See_TriSeg; % slices of GetMoive.m
        
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
    else
        Srdopt; % optimize modifiers of canonical subjects
    end
end
%% Simulation for heart failure patients
clear
load AllPatients.mat
secondslot = [17 79 206 256 288 325 352 355 360 361]; % patients with 2 3-month windows
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
% Shuntlist = [34 41 54 61 83 116 183 231 268 278 312]; % patients with shunt 
for PatID = 192 % any number between 1 and 370, example patient in paper is 192
    for ModelWin =  1
        [Windowdate,targets, inputs, mods] = targetVals_HF(patients,PatID,ModelWin);
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
% Order = [1	2 13 15	16	17	20	22	25	26	35	37	42	45	46	48	50	55	57	60	62	66	67	77	78	79	84	87	88	90	91	92	93	94	97	98	99	102	103	104	106	108	111	112	113	118	122	126	130	131	132	133	137	141	146	147	150	151	164	167	168	173	177	178	186	188	195	196	197	199	202	203	204	205	206	207	209	212	213	214	217	218	222	223	224	230	232	234	235	236	240	244	246	248	250	252	267	269	271	274	281	282	283	284	293	294	300	302	303	304	306	309	316	318	324	325	327	329	335	341	342	343	344	345	346	350	351	368];

for PATIENT_NO =46 % any number between 1 and 370, example patient in paper is 192
    try
        for ModelWin =  1
            % Cost_slot = [];
            MRI_flag = 1;
            % [Windowdate,targets, inputs, mods] = targetVals_HF(patients,PATIENT_NO,ModelWin,MRI_flag);
            % [INIparams, INIinit] = estiminiParams(targets,inputs);
            RUNOPT = 0; % 0 for simulation and 1 for optimazation
            if RUNOPT ==1
                HFopt; % optimize modifiers of HF patients, Current HFopt is still used for full model
            end
            if RUNOPT == 0
                load(sprintf('SimsUMFinal/P_NO%d.mat',PATIENT_NO)); % load results from great lake clusters
                targets = output.targets;
                inputs = output.inputs;
                params = output.params;
                init = output.init;
                % params.k_act_RV = 1.0651e5;
                % params.k_pas_RV = 1.0651e4;
                % params.R_tPA = 0.0008;
                % % params.R_p_o = 0.008;
                % params.C_PA = 10;
                % params.Vw_LV = 150.56;
                % params.Vw_SEP = 100.37;
                % output.params = params;
                % save(sprintf('SimsUMFinal/P_NO%d',PATIENT_NO), "output");
                runSimOnGL;
            end
        end
        % Print_cost = 1;% 1 for performance, other for patient's pre-condition (PHI sensitive)
        % NplotFit; % 6-panel figure for HF patients
        % filename = fullfile(save_dir, sprintf('Patient_%02d.png', PATIENT_NO));
        % exportgraphics(gcf, filename, 'Resolution', 300);  % gcf = current figure
        % % GetMovie;
        % See_TriSeg;
        [PVALV, PVARV, V0LV, V0RV,PVloop_area_LV,PVloop_area_RV] = Calculate_PVA(init, params, MRI_flag, inputs, targets);
        % MVO2_LV = (1.56e-5 * PVALV + 0.00526 * o_vals.LV_m / 100) * inputs.HR;% https://academic.oup.com/eurheartj/article/13/suppl_E/85/486078#google_vignette
        % Energy_LV = 1.333e-4*PVloop_area_LV*inputs.HR;
        % CardiacEffiency_LV =  Energy_LV/(MVO2_LV*20.1);
        % MVO2_RV = (1.56e-5 * PVARV + 0.000526 * o_vals.RV_m / 100) * inputs.HR;% https://academic.oup.com/eurheartj/article/13/suppl_E/85/486078#google_vignette
        % Energy_RV = 1.333e-4*PVloop_area_RV*inputs.HR;
        % CardiacEffiency_RV =  Energy_RV/(MVO2_RV*20.1); % Based on glucose, protein 18.6, fatty acid 19.6
        % clearvars -except Y
    catch ME
        disp(['Error: ', ME.message, 'PatID',num2str(PATIENT_NO)])
    end
end


%% Simulation for heart failure patients from UW cohort without CMR info
clear
save_dir = 'FiguresFinal';
% if ~exist(save_dir, 'dir')
%     mkdir(save_dir);  % 如果目录不存在，则创建
% end
load('UWcohort.mat');
tic
Cost = NaN(64,1);
RVcost = NaN(64,1);
PredictedRVEDV = NaN(64,1);
PredictedRVEF = NaN(64,1);
PredictedTRW = NaN(64,1);
Order = [9	14	15	17	18	19	22	25	26	30	33	34	37	44	46	48	52	53];
for PATIENT_NO = 28
    MRI_flag = 1;
    load(sprintf('SimsUWFinal/P_NO%d.mat',PATIENT_NO)); % load results from great lake clusters
    targets = output.targets;
    inputs = output.inputs;
    params = output.params;
    % params.k_act_LV = 1.0939e+05;
    % params.k_act_RV = 1.2113e+04;
    % % params.R_tSA = 0.05;
    % % params.R_tPA = 0.05;
    % % params.C_PA = 5;
    % params.C_SA = 5;
    % params.R_SA = 2;
    % params.R_PA = 1;
    % % % % % params.Vw_LV = 150.56;
    % % % % % params.Vw_SEP = 100.37;
    % output.params = params;
    % save(sprintf('SimsUWFinal/P_NO%d',PATIENT_NO), "output");
    init = output.init;
    runSimOnGL;
    Print_cost = 1;% 1 for performance, other for patient's pre-condition (PHI sensitive)
    NplotFit; % 6-panel figure for HF patients
    filename = fullfile(save_dir, sprintf('UWPatient_%02d.png', PATIENT_NO));
    exportgraphics(gcf, filename, 'Resolution', 300);  % gcf = current figure
    % Cost(PATIENT_NO) = total_cost;
    % RVcost (PATIENT_NO)  = ((o_vals.RVEDV - UWpatients.RVEDV(PATIENT_NO)) / std(UWpatients.RVEDV)).^2 + ...
    %     ((o_vals.RVESV - UWpatients.RVESV(PATIENT_NO)) / std(UWpatients.RVESV)).^2 + ...
    %     ((o_vals.Hed_RW - UWpatients.Hed_RW(PATIENT_NO)) / std(UWpatients.Hed_RW)).^2;
    % PredictedRVEDV(PATIENT_NO) = o_vals.RVEDV;
    % PredictedRVEF(PATIENT_NO) = EF_RV;
    % PredictedTRW(PATIENT_NO) = o_vals.Hed_RW;
    % GetMovie;
    % See_TriSeg;
    % pause(3);
end

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
>>>>>>> Stashed changes
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
<<<<<<< Updated upstream
=======

xline(0, '-b', 'LineWidth', 1.5);  % Admission 参考线
xlabel('Days from admission');
xlim([-5 30])
ylim([-5 105])
% ylabel('Patient index (Top 100 with most RHC+TTE)');
% title('Top 100 patients with most RHC and TTE events (within ~20-day stay)');
set(gca, 'YDir','reverse','FontName','Arial','FontSize',10);
pbaspect([1.33,2,1])
yticks([ ]);
box on;
% 设置图像尺寸（单位为英寸，1 cm ≈ 0.3937 in）
set(gcf, 'Units', 'inches', 'Position', [1 1 2.66 4]);  % 9 cm x 9 cm ≈ 3.54 in x 3.54 in
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3.54 3.54]);





%%
DTTable = readtable("Close2AdmitDTinfo.csv");
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
    windowslot = [patients(PATIENT_NO).snapshots.RHCDate;patients(PATIENT_NO).snapshots.TTEDate];
    targetstable = readtable("share_addmisionDate_VADdate.xlsx","VariableNamingRule","preserve");
    [~,idx] = ismember(patients(PATIENT_NO).snapshots(1).patKey,targetstable.patkey);
    targetsDate = targetstable.date_admit(idx);
    distance = abs(windowslot(1,:) - targetsDate) + abs(windowslot(2,:) - targetsDate);
    [~, ModelWin] = min(distance);
    if ismember(PATIENT_NO,[31 38 52 103 267 295])
        ModelWin = ModelWin+1;
    elseif ismember(PATIENT_NO,[284 307 331])
        ModelWin = ModelWin-1;
    end
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
writetable(DTTable,'Close2AdmitDTinfo.csv');
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
errorLog = {};  % 创建空cell数组保存错误信息
logIdx = 1;      % 初始化日志索引
Table = readtable("VirtualParams.csv");
MRI_flag = 1;
%%
for PATIENT_NO = 18573
    try
        FakeParamsT = Table(PATIENT_NO,1:34);
        % FakeParamsT.RAV0u = FakeParamsT.LAV0u*1.1;
        % FakeParamsT.expPeri = 0.4;
        [params, init] = ProcessParamsFromVAE(FakeParamsT);
        runSimonFakeDT;
        PlotFakeDTSim;
    catch ME
        disp({ME.message num2str(PATIENT_NO)})
        errorLog{logIdx,1} = PATIENT_NO;
        errorLog{logIdx,2} = ME.message;
        logIdx = logIdx + 1;
    end
end
%%
params = output.params;
init = output.init;
runSimonFakeDT;
PlotFakeDTSim;
%% Get Training Data for VAE
clear
RawT = readtable("inputFGclassifier.csv",'VariableNamingRule','preserve');
MRI_flag = 1;
desiredFields = RawT.Properties.VariableNames; % 只填已有的列
% UMsource = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 22 23 24 25 26 27 28 29 30 31 32 33 35 36 37 38 39 40 42 43 44 45 46 47 48 49 50 51 52 53 55 56 57 58 59 60 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 108 109 110 111 112 113 114 115 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 180 181 182 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 232 233 234 235 236 237 238 239 240 241 242 243 244 246 247 248 249 250 251 252 253 254 255 256 257 259 260 261 262 263 264 265 266 267 269 270 271 272 273 274 275 276 277 279 280 281 282 283 284 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300 301 302 303 304 305 306 307 308 309 310 311 313 314 315 316 317 318 319 320 321 322 323 324 325 326 327 328 329 330 331 332 333 334 335 336 337 338 339 340 341 342 343 344 345 346 347 348 349 350 351 353 354 357 358 360 365 367 368 369 370];
UMsource = [1	2 13 15	16	17	20	22	25	26	35	37	42	45	46	48	50	55	57	60	62	66	67	77	78	79	84	87	88	90	91	92	93	94	97	98	99	102	103	104	106	108	111	112	113	118	122	126	130	131	132	133	137	141	146	147	150	151	164	167	168	173	177	178	186	188	195	196	197	199	202	203	204	205	206	207	209	212	213	214	217	218	222	223	224	230	232	234	235	236	240	244	246	248	250	252	267	269	271	274	281	282	283	284	293	294	300	302	303	304	306	309	316	318	324	325	327	329	335	341	342	343	344	345	346	350	351	368];

load AllPatients.mat
for i = 1:length(UMsource)
    PATIENT_NO = UMsource(i);
    RawT.Patient_No(i) = PATIENT_NO;
    try
        % 加载当前病人的.mat文件
        S = load(sprintf('SimsUMFinal/P_NO%d.mat', PATIENT_NO));
        output = S.output;
        PostProcessSim;
        output.o_vals = o_vals;
        % 按结构体逐个查找
        structNames = {'inputs', 'params', 'init', 'targets','o_vals'};
        for s = 1:length(structNames)
            structName = structNames{s};
            if isfield(output, structName)
                thisStruct = output.(structName);
                fields = fieldnames(thisStruct);
                for f = 1:length(fields)
                    varName = fields{f};
                    if ismember(varName, desiredFields)
                        RawT{i, varName} = thisStruct.(varName);
                    end
                end
            end
        end

    catch ME
        fprintf('病人 %d 加载失败：%s\n', PATIENT_NO, ME.message);
    end
end

%%
RawT = readtable("inputFGclassifier.csv",'VariableNamingRule','preserve');
UWsource = [9 14	15	17	18	19	22	25	26	30	33	34 37	44	46	48	52	53];
for i = 1:length(UWsource)
    PATIENT_NO = UWsource(i);
    try
        % 加载当前病人的.mat文件
        S = load(sprintf('SimsUWFinal/P_NO%d.mat', PATIENT_NO));
        output = S.output;
        PostProcessSim;
        output.o_vals = o_vals;
        % 按结构体逐个查找
        structNames = {'inputs', 'params', 'init', 'targets','o_vals'};
        for s = 1:length(structNames)
            structName = structNames{s};
            if isfield(output, structName)
                thisStruct = output.(structName);
                fields = fieldnames(thisStruct);

                for f = 1:length(fields)
                    varName = fields{f};
                    if ismember(varName, desiredFields)
                        RawT{i+128, varName} = thisStruct.(varName);
                    end
                end
            end
        end

    catch ME
        fprintf('病人 %d 加载失败：%s\n', PATIENT_NO, ME.message);
    end
end
%%
clear
RawT = readtable("inputFGclassifier.csv",'VariableNamingRule','preserve');
T = readtable('VirtualParams.csv');
HRList = T.HR;
MRI_flag = 1;
desiredFields = RawT.Properties.VariableNames; % 只填已有的列
%%
for PATIENT_NO = 1:24531
    RawT.Patient_No(PATIENT_NO) = PATIENT_NO;
    try
        % 加载当前病人的.mat文件
        thisHR = HRList(PATIENT_NO);
        foldername = sprintf('PacingSimsVirtual/PatientNO%d', PATIENT_NO);
        Newfilename = sprintf('%s/HR%d.mat', foldername, thisHR);
        S = load(Newfilename);
        output = S.output;
        % 按结构体逐个查找
        structNames = {'params', 'init', 'o_vals'};
        for s = 1:length(structNames)
            structName = structNames{s};
            if isfield(output, structName)
                thisStruct = output.(structName);
                fields = fieldnames(thisStruct);
                for f = 1:length(fields)
                    varName = fields{f};
                    if ismember(varName, desiredFields)
                        RawT{PATIENT_NO, varName} = thisStruct.(varName);
                    end
                end
            end
        end

    catch ME
        fprintf('病人 %d 加载失败：%s\n', PATIENT_NO, ME.message);
    end
end
%%
writetable(RawT,'inputFGclassifier146.CSV')
%%
RawT.Sex = RawT.Sex-1;
%%
RawT.Height = round(RawT.Height,1);
RawT.Weight = round(RawT.Weight,1);
RawT.SBP = round(RawT.SBP);
RawT.DBP = round(RawT.DBP);
RawT.EF = round(RawT.EF,1);
RawT.Hed_SW = round(RawT.Hed_SW*10,1);
RawT.Hed_LW = round(RawT.Hed_LW*10,1);
RawT.HR = round(RawT.HR,1);
RawT.LAVmax = round(RawT.LAVmax);
RawT.EAr = round(RawT.EAr,1);
RawT.IVRT = round(RawT.IVRT*1000);
RawT.LVIDd = round(RawT.LVIDd*10,1);
RawT.LVIDs = round(RawT.LVIDs*10,1);
RawT.LVEDV = round(RawT.LVEDV);
RawT.LVESV = round(RawT.LVESV);
RawT.PASP = round(RawT.PASP);
RawT.TVr(RawT.TVr<1.5) = 0;
RawT.TVr(RawT.TVr>1.5) = 1;

%% Simulation for Scott Pacing Patient
try
    MRI_flag = 0;
    [targets, inputs, mods] = targetVals_Pacing();
    [INIparams, INIinit] = estiminiParams(targets,inputs);
    RUNOPT = 0; % 0 for simulation and 1 for optimazation
    if RUNOPT ==1
        P0opt; % optimize modifiers of HF patients, Current HFopt is still used for full model
    end
    if RUNOPT == 0
        % m = zeros(1,length(mods));
        load P0m.mat
        [params, init] = optParams(INIparams, INIinit, mods,m);
        runSimOnGL;
    end
    PlotFakeDTSim;
catch ME
    disp(['Error: ', ME.message])
end

%%
clear
RawT = readtable("inputPredictors.csv",'VariableNamingRule','preserve');
OT = readtable('VirtualParams.csv');

load pacingSUMVirtualMP.mat

[~,idx] = ismember(SummaryTable.PatientNo,(1:height(OT)));
RawT.Sex = OT.Sex(idx);
RawT.Height = OT.Height(idx);
RawT.Weight = OT.Weight(idx);

HRList = zeros(height(SummaryTable),1);
for i = 1:length(HRList)
    HRList(i) = SummaryTable.PacingResults{i}.HR(1);
end
MRI_flag = 1;
desiredFields = RawT.Properties.VariableNames; % 只填已有的列
%%
for order = 1:length(HRList)
    RawT.Patient_No(order) = SummaryTable.PatientNo(order);
    try
        % 加载当前病人的.mat文件
        thisHR = HRList(order);
        foldername = sprintf('PacingSimsVirtual/PatientNO%d', SummaryTable.PatientNo(order));
        Newfilename = sprintf('%s/HR%d.mat', foldername, thisHR);
        S = load(Newfilename);
        output = S.output;
        % 按结构体逐个查找
        structNames = {'params', 'init', 'o_vals'};
        for s = 1:length(structNames)
            structName = structNames{s};
            if isfield(output, structName)
                thisStruct = output.(structName);
                fields = fieldnames(thisStruct);
                for f = 1:length(fields)
                    varName = fields{f};
                    if ismember(varName, desiredFields)
                        RawT{order, varName} = thisStruct.(varName);
                    end
                end
            end
        end

    catch ME
        fprintf('病人 %d 加载失败：%s\n', SummaryTable.PatientNo(order), ME.message);
    end
end

%%
RawT.Sex = 1-RawT.Sex;
RawT.Height = round(RawT.Height,1);
RawT.Weight = round(RawT.Weight,1);
RawT.SBP = round(RawT.SBP);
RawT.DBP = round(RawT.DBP);
RawT.EF = round(RawT.EF,1);
RawT.Hed_SW = round(RawT.Hed_SW*10,1);
RawT.Hed_LW = round(RawT.Hed_LW*10,1);
RawT.HR = round(RawT.HR,1);
RawT.LAVmax = round(RawT.LAVmax);
RawT.EAr = round(RawT.EAr,1);
RawT.IVRT = round(RawT.IVRT*1000);
RawT.LVIDd = round(RawT.LVIDd*10,1);
RawT.LVIDs = round(RawT.LVIDs*10,1);
RawT.LVEDV = round(RawT.LVEDV);
RawT.LVESV = round(RawT.LVESV);
RawT.PASP = round(RawT.PASP);
RawT.TVr(RawT.TVr<1.5) = 0;
RawT.TVr(RawT.TVr>1.5) = 1;

%%
for order = 1:length(HRList)
   RawT.LAP_B(order) = round(SummaryTable.PacingResults{order}.LAP(1),4);
   RawT.LAP_P(order) = round(SummaryTable.PacingResults{order}.LAP(2),4);
   RawT.EF_B(order) = round(SummaryTable.PacingResults{order}.EF(1),4);
   RawT.EF_P(order) = round(SummaryTable.PacingResults{order}.EF(2),4);
   RawT.CO_B(order) = round(SummaryTable.PacingResults{order}.CO(1),4);
   RawT.CO_P(order) = round(SummaryTable.PacingResults{order}.CO(2),4);
   RawT.LAV_B(order) = round(SummaryTable.PacingResults{order}.LAV(1),4);
   RawT.LAV_P(order) = round(SummaryTable.PacingResults{order}.LAV(2),4);
   RawT.LVO2E_B(order) = round(SummaryTable.PacingResults{order}.LVO2E_all(1),4);
   RawT.LVO2E_P(order) = round(SummaryTable.PacingResults{order}.LVO2E_all(2),4);
   RawT.RVO2E_B(order) = round(SummaryTable.PacingResults{order}.RVO2E_all(1),4);
   RawT.RVO2E_P(order) = round(SummaryTable.PacingResults{order}.RVO2E_all(2),4);
   RawT.LVMVO2_B(order) = round(SummaryTable.PacingResults{order}.MVO2_LV(1),4);
   RawT.LVMVO2_P(order) = round(SummaryTable.PacingResults{order}.MVO2_LV(2),4);
   RawT.RVMVO2_B(order) = round(SummaryTable.PacingResults{order}.MVO2_RV(1),4);
   RawT.RVMVO2_P(order) = round(SummaryTable.PacingResults{order}.MVO2_RV(2),4);
   RawT.LVME_B(order) = round(SummaryTable.PacingResults{order}.CE_LV_all(1),6);
   RawT.LVME_P(order) = round(SummaryTable.PacingResults{order}.CE_LV_all(2),6);
   RawT.RVME_B(order) = round(SummaryTable.PacingResults{order}.CE_RV_all(1),6);
   RawT.RVME_P(order) = round(SummaryTable.PacingResults{order}.CE_RV_all(2),6);
end

%%
writetable(RawT,'inputPredictors.csv')
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
