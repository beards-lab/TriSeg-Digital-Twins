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
for GENDER  = 1 % 1 for male and other for female
    if GENDER == 1
        [targets, inputs, mods] = targetVals_male();
        modifiers = readmatrix('modifiers_male.csv','range','2:2');
        % modifiers = ones(1, length(mods));
    else
        [targets, inputs, mods] = targetVals_female();
        modifiers = readmatrix('modifiers_female.csv','range','2:2');
        % modifiers = ones(1, length(mods));
    end
    Sim = 1; % 1 for simulation and other for optimazation
    if Sim ==1
        [INIparams, INIinit] = estiminiParams(targets,inputs);
        [params, init] = optParams(INIparams, INIinit, mods,modifiers);
        runSim;
        NplotSrd; % 4-panel figure just for canonical subjects
        % GetMovie; % TriSeg model: displacement and stress as functions of time
        % See_TriSeg; % slices of GetMoive.m
    else
        Srdopt; % optimize modifiers of canonical subjects
    end
end

%% Simulation for heart failure patients
clear
load AllPatients.mat
secondslot = [17 79 206 256 288 325 352 355 360 361]; % patients with 2 3-month windows
% Shuntlist = [34 41 54 61 83 116 183 231 268 278 312]; % patients with shunt
for PatID = 35 % any number between 1 and 370, example patient in paper is 192
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
KactRatio = NaN(64,1);
KpasRatio = NaN(64,1);
MAP = NaN(64,1);
mPAP = NaN(64,1);
RAPmax = NaN(64,1);
PCWP = NaN(64,1);
StressLV = NaN(64,1);
StressSEP = NaN(64,1);
StressRV = NaN(64,1);
StressactLV = NaN(64,1);
StressactSEP = NaN(64,1);
StressactRV = NaN(64,1);
StresspasLV = NaN(64,1);
StresspasSEP = NaN(64,1);
StresspasRV = NaN(64,1);
VwLV = NaN(64,1);
VwRV = NaN(64,1);
VwSEP = NaN(64,1);
predictedRVEDV = NaN(64,1);
predictedRVESV = NaN(64,1);
predictedTRW = NaN(64,1);
predictedTLW = NaN(64,1);
predictedTSW = NaN(64,1);
%%
tic
for PatID =4
    % try
        MRI_flag = 0;
        [targets, inputs, mods] = targetVals_UW(UWpatients, PatID, MRI_flag);
        [INIparams, INIinit] = estiminiParams(targets,inputs);
        RUNOPT = 1; % 0 for simulation and 1 for optimazation
        if RUNOPT ==1
            UWopt; % optimize modifiers of HF patients
        end
        if RUNOPT == 0
            if MRI_flag == 1
            load(sprintf('SimsUWwithCMRLVFromTTE/P_NO%d.mat',PatID)); % load results from great lake clusters
            else
            load(sprintf('SimsUWwithoutCMR0/P_NO%d.mat',PatID)); % load results from great lake clusters
            end
            modifiers = output.m;
            % modifiers  = modifiers(1:end-1); 
            % modifiers(1) = 1.1418;
            % modifiers(2) = 1;
            % modifiers(3) = 9.2677;
            % modifiers(4) = 1.5381;
            % modifiers(5) = 1.1202;
            % modifiers(6) = 1.0108;
            % modifiers(7) = 0.9252;
            % modifiers(9) = 3.9428;
            % modifiers(10) = 5.5879;
            % modifiers(11) = 2.5829;
            % modifiers(12) = 0.5925;
            % modifiers(8) = 0.7;
            % modifiers(13) = 1.4879;
            % modifiers(14) = 0.6006;
            % modifiers(15) = 0.3965;
            % modifiers(16) = 4.1595;
            % modifiers(17) = 0.0239;
            % modifiers(18) = 0.4994;
            % modifiers(19) = 1;
            % modifiers(9) = 10;
            % modifiers = [modifiers modifiers(end)];
            % modifiers = ones(1,length(mods))*1;
            % modifiers(14) = 0.8;
            [params, init] = optParams(INIparams, INIinit, mods,modifiers,targets);
            % % modifiers(9) = 2;
            % modifiers(16) = 1.1;
            % modifiers(11) = 2;
            % modifiers = ones(1,length(mods));
            % baselineHR = params.HR;
            % NewHR = baselineHR+(0:2:40);
            % params.T = 60/params.HR;
            runSim;
            Print_cost = 1;% 1 for performance, other for patient's pre-condition (PHI sensitive)
            NplotFit; % 6-panel figure for HF patients
            % GetMovie;
            % See_TriSeg;
            % pause(3);
            KactRatio(PatID) = params.k_act_LV/params.k_act_RV;
            KpasRatio(PatID) = params.k_pas_LV/params.k_pas_RV;
            VwLV(PatID)  = params.Vw_LV;
            VwRV(PatID)  = params.Vw_RV;
            VwSEP(PatID)  = params.Vw_SEP;
            MAP(PatID) = targets.SBP;
            mPAP(PatID) = targets.PASP;
            RAPmax(PatID) = targets.RAPmax;
            PCWP(PatID) = targets.PCWP;
            StressLV(PatID) = max(sigma_LV);
            StressSEP(PatID) =max(sigma_SEP);
            StressRV(PatID) = max(sigma_RV);
            StressactLV(PatID) = max(sigma_act_LV);
            StressactSEP(PatID) =max(sigma_act_SEP);
            StressactRV(PatID) = max(sigma_act_RV);
            StresspasLV(PatID) = sigma_pas_LV(1);
            StresspasSEP(PatID) =sigma_pas_SEP(1);
            StresspasRV(PatID) = sigma_pas_RV(1);
            predictedRVEDV(PatID) = o_vals.RVEDV;
            predictedRVESV(PatID) = o_vals.RVESV;
            predictedTRW(PatID) = o_vals.Hed_RW;
            predictedTLW(PatID) = o_vals.Hed_LW;
            predictedTSW(PatID) = o_vals.Hed_SW;

        end
    % catch ME
    %     disp(['Error: ', ME.message, 'PatID',num2str(PatID)])
    % end
end
toc
% save 01272025backup.mat StresspasRV StresspasSEP StresspasLV StressactRV StressactLV StressactSEP PCWP RAPmax mPAP KpasRatio KactRatio MAP StressLV StressSEP StressRV -mat
%% 
fields = fieldnames(realparams);
result = struct();

for i = 1:length(fields)
    field = fields{i};
    result.(field) = realparams.(field) / INIparams.(field);
end
%%
CombinedRW = [UWpatients.Hed_RW;predictedTRW];
StdRW = zscore(CombinedRW);
StdrealRW = StdRW(1:64);
StdpredRW = StdRW(65:end);

CombinedRVEDV = [UWpatients.RVEDV;predictedRVEDV];
StdRVEDV = zscore(CombinedRVEDV);
StdrealRVEDV = StdRVEDV(1:64);
StdpredRVEDV = StdRVEDV(65:end);

CombinedRVESV = [UWpatients.RVESV;predictedRVESV];
StdRVESV = zscore(CombinedRVESV);
StdrealRVESV = StdRVESV(1:64);
StdpredRVESV = StdRVESV(65:end);

cost = (StdrealRW-StdpredRW).^2+...
    ...(StdrealRVESV-StdpredRVESV).^2+...
    (StdrealRVEDV-StdpredRVEDV).^2;

%%
X = UWpatients.PPAm;
Y = UWpatients.EF_RV;

% 线性回归拟合
mdl = fitlm(X, Y);
disp(mdl);
% 计算标准化残差
r = mdl.Residuals.Studentized;
threshold = 5; % 一般使用 ±2 作为标准化残差阈值

% 识别离群点
outliers = abs(r) > threshold;
% 过滤掉离群点
X3_filtered = X(~outliers, :);
Y3_filtered = Y(~outliers);

% 重新拟合回归模型
mdl_filtered = fitlm(X3_filtered, Y3_filtered,'linear','RobustOpts','on');
disp(mdl_filtered);
% 绘制原始数据和离群点
figure(66); clf;
scatter(X(~outliers), Y(~outliers),'SizeData',100,'MarkerFaceColor',[0, 0.5, 0.7]); hold on;
% xlim([0.2 1]);
% ylim([0.2 1]);
pbaspect([1,1,1]);
box on;

% xticks(0:0.3:0.9);
% yticks(0:0.3:0.9);


% 添加 Y = X 的虚线
% plot([0, ylim], [0, ylim], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
% 计算 Adjusted R
adj_r2 = mdl_filtered.Rsquared.Adjusted; % Adjusted R²
beta = mdl_filtered.Coefficients.Estimate(2); % 斜率
adj_r = sign(beta) * sqrt(adj_r2); % 计算 Adjusted R
p_value = mdl_filtered.Coefficients.pValue(2); % 斜率的 p 值

% 在图上显示 Adjusted R 和 p 值
annotation('textbox', [0.2, 0.8, 0.5, 0.1], ... % 位置：[x, y, width, height]
    'String', sprintf('Adjusted R = %.2f\np = %.3g', adj_r, p_value), ...
    'FontSize', 16, 'FontName', 'Arial', ...
    'EdgeColor', 'none', ... % 取消边框
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
% 添加坐标轴标签
xlabel('Measured RV Thickness (cm)', 'FontSize', 10, 'FontName', 'Arial');
ylabel('Predicted RV Thickness (cm)', 'FontSize', 10, 'FontName', 'Arial');
scatter(X(outliers), Y(outliers), 'ro', 'filled'); % 标出离群点
plot(mdl_filtered);
title(' ');
legend off;
% 置信区间（Confidence Interval）和预测区间（Prediction Interval）
X_pred = linspace(min(X), max(X), 100)'; % 生成一系列 X 值
[~, Y_CI] = predict(mdl_filtered, X_pred,'alpha', 5e-11);   % 置信区间
[Y_pred, Y_PI] = predict(mdl_filtered, X_pred, 'Prediction', 'observation','alpha',0.01); % 预测区间

% 绘制拟合曲线、置信区间和预测区间
figure; hold on;
scatter(X, Y, 'bo'); % 绘制数据点
plot(X_pred, Y_pred, 'r-', 'LineWidth', 2); % 绘制拟合线
fill([X_pred; flipud(X_pred)], [Y_CI(:,1); flipud(Y_CI(:,2))], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % 置信区间
fill([X_pred; flipud(X_pred)], [Y_PI(:,1); flipud(Y_PI(:,2))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % 预测区间

% legend('Data', 'Regression Line', 'Confidence Interval', 'Prediction Interval');
% title('Linear Regression with Confidence and Prediction Intervals');
hold off;
% 设置字体和字号
% set(gca, 'FontName', 'Arial', 'FontSize', 26);
% 设置图像大小（可选，调整为7英寸宽）
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 7 7]);

% % 保存为 PNG（330 dpi）
%print('TRW_Scatter.png', '-dpng', '-r330');
%%
figure(33); hold on;
scatter(StresspasRV, StresspasLV, 50,'filled'); % 50 是点的大小
colormap(jet);  % 选择 colormap，例如 'jet' 或 'parula'
colorbar;       % 添加颜色条，显示 cost 对应的颜色
clim([min(cost), max(cost)]); % 设定颜色范围
% scatter(Y2_filtered, X2_filtered, 50); % 50 是点的大小
figure(44); hold on;
scatter(MAP ./ mPAP, (StressactLV + StressactSEP) ./ StressactRV, 50,'filled'); 
colormap(jet);  % 颜色映射方案（低值蓝，高值红）
colorbar;       % 显示颜色条
clim([min(cost), max(cost)]); % 颜色范围与 cost 匹配
% scatter(X1_filtered, Y1_filtered, 50);

%%
X = [UWpatients.RVESV];
Y = [predictedRVESV];
% 计算均值和差值
meanXY = (X + Y) / 2;  
diffXY = X - Y;  

% 计算一致性界限（LoA）
meanDiff = mean(diffXY);  % 差值的均值
stdDiff = std(diffXY);    % 差值的标准差
LoA_upper = meanDiff + 1.96 * stdDiff;  % 上限
LoA_lower = meanDiff - 1.96 * stdDiff;  % 下限

% 绘制 Bland-Altman Plot
figure;
scatter(meanXY, diffXY, 'bo');  % 绘制散点图
hold on;
yline(meanDiff, 'r-', 'Mean Diff');  % 画出均值线
yline(LoA_upper, 'g--', 'Upper LoA');  % 画出上界限
yline(LoA_lower, 'g--', 'Lower LoA');  % 画出下界限

% 添加图例和标签
xlabel('Mean of X and Y');
ylabel('Difference (X - Y)');
title('Bland-Altman Plot');
grid on;
hold off;

%%
% load Kconstrain.mat
% save Kconstrain.mat X2_filtered X1_filtered Y1_filtered Y2_filtered X3_filtered Y3_filtered
%% extract modifiers from cluster
% clear
for PatID = 1:64
    if PatID < 10
        matFilesDir = sprintf('Runs/p00%d/02-11/',PatID);
    elseif PatID < 100
        matFilesDir = sprintf('Runs/p0%d/02-11/',PatID);
    else
        matFilesDir = sprintf('Runs/p%d/02-11/',PatID);
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
validIdx_RVEDV =  ~isnan(mPAP);
[r_RVEDV, p_RVEDV] = corr(UWpatients.RVEDV(validIdx_RVEDV), mPAP(validIdx_RVEDV), 'Type', 'Pearson');

% 去除 NaN 对于 RVESV vs PredictedRVESV
validIdx_RVESV = ~isnan(RAPmax);
[r_RVESV, p_RVESV] = corr(UWpatients.RVESV(validIdx_RVESV), RAPmax(validIdx_RVESV), 'Type', 'Pearson');

% 去除 NaN 对于 Hed_RW vs TRW
validIdx_RW = ~isnan(KactRatio);
[r_RW, p_RW] = corr(UWpatients.Hed_RW(validIdx_RW), KactRatio(validIdx_RW), 'Type', 'Pearson');

% 显示结果
disp(['RVEDV vs PredictedRVEDV: r = ', num2str(r_RVEDV), ', p = ', num2str(p_RVEDV)]);
disp(['RVESV vs PredictedRVESV: r = ', num2str(r_RVESV), ', p = ', num2str(p_RVESV)]);
disp(['Hed_RW vs TRW: r = ', num2str(r_RW), ', p = ', num2str(p_RW)]);

% 可选：在散点图中显示 r 和 p 值
figure(66);
scatter(UWpatients.RVEDV(validIdx_RVEDV), mPAP(validIdx_RVEDV));
title(['RVEDV: r = ', num2str(r_RVEDV), ', p = ', num2str(p_RVEDV)]);
xlabel('RVEDV'); ylabel('PredictedRVEDV');
box on;
pbaspect([1,1,1]);
% 获取 X 和 Y 数据的范围
xData = UWpatients.RVEDV(validIdx_RVEDV);
yData = mPAP(validIdx_RVEDV);
minVal = min([xData; yData]); % 取最小值
maxVal = max([xData; yData]); % 取最大值

% 设置相同的 X 和 Y 坐标范围
axis([minVal maxVal minVal maxVal]);

% 添加 x = y 辅助线
hold on;
plot([minVal maxVal], [minVal maxVal], 'k--', 'LineWidth', 1);
hold off;


figure(77); scatter(UWpatients.RVESV(validIdx_RVESV), RAPmax(validIdx_RVESV));
title(['RVESV: r = ', num2str(r_RVESV), ', p = ', num2str(p_RVESV)]);
xlabel('RVESV'); ylabel('PredictedRVESV');
box on;
pbaspect([1,1,1]);
% 获取 X 和 Y 数据的范围
xData = UWpatients.RVESV(validIdx_RVESV);
yData = RAPmax(validIdx_RVESV);
minVal = min([xData; yData]); % 取最小值
maxVal = max([xData; yData]); % 取最大值

% 设置相同的 X 和 Y 坐标范围
axis([minVal maxVal minVal maxVal]);

% 添加 x = y 辅助线
hold on;
plot([minVal maxVal], [minVal maxVal], 'k--', 'LineWidth', 1);
hold off;
figure(88); scatter(UWpatients.Hed_RW(validIdx_RW), KactRatio(validIdx_RW));
title(['Hed_RW: r = ', num2str(r_RW), ', p = ', num2str(p_RW)]);
xlabel('Hed_RW'); ylabel('TRW');
box on;
pbaspect([1,1,1]);
% 获取 X 和 Y 数据的范围
xData = UWpatients.Hed_RW(validIdx_RVEDV);
yData = KactRatio(validIdx_RVEDV);
minVal = min([xData; yData]); % 取最小值
maxVal = max([xData; yData]); % 取最大值

% 设置相同的 X 和 Y 坐标范围
axis([minVal maxVal minVal maxVal]);

% 添加 x = y 辅助线
hold on;
plot([minVal maxVal], [minVal maxVal], 'k--', 'LineWidth', 1);
hold off;
t1 = [mPAP UWpatients.RVEDV];
DI1 = mPAP-UWpatients.RVEDV;
figure;
histogram(DI1, 'Normalization', 'probability'); % 归一化为频率
xlabel('DI1 Values');
ylabel('Frequency');
title('Frequency Distribution of DI1');

%%
% 去除 NaN 对于 LVEDV vs PredictedLVEDV
CO = 1000 * UWpatients.CO./ 60;
SV = 60 * CO ./ UWpatients.HR;
KpasRatio =   SV./ (UWpatients.EF.*0.01);
MAP =   KpasRatio-SV;
% PredictedLVEDV = UWpatients.LVIDd.^3*0.7851+97.32;
% PredictedLVESV = UWpatients.LVIDs.^3*0.9185+63.92;
validIdx_LVEDV = ~isnan(KpasRatio);
[r_LVEDV, p_LVEDV] = corr(UWpatients.LVEDV(validIdx_LVEDV), KpasRatio(validIdx_LVEDV), 'Type', 'Pearson');

% 去除 NaN 对于 LVESV vs PredictedLVESV
validIdx_LVESV = ~isnan(MAP);
[r_LVESV, p_LVESV] = corr(UWpatients.LVESV(validIdx_LVESV), MAP(validIdx_LVESV), 'Type', 'Pearson');

% 显示结果
disp(['LVEDV vs PredictedLVEDV: r = ', num2str(r_LVEDV), ', p = ', num2str(p_LVEDV)]);
disp(['LVESV vs PredictedLVESV: r = ', num2str(r_LVESV), ', p = ', num2str(p_LVESV)]);

% 可选：在散点图中显示 r 和 p 值
figure(66); scatter(UWpatients.LVEDV(validIdx_LVEDV), KpasRatio(validIdx_LVEDV));
title(['LVEDV: r = ', num2str(r_LVEDV), ', p = ', num2str(p_LVEDV)]);
xlabel('LVEDV'); ylabel('PredictedLVEDV');

figure(77); scatter(UWpatients.LVESV(validIdx_LVESV), MAP(validIdx_LVESV));
title(['LVESV: r = ', num2str(r_LVESV), ', p = ', num2str(p_LVESV)]);
xlabel('LVESV'); ylabel('PredictedLVESV');


% 计算和显示 PredictedLVEDV 与 LVEDV 的差值分布
DI1 = KpasRatio - UWpatients.LVEDV;
figure;
histogram(DI1, 'Normalization', 'probability'); % 归一化为频率
xlabel('DI1 Values');
ylabel('Frequency');
title('Frequency Distribution of DI1');

T = [UWpatients.Hed_RW UWpatients.RVEDV UWpatients.RVESV UWpatients.LVEDV UWpatients.LVESV];
