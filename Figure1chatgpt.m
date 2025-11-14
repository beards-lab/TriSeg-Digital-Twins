%% ------------------------------------------------------------
%  Cleaned plotting pipeline (all exports 900 dpi; LAPCO -> LVO2ECOclear)
%  Author: Feng-ready compact version
% ------------------------------------------------------------

clear; clc;

%% ---------- Global style & paths ----------
EXPORT_DPI = 900;
FONT_NAME  = 'Arial';
FONT_SIZE  = 16;

w_small = 1.2; h_small = 1.2;  aspect_sq   = [1 1 1];
w_wide  = 2.4; h_wide  = 1.2;  aspect_wide = [2 1 1];
aspect_golden = [1.618 1 1];

baseOut = fullfile('C:', 'Users', 'fenggu', ...
    'University of Michigan Dropbox', 'Feng Gu', 'DTPacingHFpEF');

% colormap once
n = 7;  % pacing 7 HR points
start_color = [0.6350, 0.0780, 0.1840];  % dark red
end_color   = [0.5,   0.5,    0.5   ];  % gray
cmap = [linspace(start_color(1), end_color(1), n)', ...
        linspace(start_color(2), end_color(2), n)', ...
        linspace(start_color(3), end_color(3), n)'];

%% ============================================================
% Normal baseline
% ============================================================
MRI_flag = 1; %#ok<NASGU>
GENDER = 2;  % 1=male, else=female
[~, inputs, ~] = targetVals_female();
[targets, params, init] = load_gender_params(GENDER);
Sim = 1;
if Sim == 1
    runSimOnGL;
    Tnormal      = t;
    P_LV_normal  = P_LV;
    V_LV_normal  = V_LV;
    P_LA_normal  = P_LA;
    P_SA_normal  = P_SA;
else
    Srdopt;
end
TNnormal = Tnormal ./ Tnormal(round(end/2));

% LAP trace (normal)
figure(77); clf;
plot(TNnormal, P_LA_normal, 'color', [0.6350 0.0780 0.1840], 'LineWidth', 2); hold on;
plot(TNnormal, P_LV_normal, 'color', [0.3010 0.7450 0.9330], 'LineWidth', 2);
trng = [0.2 1.2];
plot(trng, [targets.SBP,  targets.SBP ], '--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
plot(trng, [targets.DBP,  targets.DBP ], '--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
plot(trng, [targets.PCWP, targets.PCWP], '--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
xlim(trng); ylim([0 20]);
style_singleY_axes(gca, FONT_NAME, FONT_SIZE);
save_fig(gcf, fullfile(baseOut,'NormalLAP.png'), w_wide, h_wide, aspect_wide, EXPORT_DPI);

% SP trace (normal)
figure(777); clf;
plot(TNnormal, P_LV_normal, 'color', [0.3010 0.7450 0.9330], 'LineWidth', 2); hold on;
plot(TNnormal, P_SA_normal, 'color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
trng = [0.8 1.8];
plot(trng, [targets.SBP,  targets.SBP ], '--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
plot(trng, [targets.DBP,  targets.DBP ], '--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
plot(trng, [targets.PCWP, targets.PCWP], '--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
xlim(trng); ylim([50 130]);
style_singleY_axes(gca, FONT_NAME, FONT_SIZE);
save_fig(gcf, fullfile(baseOut,'NormalSP.png'), w_wide, h_wide, aspect_wide, EXPORT_DPI);

%% ============================================================
% Normal pacing
% ============================================================
MRI_flag = 1; %#ok<NASGU>
P_LV_cell = cell(1,n);
V_LV_cell = cell(1,n);
P_LA_cell = cell(1,n);
P_SA_cell = cell(1,n);
T_cell    = cell(1,n);

LAPNormal     = NaN(1,n);
COnormal      = NaN(1,n);
LVO2ENormal   = zeros(1,n);
LVMVO2Normal  = zeros(1,n);
EDPLV_all     = zeros(1,n);
EDVLV_all     = zeros(1,n);
EDPRV_all     = zeros(1,n);
EDVRV_all     = zeros(1,n);

[targets, params, init] = load_gender_params(GENDER); %#ok<ASGLU>
NewHR = params.HR + (0:5:30);
fix_V0LV = 0; fix_V0RV = 0; V0LV_fixed = NaN; V0RV_fixed = NaN;

for iHR = 1:n
    if iHR > 1
        foldername = sprintf('PacingSimsC/PatientNO%d', GENDER);
        Newfilename = sprintf('%s/HR%d.mat', foldername, NewHR(iHR));
        S = load(Newfilename); params = S.output.params; init = S.output.init;
    end
    runSimOnGL;
    T_cell{iHR}    = t;
    P_LV_cell{iHR} = P_LV;  V_LV_cell{iHR} = V_LV;
    P_LA_cell{iHR} = P_LA;  P_SA_cell{iHR} = P_SA;

    COnormal(iHR)  = o_vals.CO;
    LAPNormal(iHR) = P_LA(1);
    EDPLV_all(iHR) = o_vals.LVEDP;   EDVLV_all(iHR) = V_LV(1);
    EDPRV_all(iHR) = o_vals.RVEDP;   EDVRV_all(iHR) = V_RV(1);
end

% First PVA call to get bounds
coeffRV = NaN; coeffLV = NaN;
[~, ~, V0LV, V0RV, ~, ~, ~, B_sol_ub, ~, B_sol_rv_ub] = ...
    Calculate_PVA(init, params, fix_V0LV, fix_V0RV, V0LV_fixed, V0RV_fixed, coeffLV, coeffRV);

% Re-init + fix V0
[~, params, init] = load_gender_params(GENDER);
fix_V0LV = 1; fix_V0RV = 1; V0LV_fixed = V0LV; V0RV_fixed = V0RV;
[~,~,~,~,~,~,~,B_sol_lb,~,B_sol_rv_lb] = ...
    Calculate_PVA(init, params, fix_V0LV, fix_V0RV, V0LV_fixed, V0RV_fixed, coeffLV, coeffRV);

for iHR = 1:n
    if iHR > 1
        foldername = sprintf('PacingSimsC/PatientNO%d', GENDER);
        Newfilename = sprintf('%s/HR%d.mat', foldername, NewHR(iHR));
        S = load(Newfilename); params = S.output.params; init = S.output.init;
    end
    % LV EDPVR coefficients via geometric interpolation
    LVEDP = EDPLV_all(iHR);
    LVEDV = EDVLV_all(iHR);
    t_lv  = clamp01((LVEDV - EDVLV_all(end)) / (EDVLV_all(1) - EDVLV_all(end)));
    b2_lv = (B_sol_ub)^(1 - t_lv) * (B_sol_lb)^(t_lv);
    den_lv = exp(b2_lv * (LVEDV - V0LV)) - 1;
    b1_lv  = LVEDP / max(den_lv, eps);
    coeffLV = [b1_lv, b2_lv, V0LV];

    % RV
    RVEDP = EDPRV_all(iHR);
    RVEDV = EDVRV_all(iHR);
    t_rv  = clamp01((RVEDV - EDVRV_all(end)) / (EDVRV_all(1) - EDVRV_all(end)));
    b2_rv = (B_sol_rv_ub)^(1 - t_rv) * (B_sol_rv_lb)^(t_rv);
    den_rv = exp(b2_rv * (RVEDV - V0RV)) - 1;
    b1_rv  = RVEDP / max(den_rv, eps);
    coeffRV = [b1_rv, b2_rv, V0RV];

    [PVALV, PVARV, ~, ~,PVarea_LV,PVarea_RV,~,~,~,~] = ...
        Calculate_PVA(init, params, fix_V0LV, fix_V0RV, V0LV_fixed, V0RV_fixed, coeffLV, coeffRV);

    % Energetics
    MVO2_LV = (1.56e-5 * PVALV + 0.00526 * o_vals.LV_m / 100) * params.HR;
    Energy_LV = 1.333e-4 * PVarea_LV * params.HR;
    CardiacEff_LV =  Energy_LV / (MVO2_LV * (20.1+19.6)/2); %#ok<NASGU>
    MVO2_LV_100g = MVO2_LV / o_vals.LV_m * 100;

    MVO2_RV = (1.56e-5 * PVARV + 0.000526 * o_vals.RV_m / 100) * params.HR;
    Energy_RV = 1.333e-4 * PVarea_RV * params.HR;
    CardiacEff_RV =  Energy_RV / (MVO2_RV * (20.1+19.6)/2); %#ok<NASGU>
    MVO2_RV_100g = MVO2_RV / o_vals.RV_m * 100; %#ok<NASGU>

    LVMVO2Normal(iHR) = MVO2_LV_100g;
    LVO2ENormal(iHR)  = o_vals.CO*1000 / MVO2_LV;
end
NewHRNormal = NewHR; %#ok<NASGU>

% Pacing LAP (normal)
figure(105); clf;
plot_pacing_traces(T_cell, P_LA_cell, cmap, [0.2 1.2], [0 20]);
style_singleY_axes(gca, FONT_NAME, FONT_SIZE);
save_fig(gcf, fullfile(baseOut,'NormalLAPPacing.png'), w_wide, h_wide, aspect_golden, EXPORT_DPI);

% Pacing SP (normal, aligned by systolic peaks)
figure(106); clf;
plot_pacing_SA_aligned(T_cell, P_SA_cell, cmap, [0.8 1.8], [50 130]);
style_singleY_axes(gca, FONT_NAME, FONT_SIZE);
save_fig(gcf, fullfile(baseOut,'NormalSPPacing.png'), w_wide, h_wide, aspect_golden, EXPORT_DPI);

%% PV loop in normal (with EDPVR & ESPVR overlays)
figure(202); clf;
plot_pv_with_EDPVR_ES(P_LV_cell, V_LV_cell,cmap, [0 150], [0 150]);
style_singleY_axes(gca, FONT_NAME, FONT_SIZE);
save_fig(gcf, fullfile(baseOut,'PVloopPacingNormal.png'), w_small, h_small, aspect_sq, EXPORT_DPI);

%% Normal: Δ(LVO2E) & ΔCO vs HR
figure(201); clf;
HR = NewHRNormal - NewHRNormal(1);
LVCE = LVO2ENormal - LVO2ENormal(1);
CO   = COnormal   - COnormal(1);
yyaxis left;  plot(HR, LVCE, '-o', 'LineWidth', 2, 'Color', [0.8500 0.5250 0.0980]); ylim([-50 50]);
ax = gca; ax.YColor = [0.8500 0.5250 0.0980];
yyaxis right; plot(HR, CO,   '-o', 'LineWidth', 2, 'Color', [0 0 0]);    ylim([-0.4 0.4  ]);
style_axes(gca, FONT_NAME, FONT_SIZE,'none');
save_fig(gcf, fullfile(baseOut,'PacingLVO2ECONormal.png'), w_small, h_small, aspect_sq, EXPORT_DPI);

% Normal: Δ(LV MVO2_100g)
figure(211); clf;
dM = LVMVO2Normal - LVMVO2Normal(1);
plot(HR, dM, '-o', 'LineWidth', 2, 'Color', 'k'); ylim([-1 1]);
style_singleY_axes(gca, FONT_NAME, FONT_SIZE);
save_fig(gcf, fullfile(baseOut,'PacingLVMVO2Normal.png'), w_small, h_small, aspect_sq, EXPORT_DPI);

%% ============================================================
% HFpEF baseline
% ============================================================
% NO = 19;
% S = load(sprintf('SimsUWFinal/P_NO%d.mat', NO));
NO = 274;
S = load(sprintf('SimsUMFinal/P_NO%d.mat', NO));
inputs = S.output.inputs; targets = S.output.targets; params = S.output.params; init = S.output.init; %#ok<NASGU>
% if targets.PCWP<15
%     error('too small pcwp')
% end
runSimOnGL;
THFpEF = t;
P_LV_HFpEF = P_LV; V_LV_HFpEF = V_LV; P_LA_HFpEF = P_LA; P_SA_HFpEF = P_SA;
TNHFpEF = THFpEF./THFpEF(round(end/2));

%% LAP baseline HFpEF
figure(88); clf;
plot(TNHFpEF, P_LA_HFpEF, 'color', [0.6350 0.0780 0.1840], 'LineWidth', 2); hold on;
plot(TNHFpEF, P_LV_HFpEF, 'color', [0.3010 0.7450 0.9330], 'LineWidth', 2);
trng = [0.2 1.2];
plot(trng, [targets.SBP, targets.SBP], '--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
plot(trng, [targets.DBP, targets.DBP], '--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
plot(trng, [targets.PCWP,targets.PCWP],'--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
xlim(trng); ylim([5 25]);
style_singleY_axes(gca, FONT_NAME, FONT_SIZE);
% save_fig(gcf, fullfile(baseOut,'HFpEFLAP.png'), w_wide, h_wide, aspect_wide, EXPORT_DPI);
% 
% SP baseline HFpEF
figure(888); clf;
plot(TNHFpEF, P_LV_HFpEF, 'color', [0.3010 0.7450 0.9330], 'LineWidth', 2); hold on;
plot(TNHFpEF, P_SA_HFpEF, 'color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
trng = [0.8 1.8];
plot(trng, [targets.SBP, targets.SBP], '--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
plot(trng, [targets.DBP, targets.DBP], '--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
plot(trng, [targets.PCWP,targets.PCWP],'--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
xlim(trng); ylim([50 130]);
style_singleY_axes(gca, FONT_NAME, FONT_SIZE);
%% save_fig(gcf, fullfile(baseOut,'HFpEFSP.png'), w_wide, h_wide, aspect_wide, EXPORT_DPI);

% ============================================================
% HFpEF pacing (signals & PV loop)
% ============================================================
MRI_flag = 1; %#ok<NASGU>
P_LV_H_cell = cell(1,n);
V_LV_H_cell = cell(1,n);
P_LA_H_cell = cell(1,n);
P_SA_H_cell = cell(1,n);
T_H_cell    = cell(1,n);

S = load(sprintf('SimsUMFinal/P_NO%d.mat', NO));
targets = S.output.targets; params = S.output.params; init = S.output.init; %#ok<NASGU>
NewHR = params.HR + (0:5:30);
for iHR = 1:n
    if iHR > 1
        foldername = sprintf('PacingSimsUM/PatientNO%d', NO);
        Newfilename = sprintf('%s/HR%d.mat', foldername, NewHR(iHR));
        S2 = load(Newfilename);
        targets = S2.output.targets; params = S2.output.params; init = S2.output.init; %#ok<NASGU>
    end
    runSimOnGL;
    T_H_cell{iHR}    = t;
    P_LV_H_cell{iHR} = P_LV;
    V_LV_H_cell{iHR} = V_LV;
    P_LA_H_cell{iHR} = P_LA;
    P_SA_H_cell{iHR} = P_SA;
end
%%
% 汇总表 & 示例病人
T1 = load('pacingSUM.mat');
T2 = load('pacingSUMUW.mat');
SummaryTable = [T1.SummaryTable; T2.SummaryTable];
HFpEFexample = SummaryTable(134,:);
% HFpEFexample = SummaryTable(128,:);

%% Pacing LAP (HFpEF)
figure(305); clf;
plot_pacing_traces(T_H_cell, P_LA_H_cell, cmap, [0.2 1.2], [5 25]);
style_singleY_axes(gca, FONT_NAME, FONT_SIZE);
save_fig(gcf, fullfile(baseOut,'HFpEFLAPPacing.png'), w_wide, h_wide, aspect_golden, EXPORT_DPI);

%% Pacing SP (HFpEF, aligned)
figure(306); clf;
plot_pacing_SA_aligned(T_H_cell, P_SA_H_cell, cmap, [0.8 1.8], [50 130]);
style_singleY_axes(gca, FONT_NAME, FONT_SIZE);
save_fig(gcf, fullfile(baseOut,'HFpEFSPPacing.png'), w_wide, h_wide, aspect_wide, EXPORT_DPI);

%% HFpEF PV loop (with EDPVR & ESPVR)
figure(302); clf;
V_LV = HFpEFexample.PacingResults{1,1}.PVvolume;
P_LV = HFpEFexample.PacingResults{1,1}.PVpressure;
plot_pv_with_EDPVR_ES(P_LV, V_LV, cmap,[0 230], [0 150]);
style_singleY_axes(gca, FONT_NAME, FONT_SIZE);
save_fig(gcf, fullfile(baseOut,'PVloopPacing.png'), w_small, h_small, aspect_sq, EXPORT_DPI);
%
% HFpEF: Δ(LVO2E) & ΔCO vs HR  —— 文件名改为 LVO2ECOclear
figure(401); clf;
HR = HFpEFexample.PacingResults{1,1}.HR; HR = HR - HR(1);
LVO2E = HFpEFexample.PacingResults{1,1}.LVO2E_all; LVO2E = LVO2E - LVO2E(1);
CO    = HFpEFexample.PacingResults{1,1}.CO;        CO    = CO - CO(1);
yyaxis left;  plot(HR, LVO2E, '-o', 'LineWidth', 2, 'Color', [0.8500 0.5250 0.0980]); ylim([-50 50]);
ax = gca; ax.YColor = [0.8500 0.5250 0.0980];
yyaxis right; plot(HR, CO,    '-o', 'LineWidth', 2, 'Color', [0 0 0]);                 ylim([-0.4 0.4]);
style_axes(gca, FONT_NAME, FONT_SIZE,'none');
save_fig(gcf, fullfile(baseOut,'PacingLVO2ECOclear.png'), w_small, h_small, aspect_sq, EXPORT_DPI);
%
% HFpEF: Δ(LV MVO2)
figure(411); clf;
LVMVO2 = HFpEFexample.PacingResults{1,1}.MVO2_LV;
LVMVO2 = LVMVO2 - LVMVO2(1);
plot(HR, LVMVO2, '-o', 'LineWidth', 2, 'Color', 'k'); ylim([-1 1]);
style_singleY_axes(gca, FONT_NAME, FONT_SIZE);
save_fig(gcf, fullfile(baseOut,'PacingLVMVO2HFpEF.png'), w_small, h_small, aspect_sq, EXPORT_DPI);


%% ==================== Local functions ====================

function [targets, params, init] = load_gender_params(GENDER)
    if GENDER == 1
        [targets, inputs, mods] = targetVals_male();
        modifiers = readmatrix('modifiers_male.csv','range','2:2');
    else
        [targets, inputs, mods] = targetVals_female();
        modifiers = readmatrix('modifiers_female.csv','range','2:2');
    end
    [INIparams, INIinit] = estiminiParams(targets, inputs);
    [params, init] = optParams(INIparams, INIinit, mods, modifiers);
end

function style_singleY_axes(ax, FONT_NAME, FONT_SIZE) 
set(ax, 'FontName', FONT_NAME, 'FontSize', FONT_SIZE); 
set(ax, 'XTickLabel', [], 'YTickLabel', []); 
box(ax, 'on'); 
end

function style_axes(ax, FONT_NAME, FONT_SIZE, mode)
% STYLE_AXES — 统一轴样式；支持双 y 轴
% mode: 'zero' -> 两侧 y 轴都只显示 0
%       'none' -> 两侧 y 轴都不显示刻度
    if nargin < 4, mode = 'zero'; end

    set(ax, 'FontName', FONT_NAME, 'FontSize', FONT_SIZE);
    set(ax, 'XTick', []);                 
    box(ax, 'on');

    % 判断是否存在右侧 y 轴
    hasRight = false;
    try
        hasRight = strcmp(ax.YAxis(2).Visible, 'on');
    catch
        hasRight = false;
    end

    %=== 设置左右轴颜色 ===%
    try
        yyaxis(ax, 'left');
        ax.YColor = [1, 0.8, 0];   % 左侧黄色（RGB，可调整为 [1 1 0] 等）
        if hasRight
            yyaxis(ax, 'right');
            ax.YColor = [0, 0, 0]; % 右侧黑色
            yyaxis(ax, 'left');    % 还原
        end
    catch
        % 若无双轴功能，至少设置左轴颜色
        ax.YColor = [1, 0.8, 0];
    end

    %=== 模式控制 ===%
    switch lower(mode)
        case 'none'
            yyaxis(ax, 'left');  set(ax, 'YTick', [], 'YTickLabel', []);
            if hasRight
                yyaxis(ax, 'right'); set(ax, 'YTick', [], 'YTickLabel', []);
                yyaxis(ax, 'left');
            end
        case 'zero'
            yyaxis(ax, 'left');  set(ax, 'YTick', 0, 'YTickLabel', {'0'});
            if hasRight
                yyaxis(ax, 'right'); set(ax, 'YTick', 0, 'YTickLabel', {'0'});
                yyaxis(ax, 'left');
            end
        otherwise
            error('style_axes: unknown mode "%s"', mode);
    end
end



function save_fig(figH, outpath, w_in, h_in, aspect, dpi)
    set(figH, 'Units','inches','Position',[1 1 w_in h_in]);
    set(figH, 'PaperUnits','inches','PaperPosition',[0 0 w_in h_in]);
    set(figH, 'InvertHardcopy','off');
    pbaspect(aspect);
    exportgraphics(figH, outpath, 'Resolution', dpi, ...
        'ContentType', 'image', 'BackgroundColor', 'white');
end

function plot_pacing_traces(T_cell, Y_cell, cmap, xlimv, ylimv)
    n = numel(T_cell);
    for i = 1:n
        if i == 1 || i == n
            lw = 2;  alphaVal = 1.0;
        else
            lw = 0.88; alphaVal = 0.4;
        end
        t = T_cell{i};
        x = t ./ t(round(end/2));
        y = Y_cell{i};
        plot(x, y, '-', 'Color', [cmap(i,:), alphaVal], 'LineWidth', lw); hold on;
    end
    xlim(xlimv); ylim(ylimv);
end

function plot_pacing_SA_aligned(T_cell, P_SA_cell, cmap, xlimv, ylimv)
    n = numel(T_cell);
    for i = 1:n
        if i == 1 || i == n
            lw = 2;  alphaVal = 1.0;
        else
            lw = 0.88; alphaVal = 0.4;
        end
        t = T_cell{i};
        tn = t ./ t(round(end/2));
        treg = (tn(1):0.01:tn(end));
        x = treg;
        y = interp1(tn, P_SA_cell{i}, treg);
        if i == 1
            [~, locs_ref] = findpeaks(y, "MinPeakProminence", 5);
        else
            [~, locs] = findpeaks(y, "MinPeakProminence", 5);
            change = locs(end) - locs_ref(end);
            if change ~= 0, x = x - treg(change); end
        end
        plot(x, y, '-', 'Color', [cmap(i,:), alphaVal], 'LineWidth', lw); hold on;
    end
    xlim(xlimv); ylim(ylimv);
end

function plot_pv_with_EDPVR_ES(P_cell, V_cell, cmap,xlimv, ylimv)
    % Accepts either {cell} or HFpEF table fields (cell)
    if iscell(P_cell)
        P_LV = P_cell; V_LV = V_cell;
    else
        P_LV = P_cell; V_LV = V_cell; % already cell arrays
    end

    LVEDP_points = nan(numel(V_LV), 2);
    for i = 1:numel(V_LV)
        if i == 1 || i == numel(V_LV); lw = 2; a = 1.0; else, lw = 0.88; a = 0.4; end
        x = V_LV{i}; y = P_LV{i};
        if ~isempty(x) && ~isempty(y), LVEDP_points(i,:) = [x(1), y(1)]; end
        plot(x, y, '-', 'Color', [cmap(i,:), a], 'LineWidth', lw); hold on;
    end
    xlim(xlimv); ylim(ylimv);

    ok = isfinite(LVEDP_points(:,1)) & isfinite(LVEDP_points(:,2)) & LVEDP_points(:,2) > 0;
    V = LVEDP_points(ok,1); P = LVEDP_points(ok,2);
    if numel(V) < 3, return; end

    edpvr_fun = @(b,Vi) b(1) * (exp(b(2) * (Vi - b(3))) - 1);
    V0_candidates = linspace(10, min(V)-5, 200);
    bestSSE = inf; A0=[]; B0=[]; V00=[];
    for v0 = V0_candidates
        X = V - v0; if any(X <= 0), continue; end
        y = log(P + 1);
        p = polyfit(X, y, 1);
        B_est = max(p(1), 1e-4);
        A_est = max(exp(p(2)), 1e-6);
        P_hat = A_est * (exp(B_est * X) - 1);
        sse = sum((P - P_hat).^2);
        if isfinite(sse) && sse < bestSSE, bestSSE = sse; A0=A_est; B0=B_est; V00=v0; end
    end
    if isempty(A0)
        V00 = 10; X = V - V00; y = log(P + 1);
        p = polyfit(X, y, 1);
        B0 = max(p(1), 1e-4); A0 = max(exp(p(2)), 1e-6);
    end
    init0 = [A0, B0, V00];
    lb = [1e-6, 1e-4, 10]; ub = [50, 1.0, min(V)-1];
    opts = optimset('MaxFunEvals', 5000, 'Display','off');
    fit_lv = lsqcurvefit(edpvr_fun, init0, V, P, lb, ub, opts);
    A_fit = fit_lv(1); B_fit = fit_lv(2); V0_fit = fit_lv(3);
    V0_fixed = max(10, V0_fit);

    ax = gca; hold(ax,'on'); XL = xlim(ax); YL = ylim(ax);
    V_plot = linspace(V0_fixed, 220, 200);
    P_edp = edpvr_fun([A_fit, B_fit, V0_fit], V_plot);
    P_edp(P_edp<YL(1) | P_edp>YL(2)) = NaN;
    plot(ax, V_plot, P_edp, '-', 'LineWidth', 1.6, 'Color', [0.2 0.2 0.2], 'HandleVisibility','off');

    nL = numel(V_LV);
    for i = 1:nL
        v = V_LV{i}(:); p = P_LV{i}(:);
        ok2 = isfinite(v) & isfinite(p) & (v > V0_fixed) & (p > 0);
        if nnz(ok2) < 3, continue; end
        X = v(ok2) - V0_fixed; Y = p(ok2);
        k = Y ./ X; [~, idx] = max(k); k_star = k(idx);
        P_at_X1 = k_star * (XL(1) - V0_fixed);
        P_at_X2 = k_star * (XL(2) - V0_fixed);
        Pmin_by_X = min(P_at_X1, P_at_X2); Pmax_by_X = max(P_at_X1, P_at_X2);
        Pmin = max([YL(1), Pmin_by_X, 0]); Pmax = min([YL(2), Pmax_by_X]);
        if Pmax <= Pmin, continue; end
        P_line = [Pmin, Pmax-20];
        V_line = P_line./k_star + V0_fixed;
        plot(ax, V_line, P_line, '-', 'LineWidth', 1.0, 'Color', [0 0 0 0.35], 'HandleVisibility','off');
    end
    xlim(ax, XL); ylim(ax, YL);
    yticks([0 50 100 150]);
    hold(ax, 'off');
end

function y = clamp01(x), y = max(0, min(1, x)); end
