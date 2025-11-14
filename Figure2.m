clear
load('pacingSUMVirtual1013.mat');
% T1 = load('pacingSUM.mat');
% T2 = load('pacingSUMUW.mat');
% SummaryTable = [T1.SummaryTable;T2.SummaryTable];
Patients2Delete = [6476 12287 21832 22253];
SummaryTable(Patients2Delete,:) = [];

%% remake figure in Δ-space: ΔHR vs ΔFillingP (ΔLAP)
% 规则：baseline 取该曲线第一个点的 LAP；只要有任意点 < baseline（可选容差 tol）→ 标 1，否则标 0

deltaLAP_range = zeros(height(SummaryTable),1);  % 保留"变化幅度"量（range）
PredictedLabelLAP = zeros(height(SummaryTable),1);

tol = 0;  % 如担心数值抖动，可设 tol=0.5（mmHg）等

for i = 1:height(SummaryTable)
    thisLAP = SummaryTable.PacingResults{i}.LAP;
    thisHR  = SummaryTable.PacingResults{i}.HR;


    % 基线（默认取第一个点；若你的基线并非第一个点，把 baseline_idx 改成相应索引）
    baseline_idx = 1;
    baseline_LAP = thisLAP(baseline_idx);
    baseline_HR  = thisHR(baseline_idx);

    % Δ量：相对基线
    SummaryTable.PacingResults{i}.deltaHR  = thisHR  - baseline_HR;
    SummaryTable.PacingResults{i}.deltaLAP = thisLAP - baseline_LAP;

    % 可选：保留变化幅度（原始量的 range）
    deltaLAP_range(i) = range(thisLAP);

    % 分类：是否存在任何点低于 baseline（带或不带容差）
    if min(thisLAP) < baseline_LAP - tol
        PredictedLabelLAP(i) = 1;   % 有下降 → responder
    else
        PredictedLabelLAP(i) = 0;   % 全部不低于 baseline → non-responder
    end
end
%%
% 绘图：ΔHR – ΔLAP
fig = figure('Visible','on'); clf; hold on;

for i = 1:height(SummaryTable)
    dHR  = SummaryTable.PacingResults{i}.deltaHR;
    dLAP = SummaryTable.PacingResults{i}.deltaLAP;

    if PredictedLabelLAP(i) == 1
        plot(dHR, dLAP, '-', 'Color', [0.1 0.7 0.1 0.01], 'LineWidth', 0.08, 'MarkerSize', 5);
    elseif PredictedLabelLAP(i) == 0
        plot(dHR, dLAP, '-', 'Color', [0.9 0.1 0.1 0.01], 'LineWidth', 0.08, 'MarkerSize', 5);
    else
        plot(dHR, dLAP, '-', 'Color', 'k', 'LineWidth', 1.8, 'MarkerSize', 5);
    end
end


anchorX = [10, 20, 30];   % 指定三个HR点
group1 = PredictedLabelLAP == 1;  % 绿色 responder
group0 = PredictedLabelLAP == 0;  % 红色 non-responder

% 初始化
meanLAP1 = nan(size(anchorX));  stdLAP1 = nan(size(anchorX));
meanLAP0 = nan(size(anchorX));  stdLAP0 = nan(size(anchorX));

for j = 1:numel(anchorX)
    lapVals1 = []; lapVals0 = [];
    for i = 1:height(SummaryTable)
        dHR  = SummaryTable.PacingResults{i}.deltaHR;
        dLAP = SummaryTable.PacingResults{i}.deltaLAP;

        % 查找匹配 anchor
        [~, idx] = min(abs(dHR - anchorX(j)));
        if PredictedLabelLAP(i) == 1 && ~isempty(idx)
            lapVals1(end+1) = dLAP(idx);
        elseif PredictedLabelLAP(i) == 0 && ~isempty(idx)
            lapVals0(end+1) = dLAP(idx);
        end
    end
    meanLAP1(j) = mean(lapVals1, 'omitnan');
    stdLAP1(j)  = std(lapVals1,  'omitnan');
    meanLAP0(j) = mean(lapVals0, 'omitnan');
    stdLAP0(j)  = std(lapVals0,  'omitnan');
end

errorbar(anchorX, meanLAP1, stdLAP1, 'o-', 'Color', [0.1 0.7 0.1], ...
    'LineWidth', 2.0, 'MarkerFaceColor', [0.1 0.7 0.1]);
errorbar(anchorX, meanLAP0, stdLAP0, 'o-', 'Color', [0.9 0.1 0.1], ...
    'LineWidth', 2.0, 'MarkerFaceColor', [0.9 0.1 0.1]);


box on; 
pbaspect([1 1 1]);  % 正方形比例

xlim([0 35]);
ylim([-20 10]);

set(gca, 'XTick', [10 20 30], 'XTickLabel', {'','',''}, ...
         'YTickLabel', [], 'TickLength', [0.03 0.03]);


% 统一字体设置
set(gca, 'FontName', 'Arial', 'FontSize', 16);
box on;

% 图像尺寸与导出参数
w_in = 2.4;  
h_in = 2.4; 
dpi  = 600;

set(gcf, 'Units', 'inches', 'Position', [1 1 w_in h_in]);
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 w_in h_in]);
set(gcf, 'InvertHardcopy', 'off');  % 防止导出时背景反转

pbaspect([1 1 1]); 


% 输出路径（你可按需修改文件名）
outpath = fullfile('C:', 'Users', 'fenggu', ...
    'University of Michigan Dropbox', 'Feng Gu', ...
    'DTPacingHFpEF', 'DeltaHR_vs_DeltaLAP_VP.png');

exportgraphics(gcf, outpath, 'Resolution', dpi, ...
    'ContentType', 'image', 'BackgroundColor', 'white');

%% ΔHR – ΔCO（与 ΔHR–ΔLAP 同款格式）
PredictedLabelCO = zeros(height(SummaryTable),1);
tol = 0;  % 单位同 CO（L/min），可设 0.05 之类

% 计算标签 & 写入 deltaHR / deltaCO
for i = 1:height(SummaryTable)
    thisCO = SummaryTable.PacingResults{i}.CO(:);  % 行向量
    thisHR = SummaryTable.PacingResults{i}.HR(:);  % 行向量

    baseline_CO = thisCO(1);
    baseline_HR = thisHR(1);

    % 标签：任一点 CO > baseline(+tol) → 1
    if max(thisCO) > baseline_CO + tol
        PredictedLabelCO(i) = 1;
    else
        PredictedLabelCO(i) = 0;
    end

    % Δ值（相对基线）
    SummaryTable.PacingResults{i}.deltaCO = thisCO - baseline_CO;
end

%% 绘图（个体曲线 + 误差棒）
fig = figure('Visible','on'); clf; hold on;

% 个体曲线（半透明）
for i = 1:height(SummaryTable)
    dHR = SummaryTable.PacingResults{i}.deltaHR;
    dCO = SummaryTable.PacingResults{i}.deltaCO;

    if PredictedLabelCO(i) == 1
        plot(dHR, dCO, '-', 'Color', [0.1 0.7 0.1 0.01], 'LineWidth', 0.08);
    elseif PredictedLabelCO(i) == 0
        plot(dHR, dCO, '-', 'Color', [0.9 0.1 0.1 0.01], 'LineWidth', 0.08);
    else
        plot(dHR, dCO, '-', 'Color', 'k', 'LineWidth', 1.8);
    end
end

% 误差棒（ΔHR = 10/20/30）
anchorX = [10 20 30];
meanCO1 = nan(size(anchorX));  stdCO1 = nan(size(anchorX));  % 绿色（1）
meanCO0 = nan(size(anchorX));  stdCO0 = nan(size(anchorX));  % 红色（0）

for j = 1:numel(anchorX)
    vals1 = []; vals0 = [];
    for i = 1:height(SummaryTable)
        dHR = SummaryTable.PacingResults{i}.deltaHR;
        dCO = SummaryTable.PacingResults{i}.deltaCO;
        [~, idx] = min(abs(dHR - anchorX(j)));  % 最近点
        if PredictedLabelCO(i) == 1
            vals1(end+1) = dCO(idx); 
        elseif PredictedLabelCO(i) == 0
            vals0(end+1) = dCO(idx);
        end
    end
    meanCO1(j) = mean(vals1, 'omitnan');  stdCO1(j) = std(vals1, 'omitnan');
    meanCO0(j) = mean(vals0, 'omitnan');  stdCO0(j) = std(vals0, 'omitnan');
end

% 绘制 errorbar（不透明强调）
errorbar(anchorX, meanCO1, stdCO1, 'o-', 'Color', [0.1 0.7 0.1], ...
    'LineWidth', 2.0, 'MarkerFaceColor', [0.1 0.7 0.1], 'MarkerSize', 5);
errorbar(anchorX, meanCO0, stdCO0, 'o-', 'Color', [0.9 0.1 0.1], ...
    'LineWidth', 2.0, 'MarkerFaceColor', [0.9 0.1 0.1], 'MarkerSize', 5);


xlim([0 35]);                 % 需要到 35 留空
ylim([-.5  1]);                % 你可按实际 ΔCO 调整，例如 [-1 3]
set(gca, 'XTick', [10 20 30], 'XTickLabel', {'','',''}, ...
         'YTickLabel', [], 'TickLength', [0.03 0.03]);
box on; pbaspect([1 1 1]);    % 正方形比例
set(gca, 'FontName', 'Arial', 'FontSize', 16);

w_in = 2.4;  h_in = 2.4;  dpi = 600;
set(gcf, 'Units','inches','Position',[1 1 w_in h_in]);
set(gcf, 'PaperUnits','inches','PaperPosition',[0 0 w_in h_in]);
set(gcf, 'InvertHardcopy','off');

% 防止导出带工具栏（兼容写法）
if isprop(gca,'Toolbar'); set(get(gca,'Toolbar'),'Visible','off'); end

outpath = fullfile('C:', 'Users', 'fenggu', ...
    'University of Michigan Dropbox', 'Feng Gu', ...
    'DTPacingHFpEF', 'DeltaHR_vs_DeltaCO_VP.png');

exportgraphics(gcf, outpath, 'Resolution', dpi, ...
    'ContentType','image','BackgroundColor','white');

%% ΔHR – ΔLVO2E（与前两张图同款格式）
PredictedLabelLVO2E = zeros(height(SummaryTable),1);
tol = 0;  % 单位同 LVO2E，可设 0.01 等

% 计算标签 & 写入 deltaHR / deltaLVO2E
for i = 1:height(SummaryTable)
    LVO2E = SummaryTable.PacingResults{i}.LVO2E_all(:);  % 列向量

    baseE = LVO2E(1);

    % 标签：任一点 LVO2E > baseline(+tol) → 1
    if max(LVO2E) > baseE + tol
        PredictedLabelLVO2E(i) = 1;
    else
        PredictedLabelLVO2E(i) = 0;
    end

    % Δ值（相对基线）
    SummaryTable.PacingResults{i}.deltaLVO2E = LVO2E - baseE;
end

% 绘图（个体曲线）
fig = figure('Visible','on'); clf; hold on;

for i = 1:height(SummaryTable)
    dHR   = SummaryTable.PacingResults{i}.deltaHR;
    dLVOE = SummaryTable.PacingResults{i}.deltaLVO2E;

    if PredictedLabelLVO2E(i) == 1
        plot(dHR, dLVOE, '-', 'Color', [0.1 0.7 0.1 0.01], 'LineWidth', 0.08);
    elseif PredictedLabelLVO2E(i) == 0
        plot(dHR, dLVOE, '-', 'Color', [0.9 0.1 0.1 0.01], 'LineWidth', 0.08);
    else
        plot(dHR, dLVOE, '-', 'Color', 'k', 'LineWidth', 1.8);
    end
end

% 误差棒（ΔHR = 10/20/30）
anchorX = [10 20 30];
meanE1 = nan(size(anchorX));  stdE1 = nan(size(anchorX));  
meanE0 = nan(size(anchorX));  stdE0 = nan(size(anchorX));  

for j = 1:numel(anchorX)
    vals1 = []; vals0 = [];
    for i = 1:height(SummaryTable)
        dHR   = SummaryTable.PacingResults{i}.deltaHR;
        dLVOE = SummaryTable.PacingResults{i}.deltaLVO2E;
        [~, idx] = min(abs(dHR - anchorX(j)));  
        if PredictedLabelLVO2E(i) == 1
            vals1(end+1) = dLVOE(idx);
        elseif PredictedLabelLVO2E(i) == 0
            vals0(end+1) = dLVOE(idx);
        end
    end
    meanE1(j) = mean(vals1, 'omitnan');  stdE1(j) = std(vals1, 'omitnan');
    meanE0(j) = mean(vals0, 'omitnan');  stdE0(j) = std(vals0, 'omitnan');
end

% 绘制 errorbar（不透明强调）
errorbar(anchorX, meanE1, stdE1, 'o-', 'Color', [0.1 0.7 0.1], ...
    'LineWidth', 2.0, 'MarkerFaceColor', [0.1 0.7 0.1], 'MarkerSize', 5);
errorbar(anchorX, meanE0, stdE0, 'o-', 'Color', [0.9 0.1 0.1], ...
    'LineWidth', 2.0, 'MarkerFaceColor', [0.9 0.1 0.1], 'MarkerSize', 5);

% 轴与样式（与前两张一致）
xlim([0 35]);                 % 留空到 35
% 视你数据调整 ΔLVO2E 纵轴范围（例如：[-0.15 0.25]）
ylim([-100 50]);

set(gca, 'XTick', [10 20 30], 'XTickLabel', {'','',''}, ...
         'YTickLabel', [], 'TickLength', [0.03 0.03]);
box on; pbaspect([1 1 1]);
set(gca, 'FontName', 'Arial', 'FontSize', 16);

% 导出（600dpi，白底）
w_in = 2.4;  h_in = 2.4;  dpi = 600;
set(gcf, 'Units','inches','Position',[1 1 w_in h_in]);
set(gcf, 'PaperUnits','inches','PaperPosition',[0 0 w_in h_in]);
set(gcf, 'InvertHardcopy','off');

% 兼容关闭工具栏，避免导出警告
if isprop(gca,'Toolbar'); set(get(gca,'Toolbar'),'Visible','off'); end

outpath = fullfile('C:', 'Users', 'fenggu', ...
    'University of Michigan Dropbox', 'Feng Gu', ...
    'DTPacingHFpEF', 'DeltaHR_vs_DeltaLVO2E_VP.png');

exportgraphics(gcf, outpath, 'Resolution', dpi, ...
    'ContentType','image','BackgroundColor','white');

%% ΔHR – ΔSBP（同款格式）
PredictedLabelSBP = nan(height(SummaryTable),1);
tol = 8.5;   % mmHg

for i = 1:height(SummaryTable)
    Ti = SummaryTable.PacingResults{i};
    if isempty(Ti) || ~ismember('SBP', Ti.Properties.VariableNames)
        continue
    end

    HR  = Ti.HR(:);             % 行向量
    SBP = Ti.SBP(:);

    if any(~isreal(HR)) || any(~isreal(SBP))
        error('Patient %d has invalid HR/SBP (complex).', i);
    end

    % 基线：第一个非 NaN
    bidx = find(~isnan(SBP) & ~isnan(HR), 1, 'first');
    if isempty(bidx)
        continue
    end
    baseSBP = SBP(bidx);
    baseHR  = HR(bidx);

    % Δ值（写回以便后续使用）
    SummaryTable.PacingResults{i}.deltaHR  = HR  - baseHR;
    SummaryTable.PacingResults{i}.deltaSBP = SBP - baseSBP;

    % 判定：是否任何点 < baseline - tol（忽略 NaN）
    if min(SBP,[],'omitnan') < baseSBP - tol
        PredictedLabelSBP(i) = 1;     % responder（下降）
    else
        PredictedLabelSBP(i) = 0;     % non-responder
    end
end
%%
% 绘图：个体曲线（半透明）
fig = figure('Visible','on'); clf; hold on;
for i = 1:height(SummaryTable)
    Ti = SummaryTable.PacingResults{i};
    if isempty(Ti) || ~ismember('deltaSBP', Ti.Properties.VariableNames), continue; end

    dHR  = Ti.deltaHR;
    dSBP = Ti.deltaSBP;

    if PredictedLabelSBP(i) == 1
        plot(dHR, dSBP, '-', 'Color', [0.1 0.7 0.1 0.01], 'LineWidth', 0.08);
    elseif PredictedLabelSBP(i) == 0
        plot(dHR, dSBP, '-', 'Color', [0.9 0.1 0.1 0.01], 'LineWidth', 0.08);
    else
        plot(dHR, dSBP, '-', 'Color', 'k', 'LineWidth', 1.8);
    end
end

% 误差棒：ΔHR = 10/20/30 处的 mean±SD
anchorX = [10 20 30];
meanS1 = nan(size(anchorX));  stdS1 = nan(size(anchorX));  % 绿色（1）
meanS0 = nan(size(anchorX));  stdS0 = nan(size(anchorX));  % 红色（0）

for j = 1:numel(anchorX)
    vals1 = []; vals0 = [];
    for i = 1:height(SummaryTable)
        Ti = SummaryTable.PacingResults{i};
        if isempty(Ti) || ~ismember('deltaSBP', Ti.Properties.VariableNames), continue; end
        dHR  = Ti.deltaHR;
        dSBP = Ti.deltaSBP;
        if all(isnan(dHR)) || all(isnan(dSBP)), continue; end
        [~, idx] = min(abs(dHR - anchorX(j)));   % 最近点
        if ~isnan(idx) && ~isnan(dSBP(idx))
            if PredictedLabelSBP(i) == 1
                vals1(end+1) = dSBP(idx); 
            elseif PredictedLabelSBP(i) == 0
                vals0(end+1) = dSBP(idx); 
            end
        end
    end
    meanS1(j) = mean(vals1, 'omitnan');  stdS1(j) = std(vals1, 'omitnan');
    meanS0(j) = mean(vals0, 'omitnan');  stdS0(j) = std(vals0, 'omitnan');
end

errorbar(anchorX, meanS1, stdS1, 'o-', 'Color', [0.1 0.7 0.1], ...
    'LineWidth', 2.0, 'MarkerFaceColor', [0.1 0.7 0.1], 'MarkerSize', 5);
errorbar(anchorX, meanS0, stdS0, 'o-', 'Color', [0.9 0.1 0.1], ...
    'LineWidth', 2.0, 'MarkerFaceColor', [0.9 0.1 0.1], 'MarkerSize', 5);

% 轴与样式（同款）
xlim([0 35]);
% 依据数据调整 ΔSBP 范围（示例）：
ylim([-30 10]);

set(gca, 'XTick', [10 20 30], 'XTickLabel', {'','',''}, ...
         'YTickLabel', [], 'TickLength', [0.03 0.03]);
box on; pbaspect([1 1 1]);
set(gca, 'FontName', 'Arial', 'FontSize', 16);

% 导出（600 dpi，白底）
w_in = 2.4;  h_in = 2.4;  dpi = 600;
set(gcf, 'Units','inches','Position',[1 1 w_in h_in]);
set(gcf, 'PaperUnits','inches','PaperPosition',[0 0 w_in h_in]);
set(gcf, 'InvertHardcopy','off');
if isprop(gca,'Toolbar'); set(get(gca,'Toolbar'),'Visible','off'); end

outpath = fullfile('C:', 'Users', 'fenggu', ...
    'University of Michigan Dropbox', 'Feng Gu', ...
    'DTPacingHFpEF', 'DeltaHR_vs_DeltaSBP_VP.png');

exportgraphics(gcf, outpath, 'Resolution', dpi, ...
    'ContentType','image','BackgroundColor','white');

%%
% 
% figure();
% for i = 1:7
%   SummaryTable.PacingResults{15}.LVO2E_all(i) = SummaryTable.PacingResults{15}.LVO2E_all(i)*1000;
% end