% %% Make summary table
% clear
% Order = [1	2 13 15	16	17	20	22	25	26	35	37	42	45	46	48	50	55	57	60	62	66	67	77	78	79	84	87	88	90	91	92	93	94	97	98	99	102	103	104	106	108	111	112	113	118	122	126	130	131	132	133	137	141	146	147	150	151	164	167	168	173	177	178	186	188	195	196	197	199	202	203	204	205	206	207	209	212	213	214	217	218	222	223	224	230	232	234	235	236	240	244	246	248	250	252	267	269	271	274	281	282	283	284	293	294	300	302	303	304	306	309	316	318	324	325	327	329	335	341	342	343	344	345	346	350	351	368];
% SummaryTable =  table(Order');
% SummaryTable.Properties.VariableNames = {'PatientNo'};
% for i = 1:length(Order)
%     PATIENT_NO = Order(i);
%     S = load(sprintf('PacingResultsUM/PatientNO%d/PacingTable.mat', PATIENT_NO));
%     SummaryTable.PacingResults{i} = S.TableatdifferentHR;
% end
% save pacingSUM.mat SummaryTable
% %%
% clear
% Order = [9 14	15	17	18	19	22	25	26	30	33	34	37	44	46	48	52	53];
% SummaryTable =  table(Order');
% SummaryTable.Properties.VariableNames = {'PatientNo'};
% for i = 1:length(Order)
%     PATIENT_NO = Order(i);
%     S = load(sprintf('PacingResultsUW/PatientNO%d/PacingTable.mat', PATIENT_NO));
%     SummaryTable.PacingResults{i} = S.TableatdifferentHR;
% end
% save pacingSUMUW.mat SummaryTable
%%
clear
load('pacingSUMVirtual1013.mat');
% T1 = load('pacingSUM.mat');
% T2 = load('pacingSUMUW.mat');
% SummaryTable = [T1.SummaryTable;T2.SummaryTable];
Patients2Delete = [6476 12287 21832 22253];
SummaryTable(Patients2Delete,:) = [];
n = 7;  % 假设 P_LV 是颜色映射的长度
start_color = [0.6350, 0.0780, 0.1840];  % 起点（暗红）
end_color   = [0.5, 0.5, 0.5];                % 终点（蓝）
RawT = readtable("inputFGclassifierWithLabel.CSV",'VariableNamingRule','preserve');
% RawT(Patients2Delete,:) = [];
cmap = [linspace(start_color(1), end_color(1), n)', ...
    linspace(start_color(2), end_color(2), n)', ...
    linspace(start_color(3), end_color(3), n)'];
%%
% for order = 1:height(RawT)
%    RawT.LAP_B(order) = round(SummaryTable.PacingResults{order}.LAP(1),1);
%    RawT.LAP_P(order) = round(SummaryTable.PacingResults{order}.LAP(2),1);
%    RawT.EF_B(order) = round(SummaryTable.PacingResults{order}.EF(1),1);
%    RawT.EF_P(order) = round(SummaryTable.PacingResults{order}.EF(2),1);
%    RawT.CO_B(order) = round(SummaryTable.PacingResults{order}.CO(1),1);
%    RawT.CO_P(order) = round(SummaryTable.PacingResults{order}.CO(2),1);
%    RawT.LAV_B(order) = round(SummaryTable.PacingResults{order}.LAV(1),1);
%    RawT.LAV_P(order) = round(SummaryTable.PacingResults{order}.LAV(2),1);
%    RawT.LVO2E_B(order) = round(SummaryTable.PacingResults{order}.LVO2E_all(1),2);
%    RawT.LVO2E_P(order) = round(SummaryTable.PacingResults{order}.LVO2E_all(2),2);
%    RawT.RVO2E_B(order) = round(SummaryTable.PacingResults{order}.RVO2E_all(1),2);
%    RawT.RVO2E_P(order) = round(SummaryTable.PacingResults{order}.RVO2E_all(2),2);
%    RawT.LVMVO2_B(order) = round(SummaryTable.PacingResults{order}.MVO2_LV(1),2);
%    RawT.LVMVO2_P(order) = round(SummaryTable.PacingResults{order}.MVO2_LV(2),2);
%    RawT.RVMVO2_B(order) = round(SummaryTable.PacingResults{order}.MVO2_RV(1),2);
%    RawT.RVMVO2_P(order) = round(SummaryTable.PacingResults{order}.MVO2_RV(2),2);
%    RawT.LVME_B(order) = round(SummaryTable.PacingResults{order}.CE_LV_all(1),4);
%    RawT.LVME_P(order) = round(SummaryTable.PacingResults{order}.CE_LV_all(2),4);
%    RawT.RVME_B(order) = round(SummaryTable.PacingResults{order}.CE_RV_all(1),4);
%    RawT.RVME_P(order) = round(SummaryTable.PacingResults{order}.CE_RV_all(2),4);
% end
%%
% for order = 1:height(RawT)
%     SummaryTable.PacingResults{order} = ...
%         renamevars(SummaryTable.PacingResults{order}, ...
%         {'O2E_LV_all','O2E_RV_all'}, {'LVO2E_all','RVO2E_all'});
% end


% %% remake figure based on classifiers (simplified rule)
% % 规则：baseline 取该曲线第一个点的 LAP；只要有任意点 < baseline（可选容差 tol）→ 标 1，否则标 0
% 
% deltaLAP = zeros(height(SummaryTable),1);  % 若后面不用可删
% PredictedLabelLAP = zeros(height(SummaryTable),1);
% 
% tol = 0;  % 如担心数值抖动，可设 tol=0.5（mmHg）等
% 
% for i = 1:height(SummaryTable)
%     thisLAP = SummaryTable.PacingResults{i}.LAP;
% 
%     if ~all(isreal(thisLAP))
%         error('Patient %d has invalid LAP (complex/NaN).', i);
%     end
% 
%     % 可选：保留一个"变化幅度"的量
%     deltaLAP(i) = range(thisLAP);
% 
%     % baseline 取第一个点（如你的 baseline 不在第一个点，请把 baseline_idx 改成相应索引）
%     baseline_val = thisLAP(1);
% 
%     % 判定：是否存在任何点低于 baseline（带或不带容差）
%     if min(thisLAP) < baseline_val - tol
%         PredictedLabelLAP(i) = 1;   % 有下降 → responder
%     else
%         PredictedLabelLAP(i) = 0;   % 全部不低于 baseline → non-responder
%     end
% end
% 
% 
% fig = figure('Visible','off');
% clf; hold on;
% 
% for i  = 1:height(SummaryTable)
%     if PredictedLabelLAP(i) == 1
%         plot(SummaryTable.PacingResults{i}.HR, SummaryTable.PacingResults{i}.LAP, '-o', ...
%             'Color', [0.1 0.7 0.1 0.5], 'LineWidth', 1.8)
%     elseif PredictedLabelLAP(i) == 0
%         plot(SummaryTable.PacingResults{i}.HR, SummaryTable.PacingResults{i}.LAP, '-o', ...
%             'Color', [0.9 0.1 0.1 0.5], 'LineWidth', 1.8)
%     else
%         plot(SummaryTable.PacingResults{i}.HR, SummaryTable.PacingResults{i}.LAP, '-o', ...
%             'Color', 'k', 'LineWidth', 1.8)
%     end
% end
% 
% % 重要：设置坐标标签和美化
% xlabel('HR (beats/min)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica')
% ylabel('LV filling pressure (mmHg)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica')
% 
% set(gca, ...
%     'FontSize', 13, ...
%     'FontName', 'Helvetica', ...
%     'LineWidth', 1.2, ...
%     'TickDir', 'out', ...
%     'TickLength', [0.018 0.018])
% box off;
% pbaspect([1 1 1]); % 正方形比例
% 
% % 可选：加图例，让不同类别一目了然
% h1 = plot(nan, nan, '--', 'Color', [0.1 0.7 0.1], 'LineWidth', 1.8); % 用nan, nan画假线用于legend
% h2 = plot(nan, nan, '--', 'Color', [0.9 0.1 0.1], 'LineWidth', 1.8);
% legend([h1 h2], {'Responders', 'Non-responders'}, ...
%     'Location', 'northeast', 'Box', 'off', 'FontSize', 12);
% 
% % 增加适当留白
% ax = gca;
% ax.Position = [0.15 0.17 0.78 0.76];
% 
% % 网格线（可选）
% grid on
% % 保存为 PNG，300 dpi
% output_filename = 'LAP_vs_HR.png';
% exportgraphics(fig, output_filename, 'Resolution', 300);
% close(fig);  % 关闭 figure，避免显示
% 
% RawT.LabelFillingP = PredictedLabelLAP;
% %% 2 Based on EF (simplified)
% % 规则：若该患者的任意 EF 点 < 50% → 标 0；否则（全程 ≥ 50%）→ 标 1
% % 兼容 EF 以小数(0~1)或百分比(0~100)两种刻度
% 
% PredictedLabelEF = zeros(height(SummaryTable),1);
% tol = 0;          % 可选容差；若担心数值抖动，设 tol=0.005(=0.5%) 等
% 
% for i = 1:height(SummaryTable)
%     thisEF = SummaryTable.PacingResults{i}.EF(:)';  % 行向量
% 
%     if ~all(isreal(thisEF))
%         error('Patient %d has invalid EF (complex/NaN).', i);
%     end
% 
%     % 归一化到 0~1（若数据本身已是百分比）
%     if max(thisEF) > 1
%         thisEF = thisEF / 100;
%     end
% 
%     % 判定：是否有任何点 < 0.50（带容差）
%     if min(thisEF) < (0.50 - tol)
%         PredictedLabelEF(i) = 0;   % 有下降到 50% 以下
%     else
%         PredictedLabelEF(i) = 1;   % 全程 ≥ 50%
%     end
% 
%     % 可选：若认为应当单调下降，可做一致性检查但不影响标签
%     % dEF = diff(thisEF);
%     % if any(dEF > tol*0.1)
%     %     warning('Patient %d: EF not strictly decreasing (FYI).', i);
%     % end
% end
% 
% fig = figure('Visible','off');
% clf; hold on;
% 
% for i  = 1:height(SummaryTable)
%     thisEF = SummaryTable.PacingResults{i}.EF;
%     if PredictedLabelEF(i) == 1
%         plot(SummaryTable.PacingResults{i}.HR,thisEF, '-o', ...
%             'Color', [0.1 0.7 0.1 0.5], 'LineWidth', 1.8)
%     elseif PredictedLabelEF(i) == 0
%         plot(SummaryTable.PacingResults{i}.HR, thisEF, '-o', ...
%             'Color', [0.9 0.1 0.1 0.5], 'LineWidth', 1.8)
%     else
%         plot(SummaryTable.PacingResults{i}.HR, thisEF, '-o', ...
%             'Color', 'k', 'LineWidth', 1.8)
%     end
% end
% 
% % 重要：设置坐标标签和美化
% xlabel('HR (beats/min)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica')
% ylabel('EF (%)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica')
% 
% set(gca, ...
%     'FontSize', 13, ...
%     'FontName', 'Helvetica', ...
%     'LineWidth', 1.2, ...
%     'TickDir', 'out', ...
%     'TickLength', [0.018 0.018])
% box off;
% pbaspect([1 1 1]); % 正方形比例
% 
% % 可选：加图例，让不同类别一目了然
% h1 = plot(nan, nan, '--', 'Color', [0.1 0.7 0.1], 'LineWidth', 1.8); % 用nan, nan画假线用于legend
% h2 = plot(nan, nan, '--', 'Color', [0.9 0.1 0.1], 'LineWidth', 1.8);
% legend([h1 h2], {'Responders', 'Non-responders'}, ...
%     'Location', 'northeast', 'Box', 'off', 'FontSize', 12);
% 
% % 增加适当留白
% ax = gca;
% ax.Position = [0.15 0.17 0.78 0.76];
% 
% % 网格线（可选）
% grid on
% 
% % 保存为 PNG，300 dpi
% output_filename = 'EF_vs_HR.png';
% exportgraphics(fig, output_filename, 'Resolution', 300);
% 
% % 关闭 figure
% close(fig);
% RawT.LabelEF = PredictedLabelEF;
% %% 3 Based on CO (simplified)
% % 规则：baseline 取该曲线第一个点（若你的 baseline HR 已知，见下方替代写法）
% % 只要存在 CO > baseline(+tol) → 1；否则 → 0
% 
% PredictedLabelCO = zeros(height(SummaryTable),1);
% tol = 0;  % 可选容差（单位同 CO，例如 L/min）。若担心数值抖动，可设 0.05 等
% 
% for i = 1:height(SummaryTable)
%     thisCO = SummaryTable.PacingResults{i}.CO(:)';  % 行向量
% 
%     if ~all(isreal(thisCO))
%         error('Patient %d has invalid CO (complex/NaN).', i);
%     end
% 
%     % baseline = 第一个点
%     baseline_val = thisCO(1);
% 
%     % 只要有一个点高于 baseline(+tol) 就标 1
%     if max(thisCO) > baseline_val + tol
%         PredictedLabelCO(i) = 1;
%     else
%         PredictedLabelCO(i) = 0;
%     end
% end
% 
% 
% fig = figure('Visible','off');
% clf; hold on;
% 
% for i  = 1:height(SummaryTable)
%     thisCO = SummaryTable.PacingResults{i}.CO;
%     if PredictedLabelCO(i) == 1
%         plot(SummaryTable.PacingResults{i}.HR,thisCO, '-o', ...
%             'Color', [0.1 0.7 0.1 0.5], 'LineWidth', 1.8)
%     elseif PredictedLabelCO(i) == 0
%         plot(SummaryTable.PacingResults{i}.HR, thisCO, '-o', ...
%             'Color', [0.9 0.1 0.1 0.5], 'LineWidth', 1.8)
%     else
%         plot(SummaryTable.PacingResults{i}.HR, thisCO, '-o', ...
%             'Color', 'k', 'LineWidth', 1.8)
%     end
% end
% 
% % 重要：设置坐标标签和美化
% xlabel('HR (beats/min)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica')
% ylabel('CO (L/min)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica')
% 
% set(gca, ...
%     'FontSize', 13, ...
%     'FontName', 'Helvetica', ...
%     'LineWidth', 1.2, ...
%     'TickDir', 'out', ...
%     'TickLength', [0.018 0.018])
% box off;
% pbaspect([1 1 1]); % 正方形比例
% 
% % 可选：加图例，让不同类别一目了然
% h1 = plot(nan, nan, '--', 'Color', [0.1 0.7 0.1], 'LineWidth', 1.8); % 用nan, nan画假线用于legend
% h2 = plot(nan, nan, '--', 'Color', [0.9 0.1 0.1], 'LineWidth', 1.8);
% legend([h1 h2], {'Responders', 'Non-responders'}, ...
%     'Location', 'northeast', 'Box', 'off', 'FontSize', 12);
% 
% % 增加适当留白
% ax = gca;
% ax.Position = [0.15 0.17 0.78 0.76];
% 
% % 网格线（可选）
% grid on
% output_filename = 'CO_vs_HR.png';
% exportgraphics(fig, output_filename, 'Resolution', 300);
% 
% % 关闭 figure
% close(fig);
% 
% RawT.LabelCO = PredictedLabelCO;
% 
% %% Based on LAV (simplified)
% % 规则：baseline 取曲线第一个点（若有固定 baseline HR，见下方替代写法）
% % 只要 min(LAV) < baseline(+tol) → 1；否则 → 0
% 
% PredictedLabelLAV = zeros(height(SummaryTable),1);
% tol = 0;  % 可选容差（单位同 LAV；如担心数值抖动可设成 0.5 等）
% 
% for i = 1:height(SummaryTable)
%     thisLAV = SummaryTable.PacingResults{i}.LAV(:)';  % 行向量
% 
%     if ~all(isreal(thisLAV))
%         error('Patient %d has invalid LAV (complex/NaN).', i);
%     end
% 
%     % baseline = 第一个点
%     baseline_val = thisLAV(1);
% 
%     % 只要有一个点低于 baseline(+tol) 就标 1（下降）
%     if min(thisLAV) < baseline_val - tol
%         PredictedLabelLAV(i) = 1;
%     else
%         PredictedLabelLAV(i) = 0;
%     end
% end
% 
% 
% fig = figure('Visible','off');
% clf; hold on;
% 
% for i  = 1:height(SummaryTable)
%     thisLAV = SummaryTable.PacingResults{i}.LAV;
%     if PredictedLabelLAV(i) == 1
%         plot(SummaryTable.PacingResults{i}.HR,thisLAV, '-o', ...
%             'Color', [0.1 0.7 0.1 0.5], 'LineWidth', 1.8)
%     elseif PredictedLabelLAV(i) == 0
%         plot(SummaryTable.PacingResults{i}.HR, thisLAV, '-o', ...
%             'Color', [0.9 0.1 0.1 0.5], 'LineWidth', 1.8)
%     else
%         plot(SummaryTable.PacingResults{i}.HR, thisLAV, '-o', ...
%             'Color', 'k', 'LineWidth', 1.8)
%     end
% end
% 
% % 重要：设置坐标标签和美化
% xlabel('HR (beats/min)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica')
% ylabel('LA Volume (ml)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica')
% 
% set(gca, ...
%     'FontSize', 13, ...
%     'FontName', 'Helvetica', ...
%     'LineWidth', 1.2, ...
%     'TickDir', 'out', ...
%     'TickLength', [0.018 0.018])
% box off;
% pbaspect([1 1 1]); % 正方形比例
% 
% % 可选：加图例，让不同类别一目了然
% h1 = plot(nan, nan, '--', 'Color', [0.1 0.7 0.1], 'LineWidth', 1.8); % 用nan, nan画假线用于legend
% h2 = plot(nan, nan, '--', 'Color', [0.9 0.1 0.1], 'LineWidth', 1.8);
% legend([h1 h2], {'Responders', 'Non-responders'}, ...
%     'Location', 'northeast', 'Box', 'off', 'FontSize', 12);
% 
% % 增加适当留白
% grid on
% ax = gca;
% ax.Position = [0.15 0.17 0.78 0.76];
% output_filename = 'LAV_vs_HR.png';
% exportgraphics(fig, output_filename, 'Resolution', 300);
% 
% 
% % 关闭 figure
% close(fig);
% 
% RawT.LabelLAV = PredictedLabelLAV;
% 
%% remake figure based on classifiers
% 5 Based on O2 effiency' Dan's creation
type = zeros(height(SummaryTable),1);
for i = 1:height(SummaryTable)
    PATIENT_NO = SummaryTable.PatientNo(i);
    % if i<= 128
    % load(sprintf('SimsUMFinal/P_NO%d.mat',PATIENT_NO)); % load results from great lake clusters
    % else
    %     load(sprintf('SimsUWFinal/P_NO%d.mat',PATIENT_NO)); % load results from great lake clusters
    % end
    % rho_myo = 1.055;
    % LV_m = rho_myo * (output.params.Vw_LV+output.params.Vw_SEP);
    % RV_m = rho_myo * (output.params.Vw_RV);
    % thisO2E = SummaryTable.PacingResults{i}.CO./(SummaryTable.PacingResults{i}.MVO2_LV./100*LV_m);
    thisO2E = SummaryTable.PacingResults{i}.LVO2E_all;
    if any(thisO2E(2:end)>thisO2E(1))
        type(i) = 1;
    end
end

PredictedLabelLVO2E = type;

fig = figure('Visible','off');
clf; hold on;

for i  = 1:height(SummaryTable)
    % PATIENT_NO = SummaryTable.PatientNo(i);
    % if i<= 128
    %     load(sprintf('SimsUMFinal/P_NO%d.mat',PATIENT_NO)); % load results from great lake clusters
    % else
    %     load(sprintf('SimsUWFinal/P_NO%d.mat',PATIENT_NO)); % load results from great lake clusters
    % end
    % rho_myo = 1.055;
    % LV_m = rho_myo * (output.params.Vw_LV+output.params.Vw_SEP);
    % RV_m = rho_myo * (output.params.Vw_RV);
    % thisO2E = SummaryTable.PacingResults{i}.CO*1000./(SummaryTable.PacingResults{i}.MVO2_LV./100*LV_m);
    thisO2E = SummaryTable.PacingResults{i}.LVO2E_all;
    if PredictedLabelLVO2E(i) == 1
        plot(SummaryTable.PacingResults{i}.HR,thisO2E, '-o', ...
            'Color', [0.1 0.7 0.1 0.5], 'LineWidth', 1.8)
    elseif PredictedLabelLVO2E(i) == 0
        plot(SummaryTable.PacingResults{i}.HR, thisO2E, '-o', ...
            'Color', [0.9 0.1 0.1 0.5], 'LineWidth', 1.8)
    else
        plot(SummaryTable.PacingResults{i}.HR, thisO2E, '-o', ...
            'Color', 'k', 'LineWidth', 1.8)
    end
end

% 重要：设置坐标标签和美化
xlabel('HR (beats/min)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica')
ylabel('LV O2 effiency', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica')

set(gca, ...
    'FontSize', 13, ...
    'FontName', 'Helvetica', ...
    'LineWidth', 1.2, ...
    'TickDir', 'out', ...
    'TickLength', [0.018 0.018])
box off;
pbaspect([1 1 1]); % 正方形比例

% 可选：加图例，让不同类别一目了然
h1 = plot(nan, nan, '--', 'Color', [0.1 0.7 0.1], 'LineWidth', 1.8); % 用nan, nan画假线用于legend
h2 = plot(nan, nan, '--', 'Color', [0.9 0.1 0.1], 'LineWidth', 1.8);
legend([h1 h2], {'Responders', 'Non-responders'}, ...
    'Location', 'northeast', 'Box', 'off', 'FontSize', 12);

% 增加适当留白
ax = gca;
ax.Position = [0.15 0.17 0.78 0.76];

% 网格线（可选）
grid on
output_filename = 'LVO2E_vs_HR.png';
exportgraphics(fig, output_filename, 'Resolution', 300);

% 关闭 figure
close(fig);

RawT.LabelLVO2E = PredictedLabelLVO2E;
%%
% 
% %%
% % patientNO = [];
% % for i = 1:height(SummaryTable)
% %     thisMVO2 = SummaryTable.PacingResults{i}.RVO2E_all;
% %     if any(thisMVO2(1:end)>8000)
% %        patientNO = [patientNO SummaryTable.PatientNo(i)];
% %     end
% % end
% %%
% % fig2 = figure('Visible','on');
% % clf; hold on;
% %
% % for i = 1:height(SummaryTable)
% %     HR = SummaryTable.PacingResults{i}.HR(:);
% %     thisO2E = SummaryTable.PacingResults{i}.LVO2E_all(:);
% %
% %     % Skip if missing or too short
% %     if numel(HR) < 1 || numel(thisO2E) < 1 || any(isnan(HR(1))) || any(isnan(thisO2E(1)))
% %         continue
% %     end
% %
% %     % Deltas relative to baseline point (index 1)
% %     dHR   = HR - HR(1);
% %     dO2E  = thisO2E - thisO2E(1);
% %
% %     if PredictedLabelLVO2E(i) == 1
% %         plot(dHR, dO2E, '-', 'Color', [0.1 0.7 0.1 0.5], 'LineWidth', 1.2)
% %     elseif PredictedLabelLVO2E(i) == 0
% %         plot(dHR, dO2E, '-', 'Color', [0.9 0.1 0.1 0.5], 'LineWidth', 1.2)
% %     else
% %         plot(dHR, dO2E, '-', 'Color', 'k', 'LineWidth', 1.2)
% %     end
% % end
% %
% % % Axes & styling
% % xlabel('\DeltaHR (beats/min)', 'FontSize', 12, 'FontName','Arial')
% % ylabel('\DeltaCE', 'FontSize', 12, 'FontName','Arial')
% %
% % set(gca, ...
% %     'FontSize', 12, ...
% %     'FontName', 'Arial', ...
% %     'LineWidth', 1, ...
% %     'TickDir', 'out', ...
% %     'TickLength', [0.018 0.018])
% % box off;
% % pbaspect([1 1 1]);
% % grid off
% %
% % % Legend with dummy handles
% % h1 = plot(nan, nan, '-', 'Color', [0.1 0.7 0.1], 'LineWidth', 1.2);
% % h2 = plot(nan, nan, '-', 'Color', [0.9 0.1 0.1], 'LineWidth', 1.2);
% % legend([h1 h2], {'Responders', 'Non-responders'}, ...
% %     'Location', 'northwest', 'Box', 'off', 'FontSize', 12, 'FontName','Arial');
% %
% % % Layout
% % % ax = gca;
% % % ax.Position = [0.2 0.2 0.72 0.72];  % 保持正方形比例
% %
% % % 设置导出图大小：2.5 × 2.5 inch
% % set(fig2, 'Units', 'inches', 'Position', [1 1 3 2.5])
% %
% % % Export
% % output_filename2 = 'Delta_LVO2E_vs_Delta_HR.png';
% % exportgraphics(fig2, output_filename2, 'Resolution', 300);
% 
% 
% 
% %% remake figure based on classifiers
% 6 Based on O2 effiency' Dan's creation
type = zeros(height(SummaryTable),1);
patientN = [];
for i = 1:height(SummaryTable)
    PATIENT_NO = SummaryTable.PatientNo(i);
    % if i<= 128
    % load(sprintf('SimsUMFinal/P_NO%d.mat',PATIENT_NO)); % load results from great lake clusters
    % else
    %     load(sprintf('SimsUWFinal/P_NO%d.mat',PATIENT_NO)); % load results from great lake clusters
    % end
    % rho_myo = 1.055;
    % LV_m = rho_myo * (output.params.Vw_LV+output.params.Vw_SEP);
    % RV_m = rho_myo * (output.params.Vw_RV);
    % thisO2E = SummaryTable.PacingResults{i}.CO./(SummaryTable.PacingResults{i}.MVO2_RV./100*RV_m);
    thisO2E = SummaryTable.PacingResults{i}.LVO2E_all;
    % deltaO2E = abs(diff(thisO2E(~(isnan(thisO2E)))));
    % if any(deltaO2E>7.*median(abs(deltaO2E))) && any(deltaO2E>50)
    %     patientN = [patientN;PATIENT_NO];
    % end
    if any(isnan(thisO2E))
        error('find this patient')
    end
    if any(thisO2E(2:end)>thisO2E(1))
        type(i) = 1;
    end
end
% PredictedLabelRVO2E = type;
% fig = figure('Visible','on');
% clf; hold on;
% 
% for i  = 1:height(SummaryTable)
%     % PATIENT_NO = SummaryTable.PatientNo(i);
%     % if i<= 128
%     %     load(sprintf('SimsUMFinal/P_NO%d.mat',PATIENT_NO)); % load results from great lake clusters
%     % else
%     %     load(sprintf('SimsUWFinal/P_NO%d.mat',PATIENT_NO)); % load results from great lake clusters
%     % end
%     % rho_myo = 1.055;
%     % LV_m = rho_myo * (output.params.Vw_LV+output.params.Vw_SEP);
%     % RV_m = rho_myo * (output.params.Vw_RV);
%     % thisO2E = SummaryTable.PacingResults{i}.CO*1000./(SummaryTable.PacingResults{i}.MVO2_RV./100*RV_m);
%     thisO2E = SummaryTable.PacingResults{i}.RVO2E_all;
%     if PredictedLabelRVO2E(i) == 1
%         plot(SummaryTable.PacingResults{i}.HR,thisO2E, '-o', ...
%             'Color', [0.1 0.7 0.1 0.5], 'LineWidth', 1.8)
%     elseif PredictedLabelRVO2E(i) == 0
%         plot(SummaryTable.PacingResults{i}.HR, thisO2E, '-o', ...
%             'Color', [0.9 0.1 0.1 0.5], 'LineWidth', 1.8)
%     else
%         plot(SummaryTable.PacingResults{i}.HR, thisO2E, '-o', ...
%             'Color', 'k', 'LineWidth', 1.8)
%     end
% end
% 
% % 重要：设置坐标标签和美化
% xlabel('HR (beats/min)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica')
% ylabel('RV O2 effiency', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica')
% 
% set(gca, ...
%     'FontSize', 13, ...
%     'FontName', 'Helvetica', ...
%     'LineWidth', 1.2, ...
%     'TickDir', 'out', ...
%     'TickLength', [0.018 0.018])
% box off;
% pbaspect([1 1 1]); % 正方形比例
% 
% % 可选：加图例，让不同类别一目了然
% h1 = plot(nan, nan, '--', 'Color', [0.1 0.7 0.1], 'LineWidth', 1.8); % 用nan, nan画假线用于legend
% h2 = plot(nan, nan, '--', 'Color', [0.9 0.1 0.1], 'LineWidth', 1.8);
% legend([h1 h2], {'Responders', 'Non-responders'}, ...
%     'Location', 'northeast', 'Box', 'off', 'FontSize', 12);
% 
% % 增加适当留白
% ax = gca;
% ax.Position = [0.15 0.17 0.78 0.76];
% 
% % 网格线（可选）
% grid on
% output_filename = 'RVO2E_vs_HR.png';
% exportgraphics(fig, output_filename, 'Resolution', 300);
% 
% 
% % 关闭 figure
% close(fig);
% 
% RawT.LabelRVO2E = PredictedLabelRVO2E;
% 
% 
% %% remake figure based on classifiers
% % 7 Based on MVO2
% type = zeros(height(SummaryTable),1);
% patientN = [];
% for i = 1:height(SummaryTable)
%     thisMVO2 = SummaryTable.PacingResults{i}.MVO2_LV;
%     if any(thisMVO2(2:end)<thisMVO2(1))
%         type(i) = 1;
%     end
%     % if range(thisMVO2)>10
%     %    patientN = [patientN i];
%     % end
% end
% 
% % for i = 1:height(SummaryTable)
% %     thisMVO2 = SummaryTable.PacingResults{i}.MVO2_RV;
% %     if any(thisMVO2(2:end)<thisMVO2(1))
% %         type(i) = 1;
% %     end
% %     if range(thisMVO2)>10
% %        patientN = [patientN i];
% %     end
% % end
% % %
% % for i = 1:height(SummaryTable)
% %     thisME = SummaryTable.PacingResults{i}.CE_LV_all;
% %     if any(thisME(2:end)>thisME(1))
% %         type(i) = 1;
% %     end
% %     if any(thisME>0.65) || any(thisME<0.01)
% %        patientN = [patientN i];
% %     end
% % end
% %
% % for i = 1:height(SummaryTable)
% %     thisME = SummaryTable.PacingResults{i}.CE_RV_all;
% %     if any(thisME(2:end)>thisME(1))
% %         type(i) = 1;
% %     end
% %     % if any(thisME>0.65) || any(thisME<0.1)
% %     %    patientN = [patientN SummaryTable.PatientNo(i)];
% %     %    error('maywrong')
% %     % end
% % end
% 
% %%
% PredictedLabelLVMVO2 = type;
% 
% fig = figure('Visible','on');
% clf; hold on;
% 
% for i  = 1:height(SummaryTable)
%     thisMVO2 = SummaryTable.PacingResults{i}.MVO2_LV;
%     if PredictedLabelLVMVO2(i) == 1
%         plot(SummaryTable.PacingResults{i}.HR,thisMVO2, '-o', ...
%             'Color', [0.1 0.7 0.1 0.5], 'LineWidth', 1.8)
%     elseif PredictedLabelLVMVO2(i) == 0
%         plot(SummaryTable.PacingResults{i}.HR, thisMVO2, '-o', ...
%             'Color', [0.9 0.1 0.1 0.5], 'LineWidth', 1.8)
%     else
%         plot(SummaryTable.PacingResults{i}.HR, thisMVO2, '-o', ...
%             'Color', 'k', 'LineWidth', 1.8)
%     end
% end
% 
% % 重要：设置坐标标签和美化
% xlabel('HR (beats/min)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica')
% ylabel('LV MVO2', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica')
% 
% set(gca, ...
%     'FontSize', 13, ...
%     'FontName', 'Helvetica', ...
%     'LineWidth', 1.2, ...
%     'TickDir', 'out', ...
%     'TickLength', [0.018 0.018])
% box off;
% pbaspect([1 1 1]); % 正方形比例
% 
% % 可选：加图例，让不同类别一目了然
% h1 = plot(nan, nan, '--', 'Color', [0.1 0.7 0.1], 'LineWidth', 1.8); % 用nan, nan画假线用于legend
% h2 = plot(nan, nan, '--', 'Color', [0.9 0.1 0.1], 'LineWidth', 1.8);
% legend([h1 h2], {'Responders', 'Non-responders'}, ...
%     'Location', 'northeast', 'Box', 'off', 'FontSize', 12);
% 
% % 增加适当留白
% ax = gca;
% ax.Position = [0.15 0.17 0.78 0.76];
% 
% % 网格线（可选）
% grid on
% output_filename = 'LVMOV2_vs_HR.png';
% exportgraphics(fig, output_filename, 'Resolution', 300);
% 
% % 关闭 figure
% close(fig);
% 
% RawT.LabelLVMVO2 = PredictedLabelLVMVO2;
% 
% %% remake figure based on classifiers
% % 8 Based on MVO2
% type = zeros(height(SummaryTable),1);
% for i = 1:height(SummaryTable)
%     thisMVO2 = SummaryTable.PacingResults{i}.MVO2_RV;
%     if any(thisMVO2(2:end)<thisMVO2(1))
%         type(i) = 1;
%     end
% end
% 
% PredictedLabelRVMVO2 = type;
% 
% fig = figure('Visible','off');
% clf; hold on;
% 
% for i  = 1:height(SummaryTable)
%     thisMVO2 = SummaryTable.PacingResults{i}.MVO2_RV;
%     if PredictedLabelRVMVO2(i) == 1
%         plot(SummaryTable.PacingResults{i}.HR,thisMVO2, '-o', ...
%             'Color', [0.1 0.7 0.1 0.5], 'LineWidth', 1.8)
%     elseif PredictedLabelRVMVO2(i) == 0
%         plot(SummaryTable.PacingResults{i}.HR, thisMVO2, '-o', ...
%             'Color', [0.9 0.1 0.1 0.5], 'LineWidth', 1.8)
%     else
%         plot(SummaryTable.PacingResults{i}.HR, thisMVO2, '-o', ...
%             'Color', 'k', 'LineWidth', 1.8)
%     end
% end
% 
% % 重要：设置坐标标签和美化
% xlabel('HR (beats/min)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica')
% ylabel('RV MVO2', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica')
% 
% set(gca, ...
%     'FontSize', 13, ...
%     'FontName', 'Helvetica', ...
%     'LineWidth', 1.2, ...
%     'TickDir', 'out', ...
%     'TickLength', [0.018 0.018])
% box off;
% pbaspect([1 1 1]); % 正方形比例
% 
% % 可选：加图例，让不同类别一目了然
% h1 = plot(nan, nan, '--', 'Color', [0.1 0.7 0.1], 'LineWidth', 1.8); % 用nan, nan画假线用于legend
% h2 = plot(nan, nan, '--', 'Color', [0.9 0.1 0.1], 'LineWidth', 1.8);
% legend([h1 h2], {'Responders', 'Non-responders'}, ...
%     'Location', 'northeast', 'Box', 'off', 'FontSize', 12);
% 
% % 增加适当留白
% ax = gca;
% ax.Position = [0.15 0.17 0.78 0.76];
% 
% % 网格线（可选）
% grid on
% output_filename = 'RVMVO2_vs_HR.png';
% exportgraphics(fig, output_filename, 'Resolution', 300);
% 
% 
% % 关闭 figure
% close(fig);
% 
% RawT.LabelRVMVO2 = PredictedLabelRVMVO2;
% 
% %% remake figure based on classifiers
% % 9 Based on ME
% type = zeros(height(SummaryTable),1);
% for i = 1:height(SummaryTable)
%     thisME = SummaryTable.PacingResults{i}.CE_LV_all;
%     if any(thisME(2:end)>thisME(1))
%         type(i) = 1;
%     end
% end
% 
% PredictedLabelLVME = type;
% 
% fig = figure('Visible','off');
% clf; hold on;
% 
% for i  = 1:height(SummaryTable)
%     thisME = SummaryTable.PacingResults{i}.CE_LV_all;
%     if PredictedLabelLVME(i) == 1
%         plot(SummaryTable.PacingResults{i}.HR,thisME, '-o', ...
%             'Color', [0.1 0.7 0.1 0.5], 'LineWidth', 1.8)
%     elseif PredictedLabelLVME(i) == 0
%         plot(SummaryTable.PacingResults{i}.HR, thisME, '-o', ...
%             'Color', [0.9 0.1 0.1 0.5], 'LineWidth', 1.8)
%     else
%         plot(SummaryTable.PacingResults{i}.HR, thisME, '-o', ...
%             'Color', 'k', 'LineWidth', 1.8)
%     end
% end
% 
% % 重要：设置坐标标签和美化
% xlabel('HR (beats/min)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica')
% ylabel('LV ME', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica')
% 
% set(gca, ...
%     'FontSize', 13, ...
%     'FontName', 'Helvetica', ...
%     'LineWidth', 1.2, ...
%     'TickDir', 'out', ...
%     'TickLength', [0.018 0.018])
% box off;
% pbaspect([1 1 1]); % 正方形比例
% 
% % 可选：加图例，让不同类别一目了然
% h1 = plot(nan, nan, '--', 'Color', [0.1 0.7 0.1], 'LineWidth', 1.8); % 用nan, nan画假线用于legend
% h2 = plot(nan, nan, '--', 'Color', [0.9 0.1 0.1], 'LineWidth', 1.8);
% legend([h1 h2], {'Responders', 'Non-responders'}, ...
%     'Location', 'northeast', 'Box', 'off', 'FontSize', 12);
% 
% % 增加适当留白
% ax = gca;
% ax.Position = [0.15 0.17 0.78 0.76];
% 
% % 网格线（可选）
% grid on
% output_filename = 'LVME_vs_HR.png';
% exportgraphics(fig, output_filename, 'Resolution', 300);
% 
% % 关闭 figure
% close(fig);
% 
% RawT.LabelLVME = PredictedLabelLVME;
% 
% %% remake figure based on classifiers
% % 10 Based on ME
% type = zeros(height(SummaryTable),1);
% for i = 1:height(SummaryTable)
%     thisME = SummaryTable.PacingResults{i}.CE_RV_all;
%     % if any(thisME(2:end)>1)
%     %    error('find this one')
%     % end
%     if any(thisME(2:end)>thisME(1))
%         type(i) = 1;
%     end
% end
% PredictedLabelRVME = type;
% 
% fig = figure('Visible','off');
% clf; hold on;
% 
% for i  = 1:height(SummaryTable)
%     thisME = SummaryTable.PacingResults{i}.CE_RV_all;
%     if PredictedLabelRVME(i) == 1
%         plot(SummaryTable.PacingResults{i}.HR,thisME, '-o', ...
%             'Color', [0.1 0.7 0.1 0.5], 'LineWidth', 1.8)
%     elseif PredictedLabelRVME(i) == 0
%         plot(SummaryTable.PacingResults{i}.HR, thisME, '-o', ...
%             'Color', [0.9 0.1 0.1 0.5], 'LineWidth', 1.8)
%     else
%         plot(SummaryTable.PacingResults{i}.HR, thisME, '-o', ...
%             'Color', 'k', 'LineWidth', 1.8)
%     end
% end
% 
% % 重要：设置坐标标签和美化
% xlabel('HR (beats/min)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica')
% ylabel('RV ME', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica')
% 
% set(gca, ...
%     'FontSize', 13, ...
%     'FontName', 'Helvetica', ...
%     'LineWidth', 1.2, ...
%     'TickDir', 'out', ...
%     'TickLength', [0.018 0.018])
% box off;
% pbaspect([1 1 1]); % 正方形比例
% 
% % 可选：加图例，让不同类别一目了然
% h1 = plot(nan, nan, '--', 'Color', [0.1 0.7 0.1], 'LineWidth', 1.8); % 用nan, nan画假线用于legend
% h2 = plot(nan, nan, '--', 'Color', [0.9 0.1 0.1], 'LineWidth', 1.8);
% legend([h1 h2], {'Responders', 'Non-responders'}, ...
%     'Location', 'northeast', 'Box', 'off', 'FontSize', 12);
% 
% % 增加适当留白
% ax = gca;
% ax.Position = [0.15 0.17 0.78 0.76];
% 
% % 网格线（可选）
% grid on
% 
% RawT.LabelRVME = PredictedLabelRVME;
% output_filename = 'RVME_vs_HR.png';
% exportgraphics(fig, output_filename, 'Resolution', 300);
% 
% % 关闭 figure
% close(fig);

%% 规则：baseline 取该曲线第一个点的 LVSP；只要有任意点 < baseline（可选容差 tol）→ 标 1，否则标 0

deltaLVSP = zeros(height(SummaryTable),1);  % 若后面不用可删
PredictedLabelLVSP = zeros(height(SummaryTable),1);

tol = 6;  % 如担心数值抖动，可设 tol=0.5（mmHg）等
for i = 1:height(SummaryTable)
    thisLVSP = SummaryTable.PacingResults{i}.meanLVSP;

    if ~all(isreal(thisLVSP))
        error('Patient %d has invalid LVSP (complex/NaN).', i);
    end

    % 可选：保留一个"变化幅度"的量
    deltaLVSP(i) = range(thisLVSP);
    
    % baseline 取第一个点（如你的 baseline 不在第一个点，请把 baseline_idx 改成相应索引）
    baseline_val = thisLVSP(1);

    % 判定：是否存在任何点低于 baseline（带或不带容差）
    if min(thisLVSP) < baseline_val - tol
        PredictedLabelLVSP(i) = 1;   % 有下降 → responder
    else
        PredictedLabelLVSP(i) = 0;   % 全部不低于 baseline → non-responder
    end
end

% fig = figure('Visible','on');
% clf; hold on;
% 
% for i  = 1:height(SummaryTable)
%     if PredictedLabelLVSP(i) == 1
%         plot(SummaryTable.PacingResults{i}.HR, SummaryTable.PacingResults{i}.meanLVSP, '-o', ...
%             'Color', [0.1 0.7 0.1 0.5], 'LineWidth', 1.8)
%     elseif PredictedLabelLVSP(i) == 0
%         plot(SummaryTable.PacingResults{i}.HR, SummaryTable.PacingResults{i}.meanLVSP, '-o', ...
%             'Color', [0.9 0.1 0.1 0.5], 'LineWidth', 1.8)
%     else
%         plot(SummaryTable.PacingResults{i}.HR, SummaryTable.PacingResults{i}.meanLVSP, '-o', ...
%             'Color', 'k', 'LineWidth', 1.8)
%     end
% end
% 
% % fig = figure('Visible','on'); clf;
% % 
% % % 更现代 & 好控版式：tiledlayout
% % tl = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');
% % 
% % % 预创建四个坐标轴，并加上分面标题
% % ax(1) = nexttile; title(ax(1), 'O2E↑ & ME↑');
% % ax(2) = nexttile; title(ax(2), 'O2E↑ & ME↓');
% % ax(3) = nexttile; title(ax(3), 'O2E↓ & ME↓');
% % ax(4) = nexttile; title(ax(4), 'O2E↓ & ME↑');
% % 
% % % 统一 hold on
% % for k = 1:4, hold(ax(k), 'on'); end
% % 
% % % 颜色（RGB）
% % col_g = [0.1 0.7 0.1];
% % col_r = [0.9 0.1 0.1];
% % col_k = [0 0 0];
% % col_b = [0 0 1];
% % 
% % % 按标签把曲线画到各自的子图
% % for i = 1:height(SummaryTable)
% %     % 小心防守：缺列就跳过
% %     Ti = SummaryTable.PacingResults{i};
% %     if ~all(ismember({'HR','meanLVSP'}, Ti.Properties.VariableNames)), continue; end
% % 
% %     if RawT.LabelLVO2E(i)==1 && RawT.LabelLVME(i)==1
% %         targetAx = ax(1); c = col_g;
% %     elseif RawT.LabelLVO2E(i)==1 && RawT.LabelLVME(i)==0
% %         targetAx = ax(2); c = col_r;
% %     elseif RawT.LabelLVO2E(i)==0 && RawT.LabelLVME(i)==0
% %         targetAx = ax(3); c = col_k;
% %     else % RawT.LabelLVO2E(i)==0 && RawT.LabelLVME(i)==1
% %         targetAx = ax(4); c = col_b;
% %     end
% % 
% %     plot(targetAx, Ti.HR, Ti.meanLVSP, '-o', 'Color', c, 'LineWidth', 1.8, ...
% %         'MarkerSize', 4, 'MarkerFaceColor', 'w');  % 适度强调点
% % end
% % 
% % % —— 统一美化到每个子图 ——
% % for k = 1:4
% %     set(ax(k), 'FontSize', 13, 'FontName','Helvetica', ...
% %         'LineWidth', 1.2, 'TickDir', 'out', 'TickLength', [0.018 0.018]);
% %     grid(ax(k), 'on');
% %     box(ax(k), 'off');
% %     pbaspect(ax(k), [1 1 1]);
% %     ylim(ax(k), [60 180]);
% % end
% 
% % % 全局坐标轴标签（推荐这种，不会只作用于某一个子图）
% % xlabel(tl, 'HR (beats/min)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica');
% % ylabel(tl, 'mean LV pressure during ejection (mmHg)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica');
% %
% % % 全局图例：在第一个轴上放“假线”，然后把图例停靠到布局下方
% % h1 = plot(ax(1), nan, nan, '-', 'Color', col_g, 'LineWidth', 1.8);
% % h2 = plot(ax(1), nan, nan, '-', 'Color', col_r, 'LineWidth', 1.8);
% % h3 = plot(ax(1), nan, nan, '-', 'Color', col_k, 'LineWidth', 1.8);
% % h4 = plot(ax(1), nan, nan, '-', 'Color', col_b, 'LineWidth', 1.8);
% % lgd = legend(ax(1), [h1 h2 h3 h4], ...
% %     {'all increase', 'O2E up and ME down', 'all decrease', 'O2E down and ME up'}, ...
% %     'Box','off', 'FontSize', 12);
% % lgd.Layout.Tile = 'south';  % 全局图例放在最底部
% 
% % 重要：设置坐标标签和美化
% xlabel('HR (beats/min)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica')
% ylabel('mLVPs (mmHg)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica')
% 
% set(gca, ...
%     'FontSize', 13, ...
%     'FontName', 'Helvetica', ...
%     'LineWidth', 1.2, ...
%     'TickDir', 'out', ...
%     'TickLength', [0.018 0.018])
% box off;
% pbaspect([1 1 1]); % 正方形比例
% 
% % 可选：加图例，让不同类别一目了然
% h1 = plot(nan, nan, '--', 'Color', [0.1 0.7 0.1], 'LineWidth', 1.8); % 用nan, nan画假线用于legend
% h2 = plot(nan, nan, '--', 'Color', [0.9 0.1 0.1], 'LineWidth', 1.8);
% legend([h1 h2], {'Responders', 'Non-responders'}, ...
%     'Location', 'northeast', 'Box', 'off', 'FontSize', 12);
% 
% % 增加适当留白
% ax = gca;
% ax.Position = [0.15 0.17 0.78 0.76];
% 
% % 网格线（可选）
% grid on
% 
% % 导出
% exportgraphics(fig, 'LVSP_vs_HR.png', 'Resolution', 300);
RawT.LabelLVSP = PredictedLabelLVSP;

%%
% —— 统计与打标签（基于 meanRVSP）——
deltaRVSP = zeros(height(SummaryTable),1);
PredictedLabelRVSP = zeros(height(SummaryTable),1);

tol = 4.5;   % mmHg
for i = 1:height(SummaryTable)
    % 保护：缺列就跳过
    if isempty(SummaryTable.PacingResults{i}) || ...
       ~all(ismember({'HR','meanRVSP'}, SummaryTable.PacingResults{i}.Properties.VariableNames))
        PredictedLabelRVSP(i) = NaN;
        deltaRVSP(i) = NaN;
        continue
    end

    thisRVSP = SummaryTable.PacingResults{i}.meanRVSP;

    if ~all(isreal(thisRVSP) | isnan(thisRVSP))
        error('Patient %d has invalid RVSP (complex).', i);
    end

    % 变化幅度（忽略 NaN）
    deltaRVSP(i) = range(thisRVSP,'omitnan');

    % baseline 取“第一个非 NaN 值”
    bidx = find(~isnan(thisRVSP), 1, 'first');
    if isempty(bidx)
        PredictedLabelRVSP(i) = NaN;   % 全 NaN，无从判断
        continue
    end
    baseline_val = thisRVSP(bidx);

    % 判定：是否存在任何点低于 baseline - tol（忽略 NaN）
    if min(thisRVSP,[],'omitnan') < baseline_val - tol
        PredictedLabelRVSP(i) = 1;   % responder（下降）
    else
        PredictedLabelRVSP(i) = 0;   % non-responder（未下降）
    end
end

% % —— 画图（按 responder / non-responder 分色）——
% fig = figure('Visible','on');
% clf; hold on;
% 
% for i = 1:height(SummaryTable)
%     Ti = SummaryTable.PacingResults{i};
%     if isempty(Ti) || ~all(ismember({'HR','meanRVSP'}, Ti.Properties.VariableNames))
%         continue
%     end
% 
%     if PredictedLabelRVSP(i) == 1
%         plot(Ti.HR, Ti.meanRVSP, '-o', 'Color', [0.1 0.7 0.1 0.5], 'LineWidth', 1.8)
%     elseif PredictedLabelRVSP(i) == 0
%         plot(Ti.HR, Ti.meanRVSP, '-o', 'Color', [0.9 0.1 0.1 0.5], 'LineWidth', 1.8)
%     else
%         plot(Ti.HR, Ti.meanRVSP, '-o', 'Color', 'k', 'LineWidth', 1.8)
%     end
% end
% 
% xlabel('HR (beats/min)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica');
% ylabel('mRVPs (mmHg)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica');
% 
% set(gca, 'FontSize', 13, 'FontName','Helvetica', 'LineWidth', 1.2, ...
%          'TickDir','out', 'TickLength',[0.018 0.018]);
% box off; pbaspect([1 1 1]);
% grid on
% 
% % 图例
% h1 = plot(nan, nan, '--', 'Color', [0.1 0.7 0.1], 'LineWidth', 1.8);
% h2 = plot(nan, nan, '--', 'Color', [0.9 0.1 0.1], 'LineWidth', 1.8);
% legend([h1 h2], {'Responders', 'Non-responders'}, ...
%        'Location','northeast', 'Box','off', 'FontSize',12);
% 
% ax = gca; ax.Position = [0.15 0.17 0.78 0.76];
% 
% % 导出
% exportgraphics(fig, 'RVSP_vs_HR.png', 'Resolution', 300);

% 存回标签（如需要）
RawT.LabelRVSP = PredictedLabelRVSP;
%%
% —— 统计与打标签（基于 SBP）——
deltaSBP = zeros(height(SummaryTable),1);
PredictedLabelSBP = zeros(height(SummaryTable),1);

tol = 8.5;   % mmHg
for i = 1:height(SummaryTable)
    % 保护：缺列就跳过
    if isempty(SummaryTable.PacingResults{i}) || ...
       ~ismember('SBP', SummaryTable.PacingResults{i}.Properties.VariableNames)
        PredictedLabelSBP(i) = NaN;
        deltaSBP(i) = NaN;
        continue
    end

    thisSBP = SummaryTable.PacingResults{i}.SBP;

    if ~all(isreal(thisSBP) | isnan(thisSBP))
        error('Patient %d has invalid SBP (complex).', i);
    end

    % 变化幅度（忽略 NaN）
    deltaSBP(i) = range(thisSBP,'omitnan');

    % baseline 取第一个非 NaN
    bidx = find(~isnan(thisSBP), 1, 'first');
    if isempty(bidx)
        PredictedLabelSBP(i) = NaN;
        continue
    end
    baseline_val = thisSBP(bidx);

    % 判定：是否有点 < baseline - tol
    if min(thisSBP,[],'omitnan') < baseline_val - tol
        PredictedLabelSBP(i) = 1;   % responder（下降）
    else
        PredictedLabelSBP(i) = 0;   % non-responder（未下降）
    end
end

% % —— 绘图 —— 
% fig = figure('Visible','on'); clf; hold on;
% for i = 1:height(SummaryTable)
%     Ti = SummaryTable.PacingResults{i};
%     if isempty(Ti) || ~ismember('SBP', Ti.Properties.VariableNames), continue; end
% 
%     if PredictedLabelSBP(i) == 1
%         plot(Ti.HR, Ti.SBP, '-o', 'Color', [0.1 0.7 0.1 0.5], 'LineWidth', 1.8)
%     elseif PredictedLabelSBP(i) == 0
%         plot(Ti.HR, Ti.SBP, '-o', 'Color', [0.9 0.1 0.1 0.5], 'LineWidth', 1.8)
%     else
%         plot(Ti.HR, Ti.SBP, '-o', 'Color', 'k', 'LineWidth', 1.8)
%     end
% end
% 
% xlabel('HR (beats/min)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica');
% ylabel('SBP (mmHg)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName','Helvetica');
% 
% set(gca, 'FontSize', 13, 'FontName','Helvetica', 'LineWidth', 1.2, ...
%          'TickDir','out', 'TickLength',[0.018 0.018]);
% box off; pbaspect([1 1 1]);
% grid on
% 
% % 图例
% h1 = plot(nan, nan, '--', 'Color', [0.1 0.7 0.1], 'LineWidth', 1.8);
% h2 = plot(nan, nan, '--', 'Color', [0.9 0.1 0.1], 'LineWidth', 1.8);
% legend([h1 h2], {'Responders', 'Non-responders'}, ...
%        'Location','northeast', 'Box','off', 'FontSize',12);
% 
% % 导出
% exportgraphics(fig, 'SBP_vs_HR.png', 'Resolution', 300);

% 存回标签
RawT.LabelSBP = PredictedLabelSBP;

%%
writetable(RawT,'inputFGclassifierWithLabel1013.CSV')