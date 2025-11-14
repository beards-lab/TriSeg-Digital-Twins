%% Simulation for heart failure patients from Umich cohort without CMR info
clear
Order = [1	2 13 15	16	17	20	22	25	26	35	37	42	45	46	48	50	55	57	60	62	66	67	77	78	79	84	87	88	90	91	92	93	94	97	98	99	102	103	104	106	108	111	112	113	118	122	126	130	131	132	133	137	141	146	147	150	151	164	167	168	173	177	178	186	188	195	196	197	199	202	203	204	205	206	207	209	212	213	214	217	218	222	223	224	230	232	234	235	236	240	244	246	248	250	252	267	269	271	274	281	282	283	284	293	294	300	302	303	304	306	309	316	318	324	325	327	329	335	341	342	343	344	345	346	350	351	368];
%%

MRI_flag = 1;
for GENDER  = 2 % 1 for male and other for female
    if GENDER == 1
        [targets, inputs, mods] = targetVals_male();
        modifiers = readmatrix('modifiers_male.csv','range','2:2');
        % modifiers = ones(1,length(mods));
        [INIparams, INIinit] = estiminiParams(targets,inputs);
        [params, init] = optParams(INIparams, INIinit, mods,modifiers);
    else
        [targets, inputs, mods] = targetVals_female();
        modifiers = readmatrix('modifiers_female.csv','range','2:2');
        [INIparams, INIinit] = estiminiParams(targets,inputs);
        [params, init] = optParams(INIparams, INIinit, mods,modifiers);
    end
    Sim = 1; % 1 for simulation and other for optimazation
    if Sim ==1
        runSimOnGL;
        Tnormal = t;
        P_LV_normal = P_LV;
        V_LV_normal = V_LV;
        P_LA_normal = P_LA;
        P_SA_normal = P_SA;
        % NplotSrd; % 4-panel figure just for canonical subjects
        % GetMovie; % TriSeg model: displacement and stress as functions of time
        % See_TriSeg; % slices of GetMoive.m
    else
        Srdopt; % optimize modifiers of canonical subjects
    end
end
TNnormal = Tnormal./Tnormal(round(end/2));
figure(77);clf;

% plot(THFpEF,P_LV_HFpEF,'color',[0.6350 0.0780 0.1840], 'Linewidth', 2);hold on;
% plot(THFpEF,P_SA_HFpEF,'color',[0.3010 0.7450 0.9330], 'Linewidth', 2);hold on;
% plot(THFpEF,P_LA_HFpEF,'color',[0.8500 0.5250 0.0980], 'Linewidth', 2);hold on;
plot(TNnormal,P_LV_normal,'color',[0.3010 0.7450 0.9330], 'Linewidth', 2);hold on;
plot(TNnormal,P_LA_normal,'color',[0.6350 0.0780 0.1840], 'Linewidth', 2);
% plot(TNnormal,P_LV_normal,'color',[0.6350 0.0780 0.1840], 'Linewidth', 2)
trng = [0.2 1.2];
plot(trng, [targets.SBP, targets.SBP],'color',[0.6350 0.0780 0.1840],'LineStyle','--','Linewidth', 2);
plot(trng, [targets.DBP, targets.DBP],'color',[0.6350 0.0780 0.1840],'LineStyle','--','Linewidth', 2);
plot(trng, [targets.PCWP, targets.PCWP],'color',[0.6350 0.0780 0.1840],'LineStyle','--','Linewidth', 2);
% xlabel('T (sec)', 'FontName', 'Arial', 'FontSize', 16);
xlim(trng);
ylim([0 20]);
% ylabel('Flow (ml/min)', 'FontName', 'Arial', 'FontSize', 16);

set(gca, 'XTickLabel', [], 'YTickLabel', []);

% 统一字体设置
set(gca, 'FontName', 'Arial', 'FontSize', 16);
box on;
w_in = 2.4; 
h_in = 1.2; 
dpi  = 300;

set(gcf, 'Units', 'inches', 'Position', [1 1 w_in h_in]);
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 w_in h_in]);
set(gcf, 'InvertHardcopy', 'off');  % 防止背景反转
pbaspect([2 1 1]);  % 保持2:1比例

outpath = fullfile('C:', 'Users', 'fenggu', ...
    'University of Michigan Dropbox', 'Feng Gu', ...
    'DTPacingHFpEF', 'NormalLAP.png');
exportgraphics(gcf, outpath, 'Resolution', dpi, ...
    'ContentType', 'image', 'BackgroundColor', 'white');
% TNnormal = Tnormal./Tnormal(round(end/2));
% plot(TNnormal,P_LV_normal,TNnormal,P_SV_normal);

figure(777);clf;

% plot(THFpEF,P_LV_HFpEF,'color',[0.6350 0.0780 0.1840], 'Linewidth', 2);hold on;
% plot(THFpEF,P_SA_HFpEF,'color',[0.3010 0.7450 0.9330], 'Linewidth', 2);hold on;
% plot(THFpEF,P_LA_HFpEF,'color',[0.8500 0.5250 0.0980], 'Linewidth', 2);hold on;
plot(TNnormal,P_LV_normal,'color',[0.3010 0.7450 0.9330], 'Linewidth', 2);hold on;
plot(TNnormal,P_SA_normal,'color',[0.6350 0.0780 0.1840], 'Linewidth', 2);
% plot(TNnormal,P_LV_normal,'color',[0.6350 0.0780 0.1840], 'Linewidth', 2)
trng = [0.8 1.8];
plot(trng, [targets.SBP, targets.SBP],'color',[0.6350 0.0780 0.1840],'LineStyle','--','Linewidth', 2);
plot(trng, [targets.DBP, targets.DBP],'color',[0.6350 0.0780 0.1840],'LineStyle','--','Linewidth', 2);
plot(trng, [targets.PCWP, targets.PCWP],'color',[0.6350 0.0780 0.1840],'LineStyle','--','Linewidth', 2);
% xlabel('T (sec)', 'FontName', 'Arial', 'FontSize', 16);
xlim(trng);
ylim([50 130]);
% ylabel('Flow (ml/min)', 'FontName', 'Arial', 'FontSize', 16);

set(gca, 'XTickLabel', [], 'YTickLabel', []);

% 统一字体设置
set(gca, 'FontName', 'Arial', 'FontSize', 16);
box on;
w_in = 2.4; 
h_in = 1.2; 
dpi  = 900;

set(gcf, 'Units', 'inches', 'Position', [1 1 w_in h_in]);
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 w_in h_in]);
set(gcf, 'InvertHardcopy', 'off');  % 防止背景反转
pbaspect([2 1 1]);  % 保持2:1比例

outpath = fullfile('C:', 'Users', 'fenggu', ...
    'University of Michigan Dropbox', 'Feng Gu', ...
    'DTPacingHFpEF', 'NormalSP.png');
exportgraphics(gcf, outpath, 'Resolution', dpi, ...
    'ContentType', 'image', 'BackgroundColor', 'white');
%%
% [~,ReprePatient] = find(Order==22);
% for NO = 223
for NO = 19
    % load(sprintf('SimsUMFinal/P_NO%d.mat',NO))
    load(sprintf('SimsUWFinal/P_NO%d.mat',NO))
    targets = output.targets;
    inputs = output.inputs;
    params = output.params;
    init = output.init;
    runSimOnGL;
    THFpEF = t;
    P_LV_HFpEF = P_LV;
    V_LV_HFpEF = V_LV;
    P_LA_HFpEF = P_LA;
    P_SA_HFpEF = P_SA;
end
TNHFpEF = THFpEF./THFpEF(round(end/2));

figure(88);clf;
% plot(THFpEF,P_LV_HFpEF,'color',[0.6350 0.0780 0.1840], 'Linewidth', 2);hold on;
% plot(THFpEF,P_SA_HFpEF,'color',[0.3010 0.7450 0.9330], 'Linewidth', 2);hold on;
% plot(THFpEF,P_LA_HFpEF,'color',[0.8500 0.5250 0.0980], 'Linewidth', 2);hold on;
plot(TNHFpEF,P_LV_HFpEF,'color',[0.3010 0.7450 0.9330], 'Linewidth', 2);hold on;
plot(TNHFpEF,P_LA_HFpEF,'color',[0.6350 0.0780 0.1840], 'Linewidth', 2);
% plot(TNnormal,P_LV_normal,'color',[0.6350 0.0780 0.1840], 'Linewidth', 2)
trng = [0.2 1.2];
plot(trng, [targets.SBP, targets.SBP],'color',[0.6350 0.0780 0.1840],'LineStyle','--','Linewidth', 2);
plot(trng, [targets.DBP, targets.DBP],'color',[0.6350 0.0780 0.1840],'LineStyle','--','Linewidth', 2);
plot(trng, [targets.PCWP, targets.PCWP],'color',[0.6350 0.0780 0.1840],'LineStyle','--','Linewidth', 2);
% xlabel('T (sec)', 'FontName', 'Arial', 'FontSize', 16);
xlim(trng);
ylim([5 25]);
% ylabel('Flow (ml/min)', 'FontName', 'Arial', 'FontSize', 16);
set(gca, 'XTickLabel', [], 'YTickLabel', []);

% 统一字体设置
set(gca, 'FontName', 'Arial', 'FontSize', 16);
box on;
w_in = 2.4; 
h_in = 1.2; 
dpi  = 300;

set(gcf, 'Units', 'inches', 'Position', [1 1 w_in h_in]);
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 w_in h_in]);
set(gcf, 'InvertHardcopy', 'off');  % 防止背景反转
pbaspect([2 1 1]);  % 保持2:1比例

outpath = fullfile('C:', 'Users', 'fenggu', ...
    'University of Michigan Dropbox', 'Feng Gu', ...
    'DTPacingHFpEF', 'HFpEFLAP.png');
exportgraphics(gcf, outpath, 'Resolution', dpi, ...
    'ContentType', 'image', 'BackgroundColor', 'white');

figure(888);clf;
% plot(THFpEF,P_LV_HFpEF,'color',[0.6350 0.0780 0.1840], 'Linewidth', 2);hold on;
% plot(THFpEF,P_SA_HFpEF,'color',[0.3010 0.7450 0.9330], 'Linewidth', 2);hold on;
% plot(THFpEF,P_LA_HFpEF,'color',[0.8500 0.5250 0.0980], 'Linewidth', 2);hold on;
plot(TNHFpEF,P_LV_HFpEF,'color',[0.3010 0.7450 0.9330], 'Linewidth', 2);hold on;
plot(TNHFpEF,P_SA_HFpEF,'color',[0.6350 0.0780 0.1840], 'Linewidth', 2);
% plot(TNnormal,P_LV_normal,'color',[0.6350 0.0780 0.1840], 'Linewidth', 2)
trng = [0.8 1.8];
plot(trng, [targets.SBP, targets.SBP],'color',[0.6350 0.0780 0.1840],'LineStyle','--','Linewidth', 2);
plot(trng, [targets.DBP, targets.DBP],'color',[0.6350 0.0780 0.1840],'LineStyle','--','Linewidth', 2);
plot(trng, [targets.PCWP, targets.PCWP],'color',[0.6350 0.0780 0.1840],'LineStyle','--','Linewidth', 2);
% xlabel('T (sec)', 'FontName', 'Arial', 'FontSize', 16);
xlim(trng);
ylim([50 130]);
% ylabel('Flow (ml/min)', 'FontName', 'Arial', 'FontSize', 16);

set(gca, 'XTickLabel', [], 'YTickLabel', []);

% 统一字体设置
set(gca, 'FontName', 'Arial', 'FontSize', 16);
box on;
w_in = 2.4; 
h_in = 1.2; 
dpi  = 300;

set(gcf, 'Units', 'inches', 'Position', [1 1 w_in h_in]);
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 w_in h_in]);
set(gcf, 'InvertHardcopy', 'off');  % 防止背景反转
pbaspect([2 1 1]);  % 保持2:1比例

outpath = fullfile('C:', 'Users', 'fenggu', ...
    'University of Michigan Dropbox', 'Feng Gu', ...
    'DTPacingHFpEF', 'HFpEFSP.png');
exportgraphics(gcf, outpath, 'Resolution', dpi, ...
    'ContentType', 'image', 'BackgroundColor', 'white');

% figure(99);clf;
% 
% plot(V_LV_normal(1:round(end/2)),P_LV_normal(1:round(end/2)),'color','b', 'Linewidth', 2);hold on;
% plot(V_LV_HFpEF(1:round(end/2)),P_LV_HFpEF(1:round(end/2)),'color','r', 'Linewidth', 2);hold on;
% % xlabel('Volume (ml)', 'FontName', 'Arial', 'FontSize', 16);
% % ylabel('Pressure (mmHg)', 'FontName', 'Arial', 'FontSize', 16);
% xlim([0 200]);
% ylim([0 150]);
% % 统一字体设置
% set(gca, 'FontName', 'Arial', 'FontSize', 16);
% set(gca, 'XTickLabel', [], 'YTickLabel', []);
% box on;
% pbaspect([1,1,1])
% % 设置图像尺寸（单位为英寸，1 cm ≈ 0.3937 in）
% set(gcf, 'Units', 'inches', 'Position', [1 1 2 2]);  % 9 cm x 9 cm ≈ 3.54 in x 3.54 in
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3.54 3.54]);
% 
% outpath = fullfile('C:', 'Users', 'fenggu', ...
%     'University of Michigan Dropbox', 'Feng Gu', 'DTPacingHFpEF', 'PVloopHealthHFpEF.png');
% exportgraphics(gcf, outpath, 'Resolution', 300);

%%
T1 = load('pacingSUM.mat');
T2 = load('pacingSUMUW.mat');
SummaryTable = [T1.SummaryTable;T2.SummaryTable];

figure(201); 
clf
% 取出第 NO 个病人
HFpEFexample = SummaryTable(134,:);

% 拆出数据
HR = HFpEFexample.PacingResults{1,1}.HR;
HR = HR-HR(1);
LAP = HFpEFexample.PacingResults{1,1}.LVO2E_all;
LAP = LAP-LAP(1);
CO = HFpEFexample.PacingResults{1,1}.CO;
CO = CO-CO(1);
% 画图（左轴LAP，右轴CO）
yyaxis left
hold on
plot(HR, LAP, '-o', 'LineWidth', 2, 'Color', [0.8500 0.5250 0.0980]);  % 深红色
% ylabel('LAP (mmHg)', 'FontName', 'Arial', 'FontSize', 9,'Color',[0.8500 0.5250 0.0980]);
ax = gca;
ax.YColor = [0.8500 0.5250 0.0980];
set(gca, 'XTickLabel', [], 'YTickLabel', []);
ylim([-50 50])

yyaxis right
hold on
plot(HR, CO, '-o', 'LineWidth', 2, 'Color', [0 0 0]);
% ylabel('CO (L/min)', 'FontName', 'Arial', 'FontSize', 9,'Color',[0,0,0]);
% ylim([5.5 5.63]);
ax = gca;
ax.YColor = [0,0,0];
ylim([-0.2 0.2])
% 设置 x 轴
% xlabel('HR (bpm)', 'FontName', 'Arial', 'FontSize', 9);
% xlim([45 82])
% % 图例
% legend({'LAP', 'CO'}, 'Location', 'best', 'FontName', 'Arial', 'FontSize', 10);

% 坐标轴字体
set(gca, 'FontName', 'Arial', 'FontSize', 16);
set(gca, 'XTickLabel', [], 'YTickLabel', []);

% 统一字体设置
box on;
w_in = 1.2;
h_in = 1.2;

set(gcf, 'Units', 'inches', 'Position', [1 1 w_in h_in]);
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 w_in h_in]);
set(gcf, 'InvertHardcopy', 'off');  % 防止背景反转
pbaspect([1 1 1]);  % 保持1:1比例

outpath = fullfile('C:', 'Users', 'fenggu', ...
    'University of Michigan Dropbox', 'Feng Gu', 'DTPacingHFpEF', 'PacingLAPCOHFpEF.png');
exportgraphics(gcf, outpath, 'Resolution', 900);
%%
figure(102);clf;
V_LV = HFpEFexample.PacingResults{1,1}.PVvolume;
P_LV = HFpEFexample.PacingResults{1,1}.PVpressure;
n = 7;  % 假设 P_LV 是颜色映射的长度
start_color = [0.6350, 0.0780, 0.1840];  % 起点（暗红）
end_color   = [0.5, 0.5, 0.5];                % 终点（蓝）

cmap = [linspace(start_color(1), end_color(1), n)', ...
    linspace(start_color(2), end_color(2), n)', ...
    linspace(start_color(3), end_color(3), n)'];
for i = 1:length(V_LV)
    if i == 1 || i == length(V_LV)
        lw = 2;  % 较粗
        alphaVal = 1.0;  % 不透明
    else
        lw = 0.88;  % 较细
        alphaVal = 0.4;  % 半透明
    end
    x = V_LV{i};
    y = P_LV{i};
    LVEDP_points(i,:) = [x(1), y(1)];
    plotObj = plot(x, y, '-', ...
        'Color', [cmap(i,:), alphaVal], ...
        'LineWidth', lw);
    hold on;
    % scatter(x(HFpEFexample.PacingRsults{1,1}.idxESP_all(i)),y(HFpEFexample.PacingRsults{1,1}.idxESP_all(i)));
end
xlabel('Volume (ml)', 'FontName', 'Arial', 'FontSize', 10);
ylabel('Pressure (mmHg)', 'FontName', 'Arial', 'FontSize', 10);
xlim([25 225]);
ylim([0 175]);
% 统一字体设置
set(gca, 'FontName', 'Arial', 'FontSize', 10);
box on;
pbaspect([1,1,1])
% 设置图像尺寸（单位为英寸，1 cm ≈ 0.3937 in）
set(gcf, 'Units', 'inches', 'Position', [1 1 2 2]);  % 9 cm x 9 cm ≈ 3.54 in x 3.54 in
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3.54 3.54]);

% 获取 HR 数组（共7个值）
HR_values = HFpEFexample.PacingResults{1,1}.HR;

% % 创建 colorbar 显示这些颜色
% colormap(cmap);  % 使用你刚才的自定义 cmap
% cb = colorbar;   % 添加 colorbar
% cb.Ticks = linspace(0, 1, length(HR_values));  % 设置 colorbar 的位置刻度
% cb.TickLabels = string(HR_values);  % 将 HR 显示为刻度标签
% cb.Label.String = 'HR (bpm)';
% cb.Label.FontName = 'Arial';
% cb.Label.FontSize = 10;
% cb.FontName = 'Arial';
% cb.FontSize = 9;

%%
figure(104); clf;
subplot(1,2,2)
plottedLower = false;
plottedUpper = false;
for i = 1:height(SummaryTable)
    if ~isempty(SummaryTable.PacingRsults{i})
        thisHR = SummaryTable.PacingRsults{i}.HR;
        thisCE = SummaryTable.PacingRsults{i}.CE_LV_all;

        deltaHR = thisHR - thisHR(1);
        deltaCE = thisCE - thisCE(1);
        if max(abs(deltaCE))>0.5
           error('PVA function bug')
        end
        if sign(deltaCE(2)) ~= sign(deltaCE(end))
            fprintf('patient %d weird\n', i);
        end

        if any(deltaCE > 0) && all(diff(deltaCE)>0)
            % Green line: Pacing lowers O2
            h = plot(deltaHR, deltaCE, '--', 'Color', [0.1 0.7 0.1 0.5], 'LineWidth', 1); hold on;
            if ~plottedLower
                set(h, 'DisplayName', 'Pacing lower O2');
                plottedLower = true;
            else
                set(h, 'HandleVisibility', 'off');
            end

        elseif any(deltaCE > 0)
            [~,optimalHRidx] = max(deltaCE);
            if optimalHRidx ==7
                h = plot(deltaHR, deltaCE, '--', 'Color', [0.1 0.7 0.1 0.5], 'LineWidth', 1); hold on;
            else
                h = plot(deltaHR, deltaCE, '-', 'Color', [0.1 0.7 0.1 0.5], 'LineWidth', 0.5); hold on;
                scatter(deltaHR(optimalHRidx),deltaCE(optimalHRidx),'SizeData',30,'MarkerEdgeColor', 'none','MarkerFaceColor',[0.1 0.7 0.1],'MarkerFaceAlpha',0.7)
            end

            if ~plottedLower
                set(h, 'DisplayName', 'Pacing lower O2');
                plottedLower = true;
            else
                set(h, 'HandleVisibility', 'off');
            end
        else
            % Red line: Pacing raises or maintains O2
            h = plot(deltaHR, deltaCE, '--', 'Color', [0.9 0.1 0.1 0.5], 'LineWidth', 1); hold on;
            if ~plottedUpper
                set(h, 'DisplayName', 'Pacing upper O2');
                plottedUpper = true;
            else
                set(h, 'HandleVisibility', 'off');
            end
        end
    end
end




% h3 = errorbar(nanmean(HRRes(2:end,:),2), nanmean(O2Res(2:end,:),2), nanstd(O2Res(2:end,:),0,2), nanstd(O2Res(2:end,:),0,2), ...
%     'o', 'Color', [0.1 0.7 0.1], 'CapSize', 14, 'LineWidth', 2, ...
%     'MarkerSize', 14, 'MarkerFaceColor', [0.1 0.7 0.1], 'MarkerEdgeColor', 'none','DisplayName','Responder');
% h4 = errorbar(nanmean(HRNon(2:end,:),2), nanmean(O2Non(2:end,:),2), nanstd(O2Non(2:end,:),0,2), nanstd(O2Non(2:end,:),0,2), ...
%     'o', 'Color', [0.9 0.1 0.1], 'CapSize', 14, 'LineWidth', 2, ...
%     'MarkerSize', 14, 'MarkerFaceColor', [0.9 0.1 0.1], 'MarkerEdgeColor', 'none','DisplayName','Non-responder');

% legend('Location', 'northwest',FontSize=16,Box='off');
xlabel('Δ HR (bpm)');
ylabel('Δ Cardiac Efficiency');
set(gca, 'FontName', 'Arial', 'FontSize', 10);
xlim([-5 35])
ylim([-0.05 0.05])
% 设置图形比例和边框
box on;
pbaspect([1,2,1])
% 设置图像尺寸（单位为英寸，1 cm ≈ 0.3937 in）
% set(gcf, 'Units', 'inches', 'Position', [1 1 2 4]);  % 9 cm x 9 cm ≈ 3.54 in x 3.54 in
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3.54 3.54]);
%%
subplot(1,2,1)
plottedLower = false;
plottedUpper = false;
for i = 1:height(SummaryTable)
    if ~isempty(SummaryTable.PacingRsults{i})
        thisHR = SummaryTable.PacingRsults{i}.HR;
        thisO2 = SummaryTable.PacingRsults{i}.MVO2_LV;

        deltaHR = thisHR - thisHR(1);
        deltaO2 = thisO2 - thisO2(1);
        if max(abs(deltaO2))>10
           error('PVA function bug')
        end
        if sign(deltaO2(4)) ~= sign(deltaO2(5))
            fprintf('patient %d weird\n', i);
        end
        if any(deltaO2 < 0) && all(diff(deltaO2)<0)
            % Green line: Pacing lowers O2
            h = plot(deltaHR, deltaO2, '--', 'Color', [0.1 0.7 0.1 0.5], 'LineWidth', 1); hold on;
            if ~plottedLower
                set(h, 'DisplayName', 'Pacing lower O2');
                plottedLower = true;
            else
                set(h, 'HandleVisibility', 'off');
            end

        elseif any(deltaO2 < 0)
            [~,optimalHRidx] = min(deltaO2);
            if optimalHRidx ==7
                h = plot(deltaHR, deltaO2, '--', 'Color', [0.1 0.7 0.1 0.5], 'LineWidth', 1); hold on;
            else
                h = plot(deltaHR, deltaO2, '-', 'Color', [0.1 0.7 0.1 0.5], 'LineWidth', 0.5); hold on;
                scatter(deltaHR(optimalHRidx),deltaO2(optimalHRidx),'SizeData',30,'MarkerEdgeColor', 'none','MarkerFaceColor',[0.1 0.7 0.1],'MarkerFaceAlpha',0.7)
            end

            if ~plottedLower
                set(h, 'DisplayName', 'Pacing lower O2');
                plottedLower = true;
            else
                set(h, 'HandleVisibility', 'off');
            end
        else
            % Red line: Pacing raises or maintains O2
            h = plot(deltaHR, deltaO2, '--', 'Color', [0.9 0.1 0.1 0.5], 'LineWidth', 1); hold on;
            if ~plottedUpper
                set(h, 'DisplayName', 'Pacing upper O2');
                plottedUpper = true;
            else
                set(h, 'HandleVisibility', 'off');
            end
        end
    end
end





% legend('Location', 'northwest',FontSize=16,Box='off');
xlabel('Δ HR (bpm)');
ylabel('Δ MVO_2');
set(gca, 'FontName', 'Arial', 'FontSize', 10);
xlim([-5 35])
ylim([-3 4])
% 设置图形比例和边框
box on;
pbaspect([1,2,1])
% 设置图像尺寸（单位为英寸，1 cm ≈ 0.3937 in）
set(gcf, 'Units', 'inches', 'Position', [1 1 4 4]);  % 9 cm x 9 cm ≈ 3.54 in x 3.54 in
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3.54 3.54]);

%%
figure('Units','inches','Position',[1,1,4,6]);
subplot(3,2,1);
HR = [];
LAP = [];
plottedOne = false;
for i = 1:height(SummaryTable)
    if ~isempty(SummaryTable.PacingRsults{i})
        % 计算 ΔHR 和 ΔLAP
        thisHR = SummaryTable.PacingRsults{i}.HR;
        thisLAP = SummaryTable.PacingRsults{i}.LAP;
        deltaHR = thisHR - thisHR(1);
        % deltaLAP = thisLAP - thisLAP(1);

        % 拼接进所有数据
        HR = [HR deltaHR];
        LAP = [LAP SummaryTable.PacingRsults{i}.LAP];

        % 画线
        h = plot(deltaHR,thisLAP, '-', 'Color', [0.7 0.7 0.7 0.5], 'LineWidth', 0.5); hold on;
        for ii = 1:length(deltaHR)
            scatter(deltaHR(ii), thisLAP(ii), ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', cmap(ii,:), ...
                'MarkerFaceAlpha',0.5,...
                'SizeData', 30);
        end
        if ~plottedOne
            set(h, 'DisplayName', 'Individual');
            plottedOne = true;
        else
            set(h, 'HandleVisibility', 'off');
        end
    end
end
% mean ± SD errorbar
% h2 = errorbar(mean(HR(1:end,:),2), mean(LAP(1:end,:),2), std(LAP(1:end,:),0,2), std(LAP(1:end,:),0,2), ...
%     'o', 'Color', 'b', 'CapSize', 14, 'LineWidth', 2, ...
%     'MarkerSize', 14, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none','DisplayName', 'Mean ± SD');
% figure();
x_pos = 0:5:30;  % x 轴实际位置
hold on
for i = 1:length(x_pos)
    h = boxchart( ...
        ones(size(LAP,2),1)*x_pos(i), ...   % 每组在指定 x 位置
        LAP(i,:)', ...                     % y 值（每行是一组）
        'BoxWidth', 3,...
        'notch','off');                    % 加宽 box
    h.BoxFaceColor = cmap(i,:);            % 设置颜色
    h.BoxEdgeColor = cmap(i,:);

    h.MarkerColor = 'none';
    h.LineWidth = 1.5;
end

xlim([-5 35])
ylim([-5 45])
pbaspect([1,1,1])
% legend('Location', 'southwest',FontSize=16);
xlabel('Δ HR (bpm)');
ylabel('LAP (mmHg)');
% title('Δ Pacing Response: HR vs LA Pressure', 'FontSize', 20);
% set(gca, 'FontSize', 14);  % X/Y tick label 的字体大小
% set(gca, 'LineWidth', 1.2);
% 统一字体设置
set(gca, 'FontName', 'Arial', 'FontSize', 10);
box on;
pbaspect([1,1,1])
% 设置图像尺寸（单位为英寸，1 cm ≈ 0.3937 in）
% set(gcf, 'Units', 'inches', 'Position', [1 1 2 2]);  % 9 cm x 9 cm ≈ 3.54 in x 3.54 in
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3.54 3.54]);

%%
subplot(3,2,4);
HR = [];
LAV = [];
plottedOne = false;
for i = 1:height(SummaryTable)
    if ~isempty(SummaryTable.PacingRsults{i})
        thisHR = SummaryTable.PacingRsults{i}.HR;
        thisLAV = SummaryTable.PacingRsults{i}.LAV;
        deltaHR = thisHR - thisHR(1);
        % deltaLAV = thisLAV - thisLAV(1);

        HR = [HR deltaHR];
        LAV = [LAV SummaryTable.PacingRsults{i}.LAV];

        h = plot(deltaHR, thisLAV, '-', 'Color', [0.7 0.7 0.7 0.5], 'LineWidth', 0.5); hold on;
        for ii = 1:length(deltaHR)
            scatter(deltaHR(ii), thisLAV(ii), ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', cmap(ii,:), ...
                'MarkerFaceAlpha',0.5,...
                'SizeData', 30);
        end
        if ~plottedOne
            set(h, 'DisplayName', 'Individual');
            plottedOne = true;
        else
            set(h, 'HandleVisibility', 'off');
        end
    end
end

% h2 = errorbar(mean(HR(2:end,:),2), mean(LAV(2:end,:),2), std(LAV(2:end,:),0,2), std(LAV(2:end,:),0,2), ...
%     'o', 'Color', 'b', 'CapSize', 14, 'LineWidth', 2, ...
%     'MarkerSize', 14, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none','DisplayName', 'Mean ± SD');

x_pos = 0:5:30;  % x 轴实际位置
hold on
for i = 1:length(x_pos)
    h = boxchart( ...
        ones(size(LAV,2),1)*x_pos(i), ...   % 每组在指定 x 位置
        LAV(i,:)', ...                     % y 值（每行是一组）
        'BoxWidth', 3,...
        'notch','off');                    % 加宽 box
    h.BoxFaceColor = cmap(i,:);            % 设置颜色
    h.BoxEdgeColor = cmap(i,:);

    h.MarkerColor = 'none';
    h.LineWidth = 1.5;
end

pbaspect([1,1,1])
% legend('Location', 'southwest',FontSize=16);
xlabel('Δ HR (bpm)', 'FontSize', 16);
ylabel('LAV (ml)', 'FontSize', 16);
xlim([-5 35])
ylim([10 160])
% title('Δ Pacing Response: HR vs LA Volume', 'FontSize', 20);
% set(gca, 'FontSize', 14);  % X/Y tick label 的字体大小
% set(gca, 'LineWidth', 1.2);
% 统一字体设置
set(gca, 'FontName', 'Arial', 'FontSize', 10);
box on;
pbaspect([1,1,1])
% 设置图像尺寸（单位为英寸，1 cm ≈ 0.3937 in）
% set(gcf, 'Units', 'inches', 'Position', [1 1 2 2]);  % 9 cm x 9 cm ≈ 3.54 in x 3.54 in
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3.54 3.54]);
%%
subplot(3,2,2);
HR = [];
SV = [];
plottedOne = false;
for i = 1:height(SummaryTable)
    if ~isempty(SummaryTable.PacingRsults{i})
        thisHR = SummaryTable.PacingRsults{i}.HR;
        thisSV = SummaryTable.PacingRsults{i}.SV;
        deltaHR = thisHR - thisHR(1);
        % deltaSV = thisSV - thisSV(1);

        HR = [HR deltaHR];
        SV = [SV SummaryTable.PacingRsults{i}.SV];

        h = plot(deltaHR, SummaryTable.PacingRsults{i}.SV, '-', 'Color', [0.7 0.7 0.7 0.5], 'LineWidth', 0.5); hold on;
        for ii = 1:length(deltaHR)
            scatter(deltaHR(ii), thisSV(ii), ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', cmap(ii,:), ...
                'MarkerFaceAlpha',0.5,...
                'SizeData', 30);
        end
        if ~plottedOne
            set(h, 'DisplayName', 'Individual');
            plottedOne = true;
        else
            set(h, 'HandleVisibility', 'off');
        end
    end
end

% h2 = errorbar(mean(HR(2:end,:),2), mean(SV(2:end,:),2), std(SV(2:end,:),0,2), std(SV(2:end,:),0,2), ...
%     'o', 'Color', 'b', 'CapSize', 14, 'LineWidth', 2, ...
%     'MarkerSize', 14, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none','DisplayName', 'Mean ± SD');
x_pos = 0:5:30;  % x 轴实际位置
hold on
for i = 1:length(x_pos)
    h = boxchart( ...
        ones(size(SV,2),1)*x_pos(i), ...   % 每组在指定 x 位置
        SV(i,:)', ...                     % y 值（每行是一组）
        'BoxWidth', 3,...
        'notch','off');                    % 加宽 box
    h.BoxFaceColor = cmap(i,:);            % 设置颜色
    h.BoxEdgeColor = cmap(i,:);

    h.MarkerColor = 'none';
    h.LineWidth = 1.5;
end

pbaspect([1,1,1])
% legend('Location', 'southwest','FontSize',16);
xlabel('Δ HR (bpm)', 'FontSize', 16);
ylabel('Stroke Volume (ml)', 'FontSize', 16);
xlim([-5 35])
% title('Δ Pacing Response: HR vs SV', 'FontSize', 20);
% set(gca, 'FontSize', 14);  % X/Y tick label 的字体大小
% set(gca, 'LineWidth', 1.2);
set(gca, 'FontName', 'Arial', 'FontSize', 10);
box on;
pbaspect([1,1,1])
% 设置图像尺寸（单位为英寸，1 cm ≈ 0.3937 in）
% set(gcf, 'Units', 'inches', 'Position', [1 1 2 2]);  % 9 cm x 9 cm ≈ 3.54 in x 3.54 in
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3.54 3.54]);

%%
subplot(3,2,3);
HR = [];
CO = [];
plottedOne = false;
for i = 1:height(SummaryTable)
    if ~isempty(SummaryTable.PacingRsults{i})
        thisHR = SummaryTable.PacingRsults{i}.HR;
        thisCO = SummaryTable.PacingRsults{i}.CO;
        deltaHR = thisHR - thisHR(1);
        % deltaCO = thisCO - thisCO(1);

        HR = [HR deltaHR];
        CO = [CO SummaryTable.PacingRsults{i}.CO];

        h = plot(deltaHR, thisCO, '-', 'Color', [0.7 0.7 0.7 0.5], 'LineWidth', 0.5); hold on;
        for ii = 1:length(deltaHR)
            scatter(deltaHR(ii), thisCO(ii), ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', cmap(ii,:), ...
                'MarkerFaceAlpha',0.5,...
                'SizeData', 30);
        end
        if ~plottedOne
            set(h, 'DisplayName', 'Individual');
            plottedOne = true;
        else
            set(h, 'HandleVisibility', 'off');
        end
    end
end

% h2 = errorbar(mean(HR(2:end,:),2), mean(CO(2:end,:),2), std(CO(2:end,:),0,2), std(CO(2:end,:),0,2), ...
%     'o', 'Color', 'b', 'CapSize', 14, 'LineWidth', 2, ...
%     'MarkerSize', 14, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none','DisplayName', 'Mean ± SD');

x_pos = 0:5:30;  % x 轴实际位置
hold on
for i = 1:length(x_pos)
    h = boxchart( ...
        ones(size(CO,2),1)*x_pos(i), ...   % 每组在指定 x 位置
        CO(i,:)', ...                     % y 值（每行是一组）
        'BoxWidth', 3,...
        'notch','off');                    % 加宽 box
    h.BoxFaceColor = cmap(i,:);            % 设置颜色
    h.BoxEdgeColor = cmap(i,:);

    h.MarkerColor = 'none';
    h.LineWidth = 1.5;
end

% legend('Location', 'southwest','FontSize',16);
xlabel('Δ HR (bpm)', 'FontSize', 16);
xlim([-5 35])
ylim([0.5 11])
ylabel('CO (L/min)', 'FontSize', 16);
% title('Δ Pacing Response: HR vs SV', 'FontSize', 20);
% set(gca, 'FontSize', 14);  % X/Y tick label 的字体大小
% set(gca, 'LineWidth', 1.2);

set(gca, 'FontName', 'Arial', 'FontSize', 10);
box on;
pbaspect([1,1,1])
% 设置图像尺寸（单位为英寸，1 cm ≈ 0.3937 in）
% set(gcf, 'Units', 'inches', 'Position', [1 1 2 2]);  % 9 cm x 9 cm ≈ 3.54 in x 3.54 in
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3.54 3.54]);


%%
subplot(3,2,5);
% figure(303);clf;
HR = [];
MVO2 = [];
plottedOne = false;
for i = 1:height(SummaryTable)
    if ~isempty(SummaryTable.PacingRsults{i})
        thisHR = SummaryTable.PacingRsults{i}.HR;
        thisMVO2 = SummaryTable.PacingRsults{i}.MVO2_LV;
        deltaHR = thisHR - thisHR(1);
        % deltaCO = thisCO - thisCO(1);

        HR = [HR deltaHR];
        MVO2 = [MVO2 thisMVO2];

        h = plot(deltaHR, thisMVO2, '-', 'Color', [0.7 0.7 0.7 0.5], 'LineWidth', 0.5); hold on;
        for ii = 1:length(deltaHR)
            scatter(deltaHR(ii), thisMVO2(ii), ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', cmap(ii,:), ...
                'MarkerFaceAlpha',0.5,...
                'SizeData', 30);
        end
        if ~plottedOne
            set(h, 'DisplayName', 'Individual');
            plottedOne = true;
        else
            set(h, 'HandleVisibility', 'off');
        end
    end
end

% h2 = errorbar(mean(HR(2:end,:),2), mean(CO(2:end,:),2), std(CO(2:end,:),0,2), std(CO(2:end,:),0,2), ...
%     'o', 'Color', 'b', 'CapSize', 14, 'LineWidth', 2, ...
%     'MarkerSize', 14, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none','DisplayName', 'Mean ± SD');

x_pos = 0:5:30;  % x 轴实际位置
hold on
for i = 1:length(x_pos)
    h = boxchart( ...
        ones(size(MVO2,2),1)*x_pos(i), ...   % 每组在指定 x 位置
        MVO2(i,:)', ...                     % y 值（每行是一组）
        'BoxWidth', 3,...
        'notch','off');                    % 加宽 box
    h.BoxFaceColor = cmap(i,:);            % 设置颜色
    h.BoxEdgeColor = cmap(i,:);

    h.MarkerColor = 'none';
    h.LineWidth = 1.5;
end

% legend('Location', 'southwest','FontSize',16);
xlabel('Δ HR (bpm)', 'FontSize', 16);
xlim([-5 35])
ylim([0 25])
ylabel('MVO_2 (mL/min)');
% title('Δ Pacing Response: HR vs SV', 'FontSize', 20);
% set(gca, 'FontSize', 14);  % X/Y tick label 的字体大小
% set(gca, 'LineWidth', 1.2);

set(gca, 'FontName', 'Arial', 'FontSize', 10);
box on;
pbaspect([1,1,1])
% 设置图像尺寸（单位为英寸，1 cm ≈ 0.3937 in）
% set(gcf, 'Units', 'inches', 'Position', [1 1 2 2]);  % 9 cm x 9 cm ≈ 3.54 in x 3.54 in
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3.54 3.54]);
%%
subplot(3,2,6);
HR = [];
CardiacEffiency = [];
plottedOne = false;
for i = 1:height(SummaryTable)
    if ~isempty(SummaryTable.PacingRsults{i})
        thisHR = SummaryTable.PacingRsults{i}.HR;
        thisCardiacEffiency = SummaryTable.PacingRsults{i}.MAP;
        deltaHR = thisHR - thisHR(1);
        % deltaCO = thisCO - thisCO(1);

        HR = [HR deltaHR];
        CardiacEffiency = [CardiacEffiency thisCardiacEffiency];

        h = plot(deltaHR, thisCardiacEffiency, '-', 'Color', [0.7 0.7 0.7 0.5], 'LineWidth', 0.5); hold on;
        for ii = 1:length(deltaHR)
            scatter(deltaHR(ii), thisCardiacEffiency(ii), ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', cmap(ii,:), ...
                'MarkerFaceAlpha',0.5,...
                'SizeData', 30);
        end
        if ~plottedOne
            set(h, 'DisplayName', 'Individual');
            plottedOne = true;
        else
            set(h, 'HandleVisibility', 'off');
        end
    end
end

% h2 = errorbar(mean(HR(2:end,:),2), mean(CO(2:end,:),2), std(CO(2:end,:),0,2), std(CO(2:end,:),0,2), ...
%     'o', 'Color', 'b', 'CapSize', 14, 'LineWidth', 2, ...
%     'MarkerSize', 14, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none','DisplayName', 'Mean ± SD');

x_pos = 0:5:30;  % x 轴实际位置
hold on
for i = 1:length(x_pos)
    h = boxchart( ...
        ones(size(CardiacEffiency,2),1)*x_pos(i), ...   % 每组在指定 x 位置
        CardiacEffiency(i,:)', ...                     % y 值（每行是一组）
        'BoxWidth', 3,...
        'notch','off');                    % 加宽 box
    h.BoxFaceColor = cmap(i,:);            % 设置颜色
    h.BoxEdgeColor = cmap(i,:);

    h.MarkerColor = 'none';
    h.LineWidth = 1.5;
end

% legend('Location', 'southwest','FontSize',16);
xlabel('Δ HR (bpm)', 'FontSize', 16);
xlim([-5 35])
% ylim([0.05 0.4])
ylabel('Cardiac Efficiency');
% title('Δ Pacing Response: HR vs SV', 'FontSize', 20);
% set(gca, 'FontSize', 14);  % X/Y tick label 的字体大小
% set(gca, 'LineWidth', 1.2);

set(gca, 'FontName', 'Arial', 'FontSize', 10);
box on;
pbaspect([1,1,1])
% 设置图像尺寸（单位为英寸，1 cm ≈ 0.3937 in）
% set(gcf, 'Units', 'inches', 'Position', [1 1 2 2]);  % 9 cm x 9 cm ≈ 3.54 in x 3.54 in
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3.54 3.54]);

%%
clear
% === 1. 读取数据 ===
T  = readtable("input_VAE.csv");
T  = T(:, [(1:34)]);   % 取前32列 + 最后两列
TF = readtable('generated_fake_data_correlated.csv');
TF = TF(:, [(1:34)]);  % 与 T 保持一致

% === 2. 转成矩阵，并计算相关系数 ===
X_real = table2array(T);
X_fake = table2array(TF);

R_real = corr(X_real);  % 真实数据相关性
R_fake = corr(X_fake);  % 生成数据相关性

n = size(R_real, 1);     % 变量个数，例如 34

% === 3. 初始化大矩阵 (n+1 x n+1)，用于错开上/下三角 ===
R_mix = NaN(n+1, n+1);

% 掩码
mask_lower = tril(true(n), -1);
mask_upper = triu(true(n), 1);

% === 4. 将真实数据的下三角放在左下角 ===
tmp_lower = NaN(n);  % 创建一个 n x n 全 NaN 矩阵
tmp_lower(mask_lower) = R_real(mask_lower);
R_mix(2:end, 1:end-1) = tmp_lower;

% === 5. 将生成数据的上三角放在右上角 ===
New_mask_upper = [mask_upper ; zeros(1,34)];
New_mask_upper = [zeros(1,35)' New_mask_upper];
New_mask_upper = logical(New_mask_upper);
% tmp_upper(mask_upper) = R_fake(mask_upper);
R_mix(New_mask_upper) = R_fake(mask_upper);
R_mix(isnan(R_mix)) = 0;
% === 6. 画图 ===
figure;
imagesc(R_mix);
axis square;
colormap(redbluecmap);  % 你需要有 redbluecmap 函数（见备注）
caxis([-1 1]);
colorbar;

% === 7. 设置标签 ===
% xticks(2:n+1);
% yticks(1:n);
% customLabels = arrayfun(@(x) sprintf('DT_%d', x), 1:n, 'UniformOutput', false);
% xticklabels(customLabels);
% yticklabels(customLabels);
% xtickangle(45);
colorbar('off');             % 关闭colorbar
set(gca, 'XTickLabel', []);  % 去掉x轴标签
set(gca, 'YTickLabel', []);  % 去掉y轴标签

% 美化
set(gca, 'TickLength', [0 0]);
% title('Correlation Structure: Real (lower left) vs Generated (upper right)');

set(gca, 'FontName', 'Arial', 'FontSize', 10);
axis off;

pbaspect([1,1,1])
% 设置图像尺寸（单位为英寸，1 cm ≈ 0.3937 in）
set(gcf, 'Units', 'inches', 'Position', [1 1 2 2]);  % 9 cm x 9 cm ≈ 3.54 in x 3.54 in
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3.54 3.54]);

%%
T = readtable("input_VAE.csv");
T = T(:,[(1:32) 34 35]);

TF = readtable("CholeSkyTry1.csv");
TF = TF(:,[(1:32) 34 35]);

% 转为矩阵
X_real = table2array(T);
X_fake = table2array(TF);

% 计算相关系数矩阵
R_real = corrcoef(X_real);
R_fake = corrcoef(X_fake);

% 初始化组合矩阵
R_mix = nan(size(R_real));  % 先填 NaN 是为了可视化清晰

% 填入下三角为真实数据，上三角为生成数据
R_mix(tril(true(size(R_real)), -1)) = R_real(tril(true(size(R_real)), -1));
R_mix(triu(true(size(R_fake)), 1))  = R_fake(triu(true(size(R_fake)), 1));
% R_mix = Combined_corr;
% 对角线统一设置为 1
R_mix(1:size(R_mix,1)+1:end) = 0;

% 画图
figure;
imagesc(R_mix);
colormap(redbluecmap);  % 同前面介绍的红-白-蓝 colormap
% colorbar;
caxis([-1 1]);  % 设置颜色范围

set(gca, 'TickLength', [0 0]);
% title('Correlation Structure: Real (lower left) vs Generated (upper right)');

set(gca, 'FontName', 'Arial', 'FontSize', 10);
axis off;

pbaspect([1,1,1])
% 设置图像尺寸（单位为英寸，1 cm ≈ 0.3937 in）
set(gcf, 'Units', 'inches', 'Position', [1 1 2.5 2.5]);  % 9 cm x 9 cm ≈ 3.54 in x 3.54 in
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3.54 3.54]);
%%
load UWcohort.mat
UWpatients.('TVr')(isnan(UWpatients.('TVr'))) = 1;
UWpatients.('PVr')(isnan(UWpatients.('PVr'))) = 1;
UWpatients.EF(isnan(UWpatients.('EF'))) = UWpatients.EF_LV(isnan(UWpatients.('EF')));
UWpatients.CO(isnan(UWpatients.CO)) = UWpatients.CO_TD(isnan(UWpatients.CO));
Raw_Lasso_RVEF = UWpatients.Sex * k_RVEF(1) + UWpatients.Age * k_RVEF(2) + ...
        UWpatients.SBP * k_RVEF(3) + UWpatients.PPAd * k_RVEF(4) + ...
        UWpatients.CO * k_RVEF(5) + UWpatients.EF * k_RVEF(6) + ...
        UWpatients.('TVr') * k_RVEF(7) + UWpatients.('PVr') * k_RVEF(8) + ...
        b_RVEF;
C_Lasso_RVEF = Raw_Lasso_RVEF * ck_RVEF + cb_RVEF;
C_Lasso_RVEF = C_Lasso_RVEF(1:64);
Y = UWpatients.EF_RV(1:64);
% === 原始散点图 ===
figure(1); clf;
scatter(C_Lasso_RVEF, Y, 40, 'b', 'filled'); % 蓝色实心点，40是点大小
hold on;

% === y = x 参考线 ===
plot([0 100], [0 100], 'r--', 'LineWidth', 1.5);  % 红色虚线，假设范围是 0~100

% === 线性回归线 ===
mdl = fitlm(C_Lasso_RVEF, Y);
xfit = linspace(min(C_Lasso_RVEF), max(Y), 100);
yfit = predict(mdl, xfit');
% plot(xfit, yfit, 'k-', 'LineWidth', 1.5);  % 黑色实线是拟合线

% === 计算相关系数和p值 ===
[R, P] = corr(C_Lasso_RVEF, Y);

% === 显示 R 和 P 值在图上 ===
textPosX = 0 + 0.05 * range(C_Lasso_RVEF);
textPosY = 80 - 0.1 * range(UWpatients.EF_RV);

text(textPosX, textPosY, ...
     sprintf('p = %.3g\nR = %.2f', P, R), ...
     'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'top');


% === 图像格式化 ===
xlabel('Predicted RV EF (%)');
ylabel('Measured RV EF (%)');
box on;
axis square;
xlim([0 80]);
ylim([0 80]);
grid on;

set(gcf, 'FontName', 'Arial', 'FontSize', 10);

pbaspect([1,1,1])
% 设置图像尺寸（单位为英寸，1 cm ≈ 0.3937 in）
set(gcf, 'Units', 'inches', 'Position', [1 1 2.5 2.5]);  % 9 cm x 9 cm ≈ 3.54 in x 3.54 in
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3.54 3.54]);
%%

load AllPatients.mat
C_Lasso_RVEF = [];
Y = [];
for PATIENT_NO = 1:370 % any number between 1 and 370, example patient in paper is 192
    MRI_flag = 0;
    ModelWin = 1;
    [Windowdate,targets, inputs, mods] = targetVals_HF(patients,PATIENT_NO,ModelWin,MRI_flag);
    C_Lasso_RVEF = [ C_Lasso_RVEF targets.RVEF];
    MRI_flag = 1;
    [Windowdate,targets, inputs, mods] = targetVals_HF(patients,PATIENT_NO,ModelWin,MRI_flag);
    Y = [Y (targets.RVEDV-targets.RVESV)/targets.RVEDV*100];
end
%%
% === 原始散点图 ===
figure(1); clf;
scatter(C_Lasso_RVEF, Y, 40, 'b', 'filled'); % 蓝色实心点，40是点大小
hold on;

% === y = x 参考线 ===
plot([0 100], [0 100], 'r--', 'LineWidth', 1.5);  % 红色虚线，假设范围是 0~100

% === 线性回归线 ===
mdl = fitlm(C_Lasso_RVEF, Y);
xfit = linspace(min(C_Lasso_RVEF), max(Y), 100);
yfit = predict(mdl, xfit');
% plot(xfit, yfit, 'k-', 'LineWidth', 1.5);  % 黑色实线是拟合线

% === 计算相关系数和p值 ===
[R, P] = corr(C_Lasso_RVEF', Y');

% === 显示 R 和 P 值在图上 ===
textPosX = -10 + 0.05 * range(C_Lasso_RVEF);
textPosY = 100 - 0.1 * range(UWpatients.EF_RV);

text(textPosX, textPosY, ...
     sprintf('p = %.3g\nR = %.2f', P, R), ...
     'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'top');


% === 图像格式化 ===
xlabel('Predicted RV EF (%)');
ylabel('Measured RV EF (%)');
box on;
axis square;
xlim([-10 100]);
ylim([-10 100]);
yticks([0 20 40 60 80]);
xticks([0 20 40 60 80]);
grid on;

set(gcf, 'FontName', 'Arial', 'FontSize', 10);

pbaspect([1,1,1])
% 设置图像尺寸（单位为英寸，1 cm ≈ 0.3937 in）
set(gcf, 'Units', 'inches', 'Position', [1 1 2.5 2.5]);  % 9 cm x 9 cm ≈ 3.54 in x 3.54 in
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3.54 3.54]);