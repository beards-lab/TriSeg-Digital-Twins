%% Script Summary:
% This script generates a 6-panel figure for patients, including:
% - Left pressures
% - Right pressures
% - Volumes
% - Wall thicknesses
% - Left flows
% - Right flows
% It also includes usual text output as well as target matching text output, and outputs costs in
% descending order.

% Created by Andrew Meyer and Feng Gu
% Last modified: 10/29/2024

trng = [t(1) t(end)];
figure(101); clf;
tiles = tiledlayout(2,3);
tiles.TileSpacing = 'compact';
tiles.Padding = 'tight';
set(gcf,'defaultLegendAutoUpdate','off','Position', [150, 50, 1280/1.5, 720/1.5]);

fig1 = nexttile;
xlim(fig1,trng);ylim([0 1.2*max(P_LV)]);
hold(fig1,"on");
plot(fig1,t, P_LV,'color',[0.6350 0.0780 0.1840], 'Linewidth', 2);
plot(fig1,t,P_RV, 'color', [0 0.4470 0.7410], 'Linewidth', 2);
plot(fig1,t, P_SA, 'color',[0.8500 0.5250 0.0980], 'Linewidth', 2);
plot(fig1,t, P_PA,'color',[0.3010 0.7450 0.9330], 'Linewidth', 2);
lgd = legend(fig1,'P_{LV}', 'P_{RV}', 'P_{SA}', 'P_{PA}');
lgd.FontSize = 7;
axPosition = get(fig1, 'Position');
lgdPosition = get(lgd, 'Position');
set(gca, 'FontSize', 12);
set(gca, 'FontWeight', 'bold');
lgdPosition(1) = axPosition(1) + axPosition(3) - lgdPosition(3) - 0.003;
lgdPosition(2) = axPosition(2) + axPosition(4) - lgdPosition(4) - 0.033;
set(lgd, 'Position', lgdPosition);
xlabel(fig1,'time (s)','FontSize',12); ylabel(fig1,'BP (mmHg)','FontSize',12);

fig2 = nexttile;
hold(fig2,"on");
plot(fig2,t, P_SV, 'color',[0.4660 0.6740 0.1880],'LineWidth',2);
plot(fig2,t, P_PV,'color',[0.9290 0.7940 0.1250],'LineWidth',2);
plot(fig2,t, P_RA, 'color',[0.4940 0.1840 0.5560],'LineWidth',2);
plot(fig2,t, P_LA,'color',[0.8500 0.3250 0.0980],'LineWidth',2);
xlim(fig2,trng); ylim([0.9*min([P_SV;P_PV;P_RA;P_LA]) 1.15*max([P_SV;P_PV;P_RA;P_LA])]);
plot(fig2,trng, [o_vals.PCWP, o_vals.PCWP],'color',[0.9290 0.7940 0.1250]);

lgd = legend(fig2,'P_{SV}', 'P_{PV}', 'P_{RA}', 'P_{LA}');
lgd.FontSize = 7;
axPosition = get(fig2, 'Position');
lgdPosition = get(lgd, 'Position');
set(gca, 'FontSize', 12);
set(gca, 'FontWeight', 'bold');
lgdPosition(1) = axPosition(1) + axPosition(3) - lgdPosition(3) - 0.016;
lgdPosition(2) = axPosition(2) + axPosition(4) - lgdPosition(4) - 0.03;
set(lgd, 'Position', lgdPosition);
xlabel(fig2,'time (s)','FontSize',12); ylabel(fig2,'BP (mmHg)','FontSize',12);

fig3 = nexttile;
hold(fig3,"on");
plot(fig3,t, V_LV, 'color', [0.6350 0.0780 0.1840],'Linewidth', 2);
plot(fig3,t, V_RV, 'color', [0 0.4470 0.7410],'Linewidth', 2);
plot(fig3,t, V_LA, 'color', [0.8500 0.3250 0.0980],'Linewidth', 2);
plot(fig3,t, V_RA, 'color', [0.4940 0.1840 0.5560],'Linewidth', 2);
xlim(fig3,trng); ylim([0 1.3*max([V_LV;V_RV;V_LA;V_RA])]);

lgd = legend(fig3,'LV', 'RV', 'LA', 'RA');
lgd.FontSize = 7;
axPosition = get(fig3, 'Position');
lgdPosition = get(lgd, 'Position');
set(gca, 'FontSize', 12);
set(gca, 'FontWeight', 'bold');
lgdPosition(1) = axPosition(1) + axPosition(3) - lgdPosition(3) - 0.008;
lgdPosition(2) = axPosition(2) + axPosition(4) - lgdPosition(4) - 0.023;
set(lgd, 'Position', lgdPosition);

xlabel(fig3,'time (s)','FontSize',12); ylabel(fig3,'Volume (mL)','FontSize',12);

fig4 = nexttile(4);
plot(fig4,t, Q_t*60, t, Q_p*60, 'Linewidth', 2);
xlim(fig4,trng); ylim([1.02*min([Q_t*60;Q_p*60]) 1.2*max([Q_t*60;Q_p*60])])
hold(fig4,"on");
lgd = legend(fig4,'Tricuspid Valve', 'Pulmonary Valve');
lgd.FontSize = 7;
axPosition = get(fig4, 'Position');
lgdPosition = get(lgd, 'Position');
set(gca, 'FontSize', 12);
set(gca, 'FontWeight', 'bold');
lgdPosition(1) = axPosition(1) + axPosition(3) - lgdPosition(3) - 0.005;
lgdPosition(2) = axPosition(2) + axPosition(4) - lgdPosition(4) - 0.003;
set(lgd, 'Position', lgdPosition);
xlabel(fig4,'time (s)','FontSize',12); ylabel(fig4,'Flow (ml/min)','FontSize',12);
set(gca, "YTick", [-5e4 -4e4 -3e4 -2e4 -1e4 0 1e4 2e4 3e4 4e4 5e4]);
set(gca().YAxis, 'Exponent', 4);
set(gca, 'TickLength', [0.040 0.040])
box off

fig5 = nexttile(5);
plot(fig5,t, Q_m*60, t, Q_a*60, 'Linewidth', 2);  xlim(fig5,trng);ylim([1.02*min([Q_m*60;Q_a*60]) 1.2*max([Q_m*60;Q_a*60])]);
EAr_text = "E/A: " + num2str(round(E_A_ratio,2));...
    % + " (" + targets.EAr+")";
EAr_text_position = [0.296, 0.27, 0.4, 0.195];
annotation('textbox', EAr_text_position, 'String', EAr_text, 'FitBoxToText', 'on','FontWeight','bold','FontSize',10,'BackgroundColor',"w");
lgd = legend(fig5,'Mitral Valve', 'Aortic Valve');
lgd.FontSize = 7;
axPosition = get(fig5, 'Position');
lgdPosition = get(lgd, 'Position');
set(gca, 'FontSize', 12);
set(gca, 'FontWeight', 'bold');
lgdPosition(1) = axPosition(1) + axPosition(3) - lgdPosition(3) - 0.008;
lgdPosition(2) = axPosition(2) + axPosition(4) - lgdPosition(4) - 0.018;
set(lgd, 'Position', lgdPosition);
xlabel(fig5,'time (s)','FontSize',12); ylabel(fig5,'Flow (ml/min)','FontSize',12);
set(gca, "YTick", [-5e4 -4e4 -3e4 -2e4 -1e4 0 1e4 2e4 3e4 4e4 5e4]);
set(gca().YAxis, 'Exponent', 4);
set(gca, 'TickLength', [0.040 0.040])
box off

fig6 = nexttile(6);
hold(fig6,"on");
plot(fig6,t, d_LW,'color', [0.6350 0.0780 0.1840],'LineWidth',2);
plot(fig6,t, d_SW, 'color', [0.4660 0.6740 0.1880],'LineWidth',2);
plot(fig6,t, d_RW, 'color', [0 0.4470 0.7410],'LineWidth',2);
% xlim(fig7,trng); ylim([0 1.2*max([d_SW;d_LW;d_RW;targets.Hed_LW;targets.Hed_SW])])
lgd = legend(fig6,'LV', 'Septum', 'RV');
lgd.FontSize = 7;
axPosition = get(fig6, 'Position');
lgdPosition = get(lgd, 'Position');
set(gca, 'FontSize', 12);
set(gca, 'FontWeight', 'bold');
lgdPosition(1) = axPosition(1) + axPosition(3) - lgdPosition(3) - 0.008;
lgdPosition(2) = axPosition(2) + axPosition(4) - lgdPosition(4) - 0.020;
set(lgd, 'Position', lgdPosition);
xlabel(fig6,'time (s)','FontSize',12); ylabel(fig6,'Thickness (cm)','FontSize',12);

% % 保存图像为300dpi PNG
% outputFolder = 'Virtual patients 0617';
% if ~exist(outputFolder, 'dir')
%     mkdir(outputFolder);
% end
% 
% % 构造文件名
% filename = sprintf('Virtual_Patient_%03d.png',PATIENT_NO);
% filepath = fullfile(outputFolder, filename);
% 
% % 保存为高分辨率 PNG
% print(gcf, filepath, '-dpng', '-r300');  % -r300 指定300dpi
