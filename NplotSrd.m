%% Script Summary:
% This script generates a 4-panel figure for a canonical subject, including:
% - Left pressures
% - Right pressures
% - Volumes
% - Left flows
% It also includes usual text output as well as target matching text output, and outputs costs in
% descending order.

% Created by Andrew Meyer and Feng Gu
% Last modified: 10/29/2024

figure(1);clf
trng = [t(1) t(end)];
if(isfield(targets,'CO'))
    plot_CO = targets.CO;
else
    plot_CO = NaN;
end
if(isfield(targets,'EAr'))
    plot_EA = targets.EAr;
else
    plot_EA = NaN;
end
if(~exist('total_cost','var'))
    assert(exist('struct_pXX','var'));
    total_cost = struct_pXX.output_vals.total_cost;
end

tiles = tiledlayout(2,3);
tiles.TileSpacing = 'compact';
tiles.Padding = 'tight';
set(gcf,'defaultLegendAutoUpdate','off','position',[50,50,1350,750]);

fig1 = nexttile;
xlabel(fig1,'time (s)', 'FontSize', 8, 'FontName', 'Times New Roman');
ylabel(fig1,'Pressure (mmHg)', 'FontSize', 8, 'FontName', 'Times New Roman');
xlim(fig1,trng); title(fig1,'Pressures', 'FontSize', 8, 'FontName', 'Times New Roman');
hold(fig1,"on");
plot(fig1,t, P_LV,'color',[0.6350 0.0780 0.1840], 'Linewidth', 2);
plot(fig1,t,P_RV, 'color', [0 0.4470 0.7410], 'Linewidth', 2);
plot(fig1,t, P_SA, 'color',[0.8500 0.5250 0.0980], 'Linewidth', 2);
plot(fig1,t, P_PA,'color',[0.3010 0.7450 0.9330], 'Linewidth', 2);
plot(fig1,trng, [targets.SBP, targets.SBP],'color',[0.8500 0.5250 0.0980],'LineStyle','--','Linewidth', 1);
plot(fig1,trng, [targets.DBP, targets.DBP],'color',[0.8500 0.5250 0.0980],'LineStyle','--','Linewidth', 1);
if(isfield(targets,'PASP'))
    plot(fig1,trng, [targets.PASP, targets.PASP],'color',[0.3010 0.7450 0.9330],'LineStyle','--', 'Linewidth', 1);
    plot(fig1,trng, [targets.PADP, targets.PADP],'color',[0.3010 0.7450 0.9330],'LineStyle','--', 'Linewidth', 1);
end
if(isfield(targets,'RVSP'))
    plot(fig1,trng, [targets.RVSP, targets.RVSP],'color',[0.3010 0.7450 0.9330],'LineStyle','--', 'Linewidth', 1);
end
if(isfield(targets,'LVESP'))
    plot(fig1,trng, [targets.LVESP, targets.LVESP],'color',[0.6350 0.0780 0.1840],'LineStyle','--', 'Linewidth', 1);
end
if isfield(targets, 'RVEDP')
    plot(fig1,trng, [targets.RVEDP, targets.RVEDP],'color', [0 0.2470 0.8410],'LineStyle','--');
end
if isfield(targets, 'LVEDP')
    plot(fig1,trng, [targets.LVEDP, targets.LVEDP],'color', [0.6350 0.0780 0.1840],'LineStyle','--');
end
if isfield(targets,'P_RV_min')
    plot(fig1,trng, [targets.P_RV_min, targets.P_RV_min],'color', [0 0.2470 0.8410],'LineStyle','--');
end
if isfield(targets,'P_LV_min')
    plot(fig1,trng, [targets.P_LV_min, targets.P_LV_min],'color', [0.6350 0.0780 0.1840],'LineStyle','--');
end
ylim([0 160]);

fig2 = nexttile(2);
hold(fig2,"on");
plot(fig2,t, P_SV, 'color',[0.4660 0.6740 0.1880],'LineWidth',2);
plot(fig2,t, P_PV,'color',[0.9290 0.7940 0.1250],'LineWidth',2);
plot(fig2,t, P_RA, 'color',[0.4940 0.1840 0.5560],'LineWidth',2);
plot(fig2,t, P_LA,'color',[0.8500 0.3250 0.0980],'LineWidth',2);
xlabel(fig2,'time (s)', 'FontSize', 8, 'FontName', 'Times New Roman');
ylabel(fig2,'Pressure (mmHg)', 'FontSize', 8, 'FontName', 'Times New Roman');
xlim(fig2,trng);
title(fig2,'Pressures', 'FontSize', 8, 'FontName', 'Times New Roman');
plot(fig2,trng, [o_vals.PCWP, o_vals.PCWP],'color',[0.9290 0.7940 0.1250]);
if(isfield(targets,'PCWP'))
    plot(fig2,trng, [targets.PCWP, targets.PCWP],'color',[0.9290 0.7940 0.1250],'LineStyle','--');
end
if(isfield(targets,'PCWPmax'))
    plot(fig2,trng, [targets.PCWPmax, targets.PCWPmax],'color',[0.9290 0.7940 0.1250],'LineStyle','--');
end
if(isfield(targets,'RAPmax'))
    plot(fig2,trng, [targets.RAPmax, targets.RAPmax],'color',[0.4940 0.1840 0.5560],'LineStyle','--');
end
if(isfield(targets,'RAPmean'))
    plot(fig2,trng, [targets.RAPmean, targets.RAPmean],'color',[0.4940 0.1840 0.5560],'LineStyle','--');
end
ylim([0 14]);

fig4 = nexttile(4);
hold(fig4,"on");
plot(fig4,t, V_LV, 'color', [0.6350 0.0780 0.1840],'Linewidth', 2);
plot(fig4,t, V_RV, 'color', [0 0.4470 0.7410],'Linewidth', 2);
plot(fig4,t, V_LA, 'color', [0.8500 0.3250 0.0980],'Linewidth', 2);
plot(fig4,t, V_RA, 'color', [0.4940 0.1840 0.5560],'Linewidth', 2);
xlabel(fig4,'time (s)', 'FontSize', 8, 'FontName', 'Times New Roman');
ylabel(fig4,'Volume (mL)', 'FontSize', 8, 'FontName', 'Times New Roman');
xlim(fig4,trng);
title(fig4,'Volumes', 'FontSize', 8, 'FontName', 'Times New Roman');
if(isfield(targets,'LVEDV'))
    plot(fig4,trng, [targets.LVEDV, targets.LVEDV],'color',[0.6350 0.0780 0.1840],'LineStyle','--');
end
if(isfield(targets,'RVEDV'))
    plot(fig4,trng, [targets.RVEDV, targets.RVEDV],'color',[0 0.4470 0.7410],'LineStyle','--');
end
if(isfield(targets,'LVESV'))
    plot(fig4,trng, [targets.LVESV, targets.LVESV],'color',[0.6350 0.0780 0.1840],'LineStyle','--');
end
if(isfield(targets,'RVESV'))
    plot(fig4,trng, [targets.RVESV, targets.RVESV],'color',[0 0.4470 0.7410],'LineStyle','--');
end
if(isfield(targets,'LAVmax'))
    plot(fig4,trng, [targets.LAVmax, targets.LAVmax],'color',[0.8500 0.3250 0.0980],'LineStyle','--');
end
ylim([0 220]);

fig5 = nexttile(5);
plot(fig5,t, Q_m*60, t, Q_a*60, 'Linewidth', 2);
xlabel(fig5,'time (s)', 'FontSize', 8, 'FontName', 'Times New Roman');
ylabel(fig5,'Flow (ml/min)', 'FontSize', 8, 'FontName', 'Times New Roman');
xlim(fig5,trng);
title(fig5,'Flows', 'FontSize', 8, 'FontName', 'Times New Roman');

figtext = nexttile(3, [2 1]);
targetsfn = fieldnames(targets);
inputsfn = fieldnames(inputs);
text_position = [-0.15 0.95];
if GENDER == 1
    text(figtext,0.05,text_position(2),string(sprintf('Healthy 20-year- old Male')) + newline,"FontWeight","bold","FontSize",8,"FontName","Times New Roman");
elseif GENDER == 2
    text(figtext,0.05,text_position(2),string(sprintf('Healthy 20-year- old Female')) + newline,"FontWeight","bold","FontSize",8,"FontName","Times New Roman");
end
text_spacing = -0.028;
text_position(2) = text_position(2) + text_spacing;

% Sort costed targets highest to lowest, then output as such
itemized_costs = zeros(1,length(targetsfn));
for i = 1:length(itemized_costs)
    itemized_costs(i) = (o_vals.(targetsfn{i}) - targets.(targetsfn{i}))^2 / c.(targetsfn{i})^2 * w.(targetsfn{i}) * EX;
end

temp = sort(itemized_costs,'descend');
i_decreasing_costs = zeros(1,length(temp));
for i = 1:length(temp)
    i_decreasing_costs(i) = find(itemized_costs == temp(i));
end

for i = 1:length(i_decreasing_costs)
    i_itemized_costs = i_decreasing_costs(i);
    text_position(2) = text_position(2) + text_spacing;
    text(figtext,text_position(1),text_position(2), ...
        string(sprintf('%i) %s: %1.2f (%1.2f)   %1.0f%%', ...
        i, targetsfn{i_itemized_costs}, o_vals.(targetsfn{i_itemized_costs}), ...
        targets.(targetsfn{i_itemized_costs}), o_vals.(targetsfn{i_itemized_costs})/targets.(targetsfn{i_itemized_costs})*100 - 100)) + newline,'FontSize', 8,'FontName','Times New Roman','Interpreter','none');

end
text_valves = string(sprintf(['RFm (%%) = %1.2f  (<30 | 30-39 | 40-49 | >50)\nRFa (%%) = %1.2f  (<30 | 30-39 | 40-49 | >50)\nRFt (%%) = %1.2f  (<30 | 30-39 | >= 40)\nRFp (%%) = %1.2f  (<20 | 20-40 | >40)' ...
    '\nMVmg (mmHg) = %1.2f  (<5 | 5-10 | >10)\nAVpg (mmHg) = %1.2f  (<20 | 20-40 | >40)\nTVmg (mmHg) = %1.2f  (<5 | >=5)\nPVpg (mmHg) = %1.2f  (<36 | 36-64 | >64)\n'], ...
    o_vals.RF_m, o_vals.RF_a, o_vals.RF_t, o_vals.RF_p, ...
    o_vals.MVmg, o_vals.AVpg, o_vals.TVmg, o_vals.PVpg));
text(figtext, text_position(1), text_position(2) + text_spacing -0.15, text_valves, 'FontSize', 8, 'FontName', 'Times New Roman'); % FIXME -- idk what's wrong with the text spacing between the upper and lower

axis(figtext, "off");
% Set figure size to 7 inches wide, with the height adjusted proportionally
fig = gcf;

% Set the figure units and position for both display and saving
fig.Units = 'inches';
fig.Position = [1, 1, 7, 7/1.5]; % Set position with desired aspect ratio (width = 7, height = 7/1.5)

% Set the paper size and position for saving
fig.PaperUnits = 'inches';
fig.PaperPosition = [0, 0, 7, 7/1.5]; % Set paper size for saving the figure

% Set the resolution to 330 DPI and save the figure
print(fig, 'output_figure.png', '-dpng', '-r330');
