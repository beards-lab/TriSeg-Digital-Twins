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

figure(1); clf;
tiles = tiledlayout(2,4);
tiles.TileSpacing = 'compact';
tiles.Padding = 'tight';
set(gcf,'defaultLegendAutoUpdate','off','WindowState','maximized','Position', get(0, 'Screensize'));

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
plot(fig1,trng, [targets.SBP, targets.SBP],'color',[0.8500 0.5250 0.0980 0.5],'LineStyle','--','Linewidth', 2);
plot(fig1,trng, [targets.DBP, targets.DBP],'color',[0.8500 0.5250 0.0980 0.5],'LineStyle','--','Linewidth', 2);
if(isfield(targets,'PASP'))
    plot(fig1,trng, [targets.PASP, targets.PASP],'color',[0.3010 0.7450 0.9330 0.5],'LineStyle','--', 'Linewidth', 2);
    plot(fig1,trng, [targets.PADP, targets.PADP],'color',[0.3010 0.7450 0.9330 0.5],'LineStyle','--', 'Linewidth', 2);
end
if(isfield(targets,'RVSP'))
    plot(fig1,trng, [targets.RVSP, targets.RVSP],'color',[0.3010 0.7450 0.9330 0.5],'LineStyle','--', 'Linewidth', 2);
end
if(isfield(targets,'LVESP'))
    plot(fig1,trng, [targets.LVESP, targets.LVESP],'color',[0.6350 0.0780 0.1840 0.5],'LineStyle','--', 'Linewidth', 2);
end
if isfield(targets, 'RVEDP')
    plot(fig1,trng, [targets.RVEDP, targets.RVEDP],'color', [0 0.2470 0.8410 0.5],'LineStyle','--','Linewidth', 2);
end
if isfield(targets, 'LVEDP')
    plot(fig1,trng, [targets.LVEDP, targets.LVEDP],'color', [0.6350 0.0780 0.1840 0.5],'LineStyle','--','Linewidth', 2);
end
if isfield(targets,'P_RV_min')
    plot(fig1,trng, [targets.P_RV_min, targets.P_RV_min],'color', [0 0.2470 0.8410 0.5],'LineStyle','--','Linewidth', 2);
end
if isfield(targets,'P_LV_min')
    plot(fig1,trng, [targets.P_LV_min, targets.P_LV_min],'color', [0.6350 0.0780 0.1840 0.5],'LineStyle','--','Linewidth', 2);
end

fig2 = nexttile;
hold(fig2,"on");
plot(fig2,t, P_SV, 'color',[0.4660 0.6740 0.1880],'LineWidth',2);
plot(fig2,t, P_PV,'color',[0.9290 0.7940 0.1250],'LineWidth',2);
plot(fig2,t, P_RA, 'color',[0.4940 0.1840 0.5560],'LineWidth',2);
plot(fig2,t, P_LA,'color',[0.8500 0.3250 0.0980],'LineWidth',2);
xlim(fig2,trng); ylim([0.9*min([P_SV;P_PV;P_RA;P_LA]) 1.15*max([P_SV;P_PV;P_RA;P_LA])]);
plot(fig2,trng, [o_vals.PCWP, o_vals.PCWP],'color',[0.9290 0.7940 0.1250]);
if(isfield(targets,'PCWP'))
    plot(fig2,trng, [targets.PCWP, targets.PCWP],'color',[0.9290 0.7940 0.1250 0.5],'LineStyle','--','Linewidth', 2);
end
if(isfield(targets,'PCWPmax'))
    plot(fig2,trng, [targets.PCWPmax, targets.PCWPmax],'color',[0.9290 0.7940 0.1250 0.5],'LineStyle','--','Linewidth', 2);
end
if(isfield(targets,'RAPmax'))
    plot(fig2,trng, [targets.RAPmax, targets.RAPmax],'color',[0.4940 0.1840 0.5560 0.5],'LineStyle','--','Linewidth', 2);
end
if(isfield(targets,'RAPmean'))
    plot(fig2,trng, [targets.RAPmean, targets.RAPmean],'color',[0.4940 0.1840 0.5560 0.5],'LineStyle','--','Linewidth', 2);
end
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
if(isfield(targets,'LVEDV'))
    plot(fig3,trng, [targets.LVEDV, targets.LVEDV],'color',[0.6350 0.0780 0.1840 0.5],'LineStyle','--','Linewidth', 2);
end
if(isfield(targets,'RVEDV'))
    plot(fig3,trng, [targets.RVEDV, targets.RVEDV],'color',[0 0.4470 0.7410 0.5],'LineStyle','--','Linewidth', 2);
end
if(isfield(targets,'LVESV'))
    plot(fig3,trng, [targets.LVESV, targets.LVESV],'color',[0.6350 0.0780 0.1840 0.5],'LineStyle','--','Linewidth', 2);
end
if(isfield(targets,'RVESV'))
    plot(fig3,trng, [targets.RVESV, targets.RVESV],'color',[0 0.4470 0.7410 0.5],'LineStyle','--','Linewidth', 2);
end
if(isfield(targets,'LAVmax'))
    plot(fig3,trng, [targets.LAVmax, targets.LAVmax],'color',[0.8500 0.3250 0.0980 0.5],'LineStyle','--','Linewidth', 2);
end
lgd = legend(fig3,'LV', 'RV', 'LA', 'RA');
lgd.FontSize = 7;
axPosition = get(fig3, 'Position');
lgdPosition = get(lgd, 'Position');
set(gca, 'FontSize', 12);
set(gca, 'FontWeight', 'bold');
lgdPosition(1) = axPosition(1) + axPosition(3) - lgdPosition(3) - 0.008;
lgdPosition(2) = axPosition(2) + axPosition(4) - lgdPosition(4) - 0.023;
set(lgd, 'Position', lgdPosition);
if isfield(targets,"EF")
    Vol_text = "CO: " + num2str(round(CO,2)) + " ("...
        +targets.CO +") L/min"+ newline + "EF: " + num2str(round(EF_LV * 100)+2) + "% (" + targets.EF +"%)";
elseif isfield(targets,"CO")
    Vol_text = "CO: " + num2str(round(CO,2)) + " ("...
        +targets.CO +") L/min"+ newline + "EF: " + num2str(round(EF_LV * 100)+2) + "% ";
end
Vol_text_position = [0.548, 0.78, 0.4, 0.205];
annotation('textbox', Vol_text_position, 'String', Vol_text, 'FitBoxToText', 'on','FontWeight','bold','FontSize',10,'BackgroundColor',"w");
xlabel(fig3,'time (s)','FontSize',12); ylabel(fig3,'Volume (mL)','FontSize',12);

fig5 = nexttile(5);
plot(fig5,t, Q_t*60, t, Q_p*60, 'Linewidth', 2);
xlim(fig5,trng); ylim([1.02*min([Q_t*60;Q_p*60]) 1.2*max([Q_t*60;Q_p*60])])
hold(fig5,"on");
lgd = legend(fig5,'Tricuspid Valve', 'Pulmonary Valve');
lgd.FontSize = 7;
axPosition = get(fig5, 'Position');
lgdPosition = get(lgd, 'Position');
set(gca, 'FontSize', 12);
set(gca, 'FontWeight', 'bold');
lgdPosition(1) = axPosition(1) + axPosition(3) - lgdPosition(3) - 0.005;
lgdPosition(2) = axPosition(2) + axPosition(4) - lgdPosition(4) - 0.003;
set(lgd, 'Position', lgdPosition);
xlabel(fig5,'time (s)','FontSize',12); ylabel(fig5,'Flow (ml/min)','FontSize',12);
set(gca, "YTick", [-5e4 -4e4 -3e4 -2e4 -1e4 0 1e4 2e4 3e4 4e4 5e4]);
set(gca().YAxis, 'Exponent', 4);
set(gca, 'TickLength', [0.040 0.040])
box off

fig6 = nexttile(6);
plot(fig6,t, Q_m*60, t, Q_a*60, 'Linewidth', 2);  xlim(fig6,trng);ylim([1.02*min([Q_m*60;Q_a*60]) 1.2*max([Q_m*60;Q_a*60])]);
EAr_text = "E/A: " + num2str(round(E_A_ratio,2)) + " (" + targets.EAr+")";
EAr_text_position = [0.296, 0.27, 0.4, 0.195];
annotation('textbox', EAr_text_position, 'String', EAr_text, 'FitBoxToText', 'on','FontWeight','bold','FontSize',10,'BackgroundColor',"w");
lgd = legend(fig6,'Mitral Valve', 'Aortic Valve');
lgd.FontSize = 7;
axPosition = get(fig6, 'Position');
lgdPosition = get(lgd, 'Position');
set(gca, 'FontSize', 12);
set(gca, 'FontWeight', 'bold');
lgdPosition(1) = axPosition(1) + axPosition(3) - lgdPosition(3) - 0.008;
lgdPosition(2) = axPosition(2) + axPosition(4) - lgdPosition(4) - 0.018;
set(lgd, 'Position', lgdPosition);
xlabel(fig6,'time (s)','FontSize',12); ylabel(fig6,'Flow (ml/min)','FontSize',12);
set(gca, "YTick", [-5e4 -4e4 -3e4 -2e4 -1e4 0 1e4 2e4 3e4 4e4 5e4]);
set(gca().YAxis, 'Exponent', 4);
set(gca, 'TickLength', [0.040 0.040])
box off

fig7 = nexttile(7);
hold(fig7,"on");
plot(fig7,t, d_LW,'color', [0.6350 0.0780 0.1840],'LineWidth',2);
plot(fig7,t, d_SW, 'color', [0.4660 0.6740 0.1880],'LineWidth',2);
plot(fig7,t, d_RW, 'color', [0 0.4470 0.7410],'LineWidth',2);
xlim(fig7,trng); ylim([0 1.2*max([d_SW;d_LW;d_RW;targets.Hed_LW;targets.Hed_SW])])
if(isfield(targets,'Hed_LW'))
    plot(fig7,trng, [targets.Hed_LW, targets.Hed_LW],'color',[0.6350 0.0780 0.1840 0.5],'LineStyle','--','LineWidth',2);
    plot(fig7,trng, [o_vals.Hed_LW, o_vals.Hed_LW],'color',[0.6350 0.0780 0.1840 0.5],'LineStyle','-','LineWidth',2);
end
if(isfield(targets,'Hed_SW'))
    plot(fig7,trng, [targets.Hed_SW, targets.Hed_SW], 'color',[0.4660 0.6740 0.1880 0.5],'LineStyle','--','LineWidth',2);
    plot(fig7,trng, [o_vals.Hed_SW, o_vals.Hed_SW], 'color',[0.4660 0.6740 0.1880 0.5],'LineStyle','-','LineWidth',2);
end
if(isfield(targets,'Hed_RW'))
    plot(fig7,trng, [targets.Hed_RW, targets.Hed_RW], 'color', [0 0.4470 0.7410 0.5],'LineStyle','--','LineWidth',2);
    plot(fig7,trng, [o_vals.Hed_RW, o_vals.Hed_RW], 'color', [0 0.4470 0.7410 0.5],'LineStyle','-','LineWidth',2);
end
lgd = legend(fig7,'LV', 'Septum', 'RV');
lgd.FontSize = 7;
axPosition = get(fig7, 'Position');
lgdPosition = get(lgd, 'Position');
set(gca, 'FontSize', 12);
set(gca, 'FontWeight', 'bold');
lgdPosition(1) = axPosition(1) + axPosition(3) - lgdPosition(3) - 0.008;
lgdPosition(2) = axPosition(2) + axPosition(4) - lgdPosition(4) - 0.020;
set(lgd, 'Position', lgdPosition);
xlabel(fig7,'time (s)','FontSize',12); ylabel(fig7,'Thickness (cm)','FontSize',12);

figtext = nexttile(4, [2 1]);
text_position = [-0.2 0.95];
text(figtext,0.1,text_position(2),string(sprintf('PATIENT %i  WINDOW%i',PatID,ModelWin)) + newline,"FontWeight","bold","FontSize",12);
text_spacing = -0.033;

if Print_cost == 1
    text_position(2) = text_position(2) + text_spacing;
    text(figtext, 0.1, text_position(2),string(sprintf('\t\t\t\tTOTAL COST: ¥%1.2f',total_cost)) + newline,"FontWeight","bold","FontSize",12);
    % Sort costed targets highest to lowest, then output as such
    itemized_costs = zeros(1,length(targetsfn));
    for i = 1:length(itemized_costs)
        itemized_costs(i) = (o_vals.(targetsfn{i}) - targets.(targetsfn{i}))^2 /c.(targetsfn{i})^2 *w.(targetsfn{i})*EX;
    end
    temp = sort(itemized_costs,'descend');
    i_decreasing_costs = zeros(1,length(temp));
    for i = 1:length(temp)
        i_decreasing_costs(i) = find(itemized_costs == temp(i));
    end

    for i = 1:length(i_decreasing_costs)
        i_itemized_costs = i_decreasing_costs(i);
        if  (~ismember("AVpg",targetsfn{i_itemized_costs}))&&...
                (~ismember("PVpg",targetsfn{i_itemized_costs}))&&...
                (~ismember("MVmg",targetsfn{i_itemized_costs}))&&...
                (~ismember("AVr",targetsfn{i_itemized_costs}))&&...
                (~ismember("PVr",targetsfn{i_itemized_costs}))&&...
                (~ismember("MVr",targetsfn{i_itemized_costs}))&&...
                (~ismember("TVr",targetsfn{i_itemized_costs}))&&...
                (~ismember("DNA",targetsfn{i_itemized_costs}))&&...
                (~ismember("DNP",targetsfn{i_itemized_costs}))&&...
                (~ismember("FakeLV_m",targetsfn{i_itemized_costs}))
            text_position(2) = text_position(2) + text_spacing;
            text(figtext,text_position(1),text_position(2), ...
                string(sprintf('%i) %s: %1.2f (%1.2f)   %1.0f%% ($%1.2f)', ...
                i, targetsfn{i_itemized_costs}, o_vals.(targetsfn{i_itemized_costs}), ...
                targets.(targetsfn{i_itemized_costs}), o_vals.(targetsfn{i_itemized_costs})/targets.(targetsfn{i_itemized_costs})*100 - 100, ...
                itemized_costs(i_itemized_costs))) + newline,"FontWeight","bold",'FontSize', 12,'Interpreter','none');
        end
    end
else
    if patients(PatID).snapshots(ModelWin).Sex == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),string(sprintf('%1.0f yrs old',age))...
            +(" Male") + string(sprintf(' %1.0f cm',H)) + string(sprintf(' %1.0f kg',WT))...
            + newline,"FontWeight","bold","FontSize",12);
    elseif patients(PatID).snapshots(ModelWin).Sex == 2
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),string(sprintf('%1.0f yrs old',age))...
            +(" Female") + string(sprintf(' %1.0f cm',H)) + string(sprintf(' %1.0f kg',WT))...
            + newline,"FontWeight","bold","FontSize",12);
    end
    Pre = readtable("Precondition.xlsx",'VariableNamingRule','preserve');
    text_position(2) = text_position(2) + text_spacing;
    text(figtext,text_position(1),text_position(2),string(Pre.HFbasedonexam{PatID}) + (" (")+string(Pre.HF{PatID})+(")")+ newline,"FontWeight","bold","FontSize",12);
    if ~isempty(Pre.BNP{PatID})
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),string(sprintf('BNP: %1.0f pg/ml (<100 pg/ml)',str2double(Pre.BNP{PatID})))+newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.HTN(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Hypertension")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.PHD(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Cardiopulmonary disease")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.PH(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Pulmonary hypertension")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.CAD(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Coronary artery disease")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.MI(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Myocardial infarction")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.AF(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Atrial fibrillation")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.Arrhythmia(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Arrhythmia")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.SHD(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Structural heart disease")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("Septal myomectomy")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Septal myomectomy")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.Shunt(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Shunt")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("Hypertrophic Cardiomyopathy")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Hypertrophic Cardiomyopathy")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("Noncompaction Cardiomyopathy")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Noncompaction Cardiomyopathy")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("Cardiomyopathy, restrictive")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Cardiomyopathy, restrictive")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("Cardiomyopathy, dilated")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Cardiomyopathy, dilated")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("Cardiomyopathy of undetermined type")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Cardiomyopathy of undetermined type")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("Amyloidosis")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Amyloidosis")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("Pericardium disease")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Pericardium disease")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.Pericardiectomy(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Pericardiectomy")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.HLD(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Hyperlipidemia")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.DM(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Diabetes mellitus")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.COPD(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Chronic obstructive pulmonary disease")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.CTD(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Connective tissue disease")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.LTh(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Hypothyroidism")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.HTh(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Hyperthyroidism")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("Liver failure")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Liver failure")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("Renal failure")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Renal failure")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("Respiratory failure")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Respiratory failure")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("AD/AA")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Aortic dissection/Aneurysm")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("Heart transplant")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Post heart tranplant")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("Kidney transplant")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Post kidney tranplant")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("Liver transplant")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Post liver tranplant")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("Lung transplant")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Post lung tranplant")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("Pancreas transplant")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Post pancreas tranplant")+ newline,"FontWeight","bold","FontSize",12);
    end

    if Pre.ICD(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Presence of ICD")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.LVAD(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Presence of LVAD")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.Pacemaker(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Presence of Pacemaker")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.Hypervolemia(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Hypervolemia")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.Hypovolemia(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Hypovolemia")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("HFrEF-recovered")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("EF recovered")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("Valvular disease")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("Valvular heart disease")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("MV repair")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("MV repair")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("AV repair")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("AV repair")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("TV repair")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("TV repair")+ newline,"FontWeight","bold","FontSize",12);
    end
    if Pre.("PV repair")(PatID) == 1
        text_position(2) = text_position(2) + text_spacing;
        text(figtext,text_position(1),text_position(2),("PV repair")+ newline,"FontWeight","bold","FontSize",12);
    end

    text_position(2) = text_position(2) + text_spacing;
    text(figtext,text_position(1),text_position(2),string(Pre.Notes{PatID}) ,"FontWeight","bold","FontSize",12);
end

if isfield(targets,"MVr")
    text_position(2) = text_position(2) + text_spacing;
    if length(patients(PatID).snapshots(ModelWin).MVr_str) >30
        description = patients(PatID).snapshots(ModelWin).MVr_str(1:7);
    else
        description = patients(PatID).snapshots(ModelWin).MVr_str;
    end
    text_valves = string(sprintf('RF of MV = %1.2f%% (',o_vals.RF_m))+description+")";
    text(figtext, text_position(1), text_position(2) + text_spacing, text_valves, "FontWeight","bold",'FontSize', 12);
end
if isfield(targets,"AVr")
    text_position(2) = text_position(2) + text_spacing;
    if length(patients(PatID).snapshots(ModelWin).AVr_str) >30
        description = patients(PatID).snapshots(ModelWin).AVr_str(1:7);
    else
        description = patients(PatID).snapshots(ModelWin).AVr_str;
    end
    text_valves = string(sprintf('RF of AV = %1.2f%% (',o_vals.RF_a))+description+")";
    text(figtext, text_position(1), text_position(2) + text_spacing, text_valves, "FontWeight","bold",'FontSize', 12);
end
if isfield(targets,"TVr")
    text_position(2) = text_position(2) + text_spacing;
    if length(patients(PatID).snapshots(ModelWin).TVr_str) >30
        description = patients(PatID).snapshots(ModelWin).TVr_str(1:7);
    else
        description = patients(PatID).snapshots(ModelWin).TVr_str;
    end
    text_valves = string(sprintf('RF of TV = %1.2f%% (',o_vals.RF_t))+description+")";
    text(figtext, text_position(1), text_position(2) + text_spacing, text_valves, "FontWeight","bold",'FontSize', 12);
end
if isfield(targets,"PVr")
    text_position(2) = text_position(2) + text_spacing;
    if length(patients(PatID).snapshots(ModelWin).PVr_str) >30
        description = patients(PatID).snapshots(ModelWin).PVr_str(1:7);
    else
        description = patients(PatID).snapshots(ModelWin).PVr_str;
    end
    text_valves = string(sprintf('RF of PV = %1.2f%% (',o_vals.RF_p))+description+")";
    text(figtext, text_position(1), text_position(2) + text_spacing, text_valves, "FontWeight","bold",'FontSize', 12);
end
if isfield(targets,"AVpg")
    text_position(2) = text_position(2) + text_spacing;
    text_valves = string(sprintf('AVpg = %1.2f mmHg (%1.2f mmHg)',o_vals.AVpg,targets.AVpg));
    text(figtext, text_position(1), text_position(2) + text_spacing, text_valves, "FontWeight","bold",'FontSize', 12);
end
if isfield(targets,"PVpg")
    text_position(2) = text_position(2) + text_spacing;
    text_valves = string(sprintf('PVpg = %1.2f mmHg (%1.2f mmHg)',o_vals.PVpg,targets.PVpg));
    text(figtext, text_position(1), text_position(2) + text_spacing, text_valves, "FontWeight","bold",'FontSize', 12);
end
if isfield(targets,"MVmg")
    text_position(2) = text_position(2) + text_spacing;
    text_valves = string(sprintf('MVmg = %1.2f mmHg (%1.2f mmHg)',o_vals.MVmg,targets.MVmg));
    text(figtext, text_position(1), text_position(2) + text_spacing, text_valves, "FontWeight","bold",'FontSize', 12);
end
if isfield(targets,"TVmg")
    text_position(2) = text_position(2) + text_spacing;
    text_valves = string(sprintf('TVmg = %1.2f mmHg (%1.2f mmHg)',o_vals.TVmg,targets.TVmg));
    text(figtext, text_position(1), text_position(2) + text_spacing, text_valves, "FontWeight","bold",'FontSize', 12);
end

axis(figtext,"off");
