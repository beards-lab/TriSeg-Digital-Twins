%% Script Summary:
% This script generates a 3D Volcano plot.

% Created and designed by Feng Gu
% Last modified: 10/29/2024

clear
%% Load input features
DT = readtable("LastWindowData.csv",'VariableNamingRule','preserve');
DTT = readtable("LastWindowDTinfo.csv",'VariableNamingRule','preserve');
finalslot = ~ismember(DTT.("RV Dysfunction"),{'0'}) & ~isnan(DTT.C_SA);
DT = DT(finalslot,[(5:7) (9:26) (28:35) (38:54) 59]);
DTz = array2table(zscore(table2array(DT(:,2:end)),0,'omitnan'),...
    'VariableNames', DT.Properties.VariableNames(2:end));
DT= [DT(:,1) DTz];
label = DTT(finalslot,2);
DTT = DTT(finalslot,3:83);
DTT = array2table(zscore(table2array(DTT),0,'omitnan'), ...
    'VariableNames', DTT.Properties.VariableNames);
RVoutcome = ismember(label.("RV Dysfunction"),{'Moderate'}) | ismember(label.("RV Dysfunction"),{'Severe'});
%% Perform group comparison on DT variables
varNames_DT = DT.Properties.VariableNames;
nVars_DT = width(DT);

p_values_DT = zeros(1, nVars_DT);
mean_diff_DT = zeros(1, nVars_DT);
ci_lower_DT = zeros(1, nVars_DT);
ci_upper_DT = zeros(1, nVars_DT);

for i = 1:nVars_DT
    var_data = DT{:, i};
    if strcmp(varNames_DT{i}, 'Sex')
        % Sex is categorical: perform chi-square test
        [~, p] = crosstab(var_data, RVoutcome);
        p_values_DT(i) = p;
        mean_diff_DT(i) = NaN;
        ci_lower_DT(i) = NaN;
        ci_upper_DT(i) = NaN;
    else
        % Continuous variable: two-sample t-test with CI calculation
        group0 = var_data(RVoutcome == 0);
        group1 = var_data(RVoutcome == 1);
        n0 = sum(~isnan(group0));
        n1 = sum(~isnan(group1));
        m0 = mean(group0, 'omitnan');
        m1 = mean(group1, 'omitnan');
        s0 = std(group0, 'omitnan');
        s1 = std(group1, 'omitnan');
        se = sqrt(s0^2/n0 + s1^2/n1);
        df = (s0^2/n0 + s1^2/n1)^2 / ((s0^2/n0)^2/(n0-1) + (s1^2/n1)^2/(n1-1)); % Welch-Satterthwaite
        t_crit = tinv(0.975, df);
        ci_lower_DT(i) = (m1 - m0) - t_crit * se;
        ci_upper_DT(i) = (m1 - m0) + t_crit * se;
        mean_diff_DT(i) = m1 - m0;
        [~, p] = ttest2(group1, group0);
        p_values_DT(i) = p;
    end
end

% Bonferroni correction
p_adj_DT = min(p_values_DT * nVars_DT, 1);

% Summary table
results_DT = table(varNames_DT', mean_diff_DT', ci_lower_DT', ci_upper_DT', p_values_DT', p_adj_DT', ...
    'VariableNames', {'Variable', 'MeanDifference', 'CI_Lower', 'CI_Upper', 'P_Value', 'P_Adj'});

%% Perform group comparison on DTT variables
varNames_DTT = DTT.Properties.VariableNames;
nVars_DTT = width(DTT);

p_values_DTT = zeros(1, nVars_DTT);
mean_diff_DTT = zeros(1, nVars_DTT);
ci_lower_DTT = zeros(1, nVars_DTT);
ci_upper_DTT = zeros(1, nVars_DTT);

for i = 1:nVars_DTT
    var_data = DTT{:, i};
    % Assume all variables in DTT are continuous
    group0 = var_data(RVoutcome == 0);
    group1 = var_data(RVoutcome == 1);
    n0 = sum(~isnan(group0));
    n1 = sum(~isnan(group1));
    m0 = mean(group0, 'omitnan');
    m1 = mean(group1, 'omitnan');
    s0 = std(group0, 'omitnan');
    s1 = std(group1, 'omitnan');
    se = sqrt(s0^2/n0 + s1^2/n1);
    df = (s0^2/n0 + s1^2/n1)^2 / ((s0^2/n0)^2/(n0-1) + (s1^2/n1)^2/(n1-1)); % Welch-Satterthwaite
    t_crit = tinv(0.975, df);
    ci_lower_DTT(i) = (m1 - m0) - t_crit * se;
    ci_upper_DTT(i) = (m1 - m0) + t_crit * se;
    mean_diff_DTT(i) = m1 - m0;
    [~, p] = ttest2(group1, group0);
    p_values_DTT(i) = p;
end

% Bonferroni correction
p_adj_DTT = min(p_values_DTT * nVars_DTT, 1);

% Summary table
results_DTT = table(varNames_DTT', mean_diff_DTT', ci_lower_DTT', ci_upper_DTT', p_values_DTT', p_adj_DTT', ...
    'VariableNames', {'Variable', 'MeanDifference', 'CI_Lower', 'CI_Upper', 'P_Value', 'P_Adj'});
%% Forest plot for DT comparison
% Calculate -log10(P) and assign color by direction
logP = -log10(results_DT.P_Value);
colors = repmat([0.2 0.6 0.8], height(results_DT), 1); % blue by default
colors(results_DT.MeanDifference < 0, :) = repmat([0.85 0.3 0.3], sum(results_DT.MeanDifference < 0), 1); % red for negative

% Sort by -logP
[~, sort_idx] = sort(logP, 'descend');
sorted_results = results_DT(sort_idx, :);
sorted_logP = logP(sort_idx);
sorted_colors = colors(sort_idx, :);

% Horizontal bar plot: variable name on x, -log10(P) on y
figure;
bar(sorted_logP, 'FaceColor', 'flat');
colormap([0.85 0.3 0.3; 0.2 0.6 0.8]);hold on;
for i = 1:length(sorted_logP)
    bar(i, sorted_logP(i), 'FaceColor', sorted_colors(i, :)); % Color by direction
end

xticks(1:height(sorted_results));
xticklabels(sorted_results.Variable);
xtickangle(45);
ylabel('-log_{10}(P-value)');
title('Group Comparison: DT Variables');
grid on;
%% Volcano plot
% Compute -log10(P) and prepare data
mean_diff = results_DTT.MeanDifference;
logP = -log10(results_DTT.P_Value);
p_thresh = 0.05/nVars_DTT;
logP_thresh = -log10(p_thresh);

% Color coding: gray by default
colors = repmat([0.6 0.6 0.6], height(results_DTT), 1); 

% Significant up (blue)
sig_up = (mean_diff > 0) & (results_DTT.P_Value < p_thresh);
colors(sig_up, :) = repmat([0.2 0.6 0.8], sum(sig_up), 1); 

% Significant down (red)
sig_down = (mean_diff < 0) & (results_DTT.P_Value < p_thresh);
colors(sig_down, :) = repmat([0.85 0.3 0.3], sum(sig_down), 1); 

% Plot
figure;
scatter(mean_diff, logP, 40, colors, 'filled');
yline(logP_thresh, '--k', 'Adjusted P = 0.05');
xlabel('Mean Difference (Group 1 - Group 0)');
ylabel('-log_{10}(P-value)');
title('Volcano Plot: Group Comparison (DT)');
grid on;
box on;

hold on;
for i = 1:height(results_DTT)
    if results_DTT.P_Value(i) < p_thresh
        text(mean_diff(i)-0.15, logP(i), results_DTT.Variable{i}, ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
    end
end
%% Heatmap
mean_diff = results_DT.MeanDifference;
logP = -log10(results_DT.P_Value);

% Volcano quadrant logic
quad = repmat(2, height(results_DT), 1); % default gray
quad(mean_diff > 0 & results_DT.P_Value < 0.05/nVars_DT) = 1; % red (up)
quad(mean_diff < 0 & results_DT.P_Value < 0.05/nVars_DT) = 3; % blue (down)

% Sorting: first by quadrant, then by distance from center
dist = mean_diff .* logP;
[~, sort_idx] = sortrows([quad dist], [1 -2]);  % quad ascending, dist descending
sorted_vars = results_DT.Variable(sort_idx);

% Step 2: Reorder DTT_z accordingly
DT_sorted = DT(:,sorted_vars);

% Step 3: Reorder columns: group 0 first, group 1 later
group0_idx = find(RVoutcome == 0);
group1_idx = find(RVoutcome == 1);

data_group0 = DT_sorted{group0_idx, :}';
data_group1 = DT_sorted{group1_idx, :}';

% Step 4: Draw heatmap

numG0 = length(group0_idx);
numG1 = length(group1_idx);
totalUnits = numG0 + numG1;
Widthperunit = 0.75 / totalUnits;
figure;
subplot('Position', [0.2,0.05,Widthperunit*numG0,0.9]);
hm1 = heatmap(data_group0, ...
    'Colormap', interp1(1:size(redbluecmap,1), redbluecmap, linspace(1, size(redbluecmap,1), 256)), ...
    'ColorLimits', [-2.5 2.5],'MissingDataColor', [0.8 0.8 0.8]);

hm1.YDisplayLabels = sorted_vars;  % show variable names
hm1.XDisplayLabels = repmat({''}, 1, size(data_group0 ,2)); % optional: hide column labels
hm1.FontSize = 10;
hm1.GridVisible = 'off'; % Disable the grid for a cleaner look
hm1.ColorbarVisible = 'off';

subplot('Position', [0.2 + Widthperunit*numG0 + 0.003,0.05, Widthperunit*numG1,0.9]);

hm2 = heatmap(data_group1, ...
    'Colormap', interp1(1:size(redbluecmap,1), redbluecmap, linspace(1, size(redbluecmap,1), 256)), ...
    'ColorLimits', [-2.5 2.5],'MissingDataColor', [0.8 0.8 0.8]);

hm2.YDisplayLabels = repmat({''}, 1, size(data_group1 ,1));  % show variable names
hm2.XDisplayLabels = repmat({''}, 1, size(data_group1 ,2)); % optional: hide column labels
hm2.FontSize = 10;
hm2.GridVisible = 'off'; % Disable the grid for a cleaner look
hm2.ColorbarVisible = 'off';





