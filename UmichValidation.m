PredictorsT.Gender(PredictorsT.Gender == 0) = 2; 
PredictorsT.AVr_D (isnan(PredictorsT.AVr_D)) = 1;  
PredictorsT.MVr_D (isnan(PredictorsT.MVr_D)) = 1;  
PredictorsT.TVr_D (isnan(PredictorsT.TVr_D)) = 1;
PredictorsT.PVr_D (isnan(PredictorsT.PVr_D)) = 1;
load LassoRV.mat
  Raw_Lasso_RVEDV = PredictorsT.Gender.*k_RVEDV(1) + PredictorsT.Age.*k_RVEDV(2) +...
        PredictorsT.Height.*k_RVEDV(3) + PredictorsT.Height.*PredictorsT.Weight.*k_RVEDV(4) +...
        PredictorsT.HR.*k_RVEDV(5) + PredictorsT.SBP_D.*k_RVEDV(6) +...
        PredictorsT.PASP_D.*k_RVEDV(7) + PredictorsT.PADP_D.*k_RVEDV(8) +...
        PredictorsT.CO_D.*k_RVEDV(9) + PredictorsT.AVr_D.* k_RVEDV(10) + ...
        PredictorsT.TVr_D.* k_RVEDV(11) + PredictorsT.PVr_D.* k_RVEDV(12) + ...
        b_RVEDV;
  C_Lasso_RVEDV = Raw_Lasso_RVEDV*ck_RVEDV+cb_RVEDV;

  Raw_Lasso_RVEF = PredictorsT.Gender.* k_RVEF(1) + PredictorsT.Age.* k_RVEF(2) + ...
       PredictorsT.SBP_D.* k_RVEF(3) + PredictorsT.PADP_D.* k_RVEF(4) + ...
      PredictorsT.CO_D.* k_RVEF(5) + PredictorsT.EF_D.* k_RVEF(6) + ...
      PredictorsT.TVr_D.* k_RVEF(7) + PredictorsT.PVr_D.* k_RVEF(8) + ...
      b_RVEF;
  C_Lasso_RVEF = Raw_Lasso_RVEF * ck_RVEF + cb_RVEF;
   C_Lasso_RVEF( C_Lasso_RVEF<5) = 5;
  %%
X = C_Lasso_RVEF;

Y = 100.*(PredictorsT.RVEDV_D-PredictorsT.RVESV_D)./PredictorsT.RVEDV_D;

mdl_filtered = fitlm(X, Y,'linear','RobustOpts','on');
disp(mdl_filtered);
% 绘制真实值 vs 预测值的散点图

figure(100); clf;
scatter(X, Y, 'b', 'filled');
hold on;
plot([min(Y), max(Y)], [min(Y), max(Y)], 'r--', 'LineWidth', 2); % 理想拟合线
adj_r2 = mdl_filtered.Rsquared.Ordinary; % Adjusted R²
beta = mdl_filtered.Coefficients.Estimate(2); % 斜率
adj_r = sign(beta) * sqrt(adj_r2); % 计算 Adjusted R
p_value = mdl_filtered.Coefficients.pValue(2); % 斜率的 p 值

% 在图上显示 Adjusted R 和 p 值
annotation('textbox', [0.25, 0.8, 0.5, 0.1], ... % 位置：[x, y, width, height]
    'String', sprintf('R = %.2f\np = %.3g', adj_r, p_value), ...
    'FontSize', 16, 'FontName', 'Arial', ...
    'EdgeColor', 'none', ... % 取消边框
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
hold off;
pbaspect([1,1,1])
% 添加图例和标签
xlabel('Lasso Y');
ylabel('Mesurement Y');
title('Simulation: Predicted vs True');
grid on;
