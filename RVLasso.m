load UWcohort.mat
% 处理 UWpatients.EF 中的空值，使用 UWpatients.EF_LV 代替
missing_EF = isnan(UWpatients.EF); % 找到 EF 中的空值
UWpatients.EF(missing_EF) = UWpatients.EF_LV(missing_EF); % 用 EF_LV 填充

% 处理 UWpatients.CO 中的空值，使用 UWpatients.CO_TD 代替
missing_CO = isnan(UWpatients.CO); % 找到 CO 中的空值
UWpatients.CO(missing_CO) = UWpatients.CO_TD(missing_CO); % 用 CO_TD 填充

UWpatients.HR = nanmean([UWpatients.HR UWpatients.HR_Echo UWpatients.heartRate],2);
% 检查是否还有 NaN
disp(['Remaining NaNs in EF: ', num2str(sum(isnan(UWpatients.EF)))]);
disp(['Remaining NaNs in CO: ', num2str(sum(isnan(UWpatients.CO)))]);

% 提取特征矩阵 X 和目标变量 Y
Y = UWpatients.Hed_RW;
% X = [UWpatients.Height, UWpatients.TVr, UWpatients.PVr, UWpatients.Hed_LW, ...
%      UWpatients.CO, UWpatients.EF, UWpatients.SBP, UWpatients.RAPmax,UWpatients.PPAs];
dataSubset = UWpatients(:, [5 6 (51:55) (76:87) 94 (213:215) 219 257 259 261 262]);  

% 找出完全没有 NaN 的列
validColumns = all(~isnan(table2array(dataSubset)), 1);  

% 仅保留不包含 NaN 的列
 X = dataSubset(:, validColumns);  



% 使用 MATLAB 自带的 lasso 进行 L1 正则化回归
lambdaVals = logspace(-4, 1, 100); % 生成 100 个 lambda 取值（范围从 10^-4 到 10^1）
[B, FitInfo] = lasso(X{:,:}, Y, 'Lambda', lambdaVals, 'CV', 10);

% 输出最佳 lambda 值及对应的回归系数
optimalLambda = FitInfo.LambdaMinMSE; % 交叉验证选择的最佳正则化参数
betaOptimal = B(:, FitInfo.IndexMinMSE); % 对应的 Lasso 估计系数

% 显示结果
disp(['Optimal Lambda: ', num2str(optimalLambda)]);
disp('Lasso 估计的系数:');
disp(betaOptimal);
% 显示选中的列名
disp('Selected columns for X:');
disp(X.Properties.VariableNames(~betaOptimal == 0));
% 画出 Lambda 选择路径
lassoPlot(B, FitInfo, 'PlotType', 'Lambda', 'XScale', 'log');
title('Lasso Regularization Path');
lassoPlot(B,FitInfo,PlotType="CV");
legend("show")

% 计算预测值
Y_pred = X{:,:} * betaOptimal + FitInfo.Intercept(FitInfo.IndexMinMSE); % 加上截距项
validIdx = ~isnan(Y) & ~isnan(Y_pred);
Y_clean = Y(validIdx);

Y_pred_clean = Y_pred(validIdx);
% 计算相关系数
R = corr(Y_pred_clean, Y_clean)
Y_pred_clean = (std(Y_clean)./std(Y_pred_clean)).*Y_pred_clean+mean(Y_clean)-(std(Y_clean)./std(Y_pred_clean))*mean(Y_pred_clean);

Y_pred_clean = X;
Y_clean = Y;


% 绘制真实值 vs 预测值的散点图

figure(100); clf;
scatter(Y_pred_clean, Y_clean, 'b', 'filled');
hold on;
plot([min(Y), max(Y)], [min(Y), max(Y)], 'r--', 'LineWidth', 2); % 理想拟合线
hold off;
pbaspect([1,1,1])
% 添加图例和标签
xlabel('Simulation Y');
ylabel('Mesurement Y');
title('Simulation: Predicted vs True');
grid on;

