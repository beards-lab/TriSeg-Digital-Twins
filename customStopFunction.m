% 定义自定义输出函数
function [state, options, optchanged] = customStopFunction(options, state, flag)
    optchanged = false; % 默认不修改选项
    targetValue = 50000; % 目标值
    tolerance = 1e-3; % 容差，防止浮点误差

    % 获取当前种群中最优适应度值
    currentBest = state.Best(end);

    % 动态调整停止条件
    if abs(currentBest - targetValue) < tolerance
        % 如果当前最优值等于目标值，继续运行
        disp('Fitness equals 50000, extending search...');
        options.MaxStallGenerations = 200; % 延长最大停滞代数
        optchanged = true;
    elseif currentBest < targetValue
        % 如果当前最优值小于目标值，设定较短的停滞代数
        disp('Fitness below 50000, setting short stall...');
        options.MaxStallGenerations = 8; % 设置较短的停滞代数
        optchanged = true;
    end
end