clear
% 读取数据
T1 = readtable('inputFGclassifier146WithLabel.CSV','VariableNamingRule','preserve');
T1.EF(T1.EF<50) = 50;
T2 = readtable("inputFGclassifierWithLabel.CSV",'VariableNamingRule','preserve');

RealPatient  = T1(:,2:19);     % 真实病人
VirtualPatient = T2(:,2:19);   % 虚拟病人

vars = RealPatient.Properties.VariableNames;
nVars = numel(vars);

% 布局
nCols = 4;
nRows = ceil(nVars/nCols);

figure;
for i = 1:nVars
    subplot(nRows,nCols,i);
    dataReal = RealPatient.(vars{i});
    dataVirt = VirtualPatient.(vars{i});
    
    % 统一 bin edges：根据两者的合并数据范围
    allData = [dataReal; dataVirt];
    edges = linspace(min(allData), max(allData), 30);  % 可以改成20或30，看你需要多细

    % 绘制直方图
    hold on;
    histogram(dataReal, edges, 'Normalization','probability', ...
        'FaceAlpha',0.5, 'FaceColor',[1 0.6 0.6]);
    histogram(dataVirt, edges, 'Normalization','probability', ...
        'FaceAlpha',0.5, 'FaceColor',[0.4 0.6 1]);
    hold off;

    title(vars{i}, 'Interpreter','none');
    xlabel(vars{i}, 'Interpreter','none');
    ylabel('Frequency');
end

legend({'Real','VAE'}, 'Position',[0.75 0.15 0.1 0.05]);
sgtitle('Distribution of Real vs Virtual Patients');


%%
clear
% 读取数据
T1 = readtable('input_VAE.csv','VariableNamingRule','preserve');
T1.Sex = 2-T1.Sex;
T2 = readtable("generated_fake_data_correlated.csv",'VariableNamingRule','preserve');

RealPatient  = T1(:,1:34);     % 真实病人
VirtualPatient = T2(:,1:34);   % 虚拟病人

vars = RealPatient.Properties.VariableNames;
nVars = numel(vars);

% 布局
nCols = 6;
nRows = ceil(nVars/nCols);

figure;
for i = 1:nVars
    subplot(nRows,nCols,i);
    dataReal = RealPatient.(vars{i});
    dataVirt = VirtualPatient.(vars{i});
    
    % 统一 bin edges：根据两者的合并数据范围
    allData = [dataReal];
    edges = linspace(min(allData), max(allData), 50);  % 可以改成20或30，看你需要多细

    % 绘制直方图
    hold on;
    histogram(dataReal, edges, 'Normalization','probability', ...
        'FaceAlpha',0.5, 'FaceColor',[1 0.6 0.6]);
    histogram(dataVirt, edges, 'Normalization','probability', ...
        'FaceAlpha',0.5, 'FaceColor',[0.4 0.6 1]);
    hold off;

    title(vars{i}, 'Interpreter','none');
    xlabel(vars{i}, 'Interpreter','none');
    ylabel('Frequency');
end

legend({'Real','VAE'}, 'Position',[0.75 0.15 0.1 0.05]);
sgtitle('Distribution of Real vs Virtual Patients');
