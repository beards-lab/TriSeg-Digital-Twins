slides = [RVEDP_i, Qa_pos_start, Qa_pos_end];
slides_titles = {"End of Diastole", "End of Isovolumic Contraction","End of Systole"};

% visualize the time course
plotArc = @(x0, radius,arc)  [-radius.*cos(arc) - x0;  radius.*sin(arc)];
% give the different layers different color based on the sigmaact
minData = 0;
maxData = 500;
% Custom cmap
numColors = 256;

lightBlue = [0,0.45,0.74];
lightGreen = [0, 1, 0.8];  %
yellow    = [1, 1, 0];
red = [0.84,0.08,0.18];
% Chatgpt teach me how to make a custom Jet, thank you chatgpt!
numPerSection = floor(numColors / 3);

cmap = zeros(numColors, 3);

% fill color for 3 columns
% blue to green
for i = 1:numPerSection
    cmap(i,:) = lightBlue + (lightGreen - lightBlue) * ((i-1) / (numPerSection-1));
end
% green to yellow
for i = 1:numPerSection
    cmap(numPerSection+i,:) = lightGreen + (yellow - lightGreen) * ((i-1) / (numPerSection-1));
end
% yellow to red
for i = 1:(numColors - 2 * numPerSection)
    cmap(2*numPerSection+i,:) = yellow + (red - yellow) * ((i-1) / (numPerSection-1));
end

% fill other
if numColors > 3 * numPerSection
    remainder = numColors - 3 * numPerSection;
    for i = 1:remainder
        cmap(3*numPerSection+i,:) = yellow + (red - yellow) * (i / remainder);
    end
end

Nsigma_LV = (sigma_LV - minData) / (maxData - minData); % Normalize
Nsigma_SEP = (sigma_SEP - minData) / (maxData - minData);
Nsigma_RV = (sigma_RV - minData) / (maxData - minData);
colorIndex_LV = 1 + floor(Nsigma_LV * (size(cmap, 1) - 1));
colorIndex_SEP = 1 + floor(Nsigma_SEP * (size(cmap, 1) - 1));
colorIndex_RV = 1 + floor(Nsigma_RV * (size(cmap, 1) - 1));
colorIndex_LV = min(colorIndex_LV, 256);
colorIndex_SEP = min(colorIndex_SEP, 256);
colorIndex_RV = min(colorIndex_RV, 256);
colors_LV = cmap(colorIndex_LV, :);
colors_SEP = cmap(colorIndex_SEP, :);
colors_RV = cmap(colorIndex_RV, :);
figure(102);
clf;
% set(gcf,'defaultLegendAutoUpdate','off','position',[2000,0,1750,1000]); hold on;
set(gcf,'defaultLegendAutoUpdate','off','WindowState','maximized','Position', get(0, 'Screensize'));
positions = [0.05, 0.3, 0.27, 0.4; 0.36, 0.3, 0.27, 0.4; 0.67, 0.3, 0.27, 0.4];
for i = 1:length(slides)
    xm_lv = xm_LV(slides(i));
    xm_sep = xm_SEP(slides(i));
    xm_rv = xm_RV(slides(i));
    ym_all = ym(slides(i));
    h_lv = d_LW(slides(i));
    h_rv = d_RW(slides(i));
    h_sep = d_SW(slides(i));
    Cm_lv = Cm_LV(slides(i));
    Cm_rv = Cm_RV(slides(i));
    Cm_sep = Cm_SEP(slides(i));
    R_lv = 1./Cm_lv;
    R_sep = 1./Cm_sep;
    R_rv = 1./Cm_rv;
    r_inner_LV = -r_LV(slides(i));

    r_inner_RV =  r_RV(slides(i));
    r_outer_LV =  r_inner_LV-h_lv;
    r_outer_RV = r_inner_RV+h_rv;
    if xm_sep >= 0
        r_inner_SEP =  r_SEP(slides(i));
        r_outer_SEP = r_inner_SEP+h_sep;
    else
        r_inner_SEP = - r_SEP(slides(i));
        r_outer_SEP = r_inner_SEP-h_sep;
    end
    
    subplot('Position', positions(i, :)); hold on;
    set(gca, 'Color', 'none');
    x0_lv = xm_lv - R_lv;
    x0_sep = xm_sep - R_sep;
    x0_rv = xm_rv - R_rv;

    a_lv = atan2(ym_all, x0_lv);
    if xm_sep>0
        a_sep = atan2(ym_all, -x0_sep);
    else
        a_sep = atan2(ym_all, x0_sep);
    end
    a_rv = atan2(ym_all, -x0_rv);

    a_start_lv = 0*pi/2 - a_lv;
    a_stop_lv = 0*pi/2 + a_lv;

    a_start_sep = 0*pi/2 -a_sep;
    a_stop_sep = 0*pi/2 + a_sep;

    a_start_rv = 0*pi/2 -a_rv;
    a_stop_rv = 0*pi/2 + a_rv;
    set(gca, 'FontSize', 12);
    set(gca, 'FontWeight', 'bold');

    % depends on the total width too!
    totalwidth = max(xm_rv - xm_lv);
    range = totalwidth/2*1.5;
    axis([-range range -range range])
    % axis tight
    axis equal;
    set(gca, "YTick", [-10 -5, 0, 5 10]);
    set(gca, "XTick", [-10 -5, 0, 5 10]);
    set(gca, 'TickLength', [0.040 0.040])
    axis off;
    if i == 1
        set(gca, "YTicklabel",cellstr({'-10cm','-5cm', '0cm', '5cm','10cm'}));
    else
        % no tick labels for insides
        set(gca, "YTicklabel", []);
    end
    % title({slides_titles{i}});
    set(gca, "XTicklabel", cellstr({'-10cm','-5cm', '0cm', '5cm','10cm'}));


    % % Outer LV
    arc_lv = linspace(a_start_lv, a_stop_lv, 20);
    xy_lv = plotArc(x0_lv, r_outer_LV, arc_lv);
    plot(xy_lv(1, :), xy_lv(2, :), 'LineWidth',2, 'Color', colors_LV(slides(i),:));
    % midline LV
    arc_lv = linspace(a_start_lv, a_stop_lv, 20);
    xy_lv = plotArc(x0_lv, R_lv, arc_lv);
    plot(xy_lv(1, :), xy_lv(2, :), '--', 'LineWidth',0.5, 'Color', colors_LV(slides(i),:));
    % inner LV
    arc_lv = linspace(a_start_lv, a_stop_lv, 20);
    xy_lv = plotArc(x0_lv,r_inner_LV, arc_lv);
    plot(xy_lv(1, :), xy_lv(2, :), 'LineWidth',2, 'Color', colors_LV(slides(i),:))


    % LV side
    arc_sep = linspace(a_start_sep, a_stop_sep, 20);
    xy_sep = plotArc(x0_sep, r_inner_SEP, arc_sep);
    plot(xy_sep(1, :), xy_sep(2, :), 'LineWidth',2, 'Color', colors_SEP(slides(i),:))
    % mid SEP
    arc_sep = linspace(a_start_sep, a_stop_sep, 20);
    xy_sep = plotArc(x0_sep, R_sep, arc_sep);
    plot(xy_sep(1, :), xy_sep(2, :), '--', 'Color', colors_SEP(slides(i),:))
    % RV side SEP
    arc_sep = linspace(a_start_sep, a_stop_sep, 20);
    xy_sep = plotArc(x0_sep, r_outer_SEP, arc_sep);
    plot(xy_sep(1, :), xy_sep(2, :), 'LineWidth',2, 'Color', colors_SEP(slides(i),:))

    % RV inner
    arc_rv = linspace(a_start_rv, a_stop_rv, 20);
    xy_rv = plotArc(x0_rv, r_inner_RV, arc_rv);
    plot(xy_rv(1, :), xy_rv(2, :), 'LineWidth',2, 'Color', colors_RV(slides(i),:))
    % mid
    arc_rv = linspace(a_start_rv, a_stop_rv, 20);
    xy_rv = plotArc(x0_rv, R_rv, arc_rv);
    plot(xy_rv(1, :), xy_rv(2, :), '--', 'Color', colors_RV(slides(i),:))
    % RV outer
    arc_rv = linspace(a_start_rv, a_stop_rv, 20);
    xy_rv = plotArc(x0_rv, r_outer_RV, arc_rv);
    plot(xy_rv(1, :), xy_rv(2, :), 'LineWidth',2, 'Color', colors_RV(slides(i),:))
    pbaspect([1,1,1]);
      scale_length = 5; % 设置比例尺的长度
    scale_position = [-8, -8]; % 设置比例尺的起始位置
    line([scale_position(1), scale_position(1) + scale_length], ...
         [scale_position(2), scale_position(2)], 'Color', 'k', 'LineWidth', 2);
    text(scale_position(1) + scale_length / 2, scale_position(2) - 0.5, ...
         [num2str(scale_length), ' cm'], 'HorizontalAlignment', 'center', ...
         'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
end
% 创建一个颜色条并设置刻度
colormap(cmap); % 使用自定义的颜色映射
c = colorbar('Position', [0.95, 0.3, 0.02, 0.4]); % 在右侧单独创建一个颜色条
clim([minData maxData]); % 设置颜色条的范围
c.Ticks = linspace(minData, maxData, 6); % 设置刻度为0到500，分成6段
c.TickLabels = {[],[],[],[],[],[]}; % 将刻度转换为字符串
c.Label.String = 'Data Value';
subplots = findobj(gcf, 'Type', 'Axes');

% for i = 1:3
%     subplotPos = get(subplots(i), 'Position');
%     subplotPos(1) = subplotPos(1) - 0.05; % shift
%     subplotPos(3) = subplotPos(3) + 0.05;  % increase width
%     set(subplots(i), 'Position', subplotPos);
% end