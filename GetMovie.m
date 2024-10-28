% set up a new figure and get the data for video;

LVmid = interp1(t,xm_LV,Tnew);
SEPmid = interp1(t,xm_SEP,Tnew);
RVmid = interp1(t,xm_RV,Tnew);
Ymid = interp1(t,ym,Tnew);
h_lv = interp1(t,d_LW,Tnew);
h_rv = interp1(t,d_RW,Tnew);
h_sep = interp1(t,d_SW,Tnew);
R_lv = 1./interp1(t,Cm_LV,Tnew);
R_sep = 1./interp1(t,Cm_SEP,Tnew);
R_rv = 1./interp1(t,Cm_RV,Tnew);
r_inner_LV = -interp1(t,r_LV,Tnew);
r_inner_SEP =  interp1(t,r_SEP,Tnew);
r_inner_RV =  interp1(t,r_RV,Tnew);



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

Sigma_LV = interp1(t,sigma_LV,Tnew);
Sigma_SEP = interp1(t,sigma_SEP,Tnew);
Sigma_RV = interp1(t,sigma_RV,Tnew);
Nsigma_LV = (Sigma_LV - minData) / (maxData - minData); % Normalize
Nsigma_SEP = (Sigma_SEP - minData) / (maxData - minData);
Nsigma_RV = (Sigma_RV - minData) / (maxData - minData);
colorIndex_LV = 1 + floor(Nsigma_LV * (size(cmap, 1) - 1));
colorIndex_SEP = 1 + floor(Nsigma_SEP * (size(cmap, 1) - 1));
colorIndex_RV = 1 + floor(Nsigma_RV * (size(cmap, 1) - 1));
colorIndex_LV = min(colorIndex_LV, 256);
colorIndex_SEP = min(colorIndex_SEP, 256);
colorIndex_RV = min(colorIndex_RV, 256);
colors_LV = cmap(colorIndex_LV, :);
colors_SEP = cmap(colorIndex_SEP, :);
colors_RV = cmap(colorIndex_RV, :);
% set up a new video
if exist("GENDER",'var')
    if GENDER == 1
        v = VideoWriter('TriSegMovingMale.mp4', 'MPEG-4');
    else
        v = VideoWriter('TriSegMovingFemale.mp4', 'MPEG-4');
    end
else
    v = VideoWriter(sprintf('Movies/P_NO%dWindow%d',PatID,ModelWin), 'MPEG-4');
end
v.FrameRate = 40; % per second
open(v);
totalwidth = max(RVmid - LVmid);
range = totalwidth/2*1.5;
% get every frame
for i = 1:10:length(Tnew)
    r_in_LV = r_inner_LV(i);
    r_in_RV = r_inner_RV(i);
    r_out_LV = r_inner_LV(i)-h_lv(i);
    r_out_RV = r_inner_RV(i)+h_rv(i);
    if SEPmid(i) > 0
        r_in_sep = r_inner_SEP(i);
        r_out_sep = r_in_sep+h_sep(i);
    else
        r_in_sep = - r_inner_SEP(i);
        r_out_sep = r_in_sep-h_sep(i);
    end

    x0_lv = LVmid(i) - R_lv(i);
    x0_sep = SEPmid(i) - R_sep(i);
    x0_rv = RVmid(i) - R_rv(i);

    a_lv = atan2(Ymid(i), x0_lv);
    if SEPmid(i) >= 0
        a_sep = atan2(Ymid(i), -x0_sep);
    else
        a_sep = atan2(Ymid(i), x0_sep);
    end
    a_rv = atan2(Ymid(i), -x0_rv);

    a_start_lv = 0*pi/2 - a_lv;
    a_stop_lv = 0*pi/2 + a_lv;

    a_start_sep = 0*pi/2 -a_sep;
    a_stop_sep = 0*pi/2 + a_sep;

    a_start_rv = 0*pi/2 -a_rv;
    a_stop_rv = 0*pi/2 + a_rv;
    figure(88); clf ;
    hold on;
    % Outer LV
    arc_lv = linspace(a_start_lv, a_stop_lv, 20);
    xy_lv = plotArc(x0_lv, r_out_LV, arc_lv);
    plot(xy_lv(1, :), xy_lv(2, :), 'LineWidth',2,'Color',colors_LV(i,:));
    % midline LV
    arc_lv = linspace(a_start_lv, a_stop_lv, 20);
    xy_lv = plotArc(x0_lv, R_lv(i), arc_lv);
    plot(xy_lv(1, :), xy_lv(2, :), '--', 'LineWidth',0.5,'Color',colors_LV(i,:));
    % inner LV
    arc_lv = linspace(a_start_lv, a_stop_lv, 20);
    xy_lv = plotArc(x0_lv,r_in_LV, arc_lv);
    plot(xy_lv(1, :), xy_lv(2, :), 'LineWidth',2,'Color',colors_LV(i,:))


    % mid SEP
    arc_sep = linspace(a_start_sep, a_stop_sep, 20);
    xy_sep = plotArc(x0_sep, R_sep(i), arc_sep);
    plot(xy_sep(1, :), xy_sep(2, :), '--','Color',colors_SEP(i,:))

    if abs(SEPmid(i)) > 0.1
        % LV side
        arc_sep = linspace(a_start_sep, a_stop_sep, 20);
        xy_sep = plotArc(x0_sep, r_in_sep, arc_sep);
        plot(xy_sep(1, :), xy_sep(2, :), 'LineWidth',2,'Color',colors_SEP(i,:))
        % RV side SEP
        arc_sep = linspace(a_start_sep, a_stop_sep, 20);
        xy_sep = plotArc(x0_sep, r_out_sep, arc_sep);
        plot(xy_sep(1, :), xy_sep(2, :), 'LineWidth',2,'Color',colors_SEP(i,:))
    else
        plot(xy_sep(1, :)-h_sep(i)/2, xy_sep(2, :), 'LineWidth',2,'Color',colors_SEP(i,:))
        plot(xy_sep(1, :)+h_sep(i)/2, xy_sep(2, :), 'LineWidth',2,'Color',colors_SEP(i,:))
    end

    % RV inner
    arc_rv = linspace(a_start_rv, a_stop_rv, 20);
    xy_rv = plotArc(x0_rv, r_in_RV, arc_rv);
    plot(xy_rv(1, :), xy_rv(2, :), 'LineWidth',2,'Color',colors_RV(i,:))
    % mid
    arc_rv = linspace(a_start_rv, a_stop_rv, 20);
    xy_rv = plotArc(x0_rv, R_rv(i), arc_rv);
    plot(xy_rv(1, :), xy_rv(2, :), '--','Color',colors_RV(i,:))
    % RV outer
    arc_rv = linspace(a_start_rv, a_stop_rv, 20);
    xy_rv = plotArc(x0_rv, r_out_RV, arc_rv);
    plot(xy_rv(1, :), xy_rv(2, :), 'LineWidth',2,'Color',colors_RV(i,:))
    pbaspect([1,1,1]);

    % depends on the total width too!

    axis([-range range -range range])
    set(gca, "YTick", [-10 -5, 0, 5 10]);
    set(gca, "XTick", [-10 -5, 0, 5 10]);
    set(gca, 'TickLength', [0.040 0.040])
    set(gca, "YTicklabel",cellstr({'-10cm','-5cm', '0cm', '5cm','10cm'}));
    set(gca, "XTicklabel", cellstr({'-10cm','-5cm', '0cm', '5cm','10cm'}));
    hold off;
    % save to video
    writeVideo(v, getframe(gcf));
end
close(v);

