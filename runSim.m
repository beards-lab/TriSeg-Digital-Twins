%% Script Summary:
% This script is primarily designed to solve differential equations using the function dXdT.m with
% ode15s, and to collect the corresponding simulation output.
% A cost function is embedded in the last several sections of the script .
% Created by Andrew Meyer and Feng Gu
% Last modified: 10/29/2024

%% Solve differential equations via ode solver
T = params.T;
HR = params.HR;
init_vec = cell2mat(struct2cell(init))';
M = eye(length(init_vec));
M(1,1) = 0;
M(2,2) = 0;
M(3,3) = 0;
M(4,4) = 0;
options = odeset('Mass', M, 'RelTol', 1e-7, 'AbsTol', 1e-7, 'MaxStep', T/30); % set options for ode
warning('off', 'all'); % turn off warnings message on the command window
maxTime = 10; % maximum time for odeWithTimeout function
lastwarn('', ''); % clear the warning
initialStep = 1; % Initial step size
success = false; % Flag for success status

% 11/06/2024 clean up codes
while true
    options.InitialStep = initialStep;  % Update the initial step size
    [t, y] = odeWithTimeout(@dXdT, [0, 20*T], init_vec, options, params, maxTime); % run the ODE solver
    [lastWarnMsg, lastWarnId] = lastwarn; % Get warning information
    if isempty(lastWarnMsg)
        success = true; % If no warnings, update success status
        break; % Exit the loop
    end
    lastwarn('', ''); % Clear the warning
    initialStep = initialStep / 10^(0.1); % Decrease the initial step size
    if initialStep < eps % Check if the initial step size is less than machine precision
        break; % If the initial step size is too small, exit the loop
    end
end

if ~success
    error('Unable to find a suitable initial step size. All attempts triggered a warning.'); % If the loop finishes without success, throw an error
else
    % Find the index where t is within the last two periods
    startIndex = find(t >= t(end) - 2*T, 1, 'first');
    lastTwoPeriodsT = t(startIndex:end);
    lastTwoPeriodsY = y(startIndex:end, :);
end
t = lastTwoPeriodsT-lastTwoPeriodsT(1);
y = lastTwoPeriodsY;

%% Outputs
output_no = 50;
o = zeros(output_no,length(t));
for i = 1:length(t)
    [~,o(:,i)] = dXdT(t(i),y(i,:), params);
end

xm_LV  = y(:,1);
xm_SEP = y(:,2);
xm_RV  = y(:,3);
ym     = y(:,4);
Lsc_LV  = y(:,5);
Lsc_SEP = y(:,6);
Lsc_RV  = y(:,7);
V_LV = y(:,8);
V_RV = y(:,9);
V_SA = y(:,10);
V_SV = y(:,11);
V_PA = y(:,12);
V_PV = y(:,13);
Vtot = sum(y(end,8:13)) ;
P_LV = o(1,:)';
P_SA = o(2,:)';
P_SV = o(3,:)';
P_RV = o(4,:)';
P_PA = o(5,:)';
P_PV = o(6,:)';
Vm_LV  = o(7,:)';
Vm_SEP = o(8,:)';
Vm_RV  = o(9,:)';
Am_LV  = o(19,:)';
Am_SEP = o(11,:)';
Am_RV  = o(12,:)';
Cm_LV  = o(13,:)';
Cm_SEP = o(14,:)';
Cm_RV  = o(15,:)';
eps_LV  = o(16,:)';
eps_SEP = o(17,:)';
eps_RV  = o(18,:)';
sigma_pas_LV  = o(19,:)';
sigma_pas_SEP = o(20,:)';
sigma_pas_RV  = o(21,:)';
sigma_act_LV  = o(22,:)';
sigma_act_SEP = o(23,:)';
sigma_act_RV  = o(24,:)';
sigma_LV  = o(25,:)';
sigma_SEP = o(26,:)';
sigma_RV  = o(27,:)';
Q_m = o(28,:)' ;    % Flow across mitral valve (QIN_LV)
Q_a = o(29,:)';     % Flow across aortic valve (QOUT_LV)
Q_t = o(30,:)' ;    % Flow across tricuspid valve (QIN_RV)
Q_p = o(31,:)' ;    % Flow across pulmonary valve (QOUT_RV)
Q_SA = o(32,:)' ;
Q_PA = o(33,:)' ;
Tm_LV  = o(34,:)';
Tm_SEP = o(35,:)';
Tm_RV  = o(36,:)';
Y = o(37,:)';
V_RA = o(38,:)';
V_LA = o(39,:)';
P_RA = o(40,:)';
P_LA = o(41,:)';
QIN_RA = o(42,:)';
d_LW = o(43, :)';
d_SW = o(44, :)';
d_RW = o(45, :)';
act = o(46, :)';
r_LV = o(47, :)';
r_SEP = o(48, :)';
r_RV = o(49, :)';
Pericardial = o(50, :)';
Tnew = t(1):0.001:t(end);

%% Figures
end_beat_i = find(t >= 1.02*T, 1) - 1; % Index for end of one complete cardiac cycle, sometime flow shift and a wave is in the middle of the T
% and we need to make this wider

[Qm_maxima, Qm_maxima_i,Qm_wid,Qm_prom] = findpeaks(Q_m(1:end_beat_i));
[Qt_maxima, Qt_maxima_i] = findpeaks(Q_t(1:end_beat_i));


if(length(Qm_maxima) == 2)
    if(Qm_maxima(1) < 0 || Qm_maxima(2) < 0)
        E_A_ratio = -100;
    %elseif(find(Qm_wid .* Qm_prom < 100,1))
     %   E_A_ratio = -100;
    else
        E_A_ratio = Qm_maxima(1) / Qm_maxima(2);
    end
elseif(length(Qm_maxima) > 2)
    while length(Qm_maxima) > 2
        Qm_maxima(Qm_maxima == min(Qm_maxima)) = [];
    end
    assert(length(Qm_maxima) == 2, 'EAr bug 1 runSim');
    E_A_ratio = Qm_maxima(1) / Qm_maxima(2);
else
    E_A_ratio = -100; % FIXME: AF patients doesn't have A wave. How to fix it or just exclude is
    %     error('EAr bug 2 runSim');
end


% FLOW INDICES. FIXME -- mitral and tricuspid flow profiles go to 0
% sometimes during atrial ejection
% Mitral
Qm_sign = sign(Q_m);
if(Qm_sign(1) <= 0)
    Qm_pos_start = find(Qm_sign == 1, 1);
    Qm_sign = Qm_sign(Qm_maxima_i(end):end);
    %     Qm_sign = Qm_sign(Qm_pos_start: end);
    %     Qm_pos_end = find(Qm_sign ~= 1, 1) + Qm_pos_start - 2; % reaches 0 between e and a wave?
    Qm_pos_end = find(Qm_sign ~= 1, 1) + Qm_maxima_i(end) - 2;
    Qm_neg_start = Qm_pos_end + 1;
    %     Qm_sign = Qm_sign(Qm_neg_start - Qm_pos_start + 1: end);
    Qm_sign = Qm_sign(Qm_neg_start - Qm_maxima_i(end) + 1: end);
    Qm_neg_end = find(Qm_sign == 1, 1) + Qm_neg_start - 2;
else
    Qm_neg_start = find(Qm_sign ~= 1, 1);
    Qm_sign = Qm_sign(Qm_neg_start: end);
    Qm_neg_end = find(Qm_sign == 1, 1) + Qm_neg_start - 2;
    Qm_pos_start = Qm_neg_end + 1;
    %     Qm_sign = Qm_sign(Qm_pos_start - Qm_neg_start + 1: end);
    Qm_sign = Qm_sign(Qm_maxima_i(end) - Qm_neg_start + 1:end);
    %     Qm_pos_end = find(Qm_sign ~= 1, 1) + Qm_pos_start - 2;
    Qm_pos_end = find(Qm_sign ~= 1, 1) + Qm_maxima_i(end) - 2;
end
Qm_pos = [Qm_pos_start: Qm_pos_end]; % Indices for positive mitral flow
Qm_neg = [Qm_neg_start: Qm_neg_end]; % Indices for negative mitral flow
% Aortic
Qa_sign = sign(Q_a);
if(Qa_sign(1) <= 0)
    Qa_pos_start = find(Qa_sign == 1, 1);
    Qa_sign = Qa_sign(Qa_pos_start: end);
    Qa_pos_end = find(Qa_sign ~= 1, 1) + Qa_pos_start - 2;
    Qa_neg_start = Qa_pos_end + 1;
    Qa_sign = Qa_sign(Qa_neg_start - Qa_pos_start + 1: end);
    Qa_neg_end = find(Qa_sign == 1, 1) + Qa_neg_start - 2;
else
    Qa_neg_start = find(Qa_sign ~= 1, 1);
    Qa_sign = Qa_sign(Qa_neg_start: end);
    Qa_neg_end = find(Qa_sign == 1, 1) + Qa_neg_start - 2;
    Qa_pos_start = Qa_neg_end + 1;
    Qa_sign = Qa_sign(Qa_pos_start - Qa_neg_start + 1: end);
    Qa_pos_end = find(Qa_sign ~= 1, 1) + Qa_pos_start - 2;
end
Qa_pos = [Qa_pos_start: Qa_pos_end]; % Indices for positive aortic flow
Qa_neg = [Qa_neg_start: Qa_neg_end]; % Indices for negative aortic flow
DNA = (P_SA(Qa_pos_end)+P_LV(Qa_pos_end))/2;% Dicrotic notch for pulmonary
% Tricuspid
Qt_sign = sign(Q_t);
if(Qt_sign(1) <= 0)
    Qt_pos_start = find(Qt_sign == 1, 1);
    %     Qt_sign = Qt_sign(Qt_pos_start: end);
    Qt_sign = Qt_sign(Qt_maxima_i(end): end);
    %     Qt_pos_end = find(Qt_sign ~= 1, 1) + Qt_pos_start - 2;
    Qt_pos_end = find(Qt_sign ~= 1, 1) + Qt_maxima_i(end) - 2;
    Qt_neg_start = Qt_pos_end + 1;
    %     Qt_sign = Qt_sign(Qt_neg_start - Qt_pos_start + 1: end);
    Qt_sign = Qt_sign(Qt_neg_start - Qt_maxima_i(end) + 1: end);
    Qt_neg_end = find(Qt_sign == 1, 1) + Qt_neg_start - 2;
else
    Qt_neg_start = find(Qt_sign ~= 1, 1);
    Qt_sign = Qt_sign(Qt_neg_start: end);
    Qt_neg_end = find(Qt_sign == 1, 1) + Qt_neg_start - 2;
    Qt_pos_start = Qt_neg_end + 1;
    %     Qt_sign = Qt_sign(Qt_pos_start - Qt_neg_start + 1: end);
    Qt_sign = Qt_sign(Qt_maxima_i(end) - Qt_neg_start + 1:end);
    %     Qt_pos_end = find(Qt_sign ~= 1, 1) + Qt_pos_start - 2;
    Qt_pos_end = find(Qt_sign ~= 1, 1) + Qt_maxima_i(end) - 2;

end
Qt_pos = [Qt_pos_start: Qt_pos_end]; % Indices for positive tricuspid flow
Qt_neg = [Qt_neg_start: Qt_neg_end]; % Indices for negative tricuspid flow
% Pulmonary
Qp_sign = sign(Q_p);
if(Qp_sign(1) <= 0)
    Qp_pos_start = find(Qp_sign == 1, 1);
    Qp_sign = Qp_sign(Qp_pos_start: end);
    Qp_pos_end = find(Qp_sign ~= 1, 1) + Qp_pos_start - 2;
    Qp_neg_start = Qp_pos_end + 1;
    Qp_sign = Qp_sign(Qp_neg_start - Qp_pos_start + 1: end);
    Qp_neg_end = find(Qp_sign == 1, 1) + Qp_neg_start - 2;
else
    Qp_neg_start = find(Qp_sign ~= 1, 1);
    Qp_sign = Qp_sign(Qp_neg_start: end);
    Qp_neg_end = find(Qp_sign == 1, 1) + Qp_neg_start - 2;
    Qp_pos_start = Qp_neg_end + 1;
    Qp_sign = Qp_sign(Qp_pos_start - Qp_neg_start + 1: end);
    Qp_pos_end = find(Qp_sign ~= 1, 1) + Qp_pos_start - 2;
end
Qp_pos = [Qp_pos_start: Qp_pos_end]; % Indices for positive pulmonary flow
Qp_neg = [Qp_neg_start: Qp_neg_end]; % Indices for negative pulmonary flow
DNP = (P_PA(Qp_pos_end)+P_RV(Qp_pos_end))/2;% Dicrotic notch for pulmonary
% FIXME: breaks when there is only positive flow (for some reason this
% happened when  I set the mitral resistance 100x normal)


% STENOSIS GRADING
% MS
MPG_m = trapz(t(Qm_pos), P_LA(Qm_pos) - P_LV(Qm_pos)) / (t(Qm_pos_end) - t(Qm_pos_start));

% AS
MPG_a = trapz(t(Qa_pos), P_LV(Qa_pos) - P_SA(Qa_pos)) / (t(Qa_pos_end) - t(Qa_pos_start));
Peak_AG_p = max(P_LV(Qp_pos) - P_SA(Qp_pos));

% TS
MPG_t = trapz(t(Qt_pos), P_RA(Qt_pos) - P_RV(Qt_pos)) / (t(Qt_pos_end) - t(Qt_pos_start));

% PS
MPG_p = trapz(t(Qp_pos), P_RV(Qp_pos) - P_PA(Qp_pos)) / (t(Qp_pos_end) - t(Qp_pos_start));
Peak_PG_p = max(P_RV(Qp_pos) - P_PA(Qp_pos));


% REGURGITATION GRADING
% MR
RVol_m = -trapz(t(Qm_neg), Q_m(Qm_neg)); % Regurgitant volume over period where Qm is negative
SV_LA_pos = trapz(t(Qm_pos), Q_m(Qm_pos));
SV_LA = max(V_LA) - min(V_LA); % FIXME: this was just because this wasa target in MILD_MR patient
% FIXME: Not sure how to calculate atrial stroke volume, so I'm just
% dividing by forward stroke volume
RF_m = 100 * (RVol_m / SV_LA_pos);

% AI
RVol_a = -trapz(t(Qa_neg), Q_a(Qa_neg));
SV_LV_pos = trapz(t(Qa_pos), Q_a(Qa_pos));
RF_a = 100 * (RVol_a / SV_LV_pos);

% TR
RVol_t = -trapz(t(Qt_neg), Q_t(Qt_neg));
SV_RA_pos = trapz(t(Qt_pos), Q_t(Qt_pos));
RF_t = 100 * (RVol_t / SV_RA_pos);

% PI
RVol_p = -trapz(t(Qp_neg), Q_p(Qp_neg));
SV_RV_pos = trapz(t(Qp_pos), Q_p(Qp_pos));
RF_p = 100 * (RVol_p / SV_RV_pos);





SV_LV_tot = max(V_LV) - min(V_LV); % For some reason, this is slightly different from SV_LV_pos + RVol_m
SV_RV_tot = max(V_RV) - min(V_RV);

EF_LV = SV_LV_tot / max(V_LV); % FIXME: is this based on forward SV or total SV?
EF_RV = SV_RV_tot / max(V_RV);


CO = (SV_LV_pos-RVol_a) * HR / 1000; % FIXME: Using CO with forward flow through LV (????)
CO_RV = (SV_RV_pos-RVol_p) * HR / 1000; % FIXME: Using CO with forward flow through LV (????)

MPAP = (1/3) * (max(P_PA)) + (2/3) * (min(P_PA));

[~, LVED_i] = max(V_LV); % index of end diastole. FIXME: do we need different end-diastole indices for left and right sides?
[~, RVED_i] = max(V_RV);
[~, LVES_i] = min(V_LV);
[~, RVES_i] = min(V_RV);
Y_max_i = find(abs(Y - 1) < 0.001, 1); % max of activation function doesn't quite coincide with max pressure development
RVEDP_i = find(Y > 0, 1); % Should be the first or second index
ESP_i = find(Y > 0, 1); % Should be the first or second index


LVIDs = -xm_LV(LVES_i) + xm_SEP(LVES_i)-...
    (1/Cm_SEP(LVES_i)- r_SEP(LVES_i)...% SEP miwall- inner radias
    +(-1/Cm_LV(LVES_i))-r_LV(LVES_i)); % + LV midwall- inner radias

RVIDd = r_RV(RVED_i)...% RV inner r
    + ((xm_RV(RVED_i)-1/Cm_RV(RVED_i))-xm_SEP(RVED_i))... % center of the RV - orginal - xm_sep
    - (d_SW(RVED_i)-(1/Cm_SEP(RVED_i)-r_SEP(RVED_i))); % thickness-(midwall SEP-innerwall)
LVIDd = -xm_LV(LVED_i) + xm_SEP(LVED_i)-...
    (1/Cm_SEP(LVED_i)- r_SEP(LVED_i)...% SEP miwall- inner radias
    +(-1/Cm_LV(LVED_i))- r_LV(LVED_i)); % + LV midwall- inner radias
if LVIDd >20 || LVIDd <0
    LVIDd = r_LV(LVED_i)...% LV inner r
        + ((-xm_LV(LVED_i)+1/Cm_LV(LVED_i))+xm_SEP(LVED_i))... % center of the LV - orginal - xm_sep
        - (d_SW(LVED_i)-(-1/Cm_SEP(LVED_i)-r_SEP(LVED_i))); % thickness-(midwall SEP-innerwall)
end
if RVIDd>100 || RVIDd <0
    RVIDd = xm_RV(RVED_i) - xm_SEP(RVED_i)-...% RV inner r
        (-1/Cm_SEP(RVED_i)- r_SEP(RVED_i)...% SEP miwall- inner radias
        +(1/Cm_RV(RVED_i))- r_RV(RVED_i)); % + LV midwall- inner radias
end
RVIDs = r_RV(RVES_i)...% RV inner diameters
    + ((xm_RV(RVES_i)-1/Cm_RV(RVES_i))-xm_SEP(RVES_i))... % center of the RV - orginal - xm_sep
    - (d_SW(RVES_i)-(1/Cm_SEP(RVES_i)-r_SEP(RVES_i))); % thickness-(midwall SEP-innerwall)
if RVIDs>100 || RVIDs <0
    RVIDs = xm_RV(RVES_i) - xm_SEP(RVES_i)-...% RV inner r
        (-1/Cm_SEP(RVES_i)- r_SEP(RVES_i)...% SEP miwall- inner radias
        +(1/Cm_RV(RVES_i))- r_RV(RVES_i)); % + LV midwall- inner radias
end
% Outputting masses
% thickness ED, volume ED, x
rho_myo = 1.0550;
Hed_LW = d_LW(LVED_i);
Hed_SW = d_SW(LVED_i);
Hed_RW = d_RW(RVED_i); % <- change to RVED_i?
LV_m = rho_myo * (params.Vw_LV+params.Vw_SEP);
RV_m = rho_myo * (params.Vw_RV);


% Pressure waveforms RHC a and v wave (not very well-implemented)

% hack version
act_max_atria_i = find(abs(act - 1) < 0.001, 1);
% [RAP_maxima, RAP_maxima_i] = findpeaks(P_RA(act_max_atria_i:act_max_atria_i + end_beat_i + 1), ...
%             'SortStr','descend','NPeaks',4); % should get two biggest peaks (2 beats)
P_RA_smooth = smoothdata(P_RA,"gaussian",100);
[RAP_maxima, RAP_maxima_i] = findpeaks(P_RA_smooth, ...
    'SortStr','descend','NPeaks',4); % should get two biggest peaks (2 beats)
%add max atria act i
% act max atria i + end beat i + 1 is closest to
if(length(RAP_maxima_i) == 4)
    [~,RAP_a_i] = min(abs(act_max_atria_i + end_beat_i + 1 - RAP_maxima_i));
    % RAP_a_i = find()

    % RAP_v_i = find(abs(RAP_maxima - RAP_maxima(RAP_a_i)) > 0.1, 1);
    [~,RAP_v_i] = max(abs(RAP_maxima - RAP_maxima(RAP_a_i)));

    if(isempty(RAP_v_i))
        RAP_v_i = RAP_a_i;
    end
    % if(abs(RAP_maxima_i(RAP_v_i) - RAP_maxima_i(RAP_a_i)) < 100 || abs(abs(RAP_maxima_i(RAP_v_i) - RAP_maxima_i(RAP_a_i)) - end_beat_i) < 100)
    %     RAP_a = -100;
    %     RAP_v = -100;
    % else
    RAP_a = RAP_maxima(RAP_a_i);
    RAP_v = RAP_maxima(RAP_v_i); %idk if this will work
    % end
elseif(length(RAP_maxima_i) == 3) % bad implementation
    RAP_a = RAP_maxima(1);
    RAP_v = RAP_maxima(2);
else
    RAP_a = P_RA(act_max_atria_i);
    RAP_v = P_RA(Y_max_i);
end

% PCWP waveform
[P_PV_maxima, P_PV_maxima_i] = findpeaks(P_PV, ...
    'SortStr','descend','NPeaks',4); % should get two biggest peaks (2 beats)
%add max atria act i
% act max atria i + end beat i + 1 is closest to
if(length(P_PV_maxima_i) == 4)
    [~, PCWP_a_i] = min(abs(act_max_atria_i + end_beat_i + 1 - P_PV_maxima_i));
    % RAP_a_i = find()

    % PCWP_v_i = find(abs(P_PV_maxima - P_PV_maxima(PCWP_a_i)) > 0.1, 1);
    [~,PCWP_v_i] = max(abs(P_PV_maxima - P_PV_maxima(PCWP_a_i)));
    if(isempty(PCWP_v_i))
        PCWP_v_i = PCWP_a_i;
    end
    % if(abs(P_PV_maxima_i(PCWP_v_i) - P_PV_maxima_i(PCWP_a_i)) < 100 || abs(abs(P_PV_maxima_i(PCWP_v_i) - P_PV_maxima_i(PCWP_a_i)) - end_beat_i) < 100)
    %     PCWP_a = -100;
    %     PCWP_v = -100;
    % else
    PCWP_a = P_PV_maxima(PCWP_a_i);
    PCWP_v = P_PV_maxima(PCWP_v_i); %idk if this will work
    % end

    % PCWP_a = P_PV_maxima(PCWP_a_i);
    % PCWP_v = P_PV_maxima(PCWP_v_i); %idk if this will work
else
    PCWP_a = P_PV(act_max_atria_i);
    PCWP_v = max(P_PV_maxima);
end




%% Calculating cost function metrics
o_vals = struct(); % Relevant model outputs, to be compared to target values

o_vals.SBP = max(P_SA);
o_vals.DBP = min(P_SA);
o_vals.LVEDV = max(V_LV);
o_vals.EF = EF_LV * 100; % Based on total stroke volume (with regurgitation) on the left side
o_vals.LVESV = min(V_LV);
o_vals.EAr = E_A_ratio;
o_vals.LAVmax = max(V_LA);
o_vals.LAVmin = min(V_LA);
o_vals.SV_LA = SV_LA;
o_vals.RVEDV = max(V_RV);
o_vals.RVESV = min(V_RV);
o_vals.RVEF = EF_RV * 100;
o_vals.RAVmax = max(V_RA);
o_vals.RAVmin = min(V_RA);
o_vals.RAPmax = max(P_RA);
o_vals.RAPmin = min(P_RA);
o_vals.RAPmean = trapz(t,P_RA)/(t(end)-t(1)); % 
o_vals.PASP = max(P_PA); % FIXME: should this be at the same index as RVESP is taken?
o_vals.PADP = min(P_PA);
o_vals.PCWP = trapz(t,P_PV)/(t(end)-t(1));
o_vals.PCWPmax = max(P_PV);
o_vals.PCWPmin = min(P_PV);
o_vals.CVP = trapz(t,P_SV)/(t(end)-t(1));
o_vals.CO = CO_RV;% RHC RV
o_vals.Hed_LW = d_LW(LVED_i);
o_vals.Hed_SW = d_SW(LVED_i);
o_vals.Hed_RW = d_RW(RVED_i); % <- change to RVED_i?
o_vals.RVEDP = P_RV(RVEDP_i); % Activation function beginning coincides with start of pressure development on right side. Beginning of activation function is in the first few indices.
o_vals.P_RV_min = min(P_RV);
o_vals.LVEDP = P_LV(RVEDP_i);
o_vals.P_LV_min = min(P_LV);
o_vals.RVSP = max(P_RV);
o_vals.LVESP = max(P_LV);
o_vals.MVmg = MPG_m;
o_vals.AVpg = Peak_AG_p;
o_vals.TVmg = MPG_t;
o_vals.PVpg = Peak_PG_p;
o_vals.RF_m = RF_m;
o_vals.RF_a = RF_a;
o_vals.RF_t = RF_t;
o_vals.RF_p = RF_p;

[~,a]= max(P_LV);
if t(a) < T
    o_vals.tPLVmax = t(a);
else
    o_vals.tPLVmax = t(a)-T;
end
[~,b]= min(P_LV);
if t(b) < T
    o_vals.tPLVmin = t(b);
else
    o_vals.tPLVmin = t(b)-T;
end
% Mild MR patient
o_vals.RAPbar = mean(P_RA);
o_vals.RAP_a = RAP_a;
o_vals.RAP_v = RAP_v;
o_vals.PCWP_a = PCWP_a;
o_vals.PCWP_v = PCWP_v;

% For inputs loop -
o_vals.LVIDd = LVIDd;
o_vals.LVIDs = LVIDs;
o_vals.RVIDd = RVIDd;
o_vals.RVIDs = RVIDs;
o_vals.RV_m = RV_m;
o_vals.LV_m = LV_m;
o_vals.FakeLV_m = LV_m;
o_vals.Vtot = Vtot;

% To aviod overfit the compliance
o_vals.DNA = DNA;
o_vals.DNP = DNP;

%% grading the valvular stuff and tax
Lsc_target = 2; % zero cost
% add costs from three EDSL contriubtions
Lsc_ED = [max(Lsc_LV), max(Lsc_SEP), max(Lsc_RV)];
cost_Lsc = 0;
for i = 1:length(Lsc_ED)
    cost_Lsc = cost_Lsc + 600*(Lsc_target - Lsc_ED(i))^4;
end

cost_EA_fn = @(EA)  25*(0.167066669245687*EA.^4 -1.81349025058486*EA.^3 + 7.21416047093197*EA.^2 -11.6408665836098*EA + 5.50562700516957); % to be included in tax
cost_EA = cost_EA_fn(o_vals.EAr);
if(cost_EA < 0)
    cost_EA = 0;
end
cost_Mass_fn = @(M)  25*(1.81128506110855e-08*M.^4 -1.11713451327939e-05*M.^3+0.00267922216978568*M.^2 -0.297667946824409*M + 12.9346600216876); % to be included in tax
cost_LV_m = cost_Mass_fn(o_vals.LV_m);
if(cost_LV_m < 0)
    cost_LV_m = 0;
end
cost_RV_m = cost_Mass_fn(o_vals.RV_m);
if(cost_RV_m < 0)
    cost_RV_m = 0;
end


g = [1, 2, 3, 4]; % grades: none, mild, moderate, severe.
g_t = [1, 2]; % Tricuspid grades. There is no grading convention for tricuspid stenosis, so this is supposed to mean 1 for no stenosis and 2 for stenosis.

% Mitral
gmt_MR = [0, 20, 40, 60];
% Aortic
gmt_AI = [0, 20, 40, 60];

% Tricuspid
gmt_TR = [0, 17.5, 35, 52.5];

% Pulmonary
gmt_PI = [0, 15, 30 ,45];




o_vals.MVr = interp1(gmt_MR, g, RF_m, 'linear', 'extrap');
o_vals.AVr = interp1(gmt_AI, g, RF_a, 'linear', 'extrap');
o_vals.TVr = interp1(gmt_TR, g, RF_t, 'linear', 'extrap');
o_vals.PVr = interp1(gmt_PI, g, RF_p, 'linear', 'extrap');


%% calibration based on healthy normal subject
targetsfn = fieldnames(targets);
inputsfn = fieldnames(inputs);
N = length(targetsfn);
c = struct();

for i  = 1:N
    c.(targetsfn{i}) = targets.(targetsfn{i});
end

if inputs.Sex == 1
    c.LV_m = 121;
    c.FakeLV_m = 121;
    c.RV_m = 66;
    c.SBP = 120;
    c.DBP = 80;
    c.LVEDV = 155;
    c.LVESV = 60;
    c.DNA = 100;
    c.DNP = 15;
    c.RVSP = 28.5;
    c.PASP = 22.5;
    c.PADP = 11.5;
    c.LVIDd = 6.5;
    c.LVIDs = 5.5;
    c.RVIDd = 4.5;
    c.RVIDs = 3.5;
    c.RVEDV = 166;
    c.RVESV = 73;
    c.Hed_LW = 0.93;
    c.Hed_SW = 0.92;
    c.Hed_RW = 0.35;
    c.AVpg = 10;
    c.TVmg = 5;
    c.MVmg = 5;
    c.PVpg = 10;
    c.EF = 60;
    c.CO = 5.7;
    c.RAPmax = 6;
    c.RAPmean = 4;
    c.PCWPmax = 13;
    c.PCWP = 9;
    c.PCWPmin = 5;
    c.RVEDP = 6;
    c.P_RV_min = 3;
    c.LVEDP = 6;
    c.P_LV_min = 3;
    c.EAr = 1.69;
    c.LAVmax = 72;
    c.tPLVmax = 0.35;
    c.tPLVmin = 0.55;
    c.RVEF = 60;
elseif inputs.Sex == 2
    c.LV_m = 83;
    c.FakeLV_m = 83;
    c.RV_m = 48;
    c.SBP = 116.5;
    c.DBP = 72.9;
    c.LVEDV = 123;
    c.LVESV = 43;
    c.DNA = 92.8;
    c.DNP = 15;
    c.RVSP = 22.5;
    c.PASP = 22.5;
    c.PADP = 11.5;
    c.LVIDd = 6.5;
    c.LVIDs = 5.5;
    c.RVIDd = 4.5;
    c.RVIDs = 3.5;
    c.RVEDV = 122;
    c.RVESV = 50;
    c.Hed_LW = 0.85;
    c.Hed_SW = 0.82;
    c.Hed_RW = 0.31;
    c.AVpg = 10;
    c.TVmg = 5;
    c.MVmg = 5;
    c.PVpg = 10;
    c.EF = 60;
    c.CO = 4.8;
    c.RAPmax = 6;
    c.RAPmean = 4;
    c.PCWPmax = 13;
    c.PCWPmin = 5;
    c.PCWP = 9;
    c.RVEDP = 6;
    c.P_RV_min = 3;
    c.LVEDP = 5.5;
    c.P_LV_min = 2.5;
    c.EAr = 1.72;
    c.LAVmax = 50;
    c.tPLVmax = 0.35;
    c.tPLVmin = 0.55;
    c.RVEF = 60;
end
%% weights
% cost function are normalized based on the percentage change. Will give
% the varibles with high base level higher weight.
w = struct();
wf1 = 55; % weights for 1st sub figure
w.SBP = wf1*10; w.DBP = wf1*10;
w.DNA = wf1*10;
w.PASP = wf1; w.PADP = wf1;
w.DNP = wf1/2;
w.LVESP = wf1;
if(isfield(targets,'PASP'))
    w.RVSP = wf1;
else
    w.RVSP = wf1*2; % RVSP and PASP
end
wf2 = 20; % weights for 2st sub figure
w.RAPmax = wf2; w.RAPmin = wf2/5;
if ~isfield(targets,'RAPmax') && ~isfield(targets,'RAPmin')
    w.RAPmean = wf2*2*0.1; % only info we know about RA
else
    w.RAPmean = wf2*0.1;
end
w.RVEDP = wf2*0.2;
w.LVEDP = wf2*0.2;
w.P_RV_min = wf2*0.1;
w.P_LV_min = wf2*0.1;
w.PCWP = wf2*0.3;
w.PCWPmax = wf2*0.3;
w.PCWPmin = wf2*0.3;
wf3 = 250;
w.LVEDV = wf3; w.LVESV = wf3*0.3;
w.LAVmax = wf3*0.1; w.LAVmin = wf3*0.1;
w.RAVmax = wf3*0.1; w.RAVmin = wf3*0.1;
w.RVEDV = wf3; w.RVESV = wf3*0.3;

% targets not plotted
w.EF = 6;
w.RVEF = 6;
w.CO = 88; % give CO a lucky number
if(isfield(targets,'EAr'))
    w.EAr = 88;
end
wt = 1.5; % use to adjust thickness
w.Hed_LW = wt*25*0.3;
w.Hed_SW = wt*25*0.3;
w.Hed_RW = wt*25*0.3;
w.LVIDd = wt*3;
w.LVIDs = wt*3;
w.RVIDd = wt*3;
w.RVIDs = wt*3;
wg = 4; % Weight of valve insuffiency grades
vlv_def = ["MVr",'AVr'];
for i = 1:length(vlv_def)
    if isfield(targets,vlv_def(i))
        if (targets.(vlv_def(i)) > 3.5) && (o_vals.(vlv_def(i)) > targets.(vlv_def(i)))
            w.(vlv_def(i)) = wg;
        else
            if targets.(vlv_def(i)) == 1.5
                if  o_vals.(vlv_def(i)) < targets.(vlv_def(i))
                    w.(vlv_def(i)) = wg;
                else
                    w.(vlv_def(i)) = 25 * wg;
                end
            else
                if abs(targets.(vlv_def(i))-o_vals.(vlv_def(i))) > .5
                    w.(vlv_def(i)) = 25 * wg;
                else
                    w.(vlv_def(i)) = wg;
                end
            end
        end
    end
end
if isfield(targets,"TVr")
    if (targets.TVr > 3.5) && (o_vals.TVr > targets.TVr)
        w.TVr = wg;
    else
        if targets.TVr == 4
            if abs(o_vals.TVr - targets.TVr) < 0.7143
                w.TVr = wg;
            else
                w.TVr = 25*wg;
            end
        end
        if targets.TVr == 3.5
            if abs(o_vals.TVr - targets.TVr) < 0.5
                w.TVr = wg;
            else
                w.TVr = 25*wg;
            end
        end
        if targets.TVr == 3
            if abs(o_vals.TVr - targets.TVr) < 0.2857
                w.TVr = wg;
            else
                w.TVr = 25*wg;
            end
        end
        if targets.TVr == 2.5
            if abs(o_vals.TVr - targets.TVr) < 0.3571
                w.TVr = wg;
            else
                w.TVr = 25*wg;
            end
        end
        if targets.TVr == 2
            if abs(o_vals.TVr - targets.TVr) <= 0.7143
                w.TVr = wg;
            else
                w.TVr = 25*wg;
            end
        end
        if targets.TVr == 1.5
            if o_vals.TVr < targets.TVr
                w.TVr = wg;
            else
                w.TVr = 25*wg;
            end
        end
        if targets.TVr == 1
            w.TVr = wg/25;
        end
    end
end

if isfield(targets,"PVr")
    if (targets.PVr > 3.5) && (o_vals.PVr > targets.PVr)
        w.PVr = wg;
    else
        if targets.PVr == 4
            if abs(o_vals.PVr - targets.PVr) < 0.3333
                w.PVr = wg;
            else
                w.PVr = 25*wg;
            end
        end
        if targets.PVr == 3.5
            if abs(o_vals.PVr - targets.PVr) < .5
                w.PVr = wg;
            else
                w.PVr = 25*wg;
            end
        end
        if targets.PVr == 3
            if abs(o_vals.PVr - targets.PVr) <= 0.6667
                w.PVr = wg;
            else
                w.PVr = 25*wg;
            end
        end
        if targets.PVr == 2.5
            if abs(o_vals.PVr - targets.PVr) < .5
                w.PVr = wg;
            else
                w.PVr = 25*wg;
            end
        end
        if targets.PVr == 2
            if abs(o_vals.PVr - targets.PVr) <= 0.3333
                w.PVr = wg;
            else
                w.PVr = 25*wg;
            end
        end
        if targets.PVr == 1.5
            if o_vals.PVr < targets.PVr
                w.PVr = wg;
            else
                w.PVr = 25*wg;
            end
        end
        if targets.PVr == 1
            w.PVr = wg/25;
        end
    end
end

w.MVmg = wg/5;
w.AVpg = wg/25;
w.TVmg = wg/5;
w.PVpg = wg/25;

%LV_m
if(isfield(targets,'LV_m'))
    w.LV_m = 100;
end
w.FakeLV_m = 50;
%RV_m
if(isfield(targets,'RV_m'))
    w.RV_m = 100;
end
w.tPLVmax = 150;
w.tPLVmin = 150;

%% cost function
cost = zeros(1, N);


if ~exist('print_sim','var')
    % lazy defensive programming
    print_sim = false;
end
EX = 7.2812;
for i = 1:N
    % normalized squared difference
    % normalizing e^2/target^2. BEST.
    cost(i) = (o_vals.(targetsfn{i}) - targets.(targetsfn{i}))^2 /c.(targetsfn{i})^2 *w.(targetsfn{i})*EX; % 11/7/2023 Dollars to Yuan exchange rate
    if print_sim
        fprintf('%i) %s: %1.2f (%1.2f)\t %1.0f%% (%1.2f¥)\n', i, targetsfn{i}, o_vals.(targetsfn{i}), targets.(targetsfn{i}), o_vals.(targetsfn{i})/targets.(targetsfn{i})*100 - 100, cost(i));
    end

end
if(exist('cost_EA','var'))
    tax = cost_Lsc + cost_EA; % add any extra costs that are independent of targets (i.e. constraints or something)
else
    tax = cost_Lsc;
end
total_cost = sum(cost) + tax;



%% kill out of bound params
% this part is used for kill out of bound params during GA opt
% Skip this part when you run standard simulation
load AllPatients.mat
if ~exist('GENDER', 'var')
    height = patients(PatID).snapshots(ModelWin).Height;
    weight = patients(PatID).snapshots(ModelWin).Weight;
    BSA = sqrt((height * weight) / 3600);
    NParams.C_SA = params.C_SA/BSA;
    NParams.C_SV = params.C_SV/BSA;
    NParams.C_PA = params.C_PA/BSA;
    NParams.C_PV = params.C_PV/BSA;
    NParams.G_SA = 1/(params.R_SA*BSA);
    NParams.G_PA = 1/(params.R_PA*BSA);
    NParams.G_tSA = 1/(params.R_tSA*BSA);
    NParams.G_tPA = 1/(params.R_tPA*BSA);
    NParams.k_act_LV = params.k_act_LV;
    NParams.k_act_RV = params.k_act_RV;
    NParams.k_pas_LV = params.k_pas_LV;
    NParams.k_pas_RV = params.k_pas_RV;
    NParams.Vw_LV = params.Vw_LV/BSA;
    NParams.Vw_RV = params.Vw_RV/BSA;
    NParams.Vw_SEP = params.Vw_SEP/BSA;
    NParams.RAV0u = params.RAV0u/BSA; NParams.LAV0u = params.RAV0u/BSA;
    NParams.RAV0c = params.RAV0c/BSA; NParams.LAV0c = params.RAV0c/BSA;
    NParams.G_SV = 1/(params.R_SV*BSA);
    NParams.G_PV = 1/(params.R_PV*BSA);
    NParams.Vh0 = params.Vh0/BSA;NParams.K_P = params.K_P;NParams.B_P = params.B_P;
    NParams.G_t_o = 1/(params.R_t_o*BSA);
    NParams.G_t_c = 1/(params.R_t_c/BSA);
    NParams.G_p_o = 1/(params.R_p_o/BSA);
    NParams.G_p_c = 1/(params.R_p_c/BSA);
    NParams.G_m_o = 1/(params.R_m_o/BSA);
    NParams.G_m_c = 1/(params.R_m_c/BSA);
    NParams.G_a_o = 1/(params.R_a_o/BSA);
    NParams.G_a_c = 1/(params.R_a_c/BSA);
    NParams.V_SV_s = init.V_SV_s/BSA; NParams.V_SA_s = init.V_SA_s/BSA;
    NParams.V_PV_s = init.V_PV_s/BSA; NParams.V_PA_s = init.V_PA_s/BSA;
    paramsname = fieldnames(NParams);
    for j = 1:length(paramsname)
        p = NParams.(paramsname{j});
        load("standardnormalization.mat");
        if inputs.Sex == 1
            p = p/standardparamstable.(paramsname{j})(1);
        elseif inputs.Sex ==2
            p = p/standardparamstable.(paramsname{j})(2);
        end
        RNParams.(paramsname{j}) = p;
    end

    if ~isfield(targets,'LVEDP')
        if  RNParams.k_pas_LV < 0.5
            total_cost = 50000;
        end
    end
    if ~isfield(targets,'RVEDP')
        if  RNParams.k_pas_RV < 0.5
            total_cost = 50000;
        end
    end
    
    output.mods = mods;
    output.modifiers = modifiers;
    output.params = params;
    output.init = init;
    output.targets = targets;
    output.inputs = inputs;
    save(sprintf('Sims/P_NO%dWindow%d.mat',PatID,ModelWin),"output");

    FID = fopen(sprintf('Sims/P_NO%dWindow%d.txt',PatID,ModelWin), 'w');


    targetsfn = fieldnames(targets);
    inputsfn = fieldnames(inputs);
    cost = zeros(1, length(targetsfn));
    EX = 7.2812;
    for i = 1:length(targetsfn)
        cost(i) = (o_vals.(targetsfn{i}) - targets.(targetsfn{i}))^2 /c.(targetsfn{i})^2 *w.(targetsfn{i})*EX; % 11/7/2023 Dollars to Yuan exchange rate
        fprintf(FID,'%i) %s: %1.2f (%1.2f)\t %1.0f%% ($%1.2f)\n', i, targetsfn{i}, o_vals.(targetsfn{i}), targets.(targetsfn{i}), o_vals.(targetsfn{i})/(targets.(targetsfn{i})+0.00000001)*100 - 100, cost(i));
    end
    if(exist('cost_EA','var'))
        tax = cost_Lsc + cost_EA; % add any extra costs that are independent of targets (i.e. constraints or something)
        fprintf(FID,strcat("\nTax Lsc: $",num2str(cost_Lsc),"\nTax E/A Ratio: $",num2str(cost_EA),"\n"));
    else
        tax = cost_Lsc;
        fprintf(FID,strcat("\nTax Lsc: $",num2str(cost_Lsc),"\n"));
    end
    fprintf(FID,'Total cost: $%1.2f \n\n', total_cost);

    fclose(FID);
end

if print_sim
    fprintf('Total cost: %1.2f \n\n', total_cost);
end








%% function to kill inf loop of ode15s

function [T,Y,TE,YE,IE] = odeWithTimeout(odefun, tspan, y0, options, params, maxTime)
ticID = tic; % record time

evalc('options.Events = @(t, y) eventFunction(t, y, maxTime, ticID);'); % set up event

% ode
[T, Y, TE, YE, IE] = ode15s(@(t,y) odefun(t,y,params), tspan, y0, options); % attention passing params is dangerous
end

% define time out
function [value, isterminal, direction] = eventFunction(t, y, maxTime, ticID)
elapsed = toc(ticID); % counting current time
value = elapsed < maxTime; %
isterminal = 1; % stop
direction = 0; % find all direction 0
end


