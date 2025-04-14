%% Script Summary:
% This script is primarily designed to solve differential equations using the function dXdT.m with
% ode15s, and to collect the corresponding simulation output.
% A cost function is embedded in the last several sections of the script

% 03/20: The biggest difference between this version and the previous one
% is that the current version incorporates several empirical correlations
% identified from the UW cohort to constrain RV mechanical and geometrical
% properties, helping to eliminate unrealistic conditions during optimization.

% Created by Andrew Meyer and Feng Gu
% Last modified: 03/20/2024

%% Solve the differential equations using the ODE solver
T = params.T;
HR = params.HR;
if init.LAVmin >= 150||init.RAVmin>=100
    init.LAVmin = init.LAVmin/3;
    init.RAVmin = init.RAVmin/3;
end
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
initialStep = 1; % initial step size
success = false; % flag for success status

% This loop is designed primarily for optimization. The differential equations are very stiff,
% especially for xm_sep, which may approach 0, causing the numerical solver to crash or get stuck.
% This loop can penalize this condition, although it is not perfect. Ideally, the solution should
% not depend on the initial step, but the equations are too stiff for the numerical solver to handle.
while true
    options.InitialStep = initialStep;  % update the initial step size
    [t, y] = odeWithTimeout(@dXdT, [0, 20*T], init_vec, options, params, maxTime); % run the ODE solver
    [lastWarnMsg, lastWarnId] = lastwarn; % get warning information
    if isempty(lastWarnMsg)
        success = true; % if no warnings, update success status
        break; % exit the loop
    end
    lastwarn('', ''); % clear the warning
    initialStep = initialStep / 10^(0.1); % decrease the initial step size
    if initialStep < eps % check if the initial step size is less than machine precision
        break; % if the initial step size is too small, exit the loop
    end
end

if ~success
    error('Unable to find a suitable initial step size. All attempts triggered a warning.'); % if the loop finishes without success, throw an error
else
    % Find the index where t is within the last two periods, which reflects the steady state
    startIndex = find(t >= t(end) - 2*T, 1, 'first');
    lastTwoPeriodsT = t(startIndex:end);
    lastTwoPeriodsY = y(startIndex:end, :);
end
t = lastTwoPeriodsT-lastTwoPeriodsT(1);
y = lastTwoPeriodsY; % solutions of ODE

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

% Clever way to aviod unreal uncondition
Tunit = (0:0.001:t(end));
interpLV = griddedInterpolant(t,V_LV,'pchip');
V_LVnew = interpLV(Tunit)';
new_min = 90;
new_max = 150;
currentLV_min = min(V_LVnew);
currentLV_max = max(V_LVnew);
V_LVnew_scaled = (V_LVnew - currentLV_min) * (new_max - new_min) / (currentLV_max - currentLV_min) + new_min;
dVLV = diff(V_LVnew_scaled);
LVpotentialPeriod = Tunit([false; abs(dVLV) < 0.02]);
[~,LVjump] = maxk(diff(LVpotentialPeriod),3);
LVjump = sort(LVjump);
LVperiod = [LVpotentialPeriod(1) LVpotentialPeriod(LVjump(1)); LVpotentialPeriod(LVjump(1)+1) LVpotentialPeriod(LVjump(2));...
    LVpotentialPeriod(LVjump(2)+1) LVpotentialPeriod(LVjump(3)); LVpotentialPeriod(LVjump(3)+1) LVpotentialPeriod(end)];
TisorelaxLV = mean([LVperiod(2,2)-LVperiod(2,1) LVperiod(4,2)-LVperiod(4,1)]);
TisocontractLV = LVperiod(3,2)-LVperiod(3,1);

interpRV = griddedInterpolant(t,V_RV,'pchip');
V_RVnew = interpRV(Tunit)';
currentRV_min = min(V_RVnew);
currentRV_max = max(V_RVnew);
V_RVnew_scaled = (V_RVnew - currentRV_min) * (new_max - new_min) / (currentRV_max - currentRV_min) + new_min;
dVRV = diff(V_RVnew_scaled);
RVpotentialPeriod = Tunit([false; abs(dVRV) < 0.02]);
[~,RVjump] = maxk(diff(RVpotentialPeriod),3);
RVjump = sort(RVjump);
RVperiod = [RVpotentialPeriod(1) RVpotentialPeriod(RVjump(1)); RVpotentialPeriod(RVjump(1)+1) RVpotentialPeriod(RVjump(2));...
    RVpotentialPeriod(RVjump(2)+1) RVpotentialPeriod(RVjump(3)); RVpotentialPeriod(RVjump(3)+1) RVpotentialPeriod(end)];
TisorelaxRV = mean([RVperiod(2,2)-RVperiod(2,1) RVperiod(4,2)-RVperiod(4,1)]);
TisocontractRV = RVperiod(3,2)-RVperiod(3,1);

% This is used to prevent the isocontract and isorelax phases from disappearing.
if params.R_t_c > 1e4 && params.R_p_c > 5e2
    if  TisocontractRV/T <0.025 || TisorelaxRV/T <0.025
        error("unreal condition")
    end
elseif params.R_m_c > 1e4 && params.R_a_c > 1e4
    if  TisocontractLV/T <0.025 ||TisorelaxLV/T <0.025
        error("unreal condition")
    end
end


%% Collect simulation outputs

output_no = 50;
o = zeros(output_no,length(t)); % outputs from simulation
for i = 1:length(t)
    [~,o(:,i)] = dXdT(t(i),y(i,:), params);
end

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
P_external = o(50, :)';

%% Simulation outputs requiring post-processing for cross-valve flow

end_beat_i = find(t >= 1.02*T, 1) - 1; % index for end of one complete cardiac cycle, sometime flow shift and a wave is in the middle of the T
[Qm_maxima, Qm_maxima_i,Qm_wid,Qm_prom] = findpeaks(Q_m(1:end_beat_i),'MinPeakHeight',max(Q_m)/5);
[Qt_maxima, Qt_maxima_i] = findpeaks(Q_t(1:end_beat_i),'MinPeakHeight',max(Q_m)/5);

% E/A ratio
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
    % dQm = gradient(Q_m, t);
    % [~,locsNeg] = findpeaks(-dQm,'MinPeakHeight',500);
    % [~,locsPos] = findpeaks(dQm);
    % [~,NegPeaklocs] = sort(-dQm(locsNeg),1,"descend");
    % [~,PosPeaklocs] = sort(dQm(locsPos),1,"descend");
    % if abs(sum(t(locsNeg(NegPeaklocs(1:2)))-t(locsPos(PosPeaklocs(1:2))))) > 0.4
    %     NegPeaklocs = NegPeaklocs(1:2);
    % else
    %     NegPeaklocs = NegPeaklocs(3:4);
    % end
    % lastIdx = find(t(locsPos)-max(t(locsNeg(NegPeaklocs)))<0,1,"last");
    % E_A_ratio = Qm_maxima(1)/Q_m(locsPos(lastIdx));
    % 03/21 I attempted to save more simulations that lack E/A waves,
    % but the results turned out to be unrealistic.
    % The optimization process tends to produce only the E wave.
    % Surface adjustments might provide a minimal possible solution.
    E_A_ratio = 100;
end


% Mitral Valve
Qm_sign = sign(Q_m);
if(Qm_sign(1) <= 0)
    Qm_pos_start = find(Qm_sign == 1, 1);
    Qm_sign = Qm_sign(Qm_maxima_i(end):end);
    Qm_pos_end = find(Qm_sign ~= 1, 1) + Qm_maxima_i(end) - 2;
    Qm_neg_start = Qm_pos_end + 1;
    Qm_sign = Qm_sign(Qm_neg_start - Qm_maxima_i(end) + 1: end);
    Qm_neg_end = find(Qm_sign == 1, 1) + Qm_neg_start - 2;
else
    Qm_neg_start = find(Qm_sign ~= 1, 1);
    Qm_sign = Qm_sign(Qm_neg_start: end);
    Qm_neg_end = find(Qm_sign == 1, 1) + Qm_neg_start - 2;
    Qm_pos_start = Qm_neg_end + 1;
    Qm_sign = Qm_sign(Qm_maxima_i(end) - Qm_neg_start + 1:end);
    Qm_pos_end = find(Qm_sign ~= 1, 1) + Qm_maxima_i(end) - 2;
end
Qm_pos = Qm_pos_start: Qm_pos_end; % indices for positive mitral flow
Qm_neg = Qm_neg_start: Qm_neg_end; % indices for negative mitral flow

% Aortic valve
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
Qa_pos = Qa_pos_start: Qa_pos_end; % indices for positive aortic flow
Qa_neg = Qa_neg_start: Qa_neg_end; % indices for negative aortic flow
DNA = (P_SA(Qa_pos_end)+P_LV(Qa_pos_end))/2;% dicrotic notch for systemic artery

% Tricuspid valve
Qt_sign = sign(Q_t);
if(Qt_sign(1) <= 0)
    Qt_pos_start = find(Qt_sign == 1, 1);
    Qt_sign = Qt_sign(Qt_maxima_i(end): end);
    Qt_pos_end = find(Qt_sign ~= 1, 1) + Qt_maxima_i(end) - 2;
    Qt_neg_start = Qt_pos_end + 1;
    Qt_sign = Qt_sign(Qt_neg_start - Qt_maxima_i(end) + 1: end);
    Qt_neg_end = find(Qt_sign == 1, 1) + Qt_neg_start - 2;
else
    Qt_neg_start = find(Qt_sign ~= 1, 1);
    Qt_sign = Qt_sign(Qt_neg_start: end);
    Qt_neg_end = find(Qt_sign == 1, 1) + Qt_neg_start - 2;
    Qt_pos_start = Qt_neg_end + 1;
    Qt_sign = Qt_sign(Qt_maxima_i(end) - Qt_neg_start + 1:end);
    Qt_pos_end = find(Qt_sign ~= 1, 1) + Qt_maxima_i(end) - 2;

end
Qt_pos = Qt_pos_start: Qt_pos_end; % indices for positive tricuspid flow
Qt_neg = Qt_neg_start: Qt_neg_end; % indices for negative tricuspid flow

% Pulmonary Valve
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
Qp_pos = Qp_pos_start: Qp_pos_end; % indices for positive pulmonary flow
Qp_neg = Qp_neg_start: Qp_neg_end; % indices for negative pulmonary flow
DNP = (P_PA(Qp_pos_end)+P_RV(Qp_pos_end))/2;% Dicrotic notch for pulmonary artery

% Quantify valve stenosis
% MS
MPG_m = trapz(t(Qm_pos), P_LA(Qm_pos) - P_LV(Qm_pos)) / (t(Qm_pos_end) - t(Qm_pos_start));% MVmg

% AS
MPG_a = trapz(t(Qa_pos), P_LV(Qa_pos) - P_SA(Qa_pos)) / (t(Qa_pos_end) - t(Qa_pos_start));
Peak_AG_p = max(P_LV(Qp_pos) - P_SA(Qp_pos)); % AVpg

% TS
MPG_t = trapz(t(Qt_pos), P_RA(Qt_pos) - P_RV(Qt_pos)) / (t(Qt_pos_end) - t(Qt_pos_start));% TVmg may not real measured

% PS
MPG_p = trapz(t(Qp_pos), P_RV(Qp_pos) - P_PA(Qp_pos)) / (t(Qp_pos_end) - t(Qp_pos_start));
Peak_PG_p = max(P_RV(Qp_pos) - P_PA(Qp_pos)); % PVpg


% Quantify regurgitation fraction
% MR
RVol_m = -trapz(t(Qm_neg), Q_m(Qm_neg)); % regurgitant volume over period where Qm is negative
SV_LA_pos = trapz(t(Qm_pos), Q_m(Qm_pos));
SV_LA = max(V_LA) - min(V_LA);
RF_m = 100 * (RVol_m / SV_LA_pos);

% AR
RVol_a = -trapz(t(Qa_neg), Q_a(Qa_neg));
SV_LV_pos = trapz(t(Qa_pos), Q_a(Qa_pos));
RF_a = 100 * (RVol_a / SV_LV_pos);

% TR
RVol_t = -trapz(t(Qt_neg), Q_t(Qt_neg));
SV_RA_pos = trapz(t(Qt_pos), Q_t(Qt_pos));
RF_t = 100 * (RVol_t / SV_RA_pos);

% PR
RVol_p = -trapz(t(Qp_neg), Q_p(Qp_neg));
SV_RV_pos = trapz(t(Qp_pos), Q_p(Qp_pos));
RF_p = 100 * (RVol_p / SV_RV_pos);

%% Other simulation outputs requiring post-processing

SV_LV_tot = max(V_LV) - min(V_LV);
SV_RV_tot = max(V_RV) - min(V_RV);

EF_LV = SV_LV_tot / max(V_LV);
EF_RV = SV_RV_tot / max(V_RV);


CO = (SV_LV_pos-RVol_a) * HR / 1000;
CO_RV = (SV_RV_pos-RVol_p) * HR / 1000;

MPAP = (1/3) * (max(P_PA)) + (2/3) * (min(P_PA));

[~, LVED_i] = max(V_LV); % index of end diastole.
[~, RVED_i] = max(V_RV);
[~, LVES_i] = min(V_LV);
[~, RVES_i] = min(V_RV);
Y_max_i = find(abs(Y - 1) < 0.001, 1); % max of activation function doesn't quite coincide with max pressure development
RVEDP_i = find(Y > 0, 1); % Should be the first or second index

% Calculate the inner diameters of the LV and RV using TriSeg geometry model parameters.
LVIDs = -xm_LV(LVES_i) + xm_SEP(LVES_i)-...
    (1/Cm_SEP(LVES_i)- r_SEP(LVES_i)...% SEP miwall- inner radias
    +(-1/Cm_LV(LVES_i))-r_LV(LVES_i)); % + LV midwall- inner radias

RVIDd = r_RV(RVED_i)...% RV inner r
    + ((xm_RV(RVED_i)-1/Cm_RV(RVED_i))-xm_SEP(RVED_i))... % center of the RV - orginal - xm_sep
    - (d_SW(RVED_i)-(1/Cm_SEP(RVED_i)-r_SEP(RVED_i))); % thickness-(midwall SEP-innerwall)
LVIDd = -xm_LV(LVED_i) + xm_SEP(LVED_i)-...
    (1/Cm_SEP(LVED_i)- r_SEP(LVED_i)...% SEP miwall- inner radias
    +(-1/Cm_LV(LVED_i))- r_LV(LVED_i)); % + LV midwall- inner radias

% This is because xm_sep at end-diastole could be either positive or negative, which is unknown
% before calculation. Therefore, penalize any unrealistic conditions.
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
% Calculate the masses and thicknesses of the LV and RV using TriSeg geometry model parameters.
rho_myo = 1.0550;
Hed_LW = d_LW(LVED_i);
Hed_SW = d_SW(LVED_i);
Hed_RW = d_RW(RVED_i);
LV_m = rho_myo * (params.Vw_LV+params.Vw_SEP);
RV_m = rho_myo * (params.Vw_RV);


% Pressure waveforms PCWP, RAP a and v wave (not very well-implemented)
% 11/07/2024: This feature is not utilized in version 1.0.1 but will be improved in future versions.
% hack version (V1.0.0)

%RAP
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

%% Calculate cost function metrics

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
o_vals.RAPmean = trapz(t,P_RA)/(t(end)-t(1));
o_vals.PASP = max(P_PA);
o_vals.PADP = min(P_PA);
o_vals.PCWP = trapz(t,P_PV)/(t(end)-t(1));
o_vals.PCWPmax = max(P_PV);
o_vals.PCWPmin = min(P_PV);
o_vals.CVP = trapz(t,P_SV)/(t(end)-t(1));
o_vals.CVPmax = max(P_SV);
o_vals.CVPmin = min(P_SV);
o_vals.CO = CO_RV;% CO from RHC report is RV
o_vals.Hed_LW = d_LW(LVED_i);
o_vals.Hed_SW = d_SW(LVED_i);
o_vals.Hed_RW = d_RW(RVED_i);
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
o_vals.LVIDd = LVIDd;
o_vals.LVIDs = LVIDs;
o_vals.RVIDd = RVIDd;
o_vals.RVIDs = RVIDs;
o_vals.RV_m = RV_m;
o_vals.LV_m = LV_m;
o_vals.FakeLV_m = LV_m;
o_vals.Vtot = Vtot;
o_vals.DNA = DNA;
o_vals.DNP = DNP;


% RAP, PCWP a and v wave
o_vals.RAP_a = RAP_a;
o_vals.RAP_v = RAP_v;
o_vals.PCWP_a = PCWP_a;
o_vals.PCWP_v = PCWP_v;

% Convery RF to grades
g = [1, 2, 3, 4]; % grades: none, mild, moderate, severe.

% Mitral
gmt_MR = [0, 20, 40, 60];
gmt_MS = [2.5, 3.75, 7.5, 12];
% Aortic
gmt_AR = [0, 20, 40, 60];
gmt_AS = [2.25, 10, 30, 50];
% Tricuspid
gmt_TR = [0, 17.5, 35, 52.5];
% Pulmonary
gmt_PR = [0, 15, 30 ,45];
gmt_PS = [0.67, 22, 50, 78]; % Peak pressure gradient

% Interpret RF into grades
o_vals.MVr = interp1(gmt_MR, g, RF_m, 'linear', 'extrap');
o_vals.MS = interp1(gmt_MS, g, MPG_m, 'linear', 'extrap');
o_vals.AVr = interp1(gmt_AR, g, RF_a, 'linear', 'extrap');
o_vals.AS = interp1(gmt_AS, g, MPG_a, 'linear', 'extrap');
o_vals.TVr = interp1(gmt_TR, g, RF_t, 'linear', 'extrap');
o_vals.PVr = interp1(gmt_PR, g, RF_p, 'linear', 'extrap');
o_vals.PS = interp1(gmt_PS, g, Peak_AG_p, 'linear', 'extrap');

outputno = struct2array(o_vals);
if any(~isreal(outputno))
    error('bad guessing')
end
%% Implent
% Adding additional costs to enforce synchronization:
% The goal is to ensure synchronized contraction and relaxation
% among the left ventricle (LV), right ventricle (RV), and septum (SEP).

[~,locsTSA] = max(P_SA);
[~,locsTPA] = max(P_PA);
if t(locsTSA) < T
    T2maxPSA = t(locsTSA);
else
    T2maxPSA = t(locsTSA)-T;
end

if T2maxPSA/T <= 0.2
    error('KactLV may too big')
end


if t(locsTPA) < T
    T2maxPPA = t(locsTPA);
else
    T2maxPPA = t(locsTPA)-T;
end

if T2maxPPA/T <= 0.2
    error('KactRV may too big')
end


Tint = (t(1):0.001:t(end));
NewPLV = interp1(t,P_LV,Tint);
NewPRV = interp1(t,P_RV,Tint);
Newd_LW = interp1(t,d_LW,Tint);
Newd_SW = interp1(t,d_SW,Tint);
Newd_RW = interp1(t,d_RW,Tint);
[~, locsMinPLV] = min(NewPLV);
if locsMinPLV >length(Tint)/2
    locsMinPLV = round(locsMinPLV-length(Tint)/2);
end
[~, locsMinPRV] = min(NewPRV);
if locsMinPRV >length(Tint)/2
    locsMinPRV = round(locsMinPRV-length(Tint)/2);
end

Maggicpoint = round(mean([locsMinPRV; locsMinPLV]));
diff_d_LW = diff(Newd_LW);
diff_d_SW = diff(Newd_SW);
diff_d_RW = diff(Newd_RW);

[~,locswiredpeakLW]=findpeaks(diff_d_LW,'MinPeakHeight',max(diff_d_LW)/3);
if ~isempty(locswiredpeakLW)
    if locswiredpeakLW(1) > Maggicpoint-round(length(Tint)/2*0.035) &&...
            locswiredpeakLW(1) < Maggicpoint+round(length(Tint)/2*0.015)
        error("Unreal LV Movement")
    end
end

[~,locswiredpeakSW]=findpeaks(diff_d_SW,'MinPeakHeight',max(diff_d_SW)/3);
if ~isempty(locswiredpeakSW)
    if locswiredpeakSW(1) > Maggicpoint-round(length(Tint)/2*0.035) &&...
            locswiredpeakSW(1) < Maggicpoint+round(length(Tint)/2*0.015)
        error("Unreal SEP Movement")
    end
end

[~,locswiredpeakRW]=findpeaks(diff_d_RW,'MinPeakHeight',max(diff_d_RW)/3);
if ~isempty(locswiredpeakRW)
    if locswiredpeakRW(1) > Maggicpoint-round(length(Tint)/2*0.035) &&...
            locswiredpeakRW(1) < Maggicpoint+round(length(Tint)/2*0.015)
        error("Unreal RV Movement")
    end
end
%% Function to kill inf loop of ode15s

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