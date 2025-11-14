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
%% Solve the differential equations using the ODE solver
try
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
    lastwarn('');   % 清空上一个warning
    [t, y] = odeWithTimeout1(@dXdTOnGL, [0, 30*T], init_vec, options, params, maxTime);
    [lastWarnMsg, lastWarnId] = lastwarn;  % 检查warning
    if ~isempty(lastWarnMsg)
        error(['ODE solver warning: ', lastWarnMsg]);  % 直接终止，并抛出warning信息
    end

    % Find the index where t is within the last two periods, which reflects the steady state
    startIndex = find(t >= t(end) - 2*T, 1, 'first');
    lastTwoPeriodsT = t(startIndex:end);
    lastTwoPeriodsY = y(startIndex:end, :);
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

    %% Collect simulation outputs

    output_no = 50;
    o = zeros(output_no,length(t)); % outputs from simulation
    for i = 1:length(t)
        [~,o(:,i)] = dXdTOnGL(t(i),y(i,:), params);
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
catch ME1
    disp(['runSim failed: ', ME1.message]);
    %% Solve the differential equations using the ODE solver
    T = params.T;
    HR = params.HR;
    init_vec = cell2mat(struct2cell(init))';
    iniGeo = init_vec(1:4);
    init_vec  = init_vec(5:end);
    M = eye(length(init_vec));
    options = odeset('Mass', M, 'MassSingular', 'yes', 'RelTol', 1e-7, 'AbsTol', 1e-7, 'MaxStep', params.T / 30);
    maxTime = 100;

    [t, y, ~, ~, ~, o] = odeWithTimeout2(@dXdT, [0, 22 * params.T], init_vec, options, params, iniGeo, maxTime);

    %% Collect simulation outputs from the last two cardiac cycles
    % Identify the starting index for the last two periods
    startIndex = find(t >= t(end) - 2*T, 1, 'first');
    lastTwoPeriodsT = t(startIndex:end);
    lastTwoPeriodsY = y(startIndex:end, :);
    lastTwoPeriodsO = o(:,startIndex:end);

    % Shift time to start from zero and update variables
    t = lastTwoPeriodsT - lastTwoPeriodsT(1);
    y = lastTwoPeriodsY;        % ODE solutions
    o = lastTwoPeriodsO;        % Post-processed outputs

    % Extract state variables
    Lsc_LV  = y(:,1);
    Lsc_SEP = y(:,2);
    Lsc_RV  = y(:,3);
    V_LV    = y(:,4);
    V_RV    = y(:,5);
    V_SA    = y(:,6);
    V_SV    = y(:,7);
    V_PA    = y(:,8);
    V_PV    = y(:,9);
    V_LA    = y(:,10);
    V_RA    = y(:,11);
    Vtot    = sum(y(end, 4:11));  % Total blood volume at the final time point

    % Parse each row of `o` into human-readable vectors or matrices
    % 1–4: Septal geometry
    xm_LV  = o(1,:)';
    xm_SEP = o(2,:)';
    xm_RV  = o(3,:)';
    ym     = o(4,:)';

    % 5–12: Pressures in different compartments
    P_LV = o(5,:)';
    P_SA = o(6,:)';
    P_SV = o(7,:)';
    P_RV = o(8,:)';
    P_PA = o(9,:)';
    P_PV = o(10,:)';
    P_RA = o(11,:)';
    P_LA = o(12,:)';

    % 13–15: Myocardial wall volumes
    Vm_LV  = o(13,:)';
    Vm_SEP = o(14,:)';
    Vm_RV  = o(15,:)';

    % 16–18: Myocardial wall areas
    Am_LV  = o(16,:)';
    Am_SEP = o(17,:)';
    Am_RV  = o(18,:)';

    % 19–21: Wall curvatures
    Cm_LV  = o(19,:)';
    Cm_SEP = o(20,:)';
    Cm_RV  = o(21,:)';

    % 22–24: Fiber strains
    eps_LV  = o(22,:)';
    eps_SEP = o(23,:)';
    eps_RV  = o(24,:)';

    % 25–27: Passive fiber stresses
    sigma_pas_LV  = o(25,:)';
    sigma_pas_SEP = o(26,:)';
    sigma_pas_RV  = o(27,:)';

    % 28–30: Active fiber stresses
    sigma_act_LV  = o(28,:)';
    sigma_act_SEP = o(29,:)';
    sigma_act_RV  = o(30,:)';

    % 31–33: Total fiber stresses (passive + active)
    sigma_LV  = o(31,:)';
    sigma_SEP = o(32,:)';
    sigma_RV  = o(33,:)';

    % 34–37: Valve flows
    Q_m = o(34,:)';  % Mitral valve
    Q_a = o(35,:)';  % Aortic valve
    Q_t = o(36,:)';  % Tricuspid valve
    Q_p = o(37,:)';  % Pulmonary valve

    % 38–41: Circulatory flows
    Q_SA  = o(38,:)';  % Systemic arterial flow
    Q_PA  = o(39,:)';  % Pulmonary arterial flow
    QIN_RA = o(40,:)'; % Inflow to right atrium
    QIN_LA = o(41,:)'; % Inflow to left atrium

    % 42–44: Tensions in the x-direction
    Tx_LV  = o(42,:)';
    Tx_SEP = o(43,:)';
    Tx_RV  = o(44,:)';

    % 45–47: Tensions in the y-direction
    Ty_LV  = o(45,:)';
    Ty_SEP = o(46,:)';
    Ty_RV  = o(47,:)';

    % 48–49: Activation functions
    Y   = o(48,:)';
    act = o(49,:)';

    % 50–52: Wall thickness
    d_LW = o(50,:)';
    d_SW = o(51,:)';
    d_RW = o(52,:)';

    % 53–55: Curvature radii
    r_LV  = o(53,:)';
    r_SEP = o(54,:)';
    r_RV  = o(55,:)';

    % 56: Pericardial constraint pressure
    P_external = o(56,:)';
end


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

%% Function to kill inf loop of ode15s

function [T,Y,TE,YE,IE] = odeWithTimeout1(odefun, tspan, y0, options, params, maxTime)
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

function [T, Y, TE, YE, IE, o] = odeWithTimeout2(odefun, tspan, y0, options, params, iniGeo, maxTime)
% Function to run ODE solver with timeout and collect custom outputs

% Variable to collect output
o_internal = [];

% Start timer
ticID = tic;

% Set up event function
options = odeset(options, 'Events', @(t, y) eventFunction(t, y, maxTime, ticID));

% Set custom output function (uses nested function to access o_internal)
options = odeset(options, 'OutputFcn', @(t, y, flag) myOutputFcn(t, y, flag));

% Run the solver
[T, Y, TE, YE, IE] = ode15s(@(t, y) odefun(t, y, params, iniGeo), tspan, y0, options);

% Return collected output
o = o_internal;

% --------- Nested OutputFcn ---------
    function status = myOutputFcn(t, y, flag)
        switch flag
            case 'init'
                [~, output0] = odefun(0, y0, params, iniGeo);
                o_internal = output0;

            case ''
                for i = 1:length(t)
                    [~, output_now] = odefun(t(i), y(:, i), params, iniGeo);
                    o_internal(:, end+1) = output_now;
                end

            case 'done'
                % nothing to do
        end
        status = 0;
    end
end