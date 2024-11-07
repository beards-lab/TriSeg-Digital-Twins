function [Error,params, init] = estimParams(targets, inputs, mods, modifiers)
%% Function Purpose:
% This function assigns and calculates the nominal values for all model parameters. 
% Inputs to this function are the outputs from the targetVals.m function. 
% Outputs of this function include a structure of parameters used in the model and initial
% conditions for the ODE solver (dXdT.m function).

% Created by EB Randall, modified by Filip Jezek, Andrew Meyer, and Feng Gu
% Last modified: 10/29/2024

% Inputs: 
% targets     - Patient measurements being fit to 
% inputs      - Other necessary variables to build the model 
% mods        - Cell array of names to adjust selected parameters 
% modifiers   - Vector of floats that adjust selected parameters

% Outputs: 
% params      - Structure of parameters used in the model 
% init        - Initial conditions for the ode15s solver (dXdT.m function) 
% Error       - Used for debugging the TriSeg geometry

% Related functions: 
% geom_0      - Computes the initial guess (idealized end-diastolic state) of TriSeg geometry 
% calc_xm_ym  - Computes specific geometrical initial conditions of end-diastolic and end-systolic
%               states, considering the balance of volume and sarcomere length
%% Unpack data structure
% In this function, we want to use targets and inputs interchangeably. Thus,
% creating a merged struct (inputData). The 'targets' struct is explicitly for
% the cost function of runSim.m
cat_structs = @(S1,S2) cell2struct([struct2cell(S1);struct2cell(S2)], [fieldnames(S1);fieldnames(S2)]);
inputData = cat_structs(targets, inputs);

%% Load data used for estimating parameters initial values
HR = inputData.HR; % 1/min  

% Blood pressures (mmHg)
P_SAs = inputData.SBP; 
P_SAd = inputData.DBP;

if(isfield(inputData,'PASP'))
    PP_pul = inputData.PASP - inputData.PADP;
else
    assert(isfield(inputData,'RVSP'))
    PP_pul = inputData.RVSP - inputData.PADP;
end
PP_sys = P_SAs - P_SAd;
 
CVP = inputData.CVP;

% End-diastolic and end-systolic and volumes (mL) 
LVESV = inputData.LVESV;                              
LVEDV = inputData.LVEDV;   RVEDV = inputData.RVEDV;

RAVmin = inputData.RAVmin;
if(isfield(inputData,'RAVmax'))
    RAVmax = inputData.RAVmax;
end
LAVmin = inputData.LAVmin;
if(isfield(inputData,'LAVmax'))
    LAVmax = inputData.LAVmax;
end

% Cardiac output (mL s^(-1))
if(isfield(inputData, 'CO'))
    CO    = 1000 * inputData.CO / 60;
    SV = 60 * CO / HR; 
else % this branch should only be executed for canonical male/female
    SV = LVEDV - LVESV; 
    CO = SV * HR / 60;
end

%% Load heart mechanical parameters from TriSeg papers

% Sarcomere length parameters (µm)
Lsref   = 2;    % coming from Lumen TriSeg paper
Lsc0    = 1.51; % reference length when maximally stressed
Lse_iso = 0.04; 

% Passive stress steepness parameter  
gamma = 7.5; % optimized from ex vivo model 

% Percentage of cardiac cycle, this is coming from Guyton physiology textbook
k_TS = 0.35; % Beginning of cardiac cycle to maximal systole  
k_TR = 0.15; % Relaxation time fraction 

%% Calculate patient-specific TriSeg heart geometrical parameters

% Estimate heart geometrical parameters based on an idealized TriSeg model at the end-diastolic
% state

% The following loop is for patients who don't have RV mass (RV_m) or the thickness of the RV free
% wall
% The goal is to estimate kpas_RV and backwards estimate the right ventricular wall volume (Vw_RV)
if(~isfield(inputData, 'RV_m')) 
    % Thicknesses are from echo, wall volumes are based on rvedv estimate
    if(isfield(inputData,'Hed_RW'))
        H_RW = inputData.Hed_RW;
        right = H_RW;
        use_Vw_RV = 0;
    else
        use_Vw_RV = 0;
        coeff = [14853.8592619735
            14531.7636186220
            14769.3961036171
            11440.8528760998
            -21873.4802314862
            456.030613986307
            -4913.60117297433
            156836.996351639
            -155270.878092212];
        k_RV = coeff(1)*100*(inputData.LVEDV-inputData.LVESV)/inputData.LVEDV+...
            coeff(2)*inputData.PASP+...
            coeff(3)*inputData.PADP+...
            coeff(4)*inputData.PCWP+...
            coeff(5)*inputData.CO+...
            coeff(6)*inputData.SBP+...
            coeff(7)*inputData.DBP+...
            100*coeff(8)*max([inputData.Hed_SW inputData.Hed_LW])+...
            100*coeff(9)*min([inputData.Hed_SW inputData.Hed_LW]);
        if inputData.Sex == 1
            k_RV = k_RV/7.2906; 
        else
            k_RV = k_RV/11.2290;
        end 
        % The above part is not ideal, but it works to get an estimation for kpas_RV and,
        % consequently, Vw_RV. The regression model used is based on Lasso regression applied to
        % patients who have Vw_RV measurements.
        % The regression results are not great. When applied to canonical patients, the estimates
        % are 7.2906 and 11.2290 times larger than the actual values. Unfortunately, we have no
        % choice but to calibrate based on these two scaling factors. This part definitely needs
        % improvement in the future.
        if (isfield(inputData,'RVEDP'))
            P_RV = inputData.RVEDP;
        elseif (isfield(inputData,'RAPmean'))
            P_RV = inputData.RAPmean;
        elseif (isfield(inputData,'PCWP'))
            P_RV = inputData.PASP*inputData.PCWP/inputData.SBP;
        end
        right = [P_RV k_RV];
    end
else
    rho_myo = 1.055; % g/mL
    Vw_RV = inputData.RV_m / rho_myo;
    right = Vw_RV;
    use_Vw_RV = 1;
end

if(~isfield(inputData, 'LV_m'))
    % Thicknesses are from echo, wall volumes are based on lvedv estimate
    if(isfield(inputData,'Hed_LW')) % always have both measurements
        H_LW_and_SW = (inputData.Hed_LW + inputData.Hed_SW) / 2; % This doesn't make sense if we're trying to fit both independently
    else
        assert(false); 
    end
    left = H_LW_and_SW;
    use_Vw_LV = 0;
else
    rho_myo = 1.055; % g/mL
    Vw_LV_and_SEP = inputData.LV_m / rho_myo; % for measured LV mass
    left = Vw_LV_and_SEP;
    use_Vw_LV = 1;
end

LvSepR = 2/3; % assume that at the end-diastolic state, the LV/Septum free wall ratio is 2/3

% Estimate Vw, which is assumed to be constant during a cardiac cycle.
% 'Left' could represent either Vw_LV_and_SEP (volume of the left ventricle and septum) or
% H_LV_and_SEP (average thickness of the left ventricle and septum). 
% 'Right' could represent either Vw_RV (volume of the right ventricle) or a combination of kpas_RV
% and RVEDP.
[Error,~, Vw0, ~] = geom_0(LVEDV, RVEDV, left, use_Vw_LV, right, use_Vw_RV, LvSepR,inputs);
Vw_LV = Vw0(1);
Vw_SEP = Vw0(2);
Vw_RV = Vw0(3);

% Adjustments for the LV/Septum free wall ratio
% The following loop is mainly used during optimization. With the optimized LV/Septum free wall
% ratio, the reference area (Amref), wall volume (Vw), and dimensions (Xm, Ym) are finally computed.
LvSepR_i = find(ismember(mods,{'LvSepR'}),1);
if ~isempty(LvSepR_i) % Vw is adjustable. 
    % Set boundary conditions to ensure that the LV/Septum free wall ratio remains reasonable
    if LvSepR * modifiers(LvSepR_i)<.9 && LvSepR * modifiers(LvSepR_i)>.1
        LvSepR = LvSepR * modifiers(LvSepR_i);
    else
        LvSepR = 2/3;
    end
    right = Vw_RV;
    left = Vw_LV + Vw_SEP;
    use_Vw_RV = 1;
    use_Vw_LV = 1;
    % All input for geom_0 is fixed except for LvSepR
    [~,Amref0, Vw0, dim0] = geom_0(LVEDV, RVEDV, left, use_Vw_LV, right, use_Vw_RV, LvSepR,inputs); 
end   
%% Approximations for initial displacements in end-diastole state. 
% The above assumptions about initial conditions consider the LV and septum as a single sphere and
% do not account for sarcomere length balance. Therefore, we introduced another function,
% calc_xm_ym.m, which ensures that Vw, Amref, and lumen volume remain unchanged while allowing for
% changes in the LV and septum shapes. In other words, this function ensures that all three layers
% (LV, septum, and RV) are still parts of spheres, but they will be independent and form three
% different spheres. This balances the geometry, which is reflected in geom_0.m, and maintains
% sarcomere length to 2um at end-diastole.
% The output from this section also serves as the initial conditions x(1) to x(7) for the ODE solver
% (dXdT.m)

Amref_LV = Amref0(1);
Amref_SEP = Amref0(2);
Amref_RV = Amref0(3);

Vw_LV = Vw0(1);
Vw_SEP = Vw0(2);
Vw_RV = Vw0(3);

xm_LV_d_0 = dim0(1);
xm_SEP_d_0 = dim0(2);
xm_RV_d_0 = dim0(3);
ym_d_0 = dim0(4);

dias = 1;
x0_d = [xm_LV_d_0; 
    xm_SEP_d_0; 
    xm_RV_d_0;
    ym_d_0;
    Amref_RV]; 

% If RVEDV is a target, we will use Amref_RV as an input to calc_xm_ym, and not let it be adjusted
% by fsolve
if(isfield(targets, 'RVEDV')) 
    assert(~any(ismember(mods,{'Amref_RV'})));
    fix_AmrefRV = 1; % use as an input
else
    fix_AmrefRV = 0; % use as a state
end

Vw    = [Vw_LV,Vw_SEP,Vw_RV]; 
Amref = [Amref_LV,Amref_SEP]; 
if fix_AmrefRV
    Amref = [Amref, Amref_RV];
end

% Assume end-diastolic sarcomere length 
SL_d    = 2; % µm 
opts = optimoptions('fsolve','Display','none',...
    'MaxFunctionEvaluations',2e3,'Algorithm','levenberg-marquardt'); 
[d0,~] = fsolve(@(x) calc_xm_ym(x,Lsref, Vw, Amref,SL_d,LVEDV, RVEDV,fix_AmrefRV,dias),x0_d,opts); % end-diastole. Always use RVEDV

% Outputs / Diastolic displacements
xm_LV_d  = d0(1);
xm_SEP_d = d0(2);
xm_RV_d  = d0(3);
ym_d     = d0(4);
Amref_RV = d0(5); %was either never used or  was adjusted 
% Sarcomere lengths (µm)
Lsc_LV_0 = SL_d;
Lsc_SEP_0 = SL_d;
Lsc_RV_0 = SL_d;    

%% Calculate initial condition of kpas for the LV and RV at end-diastole based on TriSeg papers 

% Midwall surface area (cm^2)
Am_LV_d = pi * (xm_LV_d^2  + ym_d^2);
Am_RV_d = pi * (xm_RV_d^2  + ym_d^2);       

% Midwall curvature (cm^(-1))
Cm_LV_d = 2 * xm_LV_d  / (xm_LV_d^2  + ym_d^2);
Cm_RV_d = -2 * xm_RV_d  / (xm_RV_d^2  + ym_d^2);

% Midwall ratio (dimensionless) 
z_LV_d = 3 * Cm_LV_d  * Vw_LV  / (2 * Am_LV_d); 
z_RV_d = 3 * Cm_RV_d  * Vw_RV  / (2 * Am_RV_d); 

% Instantaneous sarcomere length (µm) in end-diastole
Ls_LV_d = SL_d;         
Ls_RV_d = SL_d;  

% Passive stress 
sigma_pas_LV_d = (Ls_LV_d/Lsc0 - 1)^gamma; %EDPVR. Ls_LV_d >= Lsc0
sigma_pas_RV_d = (Ls_RV_d/Lsc0 - 1)^gamma; 

% Dimensionless combination function
Gamma_LV_d = -(2 / 3) * z_LV_d * (1 + (1 / 3) * z_LV_d^2 + (1 / 5) * z_LV_d^4);
Gamma_RV_d = -(2 / 3) * z_RV_d * (1 + (1 / 3) * z_RV_d^2 + (1 / 5) * z_RV_d^4);

% Passive stress scaling parameters
% k_pas = inputData.PCWP / (Gamma_LV_d * sigma_pas_LV_d); 
k_pas_LV = inputData.PCWP / (Gamma_LV_d * sigma_pas_LV_d); 
if(isfield(inputData,'RVEDP'))&&(~(inputData.RVEDP == 0))
    k_pas_RV = inputData.RVEDP / (Gamma_RV_d * sigma_pas_RV_d);
else % should only execute for canonical male/female
    assert(isfield(inputData,'CVP'))
    k_pas_RV = CVP / (Gamma_RV_d * sigma_pas_RV_d);
end

%% Approximations for initial displacements and Amref_rv in end-systole 
% Apply the same calc_xm_ym.m function but in end-systole state

% Initialize systolic displacements values (cm)
xm_LV_s_0 = dim0(1);
xm_SEP_s_0 = dim0(2);
xm_RV_s_0 = dim0(3);
ym_s_0 = dim0(4);

dias = 0;
x0_s = [xm_LV_s_0; 
    xm_SEP_s_0; 
    xm_RV_s_0;
    ym_s_0; 
    Amref_RV]; 

Amref = [Amref_LV,Amref_SEP, Amref_RV]; 
fix_AmrefRV = 1;
RVESV = inputData.RVESV; 
opts = optimoptions('fsolve','Display','none',...
    'MaxFunctionEvaluations',2e3,'Algorithm','levenberg-marquardt'); 
[s0,~] = fsolve(@(x) calc_xm_ym(x,Lsref,Vw,Amref,[],LVESV,RVESV,fix_AmrefRV,dias),x0_s,opts); 
xm_LV_s = s0(1);
xm_RV_s = s0(3);
ym_s    = s0(4);

%% Calculate inital conditions of k_act at end-systole state based on TriSeg papers

% Midwall surface area (cm^2)
Am_LV_s = pi * (xm_LV_s^2  + ym_s^2);
Am_RV_s = pi * (xm_RV_s^2  + ym_s^2);

% Midwall curvature (cm^(-1))
Cm_LV_s = 2 * xm_LV_s  / (xm_LV_s^2  + ym_s^2);
Cm_RV_s = - 2 * xm_RV_s  / (xm_RV_s^2  + ym_s^2);

% Midwall ratio (dimensionless)  
z_LV_s = 3 * Cm_LV_s  * Vw_LV  / (2 * Am_LV_s); 
z_RV_s = 3 * Cm_RV_s  * Vw_RV  / (2 * Am_RV_s); 

% Myofiber strain (dimensionless)
eps_LV_s = 0.5 * log(Am_LV_s  / Amref_LV) - (1/12) * z_LV_s^2  - 0.019 * z_LV_s^4; 
eps_RV_s = 0.5 * log(Am_RV_s  / Amref_RV) - (1/12) * z_RV_s^2  - 0.019 * z_RV_s^4; 

% Sarcomere length (µm)
Ls_LV_s  = Lsref * exp(eps_LV_s); 
Ls_RV_s  = Lsref * exp(eps_RV_s); 

% Activation function 
Y = .55; 

% Active stress 
sigma_act_LV_s = Y * (Ls_LV_s/Lsc0  - 1)*Lse_iso; 
sigma_act_RV_s = Y * (Ls_RV_s/Lsc0  - 1)*Lse_iso; 

% Dimensionless combination function 
Gamma_LV_s = - (2 / 3) * z_LV_s * (1 + (1 / 3) * z_LV_s^2 + (1 / 5) * z_LV_s^4);
Gamma_RV_s = - (2 / 3) * z_RV_s * (1 + (1 / 3) * z_RV_s^2 + (1 / 5) * z_RV_s^4);

% Active stress scaling parameters 
% k_act = P_SAs / (Gamma_LV_s * sigma_act_LV_s); 
k_act_LV = P_SAs / (Gamma_LV_s * sigma_act_LV_s);  
if(isfield(inputData,'RVSP'))
    k_act_RV = inputData.RVSP / (Gamma_RV_s * sigma_act_RV_s); 
else % should only execute in canonical male/female
    assert(isfield(inputData,'PASP'));
    k_act_RV = inputData.PASP / (Gamma_RV_s * sigma_act_RV_s); 
end

%% Estimate initial conditions x(8) to x(15) for the ODE solver (dXdT.m)

% Total blood volume (mL) 
Vtot  = inputData.TBV; 

% Snapped at end diastole - maximal ventricles, minimal atria
% Blood volume distribution values; sum total = 1.0 
d_SA = .15;              d_PA = .05; 
d_SV = .6;              d_PV = .2;
Vd = Vtot - LVEDV - RVEDV - RAVmin - LAVmin; % distributed volume (available for other compartments)

% Total compartment volumes 
V_SA_0 = d_SA*Vd;      V_PA_0 = d_PA*Vd; 
V_SV_0 = d_SV*Vd;      V_PV_0 = d_PV*Vd;

% Unstressed compartment volumes
V_SA_u = V_SA_0*0.7;   V_PA_u = V_PA_0*0.1; 
V_SV_u = V_SV_0*0.9;   V_PV_u = V_PV_0*0.9; 

% Stressed compartment volumes
V_SA_s = V_SA_0 - V_SA_u;  V_PA_s = V_PA_0 - V_PA_u;  
V_SV_s = V_SV_0 - V_SV_u;  V_PV_s = V_PV_0 - V_PV_u; 
                                            
%% Estimate initial conditions of elastances and compliances 

% Compliances (mL mmHg^(-1))       
C_SA = SV / PP_sys;
C_SV = V_SV_s/CVP; 
C_PA = SV / PP_pul;
C_PV = V_PV_s/inputData.PCWP;  

% Arteriolar resistances (mmHg s mL^(-1)) 
MAP = (P_SAs - P_SAd)/3 + P_SAd;
R_SA = (MAP - CVP)/CO;  % systemic resistance  
if(isfield(inputData,'PASP'))
    MPAP = (inputData.PASP - inputData.PADP)/3 + inputData.PADP; 
else
    assert(isfield(inputData,'RVSP'))
    MPAP = (inputData.RVSP - inputData.PADP)/3 + inputData.PADP; 
end
R_PA = (MPAP - inputData.PCWP)/CO;  % pulmonary resistance

if R_PA <0 % pat 59 MAP < PCWP
    R_PA = (MPAP*2 - inputData.PCWP)/CO;  % arbitrary assumption to ensure R_PA is positive
end

% Valve resistances (mmHg s mL^(-1)), o for open, c for closed 
if(~isfield(inputData, 'AVr') || inputData.AVr < 1.5)
    R_a_c = 65535; % 1 -> no regurgitation. 0% regurgitation fraction(RF)
elseif(inputData.AVr < 2.5)
    R_a_c = 1.8; % 2 -> mild. 20% RF
elseif(inputData.AVr < 3.5)
    R_a_c = 0.55; % 3 -> moderate. 40% RF
else 
    % assert(inputData.AVr <= 5 && inputData.AVr >= 4, 'Invalid AI Grade Input');
    R_a_c = 8e-2; % 4+ -> severe. 60% RF
end

if(~isfield(inputData, 'AVpg') || inputData.AVpg <= 4.5)
    R_a_o = 7.3e-3; % 1 -> no stenosis. 2.25 mmHg
elseif(inputData.AVpg < 20)
    R_a_o = 3.4e-2; % 2 -> mild. 10 mmHg
elseif(inputData.AVpg <= 40)
    R_a_o = 1.15e-1; % 3 -> moderate. 30 mmHg
else
    % assert(inputData.AVpg <= 60 && inputData.AVpg > 40, 'Invalid AS Grade Input');
    R_a_o = 2.25e-1; % 4+ -> severe. 50 mmHg
end

if(~isfield(inputData, 'MVr') || inputData.MVr < 1.5) % FIXME: should we have finite backflow resistance for grade 1.5?
    R_m_c = 65535; % 1 -> no regurgitation. 0% RF
elseif(inputData.MVr < 2.5)
%         R_m_c = 1.6; % 2 -> mild. 20% RF
    R_m_c = (1/0.20)*(R_SA + R_a_o);
elseif(inputData.MVr < 3.5)
%         R_m_c = .53; % 3 -> moderate. 40% RF
    R_m_c = (1/0.40)*(R_SA + R_a_o);
else
    % assert(inputData.MVr <= 5 && inputData.MVr >= 4, 'Invalid MR Grade Input');
%         R_m_c = 0.017; % 4+ -> severe. 60% RF
    R_m_c = (1/0.60)*(R_SA + R_a_o);
end

if(~isfield(inputData, 'MVmg') || inputData.MVmg <= 2.5)
    R_m_o = 1.6e-2; % 1 -> no stenosis. 2.5 mmHg.
elseif(inputData.MVmg < 5)
    R_m_o = 2.3e-2; % 2 -> mild. 3.75 mmHg
elseif(inputData.MVmg <= 10)
    R_m_o = 5.1e-2; % 3 -> moderate. 7.5 mmHg
else
    assert(inputData.MVmg <= 25 && inputData.MVmg > 10, 'Invalid MS Grade Input');
    R_m_o = 9.4e-2; % 4+ -> severe. 12 mmHg
end

if(~isfield(inputData, 'PVr') || inputData.PVr < 1.5)
    R_p_c = 5000; % 1 -> no regurgitation. 1% RF
elseif(inputData.PVr < 2.5)
    R_p_c = 0.38; % 2 -> mild. 15% RF
elseif(inputData.PVr < 3.5)
    R_p_c = 1.33e-1; % 3 -> moderate. 30% RF
else
    % assert(inputData.PVr <= 5 && inputData.PVr >= 4, 'Invalid PI Grade Input');
    R_p_c = 5e-2; % 4+ -> severe. 45% RF
end

if(~isfield(inputData, 'PVpg') || inputData.PVpg <= 5)
    R_p_o = 1.3e-3; % 1 -> no stenosis. 0.67 mmHg
elseif(inputData.PVpg < 36)
    R_p_o = 5.85e-2; % 2 -> mild. 22 mmHg
elseif(inputData.PVpg <= 64)
    R_p_o = 1.7e-1; % 3 -> moderate. 50 mmHg
else
    assert(inputData.PVpg <= 85 && inputData.PVpg > 64, 'Invalid PS Grade Input');
    R_p_o = 3.82e-1; % 4+ -> severe. 78 mmHg
end

if(~isfield(inputData, 'TVr') || inputData.TVr < 1.5)
    R_t_c = 65535; % 1 -> no regurgitation. 0% RF

elseif(inputData.TVr < 2.5) % 2 -> mild. 17.5% RF
    R_t_c = (1/0.175)*(R_PA + R_p_o);

elseif(inputData.TVr < 3.5) % 3 -> moderate. 35% RF
    R_t_c = (1/0.35)*(R_PA + R_p_o);
else % severe. 52.5% RF. 
    R_t_c = (1/0.525)*(R_PA + R_p_o);
end

% Tricuspid stenosis is mostly ignored here, as it is usually not evaluated in clinical settings.
if(~isfield(inputData, 'TVmg') || inputData.TVmg <= 1.5)
    R_t_o = 3e-3; % 1 -> no stenosis. 0.5 mmHg. 3e-3
else
    assert(inputData.TVmg >= 2, 'Invalid TS Grade Input');
    R_t_o = 6e-2; % 2+ -> stenosis. 6 mmHg
end

%% Output initial conditions of the ODE solver (init) and the initial guessed functional parameters (params)

% Initial conditions (end-diastole) for ode15s in runSim
init.xm_LV_d = xm_LV_d;
init.xm_SEP_d = xm_SEP_d;
init.xm_RV_d = xm_RV_d;
init.ym_d = ym_d;
init.Lsc_LV_0 = Lsc_LV_0;
init.Lsc_SEP_0 = Lsc_SEP_0;
init.Lsc_RV_0 = Lsc_RV_0;
init.LVEDV = LVEDV;
init.RVEDV = RVEDV;
init.V_SA_s = V_SA_s;
init.V_SV_s = V_SV_s;
init.V_PA_s = V_PA_s;
init.V_PV_s = V_PV_s;
init.LAVmin = LAVmin;
init.RAVmin = RAVmin;

% Percentage of cardiac cycle 
params.tau_TS = k_TS; % unitless time to maximal systole.
params.tau_TR = k_TR; % unitless relaxation time (maximal systole to baseline)

% Heart period (s) 
params.HR = HR; 
params.T = 60/HR;

% Compliances (mL mmHg^(-1))
params.C_SA = C_SA;  params.C_SV = C_SV;  
params.C_PA = C_PA; params.C_PV = C_PV;  
    
% Resistances (mmHg s mL^(-1))
params.R_SA  = R_SA; 
params.R_tSA = 0.08; 
params.R_PA  = R_PA; 
params.R_tPA = 0.01; 
params.R_Veins = 0.040; 
params.R_SV = params.R_Veins; 
params.R_PV = params.R_Veins; 

params.R_t_o = R_t_o;
params.R_t_c = R_t_c;
params.R_p_o = R_p_o;
params.R_p_c = R_p_c;
params.R_m_o = R_m_o; 
params.R_m_c = R_m_c;
params.R_a_o = R_a_o;
params.R_a_c = R_a_c;
    
% Force scaling factors (unitless) 
% params.k_pas = k_pas; 
params.k_pas_LV = k_pas_LV;
params.k_pas_RV = k_pas_RV;
% params.k_act = k_act;  
params.k_act_LV = k_act_LV; 
params.k_act_RV = k_act_RV; 

% Midwall reference surface area (cm^2)
params.Amref_LV  = Amref_LV; 
params.Amref_SEP = Amref_SEP ;  
params.Amref_RV  = Amref_RV; 
    
% Free wall volume (mL) 
params.Vw_LV  = Vw_LV; 
params.Vw_SEP = Vw_SEP; 
params.Vw_RV  = Vw_RV;
params.LvSepR = LvSepR;

% Atrial parameters
params.RAV0u = 0.9 * RAVmin;
params.LAV0u = 0.9 * LAVmin; 
params.V0u_coeff = 0.9;
params.LEa = 2.60; %Atrial active contraction parameter
params.REa = 2.60; %Atrial active contraction parameter
params.V0c_coeff = 1.2; % 
params.LAV0c = LAVmax;
if(~isfield(inputData,'RAVmax'))
    params.RAV0c = LAVmax;
else
    params.RAV0c = RAVmax;
end
params.LAV1c = 5;
params.RAV1c = 5; 
params.LEp = 0.050;
params.REp = params.LEp;
params.Pc = 10;

% Pericardium parameters
params.Vh0 = LVEDV + RVEDV + LAVmin + RAVmin + Vw_LV + Vw_RV + Vw_SEP; % (mL)
params.K_P = 1;
params.B_P = 1;

%% Parameter adjustment
geom_pars = {'LvSepR'}; % LvSepR has already been adjusted before
geom_pars_i = find(contains(mods,geom_pars));
for i = 1:length(mods)
    if(isfield(params, mods{i}) && ~ismember(i,geom_pars_i))  % eliminate adjustable geometry parameters from mods since we already adjusted them above
        % As of 08/20/2024, all modifiers directly related to blood volume have been removed.
        % The functional part in this loop now only involves R_Veins and Vw_LV.
        if(strcmp(mods{i},'V0c_coeff'))
            params.V0c_coeff = modifiers(i);
            if(isfield(inputData, 'LAVmax'))
                params.LAV0c = modifiers(i) * LAVmax;
                params.RAV0c = params.LAV0c;
            end
            if(isfield(inputData,'RAVmax'))
                params.RAV0c = modifiers(i) * RAVmax;
            end
        elseif(strcmp(mods{i},'V0u_coeff')) 
            params.V0u_coeff = modifiers(i);
            params.LAV0u = modifiers(i) * params.LAV0u;
            params.RAV0u = modifiers(i) * params.RAV0u;
        elseif(strcmp(mods{i},'R_Veins'))
            params.R_PV = modifiers(i) * params.R_Veins;
            params.R_SV = modifiers(i) * params.R_Veins;
        elseif(strcmp(mods{i},'Vw_LV'))  %<-- assume LV and SEP have the same modifier, make sure the change in thickness propotional to TTE measurement
            params.Vw_SEP = modifiers(i) * params.Vw_SEP*targets.Hed_SW/targets.Hed_LW*params.Amref_SEP/params.Amref_LV*...
                (LvSepR/(1-LvSepR));
            params.Vw_LV = modifiers(i) * params.Vw_LV;
        else
            params.(mods{i}) = params.(mods{i})*modifiers(i);
        end
    elseif isfield(init, mods{i})
        assert(isfield(init, mods{i})); 
        init.(mods{i}) = init.(mods{i})*modifiers(i);
    end
end
end