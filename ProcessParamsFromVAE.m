function [params, init] = ProcessParamsFromVAE(ParamsFormVAE)
%% Function Purpose:
% This function is desiged to post-process fake params from VAE.

% Created by Feng Gu
% Last modified: 04/02/2025

% Inputs: 
% ParamsFormVAE      - Structure faked parameters for the VAE  

% Outputs: 
% params      - Structure of parameters used in the model 
% init        - Initial conditions for the ode15s solver (dXdT.m function) 

% Related functions: 
% geom_0      - Computes the initial guess (idealized end-diastolic state) of TriSeg geometry 
% calc_xm_ym  - Computes specific geometrical initial conditions of end-diastolic and end-systolic
%               states, considering the balance of volume and sarcomere length

%% Load geometrical parameters first  

optVw_LV = ParamsFormVAE.Vw_LV;
optVw_RV = ParamsFormVAE.Vw_RV;
optAmref_LV = ParamsFormVAE.Amref_LV;
optAmref_RV = ParamsFormVAE.Amref_RV;
optLvSepR = ParamsFormVAE.Vw_LV/(ParamsFormVAE.Vw_LV+ParamsFormVAE.Vw_SEP);
optVw_SEP = ParamsFormVAE.Vw_SEP ;
optAmref_SEP = ParamsFormVAE.Amref_SEP ;
% Calculate LV side
optLVEDV = 1/6*sqrt((optAmref_LV + optAmref_SEP)^3/pi) - 0.5*(optVw_LV + optVw_SEP);
if optVw_SEP <= 0 
    error("unreal geometery")
end
beta = acos(2 * optLvSepR - 1);
if ~isreal(beta)
    error("unreal geometery")
end
% Calculate Xm Ym for every optimzation
optym_d = (optAmref_LV*sin(beta)^2/(2*pi*(1+cos(beta))))^(1/2);
optxm_LV = -optym_d / sin(beta) + (optAmref_SEP - optAmref_LV) * sin(beta) / (4 * optym_d * pi);
optxm_SEP = optym_d / sin(beta) + (optAmref_SEP - optAmref_LV) * sin(beta) / (4 * optym_d * pi);
optxm_RV = sqrt(optAmref_RV/pi-optym_d^2);
optRVEDV = (pi / 6) * optxm_RV * (optxm_RV^2 + 3 * optym_d^2)-0.5*(optVw_RV+optVw_SEP)- ((pi / 6) * optxm_SEP * (optxm_SEP^2 + 3 * optym_d^2));
if optRVEDV < 0 || optLVEDV < 0||optLvSepR >= 1 || optLvSepR < 0
    error("unreal geometery")
end
% Re-balance the sarcomere length
dias = 1;

x0_d = [optxm_LV; 
    optxm_SEP; 
    optxm_RV;
    optym_d;
    optAmref_RV]; 
fix_AmrefRV = 1; 

Vw    = [optVw_LV,optVw_SEP,optVw_RV]; 
Amref = [optAmref_LV,optAmref_SEP,optAmref_RV]; 

% Assume end-diastolic sarcomere length 
SL_d    = 2; % µm 
opts = optimoptions('fsolve','Display','none',...
    'MaxFunctionEvaluations',2e3,'Algorithm','levenberg-marquardt'); 
Lsref   = 2;
[d0,~] = fsolve(@(x) calc_xm_ym(x,Lsref, Vw, Amref,SL_d,optLVEDV, optRVEDV,fix_AmrefRV,dias),x0_d,opts); % end-diastole. Always use RVEDV

% Total blood volume (mL) 

if(ParamsFormVAE.Sex == 1) % Male
    Vtot  = (0.3669*(ParamsFormVAE.Height/100)^3 +0.03219*ParamsFormVAE.Weight+0.6041)*1000; % mL (Nadler's Equation Male)
else                      % Female
    Vtot  = (0.3561*(ParamsFormVAE.Height/100)^3 +0.03308*ParamsFormVAE.Weight+0.1833)*1000; % mL (Nadler's Equation Female)
end

% Snapped at end diastole - maximal ventricles, minimal atria
% Blood volume distribution values; sum total = 1.0 
d_SA = .15;              d_PA = .10; 
d_SV = .65;              d_PV = .10;

LAVmin = ParamsFormVAE.LAV0u/0.9; 
if ParamsFormVAE.Sex == 1
   RAVmin = 32/25*LAVmin; 
else % female
   RAVmin = 23/22*LAVmin;
end

Vd = Vtot - optLVEDV - optRVEDV - RAVmin - LAVmin; % distributed volume (available for other compartments)

% Total compartment volumes 
V_SA_0 = d_SA*Vd;      V_PA_0 = d_PA*Vd; 
V_SV_0 = d_SV*Vd;      V_PV_0 = d_PV*Vd;

% Unstressed compartment volumes
V_SA_u = V_SA_0*0.7;   V_PA_u = V_PA_0*0.1; 
V_SV_u = V_SV_0*0.9;   V_PV_u = V_PV_0*0.9; 

% Stressed compartment volumes
V_SA_s = V_SA_0 - V_SA_u;  V_PA_s = V_PA_0 - V_PA_u;  
V_SV_s = V_SV_0 - V_SV_u;  V_PV_s = V_PV_0 - V_PV_u; 

init.xm_LV_d = d0(1);
init.xm_SEP_d = d0(2);
init.xm_RV_d = d0(3);
init.ym_d = d0(4) ;
init.Lsc_LV_0 = 2;
init.Lsc_SEP_0 = 2;
init.Lsc_RV_0 = 2;
init.LVEDV = optLVEDV;
init.RVEDV = optRVEDV;
init.V_SA_s = V_SA_s;
init.V_SV_s = V_SV_s;
init.V_PA_s = V_PA_s;
init.V_PV_s = V_PV_s;
init.LAVmin = LAVmin;
init.RAVmin = RAVmin;


params.Amref_LV = optAmref_LV;
params.Amref_SEP = optAmref_SEP;
params.Amref_RV = optAmref_RV;

params.Vw_LV = optVw_LV;
params.Vw_SEP = optVw_SEP;
params.Vw_RV = optVw_RV;
params.LvSepR = optLvSepR;
params.Vh0 = optLVEDV + optRVEDV + optVw_LV + optVw_SEP + optVw_RV + LAVmin + RAVmin;
if ParamsFormVAE.Sex == 1
    QS2 = 545.17-2.117*ParamsFormVAE.HR;
else
    QS2 = 546.5-2.0*ParamsFormVAE.HR;
end
IVRT = 70*75/ParamsFormVAE.HR;
ActT = QS2+IVRT;
k_TS = (ActT*ParamsFormVAE.HR/6e4)*2/3; % Beginning of cardiac cycle to maximal systole  
k_TR = (ActT*ParamsFormVAE.HR/6e4)*1/3; % Relaxation time fraction 
params.tau_TR = k_TR;
params.tau_TS = k_TS;
params.HR = ParamsFormVAE.HR;
params.T = 60/ParamsFormVAE.HR;
params.C_PA = ParamsFormVAE.C_PA;
params.C_PV = ParamsFormVAE.C_PV;
params.C_SA = ParamsFormVAE.C_SA;
params.C_SV = ParamsFormVAE.C_SV;
params.R_PA = ParamsFormVAE.R_PA;
params.R_tPA = ParamsFormVAE.R_tPA;
params.R_SA = ParamsFormVAE.R_SA;
params.R_tSA = ParamsFormVAE.R_tSA;
% params.R_Veins = 0.04;
params.R_SV = ParamsFormVAE.R_SV;
params.R_PV = ParamsFormVAE.R_SV;
params.R_a_c = ParamsFormVAE.R_a_c;
params.R_a_o = ParamsFormVAE.R_a_o;
params.R_m_c = ParamsFormVAE.R_m_c;
params.R_m_o = ParamsFormVAE.R_m_o;
params.R_p_c = ParamsFormVAE.R_p_c;
params.R_p_o = ParamsFormVAE.R_p_o;
params.R_t_c = ParamsFormVAE.R_t_c;
params.R_t_o = 0.0030;
params.k_act_LV = ParamsFormVAE.k_act_LV;
params.k_act_RV = ParamsFormVAE.k_act_RV;
params.k_pas_LV = ParamsFormVAE.k_pas_LV;
params.k_pas_RV = ParamsFormVAE.k_pas_RV;
params.LAV0c = ParamsFormVAE.LAV0c;
params.RAV0c = ParamsFormVAE.LAV0c; 

params.LAV0u = ParamsFormVAE.LAV0u; % need extra Gen
if ParamsFormVAE.Sex == 1
    params.RAV0u = 32/25*params.LAV0u;
else
    params.RAV0u = 23/22*params.LAV0u;
end
params.LAV1c = 5;
params.RAV1c = 5; 
params.LEa = 0.60; %Atrial active contraction parameter
params.REa = 0.60; %Atrial active contraction parameter
params.LEp = 0.050;
params.REp = params.LEp;
params.Pc = 10;
params.V0u_coeff =  0.9000;
params.V0c_coeff = 1.2;
params.K1 = ParamsFormVAE.K1;
params.expPeri = ParamsFormVAE.expPeri;

end