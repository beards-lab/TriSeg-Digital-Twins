function [dxdt, outputs] = dXdT(t,x, params)     
%% Function Purpose:
% This function contains the differential equations that the ODE solver calls.

% Created by EB Randall, modified by Filip Jezek, Andrew Meyer, and Feng Gu
% Last modified: 03/20/2025

% 03/20: The biggest difference between the current model and the previous one is that I  
% changed the atrial model and pericardium model. The change to the pericardium model  
% is because I don't think we have enough information to identify it. The change to  
% the atrial model is because I don't think the previous version accurately reflected  
% atrial contraction.  


% Inputs: 
% t       - Time vector for simulation
% x       - Variables to be solved
% params  - Structure of parameter values

% Outputs: 
% dxdt    - Differential equations
% outputs - Simulation outputs

%% Variables 
% TriSeg geometries (cm)
xm_LV  = x(1); 
xm_SEP = x(2); 
xm_RV  = x(3);
ym     = x(4); 

% Contractile element length (um)
Lsc_LV  = x(5); 
Lsc_SEP = x(6); 
Lsc_RV  = x(7); 

% Volumes (mL) 
V_LV = x(8); % left ventricle
V_RV = x(9); % right ventricle
V_SA = x(10); % systemic arteries
V_SV = x(11); % systemic veins
V_PA = x(12); % pulmonary arteries
V_PV = x(13); % pulmonary veins
V_LA = x(14); % left atrium
V_RA = x(15); % right atrium

%% Activation function 
T = 60/params.HR; 
% Ventricular activation function 
tc_v = mod(t,T) / T; % Fraction of cardiac cycle elapsed (t mod T)
if tc_v >= 0 && tc_v < params.tau_TS 
  Y = 0.5*(1 - cos(pi*tc_v/params.tau_TS)); 
elseif tc_v >= params.tau_TS && tc_v < params.tau_TR + params.tau_TS 
  Y = 0.5*(1 + cos(pi*(tc_v - params.tau_TS)/params.tau_TR)); 
else
  Y = 0; 
end

%% Heart model
% Sarcomere length parameters (um)
Lsref   = 2;
Lsc0    = 1.51; 
Lse_iso = 0.04;
    
% Sarcomere length shortening velocity (um s^(-1))
v_max = 7; % 7 / 2
% Passive stress steepness parameter (dimensionless) 
gamma = 7.5; 

% Volume of spherical cap formed by midwall surface (mL)
Vm_LV  = (pi / 6) * xm_LV  * (xm_LV^2  + 3 * ym^2); 
Vm_SEP = (pi / 6) * xm_SEP * (xm_SEP^2 + 3 * ym^2); 
Vm_RV  = (pi / 6) * xm_RV  * (xm_RV^2  + 3 * ym^2); 

% Surface area of midwall surface (cm^2) 
Am_LV  = pi * (xm_LV^2  + ym^2);
Am_SEP = pi * (xm_SEP^2 + ym^2); 
Am_RV  = pi * (xm_RV^2  + ym^2); 

% Curvature of midwall surface (cm^(-1))
Cm_LV  = 2 * xm_LV  / (xm_LV^2  + ym^2);
Cm_SEP = 2 * xm_SEP / (xm_SEP^2 + ym^2);
Cm_RV  = 2 * xm_RV  / (xm_RV^2  + ym^2);

% Ratio of wall thickness to midwall radius of curvature (dimensionless)
z_LV  = 3 * Cm_LV  * params.Vw_LV  / (2 * Am_LV); 
z_SEP = 3 * Cm_SEP * params.Vw_SEP / (2 * Am_SEP); 
z_RV  = 3 * Cm_RV  * params.Vw_RV  / (2 * Am_RV);

% Wall thickness (cm)
theta_LV = asin(-ym*Cm_LV);
theta_SEP = asin(ym*abs(Cm_SEP)); % if septum curves into LV, radius will be negative, which will mess up calculations
theta_RV = asin(ym*Cm_RV); 
if(-1/Cm_LV > -xm_LV) % if the midwall radius is greater than xm_LV, phi bound is 0 to theta_LV
    r_LV = ((-Cm_LV)^(-3) - 3*params.Vw_LV / (4*pi*(1-cos(theta_LV))))^(1/3);
    H_LW = ((3*params.Vw_LV) / (2*pi*(1-cos(theta_LV))) + r_LV^3)^(1/3) - r_LV;
else % if the midwall radius is less than xm_LV (almost always), phi bound is 0 to pi - theta_LV
    r_LV = ((-Cm_LV)^(-3) - 3*params.Vw_LV / (4*pi*(1+cos(theta_LV))))^(1/3);
    H_LW = ((3*params.Vw_LV) / (2*pi*(1+cos(theta_LV))) + r_LV^3)^(1/3) - r_LV;
end
if(abs(1/Cm_SEP) > abs(xm_SEP)) % if the midwall radius is greater than xm_SEP, phi bound is 0 to theta_SEP (almost always)
    r_SEP = ((abs(Cm_SEP))^(-3) - 3*params.Vw_SEP / (4*pi*(1-cos(theta_SEP))))^(1/3);
    H_SW = ((3*params.Vw_SEP) / (2*pi*(1-cos(theta_SEP))) + r_SEP^3)^(1/3) - r_SEP;
else % if the midwall radius is less than xm_SEP, phi bound is 0 to pi - theta_SEP
    r_SEP = ((abs(Cm_SEP))^(-3) - 3*params.Vw_SEP / (4*pi*(1+cos(theta_SEP))))^(1/3);
    H_SW = ((3*params.Vw_SEP) / (2*pi*(1+cos(theta_SEP))) + r_SEP^3)^(1/3) - r_SEP;
end
if(1/Cm_RV > xm_RV) % if the midwall radius is greater than xm_RV, phi bound is 0 to theta_RV
    r_RV = ((Cm_RV)^(-3) - 3*params.Vw_RV / (4*pi*(1-cos(theta_RV))))^(1/3); 
    H_RW = ((3*params.Vw_RV) / (2*pi*(1-cos(theta_RV))) + r_RV^3)^(1/3) - r_RV;
else % if the midwall radius is less than xm_RV, phi bound is 0 to pi - theta_RV
    r_RV = ((Cm_RV)^(-3) - 3*params.Vw_RV / (4*pi*(1+cos(theta_RV))))^(1/3); 
    H_RW = ((3*params.Vw_RV) / (2*pi*(1+cos(theta_RV))) + r_RV^3)^(1/3) - r_RV;
end

% Myofiber strain (dimensionless)
eps_LV  = 0.5 * log( Am_LV  / params.Amref_LV  ) - (1/12) * z_LV^2  - 0.019 * z_LV^4; 
eps_SEP = 0.5 * log( Am_SEP / params.Amref_SEP ) - (1/12) * z_SEP^2 - 0.019 * z_SEP^4; 
eps_RV  = 0.5 * log( Am_RV  / params.Amref_RV  ) - (1/12) * z_RV^2  - 0.019 * z_RV^4; 

% Sarcomere length (um)
Ls_LV  = Lsref * exp(eps_LV); % From geometry of whole model
Ls_SEP = Lsref * exp(eps_SEP); 
Ls_RV  = Lsref * exp(eps_RV); 

% Passive stress  
sigma_pas_LV  =  params.k_pas_LV * (Ls_LV/Lsc0 - 1)^gamma;
sigma_pas_SEP =  params.k_pas_LV * (Ls_SEP/Lsc0 - 1)^gamma;
sigma_pas_RV  =  params.k_pas_RV * (Ls_RV/Lsc0 - 1)^gamma;
% sigma_pas_LV  =  params.k_pas * (Ls_LV/Lsc0 - 1)^gamma;
% sigma_pas_SEP =  params.k_pas * (Ls_SEP/Lsc0 - 1)^gamma;
% sigma_pas_RV  =  params.k_pas * (Ls_RV/Lsc0 - 1)^gamma;

% Active stress . 
% Cell model is state variable. If geometry model bigger than cell model, positive stress. 
sigma_act_LV  = params.k_act_LV * Y  * (Lsc_LV/Lsc0  - 1) * (Ls_LV/Lsc0  - Lsc_LV/Lsc0) ;
sigma_act_SEP = params.k_act_LV * Y  * (Lsc_SEP/Lsc0 - 1) * (Ls_SEP/Lsc0 - Lsc_SEP/Lsc0);
sigma_act_RV  = params.k_act_RV * Y  * (Lsc_RV/Lsc0  - 1) * (Ls_RV/Lsc0  - Lsc_RV/Lsc0);
% sigma_act_LV  = params.k_act * Y  * (Lsc_LV/Lsc0  - 1) * (Ls_LV/Lsc0  - Lsc_LV/Lsc0) ;
% sigma_act_SEP = params.k_act * Y  * (Lsc_SEP/Lsc0 - 1) * (Ls_SEP/Lsc0 - Lsc_SEP/Lsc0);
% sigma_act_RV  = params.k_act * Y  * (Lsc_RV/Lsc0  - 1) * (Ls_RV/Lsc0  - Lsc_RV/Lsc0);

% Total stress 
sigma_LV  = sigma_act_LV  + sigma_pas_LV; 
sigma_SEP = sigma_act_SEP + sigma_pas_SEP; 
sigma_RV  = sigma_act_RV  + sigma_pas_RV; 

% Representative midwall tension 
Tm_LV  = (params.Vw_LV  * sigma_LV  / (2 * Am_LV))  * (1 + (z_LV^2)/3  + (z_LV^4)/5); 
Tm_SEP = (params.Vw_SEP * sigma_SEP / (2 * Am_SEP)) * (1 + (z_SEP^2)/3 + (z_SEP^4)/5); 
Tm_RV  = (params.Vw_RV  * sigma_RV  / (2 * Am_RV))  * (1 + (z_RV^2)/3  + (z_RV^4)/5);

% Axial midwall tension component 
Tx_LV  = Tm_LV  * 2 * xm_LV  * ym / (xm_LV^2  + ym^2); 
Tx_SEP = Tm_SEP * 2 * xm_SEP * ym / (xm_SEP^2 + ym^2); 
Tx_RV  = Tm_RV  * 2 * xm_RV  * ym / (xm_RV^2  + ym^2); 

% Radial midwall tension component 
Ty_LV  = Tm_LV  * (xm_LV^2  - ym^2) / (xm_LV^2  + ym^2); 
Ty_SEP = Tm_SEP * (xm_SEP^2 - ym^2) / (xm_SEP^2 + ym^2); 
Ty_RV  = Tm_RV  * (xm_RV^2  - ym^2) / (xm_RV^2  + ym^2);

P_Peri = params.K1 * exp(params.expPeri*((V_LV+V_RV+V_LA+V_RA)/params.Vh0-1));% 0.4 is coming from an experiment fit but just 1 single point

% Ventricular pressure 
P_LV = -2 * Tx_LV / ym + P_Peri; 
P_RV = 2 * Tx_RV / ym + P_Peri; 

% Atria
LAV0u = params.LAV0u; % mL. Unstressed volume in atria
RAV0u = params.RAV0u; % mL. Unstressed volume in atria
LAV0c = params.LAV0c; % ml. Threshold for exponential PV relation. put a modifier on; drives LAVmax. make dif
RAV0c = params.RAV0c;
LAV1c = params.LAV1c; % mL. Volume constant dictating behavior of exponential. mkae dif
RAV1c = params.RAV1c;
Tact = -0.15; %.1
tc_a = tc_v - Tact - 1*((tc_v-Tact)>(0.5));
LEp = params.LEp ; % 0.050  passive <try this instead. make dif
REp = params.REp;
LEa = params.LEa; %active. 
REa = params.REa;
Pc = params.Pc; % 10 mmHg. Collagen
% 03/31 I don't know why the current Atria model doesn't work. I will just back
% to Dan's version
if tc_a >= 0 && tc_a < 0.15 
  act = 0.5*(1 - cos(pi*tc_a/0.15)); 
elseif tc_a >= 0.15 && tc_a < 0.3 
  act = 0.5*(1 + cos(pi*(tc_a - 0.15)/0.15)); 
else
  act = 0; 
end
% sigma_a = .0975; %1.5 * 0.065. How wide gaussian is
% act = exp(-(tc_a./sigma_a).^2 );
P_LA  = 2.0*(LEp*(V_LA-LAV0u) + LEa*act*(V_LA-LAV0u) + Pc*exp((V_LA-LAV0c)/LAV1c)) + P_Peri;
P_RA  = 1.0*(REp*(V_RA-RAV0u) + REa*act*(V_RA-RAV0u) + Pc*exp((V_RA-RAV0c)/RAV1c)) + P_Peri; 

%% Lumped circulatory model 
% Venous Pressure (mmHg)
P_SV = V_SV / params.C_SV; 
P_PV = V_PV / params.C_PV; 

% Atrial flows (mL / s)
QIN_LA = (P_PV - P_LA)/params.R_PV;
QIN_RA = (P_SV - P_RA)/params.R_SV;

% Atrioventricular flows (mL / s)
% Mitral
if(P_LA >= P_LV)
    QIN_LV = (P_LA - P_LV) / params.R_m_o;
else
    QIN_LV = (P_LA - P_LV) / params.R_m_c;
end

% Tricuspid
if(P_RA >= P_RV)
    QIN_RV = (P_RA - P_RV) / params.R_t_o;
else
    QIN_RV = (P_RA - P_RV) / params.R_t_c;
end 

% Aortic valve closed completely
QOUT_LV = 0;
Q_SA = (V_SA - params.C_SA*P_SV)/(params.C_SA*(params.R_SA + params.R_tSA));
P_SA = (params.R_SA*V_SA + params.C_SA*P_SV*params.R_tSA)/(params.C_SA*(params.R_SA + params.R_tSA));  % Pressure equal flow * resistence
% Aortic valve closed with regurgitation
if(P_LV >= P_SA) 
    % Forward Flow
    QOUT_LV = -(params.R_SA*V_SA - params.C_SA*P_LV*params.R_SA - params.C_SA*P_LV*params.R_tSA + params.C_SA*P_SV*params.R_tSA)/(params.C_SA*(params.R_SA*params.R_a_o + params.R_SA*params.R_tSA + params.R_a_o*params.R_tSA));
    Q_SA = (params.R_a_o*V_SA - params.C_SA*P_SV*params.R_a_o + params.C_SA*P_LV*params.R_tSA - params.C_SA*P_SV*params.R_tSA)/(params.C_SA*(params.R_SA*params.R_a_o + params.R_SA*params.R_tSA + params.R_a_o*params.R_tSA));
    P_SA = (params.R_SA*params.R_a_o*V_SA + params.C_SA*P_LV*params.R_SA*params.R_tSA + params.C_SA*P_SV*params.R_a_o*params.R_tSA)/(params.C_SA*(params.R_SA*params.R_a_o + params.R_SA*params.R_tSA + params.R_a_o*params.R_tSA));
elseif params.R_a_c < Inf
    % Backward Flow with regurgitation
    QOUT_LV = -(params.R_SA*V_SA - params.C_SA*P_LV*params.R_SA - params.C_SA*P_LV*params.R_tSA + params.C_SA*P_SV*params.R_tSA)/(params.C_SA*(params.R_SA*params.R_a_c + params.R_SA*params.R_tSA + params.R_a_c*params.R_tSA));
    Q_SA = (params.R_a_c*V_SA - params.C_SA*P_SV*params.R_a_c + params.C_SA*P_LV*params.R_tSA - params.C_SA*P_SV*params.R_tSA)/(params.C_SA*(params.R_SA*params.R_a_c + params.R_SA*params.R_tSA + params.R_a_c*params.R_tSA));
    P_SA = (params.R_SA*params.R_a_c*V_SA + params.C_SA*P_LV*params.R_SA*params.R_tSA + params.C_SA*P_SV*params.R_a_c*params.R_tSA)/(params.C_SA*(params.R_SA*params.R_a_c + params.R_SA*params.R_tSA + params.R_a_c*params.R_tSA));
end

% Pulmonary valve closed completely
QOUT_RV = 0; 
Q_PA = (V_PA - params.C_PA*P_PV)/(params.C_PA*(params.R_PA + params.R_tPA)); 
P_PA = (params.R_PA*V_PA + params.C_PA*P_PV*params.R_tPA)/(params.C_PA*(params.R_PA + params.R_tPA)); 
% Pulmonary valve closed with regurgitation
if (P_RV >= P_PA)
    % Forward Flow
    QOUT_RV = -(params.R_PA*V_PA - params.C_PA*P_RV*params.R_PA + params.C_PA*P_PV*params.R_tPA - params.C_PA*P_RV*params.R_tPA)/(params.C_PA*(params.R_PA*params.R_p_o + params.R_PA*params.R_tPA + params.R_p_o*params.R_tPA)); 
    Q_PA = (params.R_p_o*V_PA - params.C_PA*P_PV*params.R_p_o - params.C_PA*P_PV*params.R_tPA + params.C_PA*P_RV*params.R_tPA)/(params.C_PA*(params.R_PA*params.R_p_o + params.R_PA*params.R_tPA + params.R_p_o*params.R_tPA)); 
    P_PA = (params.R_PA*params.R_p_o*V_PA + params.C_PA*P_RV*params.R_PA*params.R_tPA + params.C_PA*P_PV*params.R_p_o*params.R_tPA)/(params.C_PA*(params.R_PA*params.R_p_o + params.R_PA*params.R_tPA + params.R_p_o*params.R_tPA));
elseif params.R_p_c < Inf
    % Backward Flow
    QOUT_RV = -(params.R_PA*V_PA - params.C_PA*P_RV*params.R_PA + params.C_PA*P_PV*params.R_tPA - params.C_PA*P_RV*params.R_tPA)/(params.C_PA*(params.R_PA*params.R_p_c + params.R_PA*params.R_tPA + params.R_p_c*params.R_tPA)); 
    Q_PA = (params.R_p_c*V_PA - params.C_PA*P_PV*params.R_p_c - params.C_PA*P_PV*params.R_tPA + params.C_PA*P_RV*params.R_tPA)/(params.C_PA*(params.R_PA*params.R_p_c + params.R_PA*params.R_tPA + params.R_p_c*params.R_tPA)); 
    P_PA = (params.R_PA*params.R_p_c*V_PA + params.C_PA*P_RV*params.R_PA*params.R_tPA + params.C_PA*P_PV*params.R_p_c*params.R_tPA)/(params.C_PA*(params.R_PA*params.R_p_c + params.R_PA*params.R_tPA + params.R_p_c*params.R_tPA));
end 

% Equations 1 - 4 (AE's)
eq1  = -V_LV - 0.5 * params.Vw_LV - 0.5 * params.Vw_SEP + Vm_SEP - Vm_LV; 
eq2 = Tx_LV + Tx_SEP + Tx_RV;
eq3  = V_RV + 0.5 * params.Vw_RV + 0.5 * params.Vw_SEP + Vm_SEP - Vm_RV;
eq4  = Ty_LV + Ty_SEP + Ty_RV; 

% Equations 5 - 7 (ODE's from TriSeg)
dLsc_LV  = ((Ls_LV  - Lsc_LV)  / Lse_iso - 1) * v_max;
dLsc_SEP = ((Ls_SEP - Lsc_SEP) / Lse_iso - 1) * v_max;
dLsc_RV  = ((Ls_RV  - Lsc_RV)  / Lse_iso - 1) * v_max;

% Equations 8 - 15 (ODE's from circulatory model)
dV_LV = QIN_LV - QOUT_LV; 
dV_SA = QOUT_LV - Q_SA; 
dV_SV = Q_SA - QIN_RA; %
dV_RV = QIN_RV - QOUT_RV; 
dV_PA = QOUT_RV - Q_PA; 
dV_PV = Q_PA - QIN_LA; %
dV_RA = QIN_RA - QIN_RV; % V_RA
dV_LA = QIN_LA - QIN_LV; % V_LA

dxdt = [eq1; eq2; eq3; eq4; dLsc_LV; dLsc_SEP; dLsc_RV; dV_LV; dV_RV; dV_SA; dV_SV; dV_PA; dV_PV; dV_LA; dV_RA]; 

outputs = [P_LV; P_SA; P_SV; P_RV; P_PA; P_PV;          % 1-6
        Vm_LV; Vm_SEP; Vm_RV;                           % 7-9
        Am_LV; Am_SEP; Am_RV;                           % 10-12
        Cm_LV; Cm_SEP; Cm_RV;                           % 13-15
        eps_LV; eps_SEP; eps_RV;                        % 16-18
        sigma_pas_LV; sigma_pas_SEP; sigma_pas_RV;      % 19-21
        sigma_act_LV; sigma_act_SEP; sigma_act_RV;      % 22-24
        sigma_LV; sigma_SEP; sigma_RV;                  % 25-27
        QIN_LV; QOUT_LV; QIN_RV; QOUT_RV;               % 28-31
        Q_SA; Q_PA;                                     % 32-33
        Tm_LV; Tm_SEP; Tm_RV;                           % 34-36
        Y;                                             % 37
        V_RA; V_LA; P_RA; P_LA; QIN_RA;                 % 38-42                                           
        H_LW; H_SW; H_RW;act;                           % 43-46
        r_LV;r_SEP;r_RV;                                 %47-49  
        P_Peri];                                     
            
end 