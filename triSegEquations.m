function F = triSegEquations(geomVar,t, V_LV, V_RV, Lsc_LV, Lsc_SEP,Lsc_RV,params)
% Sarcomere length parameters (um)
Lsref   = 2;
Lsc0    = 1.51; 
%% Activation function 
T = 60/params.HR; 
% Ventricular activation function 
tc_v = mod(t,T) / T; % Fraction of cardiac cycle elapsed (t mod T)
if tc_v >= 0 & tc_v < params.tau_TS 
  Y = 0.5*(1 - cos(pi*tc_v/params.tau_TS)); 
elseif tc_v >= params.tau_TS & tc_v < params.tau_TR + params.tau_TS 
  Y = 0.5*(1 + cos(pi*(tc_v - params.tau_TS)/params.tau_TR)); 
else
  Y = 0; 
end
% Passive stress steepness parameter (dimensionless) 
gamma = 7.5; 

% geomVar = [xm_LV; xm_SEP; xm_RV; ym]
xm_LV = geomVar(1);
xm_SEP = geomVar(2);
xm_RV = geomVar(3);
ym = geomVar(4);

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


% Active stress . 
% Cell model is state variable. If geometry model bigger than cell model, positive stress. 
sigma_act_LV  = params.k_act_LV * Y  * (Lsc_LV/Lsc0  - 1) * (Ls_LV/Lsc0  - Lsc_LV/Lsc0) ;
sigma_act_SEP = params.k_act_LV * Y  * (Lsc_SEP/Lsc0 - 1) * (Ls_SEP/Lsc0 - Lsc_SEP/Lsc0);
sigma_act_RV  = params.k_act_RV * Y  * (Lsc_RV/Lsc0  - 1) * (Ls_RV/Lsc0  - Lsc_RV/Lsc0);


% Total stress 
sigma_LV  = sigma_act_LV  + sigma_pas_LV; 
sigma_SEP = sigma_act_SEP + sigma_pas_SEP; 
sigma_RV  = sigma_act_RV  + sigma_pas_RV; 

% Representative midwall tension 
Tm_LV  = (params.Vw_LV  * sigma_LV  / (2 * Am_LV))  * (1 + (z_LV^2)/3  + (z_LV^4)/5); 
Tm_SEP = (params.Vw_SEP * sigma_SEP / (2 * Am_SEP)) * (1 + (z_SEP^2)/3 + (z_SEP^4)/5); 
Tm_RV  = (params.Vw_RV  * sigma_RV  / (2 * Am_RV))  * (1 + (z_RV^2)/3  + (z_RV^4)/5);

% Force per unit


% Axial midwall tension component 
Tx_LV  = Tm_LV  * 2 * xm_LV  * ym / (xm_LV^2  + ym^2); 
Tx_SEP = Tm_SEP * 2 * xm_SEP * ym / (xm_SEP^2 + ym^2); 
Tx_RV  = Tm_RV  * 2 * xm_RV  * ym / (xm_RV^2  + ym^2); 

% Radial midwall tension component 
Ty_LV  = Tm_LV  * (xm_LV^2  - ym^2) / (xm_LV^2  + ym^2); 
Ty_SEP = Tm_SEP * (xm_SEP^2 - ym^2) / (xm_SEP^2 + ym^2); 
Ty_RV  = Tm_RV  * (xm_RV^2  - ym^2) / (xm_RV^2  + ym^2);

eq1 = -V_LV - 0.5*params.Vw_LV - 0.5*params.Vw_SEP + Vm_SEP - Vm_LV;
eq2 =  V_RV  + 0.5*params.Vw_RV + 0.5*params.Vw_SEP + Vm_SEP - Vm_RV;
eq3 = Tx_LV + Tx_SEP + Tx_RV;
eq4 = Ty_LV + Ty_SEP + Ty_RV;

if isreal(xm_LV)
    F = [eq1; eq2; eq3; eq4];
else
    F = 1e10*[1;1;1;1];
end
