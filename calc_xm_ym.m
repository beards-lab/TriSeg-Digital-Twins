function f = calc_xm_ym(x,Lsref,Vw,Amref,SL,V_LV,V_RV,fix_AmrefRV,dias)
%% Function Purpose:
% This function computes specific geometrical initial conditions of end-diastolic and
% end-systolic states, considering the balance of volume and sarcomere length via fsolve.

% Created by EB Randall, modified by Andrew Meyer, and Feng Gu
% Last modified: 10/29/2024

% Inputs:
%   Lsref      - Reference sarcomere length parameter
%   SL         - Sarcomere length
%   Vw_i       - Wall volume parameters
%   Amref_i    - Midwall surface area reference parameters
%   V_LV       - Left ventricular volume
%   V_RV       - Right ventricular volume
%   fix_AmrefRV - Use Amref_RV as an input and not as a state
%   dias       - 1 for diastolic state, other for systolic state

% Solutions:
%   x          - xm, ym

    %% Parameters 
    % Wall volumes (mL)
    Vw_LV  = Vw(1); 
    Vw_SEP = Vw(2); 
    Vw_RV  = Vw(3); 
    
    % Midwall reference surface areas (cm^2)
    Amref_LV  = Amref(1); 
    Amref_SEP = Amref(2); 
    
    % Use assigned Amref_RV value      
    if fix_AmrefRV                       
        Amref_RV = Amref(3);        
    end 
        
    
    %% States
    % Displacements (cm) 
    xm_LV    = x(1); 
    xm_SEP   = x(2); 
    xm_RV    = x(3);
    ym       = x(4);     
    % Use Amref_RV as a state
    if ~fix_AmrefRV
        Amref_RV = x(5); 
    end 
    
    %% Equations    
    % Midwall surface volume (mL)
    Vm_LV  = (pi / 6) * xm_LV  * (xm_LV^2  + 3 * ym^2); 
    Vm_SEP = (pi / 6) * xm_SEP * (xm_SEP^2 + 3 * ym^2); 
    Vm_RV  = (pi / 6) * xm_RV  * (xm_RV^2  + 3 * ym^2); 
  
    % Midwall surface area (cm^2)
    Am_LV  = pi * (xm_LV^2  + ym^2);
    Am_SEP = pi * (xm_SEP^2  + ym^2);
    Am_RV  = pi * (xm_RV^2  + ym^2);
    
    % Midwall curvature (cm^(-1))
    Cm_LV  = -2 * xm_LV  / (xm_LV^2  + ym^2);
    Cm_SEP = -2 * xm_SEP / (xm_SEP^2  + ym^2);
    Cm_RV  = -2 * xm_RV  / (xm_RV^2  + ym^2);
    
    % Midwall ratio (dimensionless) 
    z_LV   = 3 * Cm_LV  * Vw_LV  / (2 * Am_LV); 
    z_SEP  = 3 * Cm_SEP * Vw_SEP / (2 * Am_SEP); 
    z_RV   = 3 * Cm_RV  * Vw_RV  / (2 * Am_RV); 
    
    % Strain (dimensionless) 
    eps_LV  = 0.5 * log( Am_LV  / Amref_LV  ) - (1/12) * z_LV^2  - 0.019 * z_LV^4; 
    eps_SEP = 0.5 * log( Am_SEP / Amref_SEP ) - (1/12) * z_SEP^2 - 0.019 * z_SEP^4; 
    eps_RV  = 0.5 * log( Am_RV  / Amref_RV  ) - (1/12) * z_RV^2  - 0.019 * z_RV^4; 
    
    % Instantaneous sarcomere length (um) 
    Ls_LV  = Lsref * exp(eps_LV); 
    Ls_SEP = Lsref * exp(eps_SEP); 
    Ls_RV  = Lsref * exp(eps_RV); 
    
    % Balance sarcomere lengths in end-diastole
    if dias
        residuals = [Ls_LV  - SL;  % SL and Lsref are now both 2 µm
                     Ls_SEP - SL; 
                     Ls_RV  - SL; 
                     -V_LV - 0.5 * Vw_LV - 0.5 * Vw_SEP + Vm_SEP - Vm_LV;
                     V_RV + 0.5 * Vw_RV + 0.5 * Vw_SEP + Vm_SEP - Vm_RV];
    else
        residuals = [-V_LV - 0.5 * Vw_LV - 0.5 * Vw_SEP + Vm_SEP - Vm_LV; 
                      V_RV + 0.5 * Vw_RV + 0.5 * Vw_SEP + Vm_SEP - Vm_RV];
    end 
    
    f = sum(residuals.^2); 
    end 