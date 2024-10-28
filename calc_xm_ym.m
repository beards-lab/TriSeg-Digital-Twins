function f = calc_xm_ym(x,Lsref,Vw,Amref,SL,V_LV,V_RV,fix_AmrefRV,dias)

    %{ 

x,Lsref, Vw, Amref, SL_d, LVEDV, RVEDV,use_AmrefRV


    This function is used in the calculation of the initial estimates for
    the xm and ym values.  

    Inputs: 
    x       - vector of states 
    inputs - container with 
        Lsref   - reference sarcomere length parameter 
        SL      - sarcomere length 
        Vw_i    - wall volume parameters 
        Amref_i - midwall surface area reference parameters 
        
        V_LV    - left ventricular volume 
        V_RV    - right ventricular volume 
    fix_AmrefRV - Use Amref_RV as an input and not as a state

    %}

    %% Parameters 
    % Wall volumes (mL)
    Vw_LV  = Vw(1); 
    Vw_SEP = Vw(2); 
    Vw_RV  = Vw(3); 
    
    % Midwall reference surface areas (cm^2)
    Amref_LV  = Amref(1); 
    Amref_SEP = Amref(2); 
    
    % use assigned Amref_RV value      
    if fix_AmrefRV                       
        Amref_RV = Amref(3);        
    end 
        
    
    %% States
    
%     x = exp(x); 
    
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

    % NOTHING BELOW HERE EXECUTES FOR SYSTOLE, AS WE AREN'T MATCHING
    % SARCOMERE LENGTHS
    
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
    
    %% Outputs    
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
    
    % For fmincon, we need to return a scalar value that represents the
    % objective. A common choice is the sum of squared residuals.
    f = sum(residuals.^2); 
    end 