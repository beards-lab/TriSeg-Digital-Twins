function [params, init] = optParams(params,mods,modifiers)
%% Function Purpose:
% This function is desiged to optimze the paramters based on initial guesses for all model parameters. 
% Inputs to this function are the outputs from the estiminiParams.m function. 
% Outputs of this function include a structure of parameters used in the model and initial
% conditions for the ODE solver (dXdT.m function).

% Created by Feng Gu
% Last modified: 03/20/2025

% Inputs: 
% params      - Structure containing the initial guessed parameters for the model  
% init        - Initial guess for the initial conditions used in the ode15s solver (dXdT.m function)  
% mods        - Cell array of names to adjust selected parameters 
% modifiers   - Vector of floats that adjust selected parameters

% Outputs: 
% params      - Structure of parameters used in the model 
% init        - Initial conditions for the ode15s solver (dXdT.m function) 

% Related functions: 
% estiminiParams.m - Serve as the input of current function
% geom_0      - Computes the initial guess (idealized end-diastolic state) of TriSeg geometry 
% calc_xm_ym  - Computes specific geometrical initial conditions of end-diastolic and end-systolic
%               states, considering the balance of volume and sarcomere length

%% Optimize geometrical parameters first  

optVw_LV = params.Vw_LV*modifiers(contains(mods,'Vw_LV'));
optVw_RV = params.Vw_RV*modifiers(contains(mods,'Vw_RV'));
optAmref_LV = params.Amref_LV*modifiers(contains(mods,'Amref_LV'));
optAmref_RV = params.Amref_RV*modifiers(contains(mods,'Amref_RV'));
optLvSepR = params.LvSepR*modifiers(contains(mods,'LvSepR'));

% Transfer LvSepR
optVw_SEP = optVw_LV /optLvSepR - optVw_LV ;
if optVw_SEP <= 0 
    error("unreal geometery")
end
beta = acos(2 * optLvSepR - 1);
if ~isreal(beta)
    error("unreal geometery")
end
optAmref_SEP = tan(beta/2)*optAmref_LV ;
% Calculate LVEDV
optLVEDV = 1/6*sqrt((optAmref_LV+optAmref_SEP)^3/pi) - 0.5*(optVw_LV+optVw_SEP);
% Calculate Xm Ym for every optimzation
syms y beta_sym A1 A2 pi_sym

eq = y^4 - (sin(beta_sym)^2 * (A1 + A2)/(2 * pi_sym * (1 + sin(beta_sym)^2))) * y^2 ...
     - ((A2 - A1)^2 * sin(beta_sym)^4/(16 * pi_sym^2 * (1 + sin(beta_sym)^2))) == 0;

y_solutions = solve(eq, y,"ReturnConditions",true);
beta_sym = beta;
A1 = optAmref_LV;
A2 = optAmref_SEP;
pi_sym = pi;
y_solutions_numeric = subs(y_solutions.y);
optym_d = double(y_solutions_numeric(arrayfun(@(s) isreal(s) && s > 0, y_solutions_numeric)));
optxm_LV = -optym_d / sin(beta) + (A2 - A1) * sin(beta) / (4 * optym_d * pi);
optxm_SEP = optym_d / sin(beta) + (A2 - A1) * sin(beta) / (4 * optym_d * pi);
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

init.xm_LV_d = d0(1);
init.xm_SEP_d = d0(2);
init.xm_RV_d = d0(3);
init.ym_d = d0(4) ;
init.LVEDV = optLVEDV;
init.RVEDV = optRVEDV;

params.Amref_LV = optAmref_LV;
params.Amref_SEP = optAmref_SEP;
params.Amref_RV = optAmref_RV;

params.Vw_LV = optVw_LV;
params.Vw_SEP = optVw_SEP;
params.Vw_RV = optVw_RV;
params.LvSepR = optLvSepR;
%% Optimize the remaining parameters  

geom_pars = {'Vw_LV','Vw_RV','Amref_LV','Amref_RV','LvSepR'};
geom_pars_i = find(contains(mods,geom_pars));
for i = 1:length(mods)
    if(isfield(params, mods{i}) && ~ismember(i,geom_pars_i))
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
        else
            params.(mods{i}) = params.(mods{i})*modifiers(i);
        end
    end
end

end