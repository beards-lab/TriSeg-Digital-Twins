function total_cost = evaluateModel(modifiers,data,PATIENT_NO,ModelWin,MRI_flag)
%% Function Purpose:
% This function is used to compute a cost function in runSim.m with new modifiers for patients
% during optimization.

% Created by Andrew Meyer and Feng Gu
% Last modified: 10/29/2024

try
    % Evaluates the model during optimalization
    MRI_flag = MRI_flag;
    [~,targets, inputs, mods] = targetVals_HF(data,PATIENT_NO,ModelWin,MRI_flag);
    [INIparams, INIinit] = estiminiParams(targets,inputs);
    [params, init] = optParams(INIparams, INIinit, mods,modifiers);
    if params.K1>=40
        error('impossibile pericardium')
    end
    if params.k_act_LV/params.k_act_RV  > 10
        error('wired way balancing')
    end
    params_no = struct2array(params);
    init_no = struct2array(init);
    if any(params_no <= 0)||any(init_no(2:end) <= 0)
        error('bad inital guessing')
    end
    runSimOn;           % 优先用GL
    o_vals_no = struct2array(o_vals);
    if any(o_vals_no <= -1)
        error('bad simulation')
    end
catch ME
    total_cost = inf;
    disp(['Error: ', ME.message])
end
