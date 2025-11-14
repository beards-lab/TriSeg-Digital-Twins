function total_cost = evaluateModelUmich(modifiers,data,PATIENT_NO,ModelWin,MRI_flag)
%% Function Purpose:
% This function is used to compute a cost function in runSim.m with new
% modifiers for patients from Umich cohort during optimization.

% Created by Andrew Meyer and Feng Gu
% Last modified: 03/20/2025

try
    % Evaluates the model during optimalization
    MRI_flag = MRI_flag; % use this to cohort if using CMR info or not
    [~,targets, inputs, mods] = targetVals_HF(data,PATIENT_NO,ModelWin,MRI_flag);
    [INIparams, INIinit] = estiminiParams(targets,inputs);
    [params, init] = optParams(INIparams, INIinit, mods,modifiers);
    params_no = struct2array(params);
    init_no = struct2array(init);
    if any(params_no <= 0)||any(init_no(2:end) <= 0)
        error('bad inital guessing')
    end
    runSim;
    o_vals_no = struct2array(o_vals);
    if any(o_vals_no <= -1)
        error('bad simulation')
    end
catch ME
    total_cost = inf;
    disp(['Error: ', ME.message])
end