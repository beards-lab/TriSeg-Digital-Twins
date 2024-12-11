function total_cost = evaluateModelUW(modifiers,data,PatID)
%% Function Purpose:
% This function is used to compute a cost function in runSim.m with new modifiers for patients
% during optimization.

% Created by Feng Gu
% Last modified: 10/29/2024

try
    % Evaluates the model during optimalization
    print_sim = true;
    [targets, inputs, mods] = targetVals_UW(data,PatID);
    [~,params, init] = estimParams(targets,inputs,mods,modifiers);
    params_no = struct2array(params);
    init_no = struct2array(init);
    if any(params_no <= 0)||any(init_no(2:end) <= 0)
        error('bad inital guessing')
    end
    runSim;
catch
    total_cost = inf;
end