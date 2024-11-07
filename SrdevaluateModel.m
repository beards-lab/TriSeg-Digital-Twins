function total_cost = SrdevaluateModel(modifiers,GENDER)
% This function is used to compute a cost function in runSim.m with new modifiers for canonical
% subjects during optimization.

% Created by Andrew Meyer and Feng Gu
% Last modified: 10/29/2024

try
    % Evaluates the model during optimalization
    print_sim = true;
    if GENDER ==1
        [targets, inputs, mods] = targetVals_male();
    else
        [targets, inputs, mods] = targetVals_female();
    end
    [~ ,params, init] = estimParams(targets,inputs,mods,modifiers);
    runSim;
    if total_cost<0
        total_cost = inf;
    end
catch
    total_cost = inf;
end
