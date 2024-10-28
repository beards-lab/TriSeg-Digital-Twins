function total_cost = SrdevaluateModel(modifiers,GENDER)
try
    % Evaluates the model during optimalization

    print_sim = true;

    if GENDER ==1
        [targets, inputs, mods] = targetVals_male();
    else
        [targets, inputs, mods] = targetVals_female();
    end

    % when we get the mods out, we should pass those that are geometrically
    % relevant to estimParams.
    % mods_geom = {'Amref_LV','Amref_RV','Amref_SEP','Vw_LV','Vw_RV','Vw_SEP'};

    % [targets, inputs, mods] = targetVals_female();
    [params, init] = estimParams(targets,inputs,mods,modifiers);
    runSim;
    if total_cost<0
        total_cost = inf;
    end
catch
    total_cost = inf;
end
