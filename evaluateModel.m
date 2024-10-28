function total_cost = evaluateModel(modifiers,data,PatID,ModelWin)
try
    % Evaluates the model during optimalization

    print_sim = true;

    [~,targets, inputs, mods] = targetVals_HF(data,PatID,ModelWin);

    [~,params, init] = estimParams(targets,inputs,mods,modifiers);
    params_no = struct2array(params);
    init_no = struct2array(init);
    if any(params_no <= 0)||any(init_no(2:end) <= 0)
        warning('bad inital guessing')
    end
    runSim;
catch
    total_cost = inf;
end
