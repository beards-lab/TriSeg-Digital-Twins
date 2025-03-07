function cost = evaluateModelVwRV(offsetsValue, m_initial, UWpatients, PATIENT_NO,Geo_Opt,offsetsPercentage)
    MRI_flag = 0;
    Geo_Opt = Geo_Opt;
    [targets, inputs, mods] = targetVals_UW(UWpatients,PATIENT_NO,MRI_flag);
    [INIparams, INIinit] = estiminiParams(targets,inputs);
    [params, init] = optParams(INIparams, INIinit, mods,m_initial,targets);
    NewVwRVtarget = params.Vw_RV*(1+offsetsPercentage);
    locsVwRV = ismember(mods, {'Vw_RV'});
    params.Vw_RV = INIparams.Vw_RV*(m_initial(locsVwRV)+offsetsValue);
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
    cost = (params.Vw_RV-NewVwRVtarget).^2;
end