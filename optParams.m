function initparams = optParams(initparams, mods, modifiers)
%% Parameter adjustment
geom_pars = {'Vw_LV','Vw_RV','Amref_LV','Amref_RV','LvSepR'}; % Geometerical parametes need to be adjusted first to get a optimized inital conditions for the ODE.
geom_pars_i = find(contains(mods,geom_pars));
for i = 1:length(mods)
    if(isfield(initparams, mods{i}) && ~ismember(i,geom_pars_i))  % eliminate adjustable geometry parameters from mods since we already adjusted them above
        % As of 08/20/2024, all modifiers directly related to blood volume have been removed.
        % The functional part in this loop now only involves R_Veins and Vw_LV.
        if(strcmp(mods{i},'V0c_coeff'))
            initparams.V0c_coeff = modifiers(i);
            if(isfield(inputData, 'LAVmax'))
                initparams.LAV0c = modifiers(i) * LAVmax;
                initparams.RAV0c = initparams.LAV0c;
            end
            if(isfield(inputData,'RAVmax'))
                initparams.RAV0c = modifiers(i) * RAVmax;
            end
        elseif(strcmp(mods{i},'V0u_coeff')) 
            initparams.V0u_coeff = modifiers(i);
            initparams.LAV0u = modifiers(i) * initparams.LAV0u;
            initparams.RAV0u = modifiers(i) * initparams.RAV0u;
        elseif(strcmp(mods{i},'R_Veins'))
            initparams.R_PV = modifiers(i) * initparams.R_Veins;
            initparams.R_SV = modifiers(i) * initparams.R_Veins;
        elseif(strcmp(mods{i},'Vw_LV'))  %<-- assume LV and SEP have the same modifier, make sure the change in thickness propotional to TTE measurement
            initparams.Vw_SEP = modifiers(i) * initparams.Vw_SEP*targets.Hed_SW/targets.Hed_LW*initparams.Amref_SEP/initparams.Amref_LV*...
                (LvSepR/(1-LvSepR));
            initparams.Vw_LV = modifiers(i) * initparams.Vw_LV;
        else
            initparams.(mods{i}) = initparams.(mods{i})*modifiers(i);
        end
    elseif isfield(init, mods{i})
        assert(isfield(init, mods{i})); 
        init.(mods{i}) = init.(mods{i})*modifiers(i);
    end
end
end