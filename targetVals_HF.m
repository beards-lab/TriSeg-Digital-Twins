function [Windowdate, targetVals, inputVals, mods] = targetVals_HF(data,Patient_no,Window_No)
%% Function Purpose:
% Retrieves data for input and target, with default modifiers assigned.
% No processing decisions—only data collection.
% Sets direct measurements from Echo, RHC, and Cardiac MRI as targets.

% Created by Feng Gu
% Last modified: 10/29/2024

% Inputs:
% Data        - Structure of data extracted from the narrative
% Patient_no  - Vector to locate the patient
% Window_No   - Vector to locate the time point

% Outputs:
% targetVals  - Structure of measurements that the model tries to fit
% inputVals   - Structure of other necessary variables to build the model 
%               which the model does not fit
% mods        - Cell array of names to adjust selected parameters
% Windowdate  - Date the model was built at

Data = data(Patient_no).snapshots(Window_No); 
Windowdate = mean([Data.RHCDate;Data.TTEDate;Data.MRIDate]);
%% Assign input values
% Default inputs from data
in_0 = {'Height', 'Weight', 'Sex'};
for i = 1:length(in_0)
    inputVals.(in_0{i}) = Data.(in_0{i});
end
% Default inputs from estimation
if(inputVals.Sex == 1) % Male
    inputVals.TBV = (0.3669*(inputVals.Height/100)^3 +0.03219*inputVals.Weight+0.6041)*1000; % mL (Nadler's Equation Male)
else                      % Female
    assert(inputVals.Sex == 2);
    inputVals.TBV = (0.3561*(inputVals.Height/100)^3 +0.03308*inputVals.Weight+0.1833)*1000; % mL (Nadler's Equation Female)
end
inputVals.HR = nanmean([Data.('HR_rhc') Data.('HR_vitals') Data.('MRI_HR')]);% Fixed: This probably get a lof of situation

%% Assign target values
% Target from RHC
if ~isnan(Data.('SAs'))
    %Somehow the cath results give me something which more like ventricle
    %not aorta
    targetVals.SBP = Data.('SAs'); % trust cath more than cuff
elseif ~isnan(Data.('NIBPs_vitals'))
    targetVals.SBP = Data.('NIBPs_vitals');
end

if ~isnan(Data.('SAd'))
    targetVals.DBP = Data.('SAd'); % trust cath more than cuff
elseif ~isnan(Data.('NIBPd_vitals'))
    targetVals.DBP = Data.('NIBPd_vitals');
end

if ~isnan(Data.('RAm'))
    targetVals.RAPmean = Data.('RAm');
end

if (~isnan(Data.('RAa'))) || (~isnan(Data.('RAv')))
    targetVals.RAPmax = max([Data.('RAa');Data.('RAv')]);
end

if ~isnan(Data.('PCW'))
    targetVals.PCWP = Data.('PCW');
end

if (~isnan(Data.('PCWa'))) || (~isnan(Data.('PCWv')))
    targetVals.PCWPmax = max([Data.('PCWa');Data.('PCWv')]);
end

if ~isnan(Data.('RVs'))
    targetVals.RVSP = Data.('RVs');
end

if ~isnan(Data.('RVd'))
    targetVals.RVEDP = Data.('RVd');
end

if ~isnan(Data.('RVmin'))
    targetVals.P_RV_min = Data.('RVmin');
end

if ~isnan(Data.('LVs'))
    targetVals.LVESP = Data.('LVs');
end

if ~isnan(Data.('LVd'))
    targetVals.LVEDP = Data.('LVd');
end

if ~isnan(Data.('LVmin'))
    targetVals.P_LV_min = Data.('LVmin');
end

if ~isnan(Data.('PAs'))
    targetVals.PASP = Data.('PAs');
end

if ~isnan(Data.('PAd'))
    targetVals.PADP = Data.('PAd');
end
if ~(Patient_no == 116)
    if  ~isnan(Data.('CO_td'))
        targetVals.CO = Data.('CO_td');
    elseif  ~isnan(Data.('CO_fick'))
        targetVals.CO = Data.('CO_fick');
    elseif  ~isnan(Data.('MRI_COrv'))
        targetVals.CO = Data.('MRI_COrv');
    elseif  ~isnan(Data.('MRI_CO'))
        targetVals.CO = Data.('MRI_CO');
    end
else
    targetVals.CO = (Data.('CO_td')+Data.('CO_fick'))/2;
end

% Targets from Echo
if ~isnan(Data.('LVIDd'))
    if Data.('LVIDd') > 10
        targetVals.LVIDd = Data.('LVIDd')/10; % from echo
    else
        targetVals.LVIDd = Data.('LVIDd');% making the unit to cm. most of TTE report are mm
    end
end

if ~isnan(Data.('LVIDs'))
    if Data.('LVIDs') > 10
        targetVals.LVIDs = Data.('LVIDs')/10; % from echo
    else
        targetVals.LVIDs = Data.('LVIDs');% making the unit to cm. most of TTE report are mm
    end
end

if ~isnan(Data.('LVPWd'))
    if Data.('LVPWd') > 3
        targetVals.Hed_LW = Data.('LVPWd')/10; % this is from ECHO
    else
        targetVals.Hed_LW = Data.('LVPWd');
    end
end

if ~isnan(Data.('IVSd'))
    if Data.('IVSd') > 3
        targetVals.Hed_SW = Data.('IVSd')/10; % this is from ECHO
    else
        targetVals.Hed_SW = Data.('IVSd');
    end
end

if ~isnan(Data.('LVEF_tte'))
    targetVals.EF = Data.('LVEF_tte'); % this is from ECHO
end

if ~isnan(Data.('EA'))
    targetVals.EAr = Data.('EA');
end

% Targets from MRI
if ~isnan(Data.('MRI_LVEDV'))
    targetVals.LVEDV = Data.('MRI_LVEDV');
end

if ~isnan(Data.('MRI_LVESV'))
    targetVals.LVESV = Data.('MRI_LVESV');
end

if ~isnan(Data.('MRI_RVEDV'))
    targetVals.RVEDV = Data.('MRI_RVEDV');
end

if ~isnan(Data.('MRI_RVESV'))
    targetVals.RVESV = Data.('MRI_RVESV');
end

if ~isnan(Data.('MRI_LVMass'))
    targetVals.LV_m = Data.('MRI_LVMass');
end

if ~isnan(Data.('MRI_RVMass'))
    targetVals.RV_m = Data.('MRI_RVMass');
end

% LAVmax can be assigned as either an input or a target
if inputVals.Sex == 1
    r = inputVals.TBV / 5700; % ratio to normal
elseif inputVals.Sex == 2
    r = inputVals.TBV / 4300; % ratio to normal
end
if ~isnan(Data.('VLA'))
    targetVals.LAVmax = Data.('VLA'); % this is form ECHO Simpson
    if targetVals.LAVmax <= 0.6*(targetVals.LVEDV-targetVals.LVESV)
        if(inputVals.Sex == 1) % male
            inputVals.LAVmin = 25*r;
        else % female
            inputVals.LAVmin = 22*r;
        end
    else
        inputVals.LAVmin = targetVals.LAVmax - 0.6*(targetVals.LVEDV-targetVals.LVESV);
    end
else % as an input
    if ~isnan(Data.('LAd'))
        if Data.('LAd') >= 10
            inputVals.LAVmax = 25.17*Data.('LAd')/10-54.92; % from . (J Am Soc Echocardiogr 2017;30:262-9.)
        else
            inputVals.LAVmax = 25.17*Data.('LAd')-54.92; % from . (J Am Soc Echocardiogr 2017;30:262-9.)
        end
        if inputVals.LAVmax <= 0.6*(targetVals.LVEDV-targetVals.LVESV)
            if(inputVals.Sex == 1) % male
                LA_LV_r = 72 / 155*0.8; % from baseline resting, Max LA: LVEDV
                inputVals.LAVmax = targetVals.LVEDV * LA_LV_r;
                inputVals.LAVmin = 25*r;
            else % female
                LA_LV_r = 64 / 123*0.8; % from baseline resting, Max LA: LVEDV
                inputVals.LAVmax = targetVals.LVEDV * LA_LV_r;
                inputVals.LAVmin = 22*r;
            end
        else
            inputVals.LAVmin = inputVals.LAVmax - 0.6*(targetVals.LVEDV-targetVals.LVESV);
        end
    end
end

if(inputVals.Sex == 1) % male
    inputVals.RAVmin = 32/25*inputVals.LAVmin; % from Baseline_resting numbers. need these for diastolic displacement for model intial conditions
else % female
    inputVals.RAVmin = 23/22*inputVals.LAVmin; % from Baseline_resting numbers
end

% Add additional targets to enhance physiological relevance of the model
targetVals.DNA = 1.32*((targetVals.SBP-targetVals.DBP)/3+targetVals.DBP)-22.6; % dicrotic Notch Aorta
targetVals.DNP = 1.004*((targetVals.PASP-targetVals.PADP)/3+targetVals.PADP)-0.5; % dicrotic Notch Pulmonary
if ~isnan(Data.('RAm'))
    if targetVals.RAPmean < 18
        inputVals.CVP = targetVals.RAPmean + 2; % assuming CVP basing on RAP
    else
        inputVals.CVP = 20;
    end
else
    inputVals.CVP = 4;
end

%% Assign synthetic target values for patients without MRILVmass
% Uses a linear regression model from 146 patients.
% Thickness predictions cluster into two groups based on presence of LV mass target.
% For echo-only patients, wall thickness appears greater; this adjustment mitigates that bias.
if isnan(Data.('MRI_LVMass'))
    FakeHed_LW = 0.4235*targetVals.Hed_LW + 0.2021;% coming from 146 patients who have LV mass r2 = 0.50
    FakeHed_SW = 0.6342*targetVals.Hed_SW - 0.003236;% r2 = 0.79
    left = (FakeHed_LW + FakeHed_SW) / 2;
    use_Vw_LV = 0;
    LvSepR = 2/3;
    if ~isnan(Data.('MRI_RVMass'))
        right = targetVals.RV_m;
        use_Vw_RV = 1;
    elseif isfield(Data,"Hed_RV")&&~isnan(Data.('Hed_RV'))
        right = targetVals.Hed_RV;
        use_Vw_RV = 0;
    else
        use_Vw_RV = 0;
        if ~isnan(Data.('RVd'))
            P_RV = targetVals.RVEDP;
        elseif ~isnan(Data.('RAm'))
            P_RV = targetVals.RAPmean;
        elseif isfield(targetVals,'PCWP')
            P_RV = targetVals.PASP*targetVals.PCWP/targetVals.SBP;
        end
        coeff = [14853.8592619735
            14531.7636186220
            14769.3961036171
            11440.8528760998
            -21873.4802314862
            456.030613986307
            -4913.60117297433
            156836.996351639
            -155270.878092212];
        k_RV = coeff(1)*100*(targetVals.LVEDV-targetVals.LVESV)/targetVals.LVEDV+...
            coeff(2)*targetVals.PASP+...
            coeff(3)*targetVals.PADP+...
            coeff(4)*targetVals.PCWP+...
            coeff(5)*targetVals.CO+...
            coeff(6)*targetVals.SBP+...
            coeff(7)*targetVals.DBP+...
            100*coeff(8)*max([targetVals.Hed_SW targetVals.Hed_LW])+...
            100*coeff(9)*min([targetVals.Hed_SW targetVals.Hed_LW]);
        if inputVals.Sex == 1
            k_RV = k_RV/7.2906; % calibrate using coefficients estimated from a canonical subject
        else
            k_RV = k_RV/11.2290; % calibrate using coefficients estimated from a canonical subject
        end
         right = [P_RV k_RV];
    end
    [~,~, Vw0, ~] = geom_0(targetVals.LVEDV, targetVals.RVEDV, left, use_Vw_LV, right, use_Vw_RV, LvSepR,inputVals);
    Vw_LV = Vw0(1);
    Vw_SEP = Vw0(2);
    targetVals.FakeLV_m = (Vw_LV + Vw_SEP) * 1.055;
end

%% Parameters requiring modification (mods), used in the function estimParams.m
% Default parameters should be always optimized
mods_0 = {'k_pas_LV','k_pas_RV','k_act_LV','k_act_RV',...
    ... A minor naming issue: initially, I mistakenly labeled R_Veins as R_atria and sent it to the
    ... cluster for computation. Later, I realized R_Veins was a more accurate name, but the cluster
    ... results still use R_atria. To keep consistent with the cluster outputs, I’ll keep R_atria for
    ... now and change it to R_Veins after applying the function unique.
    'C_SA','C_PA','R_SA','R_PA','R_atria','R_m_o',...
    'Vw_LV','LvSepR',... % LV Septum ratio provides information on both LV and RV, so it replaces Vw_SEP
    ... 'V_SV_s','V0u_coeff','V0c_coeff', % fixed blood volume 
    'K_P','B_P',...
    'R_tPA','R_tSA'};

%% Add additional mods based on the availability of the patient's data.
mods_pN = cellstr(string(cell(0)));
tg_pN = cellstr(string(cell(0)));
% Evaluate PV stenosis based on RVSP and PASP
if ~isnan(Data.('RVs'))
    if(targetVals.PASP <= 0.925 * targetVals.RVSP) % if PPAs is 7.5% less than RVSP
        mods_pN{end + 1} = 'R_p_o';
    end
end
if ~isnan(Data.('MRI_RVMass'))
        mods_pN{end + 1} = 'Vw_RV'; % only opt when having targets
end
if ~isnan(Data.('LVs'))
    if(targetVals.SBP <= 0.925 * targetVals.LVESP) % if SBP is 7.5% less than LVESP
        mods_pN{end + 1} = 'R_a_o';
    end
end

vlv_def = dictionary('MVr','R_m_c', ...
    'MVmg','R_m_o', ...
    'AVr','R_a_c', ...
    'AVpg','R_a_o', ...
    'TVr','R_t_c', ...
    'PVr','R_p_c', ...
    'PVpg','R_p_o');
for k = keys(vlv_def)'
    if(Data.(k) >= 1.2)
        mods_pN(end + 1) = cellstr(vlv_def(k));
        tg_pN(end + 1) = cellstr(k);
    end
end

mods = [mods_0 mods_pN];
mods = unique(mods);

mods{strcmp(mods,'R_atria')} = 'R_Veins'; 
for i = 1:length(tg_pN)
    targetVals.(tg_pN{i}) = Data.(tg_pN{i});
end

%% Data validation (using empirical bounds to safeguard against unit errors)

% Validate targets and inputs
tg_bounds = dictionary('SBP',{[50 300]}, ...
    'DBP',{[15 170]}, ...
    'DNP',{[0 100]}, ...
    'DNA',{[20 300]},...
    'tPLVmax',{[0 1]},...
    'tPLVmin',{[0 1]},...
    'RAPmean',{[0 50]}, ...
    'RAPmax',{[0 50]}, ...
    'PCWPmax',{[0 75]}, ...
    'RVSP',{[5 150]}, ...
    'RVEDP',{[0 50]}, ...
    'P_RV_min',{[0 50]}, ...
    'LVEDP',{[0 100]}, ...
    'LVESP',{[5 250]}, ...
    'P_LV_min',{[0 50]}, ...
    'PASP',{[0 140]}, ...% 282 has 130mmHg
    'PADP',{[0 75]}, ...
    'PCWP',{[3 50]}, ...
    'CO',{[1 15]}, ...
    'EF',{[5 85]}, ...
    'Hed_LW',{[0.2 3]}, ...
    'Hed_SW',{[0.2 3]}, ...
    'Hed_RW',{[0.2 3]}, ...
    'EAr',{[0.2 8]}, ...
    'LAVmax',{[15 400]}, ... % 263 has LAVmax 377
    'MVr',{[1 4]}, ...
    'MVmg',{[0 100]}, ...
    'AVr',{[1 4]}, ...
    'AVpg',{[0 150]}, ...
    'TVr',{[1 4]}, ...
    'PVr',{[1 4]}, ...
    'PVpg',{[0 100]}, ...
    'RVEDV',{[10 700]}, ...
    'RVESV',{[10 650]}, ...
    'LVIDd',{[2 15]}, ...
    'LVIDs',{[1 10]}, ...
    'LVEDV',{[10 1050]}, ...
    'LVESV',{[10 800]}, ...
    'RV_m',{[15 400]},...% 495 only 23
    'LV_m',{[35 400]},...
    'FakeLV_m',{[35 400]});

tg_fn = fieldnames(targetVals);
for i = 1:length(tg_fn)
    bounds_i = cell2mat(tg_bounds(tg_fn{i}));
    assert(targetVals.(tg_fn{i}) >= bounds_i(1) ...
        && targetVals.(tg_fn{i}) <= bounds_i(2),sprintf(tg_fn{i}));
end

in_bounds = dictionary('Height',{[90 230]}, ...
    'Weight',{[30 250]}, ...% 267 over 200kg
    'Sex',{[1 2]}, ...
    'HR',{[30 200]}, ...
    'Hed_RW',{[0.2 1.5]}, ...
    'LAVmax',{[5 250]},...
    'LVEDV',{[25 250]}, ...
    'LVESV',{[10 200]}, ...
    'RVEDV',{[75 325]}, ...
    'CVP',{[2 20]}); % TBV, LAVmin, RAVmin, are all derived.

in_fn = fieldnames(inputVals);
in_keys = isKey(in_bounds,string(in_fn));

for i = 1:length(in_keys)
    if(in_keys(i)) % for those inputs that aren't derived
        bounds_i = cell2mat(in_bounds(in_fn{i}));
        assert(inputVals.(in_fn{i}) >= bounds_i(1) ...
            && inputVals.(in_fn{i}) <= bounds_i(2),sprintf(in_fn{i}));
    end
end


