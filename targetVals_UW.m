function [targetVals, inputVals, mods] = targetVals_UW(data,Patient_no,MRI_flag)
%% Function Purpose:
% Retrieves data for input and target, with default modifiers assigned.
% No processing decisions—only data collection.
% Sets direct measurements from Echo and RHC as targets.

% Created by Feng Gu
% Last modified: 03/20/2025

% Inputs:
% Data        - Table of data extracted from the narrative
% Patient_no  - Vector to locate the patient
% MRI_flag    - 1 stands for reading info from CMR, other is not reading

% Outputs:
% targetVals  - Structure of measurements that the model tries to fit
% inputVals   - Structure of other necessary variables to build the model
%               which the model does not fit
% mods        - Cell array of names to adjust selected parameters

Data = data(Patient_no,:);
%% Assign input values
% Default inputs from data
in_0 = {'Height', 'Weight', 'Sex','Age'};
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
inputVals.HR = nanmean([Data.('HR') Data.('HR_Echo') Data.('heartRate')]);% calculate the average HR of echo, RHC, and initial assessment values
% inputVals.HR = 70;% calculate the average HR of echo, RHC, and initial assessment values
%% Assign target values
% Target from RHC
if ~isnan(Data.('SBP'))
    targetVals.SBP = Data.('SBP');
end

if ~isnan(Data.('DBP'))
    targetVals.DBP = Data.('DBP'); % trust cath more than cuff
end

if ~isnan(Data.('RAPmax'))
    targetVals.RAPmax = Data.('RAPmax');
end

if ~isnan(Data.('PCWP'))
    targetVals.PCWP = Data.('PCWP');
end

if ~isnan(Data.('RVSP'))
    targetVals.RVSP = Data.('RVSP');
end

if ~isnan(Data.('RVEDP'))
    targetVals.RVEDP = Data.('RVEDP');
end

if ~isnan(Data.('P_RV_min'))
    targetVals.P_RV_min = Data.('P_RV_min');
end

if ~isnan(Data.('PPAs'))
    targetVals.PASP = Data.('PPAs');
end

if ~isnan(Data.('PPAd'))
    targetVals.PADP = Data.('PPAd');
end

if  ~isnan(Data.('CO'))
    targetVals.CO = Data.('CO');
elseif  ~isnan(Data.('CO_TD'))
    targetVals.CO = Data.('CO_TD');
end


% Targets from Echo
if ~isnan(Data.('LVIDd'))
    targetVals.LVIDd = Data.('LVIDd');
end

if ~isnan(Data.('LVIDs'))
    targetVals.LVIDs = Data.('LVIDs');
end

if ~isnan(Data.('Hed_LW'))
    targetVals.Hed_LW = Data.('Hed_LW');
else
    targetVals.LV_m = Data.('LV_m');
end

if ~isnan(Data.('Hed_SW'))
    targetVals.Hed_SW = Data.('Hed_SW');
end

if MRI_flag ==1
    if ~isnan(Data.('Hed_RW'))
        targetVals.Hed_RW = Data.('Hed_RW');
    end
end

if ~isnan(Data.('EF'))
    targetVals.EF = Data.('EF'); % this is from ECHO
else
    targetVals.EF = Data.('EF_LV'); % this is from ECHO
end

if ~isnan(Data.('EAr'))
    targetVals.EAr = Data.('EAr');
end


% Targets from MRI
% There is a strange aspect here: I am not using LV information from CMR,
% regardless of the MRI_flag, since the goal is to accommodate the condition
% when CMR is not available. Instead, I use LV information from echo to identify
% the relationship between LV and RV, with RV data coming from CMR.

if MRI_flag == 1
    if ~isnan(Data.('RVEDV'))
        targetVals.RVEDV = Data.('RVEDV');
    end

    if ~isnan(Data.('RVESV'))
        targetVals.RVESV = Data.('RVESV');
    end

    if ~isnan(Data.('LVEDV'))
        targetVals.LVEDV = Data.('LVEDV');
    end

    if ~isnan(Data.('LVESV'))
        targetVals.LVESV = Data.('LVESV');
    end

    if ~isnan(Data.('LV_m'))
        targetVals.LV_m = Data.('LV_m');
    end

end


% LAVmax can be assigned as either an input or a target
if inputVals.Sex == 1
    r = inputVals.TBV / 5700; % ratio to normal
elseif inputVals.Sex == 2
    r = inputVals.TBV / 4300; % ratio to normal
end
if ~isnan(Data.('LAVmax'))
    targetVals.LAVmax = Data.('LAVmax'); % this is form ECHO Simpson
    if targetVals.LAVmax <= 0.6*(targetVals.CO/inputVals.HR*1e3)
        if(inputVals.Sex == 1) % male
            inputVals.LAVmin = 25*r;
        else % female
            inputVals.LAVmin = 22*r;
        end
    else
        inputVals.LAVmin = targetVals.LAVmax - 0.6*(targetVals.CO/inputVals.HR*1e3);
    end
else % as an input
    if ~isnan(Data.('LADiacm'))
        inputVals.LAVmax = 25.17*Data.('LADiacm')-54.92; % from . (J Am Soc Echocardiogr 2017;30:262-9.)
        if inputVals.LAVmax <= 0.6*(targetVals.CO/inputVals.HR*1e3)
            if(inputVals.Sex == 1) % male
                inputVals.LAVmax = 72*r;
                inputVals.LAVmin = 25*r;
            else % female
                inputVals.LAVmax = 64*r;
                inputVals.LAVmin = 22*r;
            end
        else
            inputVals.LAVmin = inputVals.LAVmax - 0.6*(targetVals.CO/inputVals.HR*1e3);
        end
    else
        if(inputVals.Sex == 1) % male
            inputVals.LAVmax = 72*r;
            inputVals.LAVmin = 25*r;
        else % female
            inputVals.LAVmax = 64*r;
            inputVals.LAVmin = 22*r;
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
if ~isnan(Data.('RAPmax'))
    if targetVals.RAPmax < 22
        inputVals.CVP = targetVals.RAPmax - 2; % assuming CVP basing on RAP
    else
        inputVals.CVP = 20;
    end
else
    inputVals.CVP = 4;
end

if ~MRI_flag == 1
    if (isfield(targetVals,'EF'))
        CO = 1000 * targetVals.CO / 60;
        SV = 60 * CO / inputVals.HR;
        inputVals.LVEDV =   SV/ (targetVals.EF*0.01);
        inputVals.LVESV =   inputVals.LVEDV-SV;
    elseif isfield(targetVals,'LVIDd')
        inputVals.LVEDV = targetVals.LVIDd^3*0.7851+97.32;
        inputVals.LVESV = targetVals.LVIDs^3*0.9185+63.92;
    end
end


%% Right side assumptions based on Lasso regression trained on the UW cohort
if ~MRI_flag == 1
    load LassoRV.mat
    Raw_Lasso_RVEDV = inputVals.Sex*k_RVEDV(1) + inputVals.Age*k_RVEDV(2) +...
        inputVals.Height*k_RVEDV(3) + inputVals.Height*inputVals.Weight*k_RVEDV(4) +...
        inputVals.HR*k_RVEDV(5) + targetVals.SBP*k_RVEDV(6) +...
        targetVals.PASP*k_RVEDV(7) + targetVals.PADP*k_RVEDV(8) +...
        targetVals.CO*k_RVEDV(9) + Data.('AVr') * k_RVEDV(10) + ...
        Data.('TVr') * k_RVEDV(11) + Data.('PVr') * k_RVEDV(12) + ...
        b_RVEDV;
    C_Lasso_RVEDV = Raw_Lasso_RVEDV*ck_RVEDV+cb_RVEDV;
    if  C_Lasso_RVEDV > 400
        C_Lasso_RVEDV = 400;
    end
    targetVals.RVEDV = C_Lasso_RVEDV;
    Raw_Lasso_RVEF = inputVals.Sex * k_RVEF(1) + inputVals.Age * k_RVEF(2) + ...
        targetVals.SBP * k_RVEF(3) + targetVals.PADP * k_RVEF(4) + ...
        targetVals.CO * k_RVEF(5) + targetVals.EF * k_RVEF(6) + ...
        Data.('TVr') * k_RVEF(7) + Data.('PVr') * k_RVEF(8) + ...
        b_RVEF;
    C_Lasso_RVEF = Raw_Lasso_RVEF * ck_RVEF + cb_RVEF;
    targetVals.RVEF = C_Lasso_RVEF;
    inputVals.RVESV = targetVals.RVEDV * (100-targetVals.RVEF) * 0.01;
    
end
if ~isfield(targetVals,'Hed_RW')
    load LassoRV.mat
    Raw_Lasso_TRW = inputVals.Height * k_TRW(1) + targetVals.SBP * k_TRW(2) + ...
        targetVals.RAPmax * k_TRW(3) + targetVals.RVSP * k_TRW(4) + ...
        targetVals.PADP * k_TRW(5) + Data.('TVr') * k_TRW(6) + ...
        b_TRW;
    C_Lasso_TRW = Raw_Lasso_TRW * ck_TRW + cb_TRW;
    if  C_Lasso_TRW > 0.95
        C_Lasso_TRW = 0.95;
    end
    targetVals.Hed_RW = C_Lasso_TRW;
end

%% Parameters requiring modification (mods), used in the function estimParams.m
% Default parameters should be always optimized
mods_0 = {'k_pas_LV','k_pas_RV','k_act_LV','k_act_RV',...
    'C_SA','C_PA','R_SA','R_PA','R_Veins','R_m_o',...
    'Vw_LV','LvSepR','Vw_RV',... % LV Septum ratio provides information on both LV and RV, so it replaces Vw_SEP
    'Amref_LV','Amref_RV',...
    ... 'V_SV_s','V0u_coeff','V0c_coeff', % fixed blood volume
    ... 'K1','expPeri',...
    'R_tPA','R_tSA'};

%% Add additional mods based on the availability of the patient's data.
mods_pN = cellstr(string(cell(0)));
tg_pN = cellstr(string(cell(0)));
% Evaluate PV stenosis based on RVSP and PASP
if ~isnan(Data.('RVSP'))
    if(targetVals.PASP <= 0.925 * targetVals.RVSP) % if PPAs is 7.5% less than RVSP
        mods_pN{end + 1} = 'R_p_o';
    end
end


if ~isnan(Data.('PCWP')) && ~isnan(Data.('RAPmax'))
    if abs(Data.('PCWP') - Data.('RAPmax') + 2) <= 5
        mods_pN{end + 1} = 'expPeri';
    else
        mods_pN{end + 1} = 'K1';
    end
end

vlv_def = dictionary('MVr','R_m_c', ...
    'MS','R_m_o', ...
    'AVr','R_a_c', ...
    'AS','R_a_o', ...
    'TVr','R_t_c', ...
    'PVr','R_p_c', ...
    'PS','R_p_o');
for k = keys(vlv_def)'
    if(Data.(k) >= 1.2)
        mods_pN(end + 1) = cellstr(vlv_def(k));
        tg_pN(end + 1) = cellstr(k);
    end
end

mods = [mods_0 mods_pN];
mods = unique(mods);

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
    'RVEF',{[5 85]}, ...
    'Hed_LW',{[0.2 3]}, ...
    'Hed_SW',{[0.2 3]}, ...
    'Hed_RW',{[0.2 3]}, ...
    'EAr',{[0.2 8]}, ...
    'LAVmax',{[15 400]}, ... % 263 has LAVmax 377
    'MVr',{[1 4]}, ...
    'MVr',{[1 4]}, ...
    'MS',{[1 4]}, ...
    'AVr',{[1 4]}, ...
    'AVpg',{[0 150]}, ...
    'AS',{[1 4]}, ...
    'TVr',{[1 4]}, ...
    'PVr',{[1 4]}, ...
    'PVpg',{[0 100]}, ...
    'PS',{[1 4]}, ...
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
    'LVEDV',{[10 1050]}, ...
    'LVESV',{[10 800]}, ...
    'RVEDV',{[10 700]}, ...
    'RVESV',{[10 650]}, ...
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
end

