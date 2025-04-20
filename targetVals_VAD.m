function [Windowdate, targetVals, inputVals, mods] = targetVals_VAD(data,Patient_no,ModelWin)
%% Function Purpose:
% Retrieves data for input and target, with default modifiers assigned.
% No processing decisions—only data collection.
% Sets direct measurements from Echo, RHC, and Cardiac MRI as targets.

% Created by Feng Gu
% Last modified: 3/31/2025

% Inputs:
% Data        - Structure of data extracted from the narrative
% Patient_no  - Vector to locate the patient

% Outputs:
% targetVals  - Structure of measurements that the model tries to fit
% inputVals   - Structure of other necessary variables to build the model 
%               which the model does not fit
% mods        - Cell array of names to adjust selected parameters
% Windowdate  - Date the model was built at
Data = data(Patient_no).snapshots(ModelWin);
Windowdate = mean([Data.RHCDate;Data.TTEDate]);
%% Assign input values
% Default inputs from data
in_0 = {'Height', 'Weight', 'Sex'};
for i = 1:length(in_0)
    inputVals.(in_0{i}) = Data.(in_0{i});
end
inputVals.Age =  round(years(Windowdate-Data.Birthday));
% Default inputs from estimation
if(inputVals.Sex == 1) % Male
    inputVals.TBV = (0.3669*(inputVals.Height/100)^3 +0.03219*inputVals.Weight+0.6041)*1000; % mL (Nadler's Equation Male)
else                      % Female
    assert(inputVals.Sex == 2);
    inputVals.TBV = (0.3561*(inputVals.Height/100)^3 +0.03308*inputVals.Weight+0.1833)*1000; % mL (Nadler's Equation Female)
end
inputVals.HR = nanmean([Data.('HR_rhc') Data.('HR_vitals')]);% Fixed: This probably get a lof of situation

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

if ~isnan(Data.('IVSd'))
    if Data.('IVSd') > 3
        targetVals.Hed_SW = Data.('IVSd')/10; % this is from ECHO
    else
        targetVals.Hed_SW = Data.('IVSd');
    end
elseif ~isnan(Data.('LVPWd'))
    if Data.('LVPWd') > 3
        targetVals.Hed_SW = Data.('LVPWd')/10; % this is from ECHO
    else
        targetVals.Hed_SW = Data.('LVPWd');
    end
end

if ~isnan(Data.('LVPWd'))
    if Data.('LVPWd') > 3
        targetVals.Hed_LW = Data.('LVPWd')/10; % this is from ECHO
    else
        targetVals.Hed_LW = Data.('LVPWd');
    end
elseif ~isnan(Data.('IVSd'))
    if Data.('IVSd') > 3
        targetVals.Hed_LW = Data.('IVSd')/10; % this is from ECHO
    else
        targetVals.Hed_LW = Data.('IVSd');
    end
end



if ~isnan(Data.('LVEF_tte'))
    targetVals.EF = Data.('LVEF_tte'); % this is from ECHO
end

if ~isnan(Data.('EA'))
    targetVals.EAr = Data.('EA');
end


% LAVmax can be assigned as either an input or a target
if inputVals.Sex == 1
    r = inputVals.TBV / 5700; % ratio to normal
elseif inputVals.Sex == 2
    r = inputVals.TBV / 4300; % ratio to normal
end
if ~isnan(Data.('VLA'))
    targetVals.LAVmax = Data.('VLA'); % this is form ECHO Simpson
    if  targetVals.LAVmax <= 0.6*(targetVals.CO/inputVals.HR*1e3)
        if(inputVals.Sex == 1) % male
            inputVals.LAVmin = 25*r;
        else % female
            inputVals.LAVmin = 22*r;
        end
    else
        inputVals.LAVmin = targetVals.LAVmax - 0.6*(targetVals.CO/inputVals.HR*1e3);
    end
else % as an input
    if ~isnan(Data.('LAd'))
        if Data.('LAd') >= 10
            inputVals.LAVmax = 25.17*Data.('LAd')/10-54.92; % from . (J Am Soc Echocardiogr 2017;30:262-9.)
        else
            inputVals.LAVmax = 25.17*Data.('LAd')-54.92; % from . (J Am Soc Echocardiogr 2017;30:262-9.)
        end
        if inputVals.LAVmax <= 0.6*(targetVals.CO/inputVals.HR*1e3)
            if(inputVals.Sex == 1) % male
                inputVals.LAVmax = 72*r;
                inputVals.LAVmin = 25*r;
            else % fema
                inputVals.LAVmax = 64*r;
                inputVals.LAVmin = 22*r;
            end
        else
            inputVals.LAVmin = inputVals.LAVmax -  0.6*(targetVals.CO/inputVals.HR*1e3);
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
if ~isnan(Data.('RAm'))
    if targetVals.RAPmean < 18
        inputVals.CVP = targetVals.RAPmean + 2; % assuming CVP basing on RAP
    else
        inputVals.CVP = 20;
    end
else
    inputVals.CVP = 4;
end

if (isfield(targetVals,'EF'))
        CO = 1000 * targetVals.CO / 60;
        SV = 60 * CO / inputVals.HR;
        inputVals.LVEDV =   SV/ (targetVals.EF*0.01);
        inputVals.LVESV =   inputVals.LVEDV-SV;
elseif isfield(targetVals,'LVIDd')
        inputVals.LVEDV = targetVals.LVIDd^3*0.7851+97.32;
        inputVals.LVESV = targetVals.LVIDs^3*0.9185+63.92;
        targetVals.EF = (inputVals.LVEDV-inputVals.LVESV)/inputVals.LVEDV*100; % this is from ECHO
end

%% Right Side assumption

load LassoRV.mat
Raw_Lasso_RVEDV = inputVals.Sex*k_RVEDV(1) + inputVals.Age*k_RVEDV(2) +...
    inputVals.Height*k_RVEDV(3) + inputVals.Height*inputVals.Weight*k_RVEDV(4) +...
    inputVals.HR*k_RVEDV(5) + targetVals.SBP*k_RVEDV(6) +...
    targetVals.PASP*k_RVEDV(7) + targetVals.PADP*k_RVEDV(8) +...
    targetVals.CO*k_RVEDV(9) + Data.('AVr') * k_RVEDV(10) + ...
    Data.('TVr') * k_RVEDV(11) + Data.('PVr') * k_RVEDV(12) + ...
    b_RVEDV;
C_Lasso_RVEDV = Raw_Lasso_RVEDV*ck_RVEDV+cb_RVEDV;
if  C_Lasso_RVEDV > 600
    C_Lasso_RVEDV = 600;
end
if  C_Lasso_RVEDV <50
    C_Lasso_RVEDV = 50;
end
targetVals.RVEDV = C_Lasso_RVEDV;
Raw_Lasso_RVEF = inputVals.Sex * k_RVEF(1) + inputVals.Age * k_RVEF(2) + ...
    targetVals.SBP * k_RVEF(3) + targetVals.PADP * k_RVEF(4) + ...
    targetVals.CO * k_RVEF(5) + targetVals.EF * k_RVEF(6) + ...
    Data.('TVr') * k_RVEF(7) + Data.('PVr') * k_RVEF(8) + ...
    b_RVEF;
C_Lasso_RVEF = Raw_Lasso_RVEF * ck_RVEF + cb_RVEF;
if  C_Lasso_RVEF > 75
    C_Lasso_RVEF = 75;
end
if  C_Lasso_RVEF < 10
    C_Lasso_RVEF = 10;
end
targetVals.RVEF = C_Lasso_RVEF;
inputVals.RVESV = targetVals.RVEDV * (100-targetVals.RVEF) * 0.01;

if  (~isnan(Data.('RAa'))) || (~isnan(Data.('RAv')))
    coeff3= targetVals.RAPmax;
elseif ~isnan(Data.('RAm'))
    coeff3 = targetVals.RAPmean+3.5;
else
    coeff3 = 15.11;
end


if  ~isnan(Data.('RVs'))
    coeff4 = targetVals.RVSP;

else
    coeff4 = targetVals.PASP;
end

Raw_Lasso_TRW = inputVals.Height * k_TRW(1) + targetVals.SBP * k_TRW(2) + ...
    coeff3* k_TRW(3) + coeff4 * k_TRW(4) + ...
    targetVals.PADP * k_TRW(5) + Data.('TVr') * k_TRW(6) + ...
    b_TRW;
C_Lasso_TRW = Raw_Lasso_TRW * ck_TRW + cb_TRW;
if  C_Lasso_TRW > 1
    C_Lasso_TRW = 1;
end
targetVals.Hed_RW = C_Lasso_TRW;


%% Parameters requiring modification (mods), used in the function estimParams.m
% Default parameters should be always optimized
mods_0 = {'k_pas_LV','k_pas_RV','k_act_LV','k_act_RV',...
    'C_SA','C_PA','R_SA','R_PA','R_Veins','R_m_o',...
    'Vw_LV','LvSepR','Vw_RV',... % LV Septum ratio provides information on both LV and RV, so it replaces Vw_SEP
    'Amref_LV','Amref_RV',...
    ... 'V_SV_s','V0u_coeff','V0c_coeff', % fixed blood volume
    'expPeri','K1',...
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
    'PCWPmax',{[0 100]}, ...
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
    'MS',{[1 4]}, ...
    'AS',{[1 4]}, ...
    'PS',{[1 4]}, ...
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
    'LVIDs',{[1 12]}, ...% VAD patient 162
    'LVEDV',{[10 1050]}, ...
    'LVESV',{[10 900]}, ...% VAD patient 149
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
    'LVEDV',{[10 1550]}, ...% VAD patient 162
    'LVESV',{[10 1350]}, ...% VAD patient 149
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
