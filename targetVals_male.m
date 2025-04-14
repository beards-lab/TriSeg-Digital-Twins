function [targetVals, inputVals, mods] = targetVals_male()
%% Function Purpose:
% Generate target and input data for a healthy 20-year-old male based on literature.

% Created by Andrew Meyer
% Last modified: 10/29/2024

%% General Information
inputVals.Sex = 1;
inputVals.HR = 60; 
inputVals.TBV = 5700; % 200 pounds 6 feets
targetVals.SBP = 123.5;
targetVals.DBP = 75.7;

%% Echo measurements
targetVals.Hed_LW = 0.93; 
targetVals.Hed_SW = 0.92; 
targetVals.EAr = 1.36; 
inputVals.Hed_RW = 0.35;

%% CMR measurements
targetVals.LVEDV = 145; % b-SSFP, papillary muscles in LV lumen volume
targetVals.LVESV = 53; % b-SSFP, papillary muscles in LV lumen volume
% inputVals.LVESV = 55; % b-SSFP, exclude papillary muscles

targetVals.RVEDV = 166; % b-SSFP, papillary muscles RV lumen volume
targetVals.RVESV = 73; % b-SSFP, papillary muscles RV lumen volume
% inputVals.RVEDV = 150; % b-SSFP, papillary muscles RV lumen volume
% inputVals.RVESV = 60; % b-SSFP, papillary muscles RV lumen volume
% inputVals.RVESV = 73; % b-SSFP, exclude papillary muscles
 
targetVals.LAVmax = 72; % b-SSFP Biplane area-length
targetVals.LAVmin = 25; % b-SSFP Biplane area-length
targetVals.RAVmax = 65; % b-SSFP Biplane area-length
targetVals.RAVmin = 32; % b-SSFP Biplane area-length
 
targetVals.LV_m = 121; % b-SSFP, include papillary muscles. We should constrain septal:left thickness ratio so that it matches the ratio from echo
targetVals.RV_m = 66; % b-SSFP, include papillary muscles

%% RHC measurements
targetVals.RAPmax = 6; % not sex-specific
targetVals.RAPmin = 2; % not sex-specific
targetVals.RVEDP = 3; % not sex-specific
targetVals.RVSP = 20.8; % sPAP. don't give sPAP target, use this ...
targetVals.PASP = 20; 
targetVals.PADP = 8.8; % not sex-specific
targetVals.PCWP = 8; % not sex-specific 
inputVals.CVP = 4; %
% targetVals.CO = 5.76; 

% Feng dichrotic notch relationships
targetVals.DNA = 1.32*((targetVals.SBP-targetVals.DBP)/3+targetVals.DBP)-22.6; % dicrotic notch aorta
targetVals.DNP = 1.004*((targetVals.RVSP - targetVals.PADP) / 3 + targetVals.PADP)-0.5; % dicrotic notch pulmonary

inputVals.AVr = 1;
inputVals.MVmg = 1;
inputVals.MVr = 1;
inputVals.AVpg = 1; 
inputVals.TVr = 1;
inputVals.PVr = 1;
inputVals.Pvpg = 0.5;
% mods = {'k_pas_LV','k_pas_RV','k_act_LV','k_act_RV','C_SA','C_PA','R_SA','R_PA','R_atria','R_m_o','Vw_LV','Vw_RV','LvSepR','R_tPA','R_tSA','K_P','B_P'};
mods = {'k_pas_LV','k_pas_RV','k_act_LV','k_act_RV',...
    'C_SA','C_PA','R_SA','R_PA','R_Veins','R_m_o',...
    'Vw_LV','LvSepR','Vw_RV',... % LV Septum ratio provides information on both LV and RV, so it replaces Vw_SEP
    'Amref_LV','Amref_RV',...
    'K1','expPeri',... 
    'R_tPA','R_tSA'};

