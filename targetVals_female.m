function [targetVals, inputVals, mods] = targetVals_female()
%% Function Purpose:
% Generate target and input data for a healthy 20-year-old female based on literature.

% Created by Andrew Meyer
% Last modified: 10/29/2024

%% General Information
inputVals.Sex = 2;
inputVals.HR = 60; 
inputVals.TBV = 4300; % 165 pounds 5'5 Butterworth JF, Mackey DC, Wasnick JD. Morgan and Mikhail3s Clinical Anesthesiology, 7th edition. McGraw-Hill Education. 2022.
targetVals.SBP = 116.5;
targetVals.DBP = 72.9;

%% Echo measurements
targetVals.Hed_LW = 0.85; % Posterior LW thickness, end-diastole (cm)
targetVals.Hed_SW = 0.82; % SW thickness, end-diastole (cm)
targetVals.EAr = 1.38;

%% CMR measurements
targetVals.LVEDV = 112; % b-SFFP, papillary muscles included in left ventricular volume
targetVals.LVESV = 39; % b-SFFP, papillary muscles included in left ventricular volume
% inputVals.LVESV = 39; % b-SFFP, exclude papillary muscles
targetVals.RVEDV = 122; % b-SFFP, papillary muscles included in right ventricular volume
targetVals.RVESV = 50; % b-SFFP, papillary muscles included in right ventricular volume
% inputVals.RVESV = 50; % b-SFFP, papillary muscles included in right ventricular volume
targetVals.LAVmax = 64; % SFFP, Simpson6s method; LA appendage excluded
targetVals.LAVmin = 22; % SFFP, Simpsonis method; LA appendage excluded
targetVals.RAVmax = 53; % SFFP, Simpson's method; RA appendage excluded
targetVals.RAVmin = 23; % SFFP, Simpson's method; RA appendage excluded 
targetVals.LV_m = 83; % b-SFFP, include papillary muscles. We should constrain septal:left thickness ratio so that it matches the ratio from echo
targetVals.RV_m = 48; % b-SFFP, include papillary muscles

%% RHC measurements
targetVals.RAPmax = 6; % not sex-specific 
targetVals.RAPmin = 2; % not sex-specific
targetVals.RVEDP = 3; % not sex-specific
targetVals.RVSP = 20.8; % sPAP. don't give sPAP target, use this ...
targetVals.PADP = 8.8; % not sex-specific
targetVals.PCWP = 8; % not sex-specific
inputVals.CVP = 4; 
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

mods = {'k_pas_LV','k_pas_RV','k_act_LV','k_act_RV','C_SA','C_PA','R_SA','R_PA','R_atria','R_m_o','Vw_LV','Vw_RV','LvSepR','R_tPA','R_tSA','K_P','B_P'};
