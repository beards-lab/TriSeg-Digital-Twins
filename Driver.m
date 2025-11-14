%% Script Summary:
% This script calls the functions: targetVal.m and estimParams.m, as well as the scripts: XXopt.m
% and runSim.m.
% The first section generates simulations for canonical male and female subjects, and the
% second section generates simulations for heart failure patients.

% Created by Feng Gu
% 11/14/2025 current script can fit patients with/without CMR

%% Simulation for Canonical Subjects
clear
MRI_flag = 1; % Canonical can only have with CMR
for GENDER  = 2 % 1 for male and other for female
    if GENDER == 1
        [targets, inputs, mods] = targetVals_male();
        modifiers = readmatrix('modifiers_male.csv','range','2:2');
        [INIparams, INIinit] = estiminiParams(targets,inputs);
    else
        [targets, inputs, mods] = targetVals_female();
        modifiers = readmatrix('modifiers_female.csv','range','2:2');
        [INIparams, INIinit] = estiminiParams(targets,inputs);
    end
    Sim = 1; % 1 for simulation and other for optimazation
    if Sim ==1
        [params, init] = optParams(INIparams, INIinit, mods,modifiers);
        runSimOnGL;
        NplotSrd; % 4-panel figure just for canonical subjects
        GetMovie; % TriSeg model: displacement and stress as functions of time
        See_TriSeg; % slices of GetMoive.m       
    else
        Srdopt; % optimize modifiers of canonical subjects
    end
end


%% Simulation for heart failure patients from Umich cohort with/without CMR info

load AllPatients.mat
secondslot = [17 79 206 256 288 325 352 355 360 361]; % patients with 2 3-month windows
% HFpEFOrder = [1	2 13 15	16	17	20	22	25	26	35	37	42	45	46	48	50	55	57	60	62	66	67	77	78	79	84	87	88	90	91	92	93	94	97	98	99	102	103	104	106	108	111	112	113	118	122	126	130	131	132	133	137	141	146	147	150	151	164	167	168	173	177	178	186	188	195	196	197	199	202	203	204	205	206	207	209	212	213	214	217	218	222	223	224	230	232	234	235	236	240	244	246	248	250	252	267	269	271	274	281	282	283	284	293	294	300	302	303	304	306	309	316	318	324	325	327	329	335	341	342	343	344	345	346	350	351	368];

for PATIENT_NO =192 % any number between 1 and 370, example patient in paper is 192
    try
        for ModelWin =  1
            MRI_flag = 1; % 1:with, 0:without
            [Windowdate,targets, inputs, mods] = targetVals_HF(patients,PATIENT_NO,ModelWin,MRI_flag);
            [INIparams, INIinit] = estiminiParams(targets,inputs);
            RUNOPT = 0; % 0 for simulation and 1 for optimazation
            if RUNOPT ==1
                HFopt; % optimize modifiers of HF patients, Current HFopt is still used for full model
            end
            if RUNOPT == 0
                load(sprintf('SimsUMFinal/P_NO%d.mat',PATIENT_NO)); % load results from great lake clusters
                targets = output.targets;
                inputs = output.inputs;
                params = output.params;
                init = output.init;
                % if you want to run this function yourself... otherwise I
                % have saved everything
                % modifiers = output.m; 
                % [params, init] = optParams(INIparams, INIinit, mods,modifiers);
                runSim;
            end
        end
        Print_cost = 1;% 1 for performance, other for patient's pre-condition (PHI sensitive)
        NplotFit; % 6-panel figure for HF patients
        % filename = fullfile(save_dir, sprintf('Patient_%02d.png', PATIENT_NO));
        % exportgraphics(gcf, filename, 'Resolution', 300);  % gcf = current figure
        % GetMovie;
        See_TriSeg;
    catch ME
        disp(['Error: ', ME.message, 'PatID',num2str(PATIENT_NO)])
    end
end
