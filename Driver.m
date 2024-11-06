%% Script Summary:
% This script calls the functions: targetVal.m and estimParams.m, as well as the scripts: XXopt.m
% and runSim.m.
% The first section generates simulations for canonical male and female subjects, and the
% second section generates simulations for heart failure patients.
% Created by Feng Gu
% Last modified: 10/29/2024

%% Simulation for Canonical Subjects
clear
for GENDER  = 1 % 1 for male and other for female
    if GENDER == 1
        [targets, inputs, mods] = targetVals_male();
        modifiers = readmatrix('modifiers_male.csv','range','2:2');
        [~,params, init] = estimParams(targets,inputs,mods,modifiers);
    else
        [targets, inputs, mods] = targetVals_female();
        modifiers = readmatrix('modifiers_female.csv','range','2:2');
        [~,params, init] = estimParams(targets,inputs,mods,modifiers);
    end
    Sim = 1; % 1 for simulation and other for optimazation
    if Sim ==1
        runSim; 
        NplotSrd; % 4-panel figure just for canonical subjects
        GetMovie; % TriSeg model: displacement and stress as functions of time
        See_TriSeg; % slices of GetMoive.m
    else
        Srdopt; % optimize modifiers of canonical subjects
    end
end
%% Simulation for heart failure patients
clear
load AllPatients.mat
secondslot = [17 79 206 256 288 325 352 355 360 361]; % patients with 2 3-month windows
% Shuntlist = [34 41 54 61 83 116 183 231 268 278 312]; % patients with shunt 
for PatID = 192 % any number between 1 and 370, example patient in paper is 192
    for ModelWin =  1
        [Windowdate,targets, inputs, mods] = targetVals_HF(patients,PatID,ModelWin);
        RUNOPT = 0; % 0 for simulation and 1 for optimazation
        if RUNOPT ==1
            HFopt; % optimize modifiers of HF patients
        end
        if RUNOPT == 0
            load(sprintf('Sims/P_NO%dWindow%d.mat',PatID,ModelWin)); % load results from great lake clusters
            modifiers = output.modifiers;
            [Error, params, init] = estimParams(targets,inputs,mods, modifiers);
            runSim;
            Print_cost = 1;% 1 for performance, other for patient's pre-condition (PHI sensitive)
            NplotFit; % 6-panel figure for HF patients
            GetMovie;
            See_TriSeg;
        end
    end
end
