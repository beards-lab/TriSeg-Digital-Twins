%% canonical subject simulation
clear
for GENDER  = 1 % 1 for male, 2 for female
    if GENDER == 1
        [targets, inputs, mods] = targetVals_male();
        modifiers = readmatrix('modifiers_male.csv','range','2:2');
        [~,params, init] = estimParams(targets,inputs,mods,modifiers);
    else
        [targets, inputs, mods] = targetVals_female();
        modifiers = readmatrix('modifiers_female.csv','range','2:2');
        [~,params, init] = estimParams(targets,inputs,mods,modifiers);
    end
    Sim = 1; % 1 for run sim 0 for opt
    if Sim ==1
        runSim;
        NplotSrd;
        % GetMovie;
        % See_TriSeg;
    else
        Srdopt;
    end
end
%% Show patients' example
clear
load AllPatients.mat
secondslot = [17 79 206 256 288 325 352 355 360 361]; % patients with 2 3-month windows
% Shuntlist = [34 41 54 61 83 116 183 231 268 278 312]; % patients with shunt 
for PatID = 3 % any patient NO between 1 and 370
    for ModelWin =  1
        % try
        [Windowdate,targets, inputs, mods] = targetVals_HF(patients,PatID,ModelWin);
        RUNOPT = 0;% 1 for opt, 0 for Sim
        if RUNOPT ==1
            HFopt;
        end
        if RUNOPT == 0
            load(sprintf('Sims/P_NO%dWindow%d.mat',PatID,ModelWin));
            modifiers = output.modifiers;
            % Calculate BSA using Mosteller formula
            [Error, params, init] = estimParams(targets,inputs,mods, modifiers);
            runSim;
            Print_cost = 1;% 1 for performance, other for patient's pre-condition
            NplotFit;
            GetMovie;
            See_TriSeg;
        end
    end
end
