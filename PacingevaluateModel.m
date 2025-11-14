function Pacing_cost = PacingevaluateModel(mnew,NewHR,PATIENT_NO,MRI_flag,NewMAPtarget,params,init,ActT,HRindex)
%% Function Purpose:
% This function is used to compute a cost function in runSim.m with new modifiers for patients
% during optimization.

% Created by Andrew Meyer and Feng Gu
% Last modified: 10/29/2024

try
    % Evaluates the model during optimalization
    params = estimPacingParams(params,mnew,NewHR,ActT);
    params_no = struct2array(params);
    init_no = struct2array(init);
    if any(params_no <= 0)||any(init_no(2:end) <= 0)
        error('bad inital guessing')
    end
    MRI_flag = MRI_flag;
    runSimonFakeDT;
    Pacing_cost = abs(o_vals.MAP - NewMAPtarget);
    fprintf('HR = %.2f, MAP = %.2f, Target = %.2f, Cost = %.4f\n', NewHR, o_vals.MAP, NewMAPtarget, Pacing_cost);
catch
   Pacing_cost = inf;
end
end