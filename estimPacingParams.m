function [params] = estimPacingParams(params,mnew,NewHR,ActT)

params.C_SV = params.C_SV*mnew;
params.HR = NewHR;
params.T = 60/params.HR;
params.tau_TS = (ActT*params.HR/6e4)*2/3;
params.tau_TR = (ActT*params.HR/6e4)*1/3;

end