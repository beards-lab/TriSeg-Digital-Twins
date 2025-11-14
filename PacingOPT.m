%% New Optimization Block
% If HRindex > 1, perform pacing parameter optimization.
% You may later implement: if cost does not decrease after ~100 runs, stop optimization.
if HRindex > 1
    % Evaluate baseline pacing cost at current HR
    Pacing_cost = PacingevaluateModel( ...
        1, NewHR(HRindex), GENDER, MRI_flag, NewMAPtarget, ...
        paramsbaseline, initbaseline, ActT, HRindex);

    %% fminsearch options (local optimizer)
    options = optimset( ...
        'Display', 'iter', ...         % Show iteration details
        'PlotFcns', @optimplotfval, ...% Plot objective value
        'TolFun', 1e-4, ...            % Function tolerance
        'TolX', 1e-3, ...              % Parameter tolerance
        'MaxIter', 50);                % Max iterations (adjust as needed)

    %% Genetic Algorithm (GA) options — global search
    % Use GA to obtain a good initial guess for fminsearch.
    initPop = 0.1;                     % Initial seed for population
    popSize = 50;                      % Population size
    InitialPopulationMatrix = initPop + rand(popSize, 1) * 2;
    InitialPopulationMatrix(1) = 1;    % Ensure first candidate = 1

    gaOptions = optimoptions('ga', ...
        'Display', 'iter', ...         % Show GA progress
        'MaxGenerations', 20, ...      % Number of GA generations
        'PopulationSize', popSize, ... % Population size
        'InitialPopulationMatrix', InitialPopulationMatrix, ...
        'UseParallel', true, ...       % Enable parallel computing
        'MaxStallGenerations', 5);     % Early stopping: terminate if no improvement for 5 generations

    nVars = 1;  % Number of optimization variables (update if model dimension changes)

    %% Global optimization using GA
    [mnew_ga, fval_ga] = ga( ...
        @(m) PacingevaluateModel(m, NewHR(HRindex), GENDER, MRI_flag, ...
                                 NewMAPtarget, paramsbaseline, initbaseline, ActT, HRindex), ...
        nVars, [], [], [], [], 0.1, 10, [], gaOptions);

    %% Local refinement using fminsearch
    [mnew_final, fval_final] = fminsearch( ...
        @(m) PacingevaluateModel(m, NewHR(HRindex), GENDER, MRI_flag, ...
                                 NewMAPtarget, paramsbaseline, initbaseline, ActT, HRindex), ...
        mnew_ga, options);

    %% Update parameters and run simulation
    params = estimPacingParams(paramsbaseline, mnew_final, NewHR(HRindex), ActT);
    init = initbaseline;
    runSimonFakeDT;
end

%% Output structure
output.params = params;
output.init   = init;
output.o_vals = o_vals;

%% Save results
foldername = sprintf('PacingSimsC/PatientNO%d', GENDER);
Newfilename = sprintf('%s/HR%d.mat', foldername, NewHR(HRindex));

% Create folder if it does not exist
if ~exist(foldername, 'dir')
    mkdir(foldername);
end

% Save output
save(Newfilename, "output");
