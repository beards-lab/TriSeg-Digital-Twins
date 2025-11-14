# TriSeg-Digital-Twins
This repository supports the manuscript "Identification of Digital Twins to Guide Interpretable AI for Diagnosis and Prognosis in Heart Failure" by Feng et al. It contains the identification framework and necessary files for creating digital twins for heart failure patients.

## 1. File Descriptions

### i) Source Data
- **AllPatients.mat**: A struct containing de-identified clinical measurements of example patients.
- **Sims folder**: Stores optimization results from supercomputer nodes.

### ii) Functions
- **Driver.m**: Runs a patient simulation. Several boolean options control execution, such as running pre-saved simulations, assigning modifier values, and printing/plotting options.
- **targetVals____.m**: A series of "targetVals" functions that return target, input structs, and an array of modifier values for adjustable parameters. Non-void versions (e.g., `targetVals_HF.m`) use patient numbers to assign targets and parameters and validate patient clinical values.
- **estimParams.m**: Accepts targets and inputs for a patient, returning a parameter struct and initial conditions for the model’s ODE solver, estimating parameters to replicate clinical measurements.
- **calc_xm_ym.m** and **gemo_0.m**: Functions that estimate heart geometry metrics. `gemo_0.m` calculates wall and lumen volumes, while `calc_xm_ym.m` computes sarcomere reference length and solves for xm and ym values for consistent initial conditions.
- **dXdTDAE.m and dXdTode.m**: Defines the DAE or ODE system governing the TriSeg model. Inputs include time, state variables, and model parameters, and outputs constraints and state variable differentials for the MATLAB ODE solver (`ode15s`).
- **runSim.m**: Runs the model using estimated parameters and initial conditions. It calls the ODE solver to simulate the cardiac cycle, adjusting ODE options as necessary for performance.
- **HF_opt.m** and **Srd_opt.m**: Used for parameter estimation, incorporating Genetic Algorithm (GA), Pattern Search, and `fminsearch`.
- **evaluateModel.m**: Accepts a patient number and modifier vector, returning the cost of one simulation. Used by optimizers as the objective function.
- **NplotFit.m**, **GetMovie.m**, and **See_TriSeg.m**: For graphical outputs of simulations and TriSeg geometry after running `runSim`.

## 2. User Guide

### i) Running a Patient
1. Choose a patient number and pass it to `targetVals_HF.m` to obtain target, input, and adjustable parameter data.
2. Feed targets and inputs into `estimParams.m` to get model parameters and initial conditions.
3. Optionally, adjust the parameter structure for optimization purposes.
4. Execute `runSim.m` with `print_sim` enabled for target comparison.
5. For graphical outputs, run `NplotFit`, `GetMovie`, or `See_TriSeg`.

### ii) Optimizing a Patient
1. To fit a patient’s targets, use **HF_opt.m** and choose an optimization algorithm (e.g., global optimizers like GA followed by gradient-based methods like `fminsearch`).
2. Carefully assign adjustable parameters in `targetVals_HF.m`.
3. Use previous modifier vectors as starting points when optimizing.
4. Scripts are preconfigured for optimization runs with minimal setup.

