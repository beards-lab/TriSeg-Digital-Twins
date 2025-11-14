%% Simulation for heart failure patients from Umich cohort without CMR info
clear
load AllPatients.mat
Order = [1	2 13 15	16	17	20	22	25	26	35	37	42	45	46	48	50	55	57	60	62	66	67	77	78	79	84	87	88	90	91	92	93	94	97	98	99	102	103	104	106	108	111	112	113	118	122	126	130	131	132	133	137	141	146	147	150	151	164	167	168	173	177	178	186	188	195	196	197	199	202	203	204	205	206	207	209	212	213	214	217	218	222	223	224	230	232	234	235	236	240	244	246	248	250	252	267	269	271	274	281	282	283	284	293	294	300	302	303	304	306	309	316	318	324	325	327	329	335	341	342	343	344	345	346	350	351	368];
SummaryTable =  table(Order');
SummaryTable.Properties.VariableNames = {'PatientNo'};

for NO = 1:length(Order)
    PATIENT_NO = Order(NO);
    foldername = sprintf('PacingResultsUM/PatientNO%d', PATIENT_NO);
    Newfilename = sprintf('%s/PacingTable.mat', foldername);
    load(Newfilename)
    SummaryTable.PacingResults{NO} = TableatdifferentHR;
end

%%
save pacingSUM.mat SummaryTable

%% Simulation for heart failure patients from UW cohort without CMR info
clear
load('UWcohort.mat');
Order = [9 14	15	17	18	19	22	25	26	30	33	34	37	44	46	48	52	53];
SummaryTable =  table(Order');
SummaryTable.Properties.VariableNames = {'PatientNo'};
%%
for NO = 1:length(Order)
    PATIENT_NO = Order(NO);
    foldername = sprintf('PacingResultsUW/PatientNO%d', PATIENT_NO);
    Newfilename = sprintf('%s/PacingTable.mat', foldername);
    load(Newfilename)
    SummaryTable.PacingResults{NO} = TableatdifferentHR;
end
save pacingSUMUW.mat SummaryTable

%%
clear
MRI_flag = 1;
for GENDER  = 1 % 1 for male and other for female
    if GENDER == 1
        [targets, inputs, mods] = targetVals_male();
        modifiers = readmatrix('modifiers_male.csv','range','2:2');
        % modifiers = ones(1,length(mods));
        [INIparams, INIinit] = estiminiParams(targets,inputs);
        [params, init] = optParams(INIparams, INIinit, mods,modifiers);
        QS2 = 545.17-2.117*params.HR;
    else
        [targets, inputs, mods] = targetVals_female();
        modifiers = readmatrix('modifiers_female.csv','range','2:2');
        [INIparams, INIinit] = estiminiParams(targets,inputs);
        [params, init] = optParams(INIparams, INIinit, mods,modifiers);
        QS2 = 546.5-2.0*params.HR;
    end
    NewHR = params.HR+(0:5:30);% I am just using this as an example
    MRI_flag = 1;
    IVRT = 70*75/params.HR;
    ActT = QS2+IVRT;  
    runSimonFakeDT;
    NewMAPtarget = o_vals.MAP;
    paramsbaseline = params;
    initbaseline = init;
    for HRindex = 1:length(NewHR)       
        PacingOPT;
    end
end