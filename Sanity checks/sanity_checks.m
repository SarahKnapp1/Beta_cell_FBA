%Sarah Knapp
%last edited at december 2022
close all
%sanity checks 
%% Content:
% The tests include:
% 
% * leak test
% * production of protons from nothing as well as from water, and/or oxygen 
% alone
% * production of matter when atp hydrolysis reaction is allowed to work but 
% all uptakes are closed
% % * duplicated reactions
% * empty colmuns in the model.rxnGeneMat
% * the single gene deletion analysis runs smoothly
% * ATP yield from different carbon sources
% * metabolic objective functions
% * demand reactions with negative lower bound (should not occur based on definition 
% of demand reactions)
% * Test if ATP is product from glucose

%%


%if necessary, initialize the cobra toolbox:
%initCobraToolbox(false) %false, as we don't want to update
%% 
% For solving linear programming problems in FBA analysis, certain solvers 
% are required:
%cplex interface added to MATLAB path.

changeCobraSolver('ibm_cplex')

%load model
model = load('RedHumanSarah.mat');

model=(model.RedModel);
%%

%EX_glc(e)- D-glucose exchange - Two monosaccharide isomeric enantiomers are called glucose,
%and only one of them, the D-glucose enantiomer, is biologically active.
%1) defining the allowed carbon sources
% %2) defining the oxygen uptake level,allow unlimited oxygen
%%
% Replace reaction abbreviation for the ATP hydrolysis (DM_atp_c_) and Biomass 
% reaction used differently in various models.
model.rxns(find(ismember(model.rxns,'ATPM')))={'DM_atp_c_'};
model.rxns(find(ismember(model.rxns,'ATPhyd')))={'DM_atp_c_'};
model.rxns(find(ismember(model.rxns,'DM_atp(c)')))={'DM_atp_c_'};
model.rxns(find(ismember(model.rxns,'EX_biomass_reaction')))={'biomass_reaction'};
model.rxns(find(ismember(model.rxns,'EX_biomass_maintenance')))={'biomass_maintenance'};
model.rxns(find(ismember(model.rxns,'EX_biomass_maintenance_noTrTr')))={'biomass_maintenance_noTrTr'};
%%
%block all Exchange/demand/sink reactions/biomass LB=o,ub=defult
[model,modelexchanges]=block_uptake_reactions(model);
modelClosed=model;
%set upper bround to 100
%modelClosed.ub(selExc) =0;% 100;
%modelClosed.ub(find(ismember(modelClosed.rxns, modelClosed.rxns(modelexchanges))))=0;

modelClosedOri = modelClosed;
%% Start with tests.
% Define some parameters that we will need.

cnt = 1;
tol = 1e-6;
%%
%Perform leak test, i.e., whether the closed modelcan produce any 
%exchanged metabolite, as defined in the model, from nothing.

modelClosed = modelClosedOri;
[LeakRxns,modelTested,LeakRxnsFluxVector] = fastLeakTest(modelClosed,modelClosed.rxns(modelexchanges),'false');
TableChecks{cnt,1} = 'fastLeakTest 1';
if length(LeakRxns)>0
warning('model leaks metabolites!')
TableChecks{cnt,2} = 'Model leaks metabolites!';
else
TableChecks{cnt,2} = 'Leak free!';
end
cnt = cnt + 1;
%% 
%if something leaks when demand reactions for each metabolite in the model are added.
modelClosed = modelClosedOri;
%As these measures may not identify all exchange and sink reactions in a model,so grab all reactions based on stoichiomettry.
%identify all reactions that contain only one non-zero entry in the S matrix (column).
selExc = (find(full((sum(abs(modelClosed.S)==1, 1)==1) & (sum(modelClosed.S~=0) == 1))))';
[LeakRxnsDM,modelTestedDM,LeakRxnsFluxVectorDM] = fastLeakTest(modelClosed,modelClosed.rxns(selExc),'true');
TableChecks{cnt,1} = 'fastLeakTest 2 - add demand reactions for each metabolite in the model';
if length(LeakRxnsDM)>0
TableChecks{cnt,2} = 'Model leaks metabolites when demand reactions are added!';
else
TableChecks{cnt,2} = 'Leak free when demand reactions are added!';
end
cnt = cnt + 1;
%%
%Test if the model produces energy from water!
modelClosed = modelClosedOri;
tol = 1e-6;
modelClosedATP = changeObjective(modelClosed,'DM_atp_c_');
modelClosedATP = changeRxnBounds(modelClosedATP,'DM_atp_c_',0,'l'); 
modelClosedATP = changeRxnBounds(modelClosedATP,'EX_h2o_e',-1,'l');
FBA3=optimizeCbModel(modelClosedATP,'max',0);
TableChecks{cnt,1} = 'Exchanges, sinks, and demands have lb = 0, except h2o';
if abs(FBA3.f) > 1e-6
TableChecks{cnt,2} = 'model produces energy from water!';
else
TableChecks{cnt,2} = 'model DOES NOT produce energy from water!';
end
cnt = cnt + 1;
%%
%Test if the model produces energy from water and oxygen!
modelClosed = modelClosedOri;
modelClosedATP = changeObjective(modelClosed,'DM_atp_c_');
modelClosedATP = changeRxnBounds(modelClosedATP,'DM_atp_c_',0,'l');
modelClosedATP = changeRxnBounds(modelClosedATP,'EX_h2o_e',-50,'l');
modelClosedATP = changeRxnBounds(modelClosedATP,'EX_o2_e',-50,'l');
FBA6=optimizeCbModel(modelClosedATP,'max',0);
TableChecks{cnt,1} = 'Exchanges, sinks, and demands have lb = 0, except h2o and o2';
if abs(FBA6.f) > 1e-6
TableChecks{cnt,2} = 'model produces energy from water and oxygen!';
else
TableChecks{cnt,2} = 'model DOES NOT produce energy from water and oxygen!';
end
cnt = cnt + 1;
%%
% Test if the model produces matter when atp demand is reversed!

modelClosed = modelClosedOri;
modelClosed = changeObjective(modelClosed,'DM_atp_c_');
modelClosed.lb(find(ismember(modelClosed.rxns,'DM_atp_c_'))) = -100;
FBA = optimizeCbModel(modelClosed);
TableChecks{cnt,1} = 'Exchanges, sinks, and demands have  lb = 0, allow DM_atp_c_ to be reversible';
if abs(FBA.f) > 1e-6
    TableChecks{cnt,2} = 'model produces matter when atp demand is reversed!';
else
    TableChecks{cnt,2} = 'model DOES NOT produce matter when atp demand is reversed!';
end
cnt = cnt + 1;
%% 
% Test if the model has flux through h[m] demand !

modelClosed = modelClosedOri;
modelClosed = addDemandReaction(modelClosed,'h_m');
modelClosed = changeObjective(modelClosed,'DM_h_m');
modelClosed.ub(find(ismember(modelClosed.rxns,'DM_h_m'))) = 100;
FBA = optimizeCbModel(modelClosed,'max');
TableChecks{cnt,1} = 'Exchanges, sinks, and demands have  lb = 0, test flux through DM_h[m] (max)';
if abs(FBA.f) > 1e-6
    TableChecks{cnt,2} = 'model has flux through h[m] demand (max)!';
else
    TableChecks{cnt,2} = 'model has NO flux through h[m] demand (max)!';
end
cnt = cnt + 1;
%% 
% Test if the  model has flux through h[c] demand !

modelClosed = modelClosedOri;
modelClosed = addDemandReaction(modelClosed,'h_c');
modelClosed = changeObjective(modelClosed,'DM_h_c');
modelClosed.ub(find(ismember(modelClosed.rxns,'DM_h_c'))) = 100;
FBA = optimizeCbModel(modelClosed,'max');
TableChecks{cnt,1} = 'Exchanges, sinks, and demands have  lb = 0, test flux through DM_h[c] (max)';
if abs(FBA.f) > 1e-6
    TableChecks{cnt,2} = 'model has flux through h[c] demand (max)!';
else
    TableChecks{cnt,2} = 'model has NO flux through h[c] demand (max)!';
end
cnt = cnt + 1;

%%
% Compute ATP yield. This test is identical to the material covered in the 
% tutorial testModelATPYield. to test if the correct ATP yield from different carbon sources can be realized by the model.

TableChecks{cnt,1} = 'Compute ATP yield';
if 1 % test ATP yield
    [Table_csources, TestedRxns, Perc] = testATPYieldFromCsources(model);
    TableChecks{cnt,2} = 'Done. See variable Table_csources for results.';
else
    TableChecks{cnt,2} = 'Not performed.';
end
cnt = cnt + 1;
%% 
% Check for duplicated reactions in the model.

TableChecks{cnt,1} = 'Check duplicated reactions';
method='S';
removeFlag=0;
[modelOut,removedRxnInd, keptRxnInd] = checkDuplicateRxn(model,method,removeFlag,0);
if isempty(removedRxnInd)
    TableChecks{cnt,2} = 'No duplicated reactions in model.';
else
    TableChecks{cnt,2} = 'Duplicated reactions in model.';
end
cnt = cnt + 1;


%%
% Check empty columns in 'model.rxnGeneMat'.
TableChecks{cnt,1} = 'Check empty columns in rxnGeneMat';
E = find(sum(model.rxnGeneMat)==0);
if isempty(E)
    TableChecks{cnt,2} = 'No empty columns in rxnGeneMat.';
else
    TableChecks{cnt,2} = 'Empty columns in rxnGeneMat.';
end
cnt = cnt + 1;

%% 
% Check that demand reactions have a lb >= 0.

TableChecks{cnt,1} = 'Check that demand reactions have a lb >= 0';
DMlb = find(model.lb(strmatch('DM_',model.rxns))<0);
if isempty(DMlb)
    TableChecks{cnt,2} = 'No demand reaction can have flux in backward direction.';
else
    TableChecks{cnt,2} = 'Demand reaction can have flux in backward direction.';
end
cnt = cnt + 1;
%% 
% Check whether singleGeneDeletion runs smoothly.

TableChecks{cnt,1} = 'Check whether singleGeneDeletion runs smoothly';
try
    [grRatio,grRateKO,grRateWT,hasEffect,delRxns,fluxSolution] = singleGeneDeletion(model);
    TableChecks{cnt,2} = 'singleGeneDeletion finished without problems.';
catch
    TableChecks{cnt,2} = 'There are problems with singleGeneDeletion.';
end
cnt = cnt + 1;

%% 
% Display all results.

TableChecks
%% 
%% 
% Save all results.

resultsFileName = 'TestResults';
save(strcat(resultsFileName,'.mat'));
%%

model=block_uptake_reactions(model);

%%
%Test if ATP is product from glucose
% allow water uptake
model.lb(find(ismember(model.rxns,'EX_h2o_e'))) =-100;

model= changeObjective(model,'ATPS4m');
FBAsolution = optimizeCbModel(model,'max',0);

if FBAsolution.stat == 1
disp('Relaxed model is feasible');
end
% set both the upper and lower bounds of the glucose exchange reaction to values between 0 and -20 mmol gDW-1hr-1
% use optimize the Model to perform FBA with each uptake rate.
objFlux = zeros(201,1);
flux_EX_glc_e=zeros(201,1);
for i=0:200
    N=i*0.05
    model = changeRxnBounds(model,'EX_glc_e',-N,'l');% 
    FBAsolution = optimizeCbModel(model,'max',0);%optimize model

    flux_EX_glc_e(i+1)=-FBAsolution.x(find(ismember(model.rxns,'EX_glc_e')));


    objFlux(i+1) = FBAsolution.x(find(ismember(model.rxns,'ATPS4m')));%FBAsolution.f;

end

    figure
    controlFlux=[0:0.1:20];
    objFlux = objFlux';
    length(objFlux);
    plot(flux_EX_glc_e,objFlux,'c','LineWidth',5)
    xlabel(['Glucose flux itself [mmol/(gDW*h)]']);
    axis tight;
    ylabel(['ATP synthesis[mmol/(gDW*h)]']);



%%
% %%
% [P, C, vP, vC] = computeFluxSplits(model, {'adp_c'}, FBAsolution.x);
% total_adp_flux = sum(vP)
% 
% [P, C, vP, vC] = computeFluxSplits(model, {'atp_c'}, FBAsolution.x);
% total_aTp_flux = sum(vP)

