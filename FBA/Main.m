%Sarah knapp last edited at december 2022

%adjusted the model to better describe the beta cell in terms of reactions and objective function.
%GEM simulations. 

clear 
close all

%CPLEX interface added to MATLAB path.
addpath('C:\Program Files\IBM\ILOG\CPLEX_Enterprise_Server1210\CPLEX_Studio\cplex\matlab\x64_win64')
savepath
changeCobraSolver('ibm_cplex','LP')
changeCobraSolver ('ibm_cplex', 'QP');  
%changeCobraSolver('ibm_cplex','all')

%%
%load model
model = load('RedHumanSarah.mat');


%%
model=(model.RedModel);

model.rxns=(strrep(model.rxns,'(','_'));
model.rxns=(strrep(model.rxns,')',''));
model.mets=(strrep(model.mets,'[','_'));
model.mets=(strrep(model.mets,']',''));
%%

modelgg=model;

%extract a model which corresponds to the genes in beta cell by removing
%the reactions with genes not expperes inbeta cells
% reactions_remove = load('gene_remove_beta_cell.mat');
% reactions_remove=reactions_remove.RXDremove;
%model= removeRxns(model,reactions_remove);
% 
%model = changeRxnBounds(model,reactions_remove,0,'b');

%%
%Adding insulin synthesis reaction 

%insulin formula
proINS_formula='C410H644N114O127S6';
%add metabolite of pre pro insulin
model=addMetabolite(model,'INS_pro_c','metFormula',proINS_formula,'KeggId','-');


% Exchange/demand reaction of INS_pro[c]
model=addReaction(model,'DM_INS_pro','reactionFormula','INS_pro_c <=>');
%The lower bound of the proinsulin exchange reaction is now 0,
%so proinsulin may not enter the system.
model = changeRxnBounds(model,'DM_INS_pro',0,'l');
model = changeRxnBounds(model,'DM_INS_pro',100,'u');
%Insulin synthesis
rxnID=['4 ala_L_c + 4 arg_L_c + 3 asn_L_c + 1 asp_L_c + 6 cys_L_c + 7 gln_L_c + 8 glu_L_c + 11 gly_c + 2 his_L_c + '...
     '2 ile_L_c + 12 leu_L_c + 2 lys_L_c + 3 phe_L_c + 3 pro_L_c + 5 ser_L_c + 3 thr_L_c'... 
     ' + 4 tyr_L_c + 6 val_L_c + 370.316 atp_c + 285.316 h2o_c <=> 370.316 adp_c + 370.316 pi_c + 6 h_c + INS_pro_c'];

 
%add reaction to model
model =addReaction(model,'INSsyn','reactionFormula',rxnID);
model = changeRxnBounds(model,'INSsyn',0,'l');
model = changeRxnBounds(model,'INSsyn',100,'u');
model.subSystems(end-1:end,:)=[];


   
% % Exchange/demand reaction of INS_sec
% model=addReaction(model,'DM_INS_sec','reactionFormula','M03169[c] <=>');
% %The lower bound of the insulin secration exchange reaction is now 0,
% % so insulin may not enter the system.
% model = changeRxnBounds(model,'DM_INS_sec',0,'l');
% 
% %Insulin secretion
% rxnID_sec='INS_pro[c] -> M03169[c]';
% model=addReaction(model,'INSsec','reactionFormula',rxnID_sec);

%%
%Adjusting the uptake medium of the model 

%block all Exchange/demand/sink reactions lb=o,ub=defult
model=block_uptake_reactions(model);

%%
% %1) defining the oxygen uptake level https://www.nature.com/articles/nbt0201_125
model.lb(find(ismember(model.rxns,'EX_o2_e'))) =-2;
% allow water uptake
model.lb(find(ismember(model.rxns,'EX_h2o_e'))) =-50;

% 
% model.lb(find(ismember(model.rxns,'EX_nh4_e'))) = -100; %
% model.lb(find(ismember(model.rxns,'EX_so4_e'))) = -70;%Exchange of Sulfate
% model.lb(find(ismember(model.rxns,'EX_pi_e'))) = -100;


%%
%2) defining the allowed carbon sources Glucose.
model= changeRxnBounds(model,'EX_glc_e',-10,'l');
%%
N=0;

%with out the essential amino acids there is not flux in inssyn reaction even if one
%of them are missing.

%change bounds of essential amino acid 
model = changeRxnBounds_essential_AA(model,N);

%model= changeRxnBounds(model,'SUCOAS1m',0,'b');
% model= changeRxnBounds(model,'PYRt2r',0,'b');
%model= changeRxnBounds(model,'r0062',0,'b');
% model= changeRxnBounds(model,'ME2m',0,'l');
% model= changeRxnBounds(model,'ME1m',0,'l');
% model= changeRxnBounds(model,'PEPCKm',0,'l');


%nonessential amino acids 
%alanine, arginine, asparagine, aspartic acid, cysteine, glutamic acid, glutamine, glycine, proline, serine, and tyrosine.
%model= changeRxnBounds(model,'EX_ala_L_e',N,'l');
%model= changeRxnBounds(model,'EX_arg_L_e',N,'l');
% model= changeRxnBounds(model,'EX_asp_L_e',N,'l');
% model= changeRxnBounds(model,'EX_asn_L_e',N,'l');
% model= changeRxnBounds(model,'EX_cys_L_e',N,'l');
% model= changeRxnBounds(model,'EX_glu_L_e',N,'l');
% model= changeRxnBounds(model,'EX_gln_L_e',N,'l');
% model= changeRxnBounds(model,'EX_gly_e',N,'l');
% model= changeRxnBounds(model,'EX_pro_L_e',N,'l');
% model= changeRxnBounds(model,'EX_ser_L_e',N,'l');
% model= changeRxnBounds(model,'EX_tyr_L_e',N,'l');
%%
%limited capacity to generate lactate.

model= changeRxnBounds(model,'LDH_L',5,'u');
%model= changeRxnBounds(model,'LDH_Lm',5,'u');
model= changeRxnBounds(model,'r0173',5,'u');
model= changeRxnBounds(model,'L_LACDcm',5,'u');


%%
% Test if the model can preduce insulin and ATP.

% change Objective function
model = changeObjective(model,{'ATPS4m','INSsyn'},[0.001,0.999]);

%optimizetion
FBAsolution = optimizeCbModel(model,'max',0);

%checking the fluxes through interesting reactions 
flux_PCm=FBAsolution.x(find(ismember(model.rxns,'PCm')))
flux_PDH=FBAsolution.x(find(ismember(model.rxns,'PDHm')))
objFlux =FBAsolution.x(find(ismember(model.rxns,'INSsyn')))
flux_ATP=FBAsolution.x(find(ismember(model.rxns,'ATPS4m')))

%%

%change bounds of essential amino acid 
model = changeRxnBounds_essential_AA(model,-1);

%robustness analysis on the glucose exchange reaction
plot_objective_as_function_glu(model,'INSsyn','EX_glc_e')


%%
% The function analyzies reactions in the TCA from insulin synthesis.
% and modeles hyperinsulinemia by increasing the uptake of EAA.
plot_objective_as_function_EAA_interesting_reactions_flux(model);

%%
% The function checks if when the demand for insulin synthesis increases
% a decrease in mitochondrial ATP production occurs in differner glucose
% input fluxs.
hyperglycemia = hyperglycemia_ATP(model)
%%
% Looking at the change in extracellular glucose levels when insulin synthesis increases
hyperglycemia=ins_as_function_of_glucose_fasting(model)



%%
%find the fluxs of the interested reactions when optimizing the 
%objective function in norminsulinmia and hyperinsulinemia.
figure
mone=0;
met_flux=[];
model = changeRxnBounds_essential_AA(model,-0.3)
model= changeRxnBounds(model,'EX_glc_e',-10,'l');
for index=1:2
    mone=mone+1;
    INS_syn_opt_matrix=flux_met(model);
    X = categorical(INS_syn_opt_matrix(1,:));
    X = reordercats(X,INS_syn_opt_matrix(1,:));
    met_flux=[met_flux;cell2mat([INS_syn_opt_matrix(2,:)])]
    %change essential AA bounds
    model = changeRxnBounds_essential_AA(model,-2)
    % side by  side
end

hb=bar(X,met_flux');%, legend( 'Norminsulinemia' ,'Hyperinsulinemia');

% just for the ratio between PCm/PDHm
% X = categorical({'PCm/PDHm'});
% X = reordercats(X,{'PCm/PDHm'});
% met_flux=[met_flux(1,1)/met_flux(1,2),met_flux(2,1)/met_flux(2,2)]; 
% hb=bar(X,met_flux);

hb(1).FaceColor = [98/255,197/255,218/255];
hb(2).FaceColor = [254/255,125/255,104/255];            
%bar([Main_Result', Simulated_Result']), legend(X);
ylabel('Flux [mmol/gDW*h]')
set(gcf,'color','w');
axis off




%%
%perform FVA and FBA analysis on your intrest reaction and ploting
%the variations in fluxes in the intreat reaction.

rxns = {'PDHm'};%;'ATPS4m';'PCm';'PDHm';'LDH_L';'CSm';'MDHm';'r0062';'INSsyn'};
FVA_FBA_reaction(model,rxns);

%%
%The function performes FVA and FBA on the reactions within the TCA
%cycle in norminsulinemia and hyperinsulinemia
FVA_TCA(model)


%%   
%Checks if carbon come from glucose or EAA or both in order to
%synthezie insulin 
Max_Fraction_EAA_ins=fraction_of_carbon_EAA_to_ins_syn(model) 

%%

%The function checked where glucose causes insulin synthesis when there is a change in EAA
%and showing the results as an heatmap.
matrix_INSsyn=Glucose_ins_syn_change(model)

%%

%generation of a heatmap of the fluxes of the reactions the produce the nonessential
%amino acids from the TCA cycle.
matric_AA_reactions=amino_acid_from_TCA(model);


%%
%check the flux which enters and leaves the TCA cycle  
Flux_in_out_TCA=flux_TCA(model)

%%
%Checks in which reaction the flux changes in hyperinsulinemia as opposed to norm-insulinom
changed_reactions=flux_change_in_hyperinsulinmia(model);



%%%%visualization
%%
%change objective function
model=changeObjective(model,{'ATPS4m','INSsyn'},[0.001,0.999]);

%optimize model
FBAsolution=optimizeCbModel(model,'max',0);

% The function prints an overlay in the minerva map based on the fluxs of
% the optimization
minerva_visualization(model,FBAsolution)

%%
%visualization of the flux in the TCA cycle by the CellDesigner
viewflux(model)















