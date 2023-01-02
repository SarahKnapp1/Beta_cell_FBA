%input:model
%output:new model
%block all Exchange/demand/sink reactions LB=o,ub=defult
function  [model,modelexchanges]=block_uptake_reactions(modelClosed)
    %set the lower bound of all exchange and sink reactions to ensure that only those metabolites that are 
    %supposed to be taken up are indded supplied to the model.
    %find all reactions based on their abbreviation
    modelexchanges1 = strmatch('Ex_', modelClosed.rxns);%Exchange reaction
    modelexchanges2 = strmatch('EX_', modelClosed.rxns);
    modelexchanges3 = strmatch('DM_', modelClosed.rxns);%Demand reaction
    modelexchanges4 = strmatch('sink_', modelClosed.rxns);%Exchange/demand reaction in the cytosol like: Sink_his_L[c]
    %Grab the biomass reaction(s) based on the reaction abbreviation.
    BM= (find(~cellfun(@isempty,strfind(lower(modelClosed.rxns),'biom'))));

    %As these measures may not identify all exchange and sink reactions in a model,so grab all reactions based on stoichiomettry.
    %identify all reactions that contain only one non-zero entry in the S matrix (column).
    selExc = (find(full((sum(abs(modelClosed.S)==1, 1)==1) & (sum(modelClosed.S~=0) == 1))))';

    %put all these identified reactions together into one variable 'modelexchanges' and set the lower bound for these reactions to 0.
    modelexchanges = unique([modelexchanges1;modelexchanges2;modelexchanges4; selExc; modelexchanges3;BM]);% modelexchanges4; selExc
    modelClosed.lb(find(ismember(modelClosed.rxns, modelClosed.rxns(modelexchanges))))=0;
    %set upper bround to 1000
    modelClosed.ub(selExc) =100;
    %modelClosed.ub(find(ismember(modelClosed.rxns, modelClosed.rxns(modelexchanges))))=100;
    model=modelClosed;
end