%% Sarah Knapp

function FVA_TCA(model)

%   input: model - GEM model
%   output: Barplot x-axis: TCA reactions, y-axis: Flux range in
%   norminsulinemia and hyperinsulinemia
%
%   The function performes FVA and FBA on the reactions within the TCA
%   cycle in norminsulinemia and hyperinsulinemia

    TCA_reactions ={'PCm','PDHm','LDH_L','CSm','ACONTm','ICDHxm','ICDHyrm','AKGDm','SUCOAS1m','SUCOASm','SUCD1m','FUMm','MDHm'};
   

    %essential amino acid change bounds
    model = changeRxnBounds_essential_AA(model,-0.2)

    model =changeObjective(model,{'INSsyn','ATPS4m'},[0.999,0.001]);
    FBAsolution=optimizeCbModel(model,'max',0);%optimize model

    INSsyn=FBAsolution.x(find(ismember(model.rxns,'INSsyn')))
    ATPS4m=FBAsolution.x(find(ismember(model.rxns,'ATPS4m')))

    
    TCA_reactions_flux=[];
    
    for i=1:length(TCA_reactions)
        flux_rec=FBAsolution.x(find(ismember(model.rxns,TCA_reactions(i))));
        TCA_reactions_flux=[TCA_reactions_flux;flux_rec];
        
    end
        
    
    model1=model;
    
    %Constraints in the model
    model1=changeRxnBounds(model1,'INSsyn',INSsyn,'u');
    model1=changeRxnBounds(model1,'INSsyn',INSsyn,'l');
        
    model1=changeRxnBounds(model1,'ATPS4m',ATPS4m,'u');
    model1=changeRxnBounds(model1,'ATPS4m',ATPS4m,'l');
        
    % Run FVA analysis for the model with the constraints that simulates beta cell conditions
    [minFlux, maxFlux, Vmin, Vmax] = fluxVariability(model1,100,'max', TCA_reactions);
    %[minFlux, maxFlux, optsol, ret, fbasol, fvamin, fvamax, statussolmin, statussolmax] = fastFVA(model1,100,'max','ibm_cplex', rxnsList);
    
    
    
    days = 1:2:26;    
    length(days)
    length(minFlux)
    % Create Data
%     figure
%     plot([days; days], [minFlux'; maxFlux'],'color', 'g', 'LineWidth',5)
%     hold all
%     
%     plot(days,TCA_reactions_flux ,'.')
    %essential amino acid change bounds
    model = changeRxnBounds_essential_AA(model,-3)

    model =changeObjective(model,{'INSsyn','ATPS4m'},[0.999,0.001]);
    FBAsolution=optimizeCbModel(model,'max',0);%optimize model

    INSsyn=FBAsolution.x(find(ismember(model.rxns,'INSsyn')))
    ATPS4m=FBAsolution.x(find(ismember(model.rxns,'ATPS4m')))
    
    TCA_reactions_flux=[];
    
    for w=1:length(TCA_reactions)
        flux_rec=FBAsolution.x(find(ismember(model.rxns,TCA_reactions(w))));
        TCA_reactions_flux=[TCA_reactions_flux;flux_rec];
        
    end

   
    
    %Constraints in the model
    model=changeRxnBounds(model,'INSsyn',INSsyn,'u');
    model=changeRxnBounds(model,'INSsyn',INSsyn,'l');
        
    model=changeRxnBounds(model,'ATPS4m',ATPS4m,'u');
    model=changeRxnBounds(model,'ATPS4m',ATPS4m,'l');
        
    % Run FVA analysis for the model with the constraints that simulates beta cell conditions
    [minFlux, maxFlux, Vmin, Vmax] = fluxVariability(model,100,'max', TCA_reactions);
    
    
     days = 2:2:26;    
    
    plot([days; days], [minFlux'; maxFlux'], 'color',[0.4 0.6 0],'LineWidth',5)
    
    
    plot(days,TCA_reactions_flux ,'.')
     

 
    set(gca, 'XTick', 1:2:26, 'XTickLabel', TCA_reactions);
    grid
    xlabel('TCA reactions')
    ylabel('Flux range [mmol/(gDW*h)]')
    legend('Norminsulinemia' ,'Hyperinsulinemia')
%     legend('Norminsulinemia' ,'Hyperinsulinemia')
%     xlabel('Reaction')
%     ylabel('Flux range [mmol/(gDW*h)]')


end