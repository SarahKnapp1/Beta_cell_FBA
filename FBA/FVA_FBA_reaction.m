%% Sarah Knapp 

function FVA_FBA_reaction(model,rxns)
%   input: model- GEM model, rxns - reactions of intreat.
%   output: A plot of the variations in fluxes in the intreat reaction,
%   including FVA and FBA
%
%   The function perform FVA and FBA analysis on your intrest reaction
%   and plots a figure of fluxes of the intreast reaction as a function of insulin synthesis.
%   Flux variablity analysis optimized for the CPLEX solver.
%   Solves LPs probleams.
    N=0;
    mone=0;
    x=101;
    objFlux_ins = zeros(x,1);
    flux_reaction= zeros(x,1);
    minFlux_FVA=zeros(x,1);
    maxFlux_FVA=zeros(x,1);
    
    for i=0:100
        mone=mone+1;
        N=-i*0.05;
        %essential amino acid change bounds
        model = changeRxnBounds_essential_AA(model,N)
        
        
        model =changeObjective(model,{'INSsyn','ATPS4m'},[0.999,0.001]);
        FBAsolution=optimizeCbModel(model,'max',0);%optimize model
        
        objFlux_ins(i+1) = FBAsolution.x(find(ismember(model.rxns,'INSsyn')));
        flux_reaction(i+1)=FBAsolution.x(find(ismember(model.rxns,rxnsList)))
        
        INSsyn=FBAsolution.x(find(ismember(model.rxns,'INSsyn')))
        ATPS4m=FBAsolution.x(find(ismember(model.rxns,'ATPS4m')))
        
        %Constraints in the model
        model1=changeRxnBounds(model,'INSsyn',INSsyn,'u');
        model1=changeRxnBounds(model,'INSsyn',INSsyn,'l');
        
        model1=changeRxnBounds(model,'ATPS4m',ATPS4m,'u');
        model1=changeRxnBounds(model,'ATPS4m',ATPS4m,'l');
        
        % Run FVA analysis for the model with the constraints that simulates beta cell conditions
        [minFlux, maxFlux, Vmin, Vmax] = fluxVariability(model1,100,'max',rxnsList);
        %[minFlux, maxFlux, optsol, ret, fbasol, fvamin, fvamax, statussolmin, statussolmax] = fastFVA(model1,100,'max','ibm_cplex', rxnsList);
     
         minFlux_FVA(i+1)=minFlux;
         maxFlux_FVA(i+1)=maxFlux;
    end

    figure
    controlFlux=[0:0.05:5];
    length(controlFlux)
    objFlux_ins = objFlux_ins;
    plot(objFlux_ins,flux_reaction,'-o','MarkerSize',1,'color',[0.4 0.6 0],'LineWidth',2)
        %fill
    %patch([objFlux_ins flipud(objFlux_ins)],[minFlux_FVA maxFlux_FVA], [0.5 0.5 0], 'edgecolor', 'none', 'facealpha', 0.1);
    hold all
    plot(objFlux_ins,maxFlux_FVA,'-o','MarkerSize',2,'color',[0.5 0.5 0],'LineWidth',3)
    
    
    
    plot(objFlux_ins,minFlux_FVA,'-o','MarkerSize',2,'color',[0.5 0.5 0],'LineWidth',3)

    
    
    axis tight;
    xlabel(['Insulin synthesis [mmol/(gDW*h)]']);
    ylabel('Fluxes');
    title(['Variations in fluxes in the ',num2str(cell2mat(rxnsList)),' reaction under beta cell conditions']);
    legend({'FBA', 'FVA'}, 'Location', 'southwest')
    set(gcf,'color','w');
    
    hold off

% Flux variablity analysis optimized for the CPLEX solver.
% Solves LPs of the form:
% [minFlux, maxFlux, optsol, ret, fbasol, fvamin, fvamax, statussolmin, statussolmax] = fastFVA(model,100, 'max','ibm_cplex');

end