%% Sarah Knapp



function hyperglycemia=ins_as_function_of_glucose_fasting(model)

%   input: model
%   output: A figure of fasting extracellular glucose as a function of insulin synthesis. 
%
%   The function shoes that higher levels of extracellular glucose are required to reach
%   sufficient ATP for insulin secretion leading to hyperglycemia.

    model = changeObjective(model,{'INSsyn','ATPS4m'},[0.999,0.001]);%'ATPS4mi''biomass_reaction'
    FBAsolution = optimizeCbModel(model,'max',0);%,0);
    
    %flux_ATP=zeros(30,1);
    if FBAsolution.stat == 1
    disp('Relaxed model is feasible');
    end
    
    figure
    x=61;
    N=0;
    hyperglycemia=[];
    ins_syn=[];
    for k=0:0.05:4
        N=-k;
        %essential amino acids 
        model = changeRxnBounds_essential_AA(model,N)
        model = changeObjective(model,{'INSsyn','ATPS4m'},[0.999,0.001]);
        FBAsolution = optimizeCbModel(model,'max',0);%,0);%optimize model
        objFlux = FBAsolution.x(find(ismember(model.rxns,'INSsyn')));%FBAsolution.f;

        for i=0:0.01:30
            lb=-i;
            % defining the allowed carbon sources Glucose.
            model= changeRxnBounds(model,'EX_glc_e',lb,'l');
            model = changeObjective(model,{'INSsyn','ATPS4m'},[0.999,0.001]);
            FBAsolution = optimizeCbModel(model,'max',0);
            ATP_flux=FBAsolution.x(find(ismember(model.rxns,'ATPS4m')));
            %flux_ATP(i+1)=ATP_flux;
            ATP_flux=round(ATP_flux,1);
            if ATP_flux==10.2
                hyperglycemia=[hyperglycemia;(i/2)];
                ins_syn=[ins_syn;objFlux];
            end
        end
    end
    
    h=plot(ins_syn,hyperglycemia,'-o','MarkerSize',2,'color',[0 0.5 0],'LineWidth',4);
    xlabel(['Basal insulin synthesis flux [mmol/(gDW*h)]']);
    axis tight;
    ylabel('Glucose levels required for secretion [mmol/(gDW*h)]');
    hold on
    plot(ins_syn(19),hyperglycemia(19),'o-','MarkerSize',6,'MarkerFaceColor',[98/255,197/255,218/255],'MarkerEdgeColor',[98/255,197/255,218/255]) 
    plot(ins_syn(88),hyperglycemia(88),'o-','MarkerSize',6,'MarkerFaceColor',[254/255,125/255,104/255],'MarkerEdgeColor',[254/255,125/255,104/255]) 
    set(gcf,'color','w');
    legend(h(1),'ATP production rate = 10');
    set(gca,'FontSize',11)

              
end
