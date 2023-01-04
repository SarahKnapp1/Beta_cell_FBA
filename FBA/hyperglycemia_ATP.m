%% Sarah Knapp


function objFlux=hyperglycemia_ATP(model)

%   input: model
%   output: figure of:
%   ATP production as a function of Insulin synthesis under different glucose
%   input flux. 

%   The function checks when the demand for insulin synthesis is high as in hyperinsulinemia,
%   a decrease in mitochondrial ATP production occurs in differner glucose
%   input fluxs.

    model = changeObjective(model,{'INSsyn','ATPS4m'},[0.999,0.001]);%'ATPS4mi''biomass_reaction'
    FBAsolution = optimizeCbModel(model,'max',0);%,0);
    
    if FBAsolution.stat == 1
    disp('Relaxed model is feasible');
    end
    figure
    model = changeRxnBounds_essential_AA(model,-1);
    hyperglycemia=[];
    
    % defining the allowed carbon sources Glucose.
    model= changeRxnBounds(model,'EX_glc_e',-10,'l');
    
    for i=0:1


        model = changeObjective(model,{'INSsyn','ATPS4m'},[0.999,0.001]);
        FBAsolution = optimizeCbModel(model,'max',0);
        ATP_flux=FBAsolution.x(find(ismember(model.rxns,'ATPS4m')));

        if i==0
            [h1,obj]=help_plot_hyperglycemia(model,[98/255,197/255,218/255],1)

           
        else
            [h2,objFlux]=help_plot_hyperglycemia(model,[254/255,125/255,104/255],0)
        end
        model= changeRxnBounds(model,'EX_glc_e',-20,'l');

            
    end
    
     legend([h1(1) h2(1)],'Glucose input flux = 5','Glucose input flux = 10');
      
     set(gcf,'color','w');
     set(gca,'FontSize',11)
      
     hold off

    
end
