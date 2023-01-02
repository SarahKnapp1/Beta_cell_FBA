%%Sarah Knapp

% input: model model,color,flag
% output: h - hist, objective Flux

% The function helps the hyperglycemia_ATP function.
% input fluxs.

function [h,objFlux]=help_plot_hyperglycemia(model,color,flag)
    x=121;
    objFlux = zeros(x,1);
    flux_ATP=zeros(x,1);
    N=0;
    mone=0
    for k=0:0.01:4
        mone=mone+1
        N=-k;
        %essential amino acids 
        model = changeRxnBounds_essential_AA(model,N);

        FBAsolution = optimizeCbModel(model,'max',0);%,0);%optimize model
        if ~(FBAsolution.stat == 1)
            disp('Relaxed model is not feasible');
        end
        objFlux(mone) = FBAsolution.x(find(ismember(model.rxns,'INSsyn')));%FBAsolution.f;
        
        flux_ATP(mone)=FBAsolution.x(find(ismember(model.rxns,'ATPS4m')));
    end
    objFlux = objFlux'
    % change color, markerstyle, x-position, etc...


    h=plot(objFlux,flux_ATP,'-o','MarkerSize',2,'color',color,'LineWidth',4)
    hold on 
    if flag==1
        plot(objFlux(19),flux_ATP(19),'o-','MarkerSize',6,'MarkerFaceColor',[0 0.5 0]','MarkerEdgeColor',[0 0.5 0]) 
    elseif flag==0
        plot(objFlux(95),flux_ATP(95),'o-','MarkerSize',6,'MarkerFaceColor',[0 0.5 0]','MarkerEdgeColor',[0 0.5 0]) 
    end
    xlabel(['Basal insulin synthesis rate [mmol/(gDW*h)]']);
    axis tight;
    ylabel('ATP production rate [mmol/(gDW*h)]');
    
end