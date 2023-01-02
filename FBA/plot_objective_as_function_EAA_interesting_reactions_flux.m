%%Sarah Knapp

% input: model
% output: figures of:  
% 1. Basal insulin synthesis rate a function of the uptake of EAA. 
% 2. a reaction as a function of insulin synthesis 
% including these reactions: pyruvate carboxylase (PC) and pyruvate dehydrogenase (PDH), leucine transaminase (LEUTA),
% mitochondrial CO2 production from the TCA, aspartate transaminase (ASPTA), citrate synthase (CSm),
% valine transaminase (VALTA), ATP production.

% The function analyzies reactions in the TCA from insulin synthesis.
% and modeles hyperinsulinemia by increasing the uptake of EAA 
%and a constant glucose uptake 

function plot_objective_as_function_EAA_interesting_reactions_flux(model)

    model = changeObjective(model,{'INSsyn','ATPS4m'},[0.999,0.001]);
    FBAsolution = optimizeCbModel(model,'max',0);%,0);
    if FBAsolution.stat == 1
    disp('Relaxed model is feasible');
    end
    % set the uptake of the essential amino acids to values between 0.7 to 4.7 mmol gDW-1hr-1
    % use optimizeCbModel to perform FBA with each uptake rate.
    x=81;
    objFlux = zeros(x,1);
    flux_PCm= zeros(x,1);
    flux_PDH= zeros(x,1);
    flux_CSm=zeros(x,1);
    flux_ATP_TCA=zeros(x,1);
    flux_Glc=zeros(x,1);
    flux_ATP=zeros(x,1);
    
    flux_ASPTA=zeros(x,1);
    flux_MDHm=zeros(x,1);
    flux_VALTA=zeros(x,1);
    flux_LEUTA =zeros(x,1);
    flux_co2=zeros(x,1);
        
    N=0;
    mone=0
    for i=0:80
        mone=mone+1;
        N=-i*0.05;
        %essential amino acid change bounds
        model = changeRxnBounds_essential_AA(model,N)
        
        
        FBAsolution = optimizeCbModel(model,'max',0);%,0);%optimize model
        objFlux(i+1) = FBAsolution.x(find(ismember(model.rxns,'INSsyn')));%FBAsolution.f;
        
%         PC=FBAsolution.x(find(ismember(model.rxns,'PCm')));
%         PDH=FBAsolution.x(find(ismember(model.rxns,'PDHm')));
%         flux_PCm(i+1)=PC+PDH;
        
        flux_PCm(i+1)=FBAsolution.x(find(ismember(model.rxns,'PCm')));
        flux_PDH(i+1)=FBAsolution.x(find(ismember(model.rxns,'PDHm')));
        flux_CSm(i+1)=FBAsolution.x(find(ismember(model.rxns,'CSm'))); 
        
        flux_ASPTA(i+1)=-(FBAsolution.x(find(ismember(model.rxns,'ASPTA')))); 
        flux_MDHm(i+1)=FBAsolution.x(find(ismember(model.rxns,'MDHm'))); 
        flux_VALTA(i+1)=FBAsolution.x(find(ismember(model.rxns,'VALTA'))); 
        flux_LEUTA(i+1)=FBAsolution.x(find(ismember(model.rxns,'LEUTA'))); 
        CO2_flux_production=FBAsolution.x(find(ismember(model.rxns,{'ICDHxm','ICDHyrm','AKGDm'}))); 
        flux_co2(i+1)=sum(CO2_flux_production)
        
        flux_Glc(i+1)=-FBAsolution.x(find(ismember(model.rxns,'EX_glc_e')));%L_LACDcm
        
        %Reaction in the TCA cycle that produce energy- ATP.
        %SUCOASm and SUCOAS1m when flux is negtive they produce ATP.
        Succinate_CoA=-FBAsolution.x(find(ismember(model.rxns,{'SUCOASm','SUCOAS1m'})));
        Dehydrogenase=FBAsolution.x(find(ismember(model.rxns,{'ICDHxm','ICDHyrm','AKGDm','SUCD1m','MDHm','MDH'})));
        flux_ATP_TCA(i+1)=(sum(Succinate_CoA)+sum(Dehydrogenase)) 
        
     
        flux_ATP(i+1)=FBAsolution.x(find(ismember(model.rxns,'ATPS4m')));
        

    end
    figure
    controlFlux=[0.7:0.05:4.7];
    length(controlFlux)
    objFlux = objFlux';
    length(objFlux)
    plot(controlFlux,objFlux,'-o','MarkerSize',2,'color',[200 100 16]/256,'LineWidth',4)
    
    axis tight;
    ylabel('Basal insulin synthesis rate [mmol/(gDW*h)]');
    xlabel('Uptake of essential amino acids [mmol/gDW*h]'); 
    set(gcf,'color','w');
    set(gca,'FontSize',11)
    
    
%     figure
%     plot(objFlux,flux_ATP,'color',[0.4940 0.1840 0.5560],'LineWidth',3)
%     xlabel(['Insulin synthesis [mmol/(gDW*h)]']);
%     axis tight;
%     ylabel('TCA reactions which produce ATP');
%     
    figure
    plot(objFlux,flux_ATP,'-o','MarkerSize',2,'color',[98/255,197/255,218/255],'LineWidth',4)
    xlabel('Basal insulin synthesis rate [mmol/(gDW*h)]');
    axis tight;
    ylabel('ATP production rate [mmol/(gDW*h)]');
    set(gcf,'color','w');
    set(gca,'FontSize',11)
   
    
    figure
    plot(objFlux,flux_PCm,'-o','MarkerSize',2,'Color','[0.6350, 0.0780, 0.1840]','LineWidth',4)
    %xlabel('Insulin synthesis[mmol/(gDW*h)]');
    %axis tight;
    %ylabel(['Pyruvate Carboxylase[mmol/(gDW*h)]']);
    
    hold on
    plot(objFlux,flux_PDH,'-o','MarkerSize',2,'Color',[98/255,197/255,218/255],'LineWidth',4)
    %xlabel('Insulin synthesis[mmol/(gDW*h)]');
    %axis tight;
    %ylabel(['Pyruvate Dehydrogenase[mmol/(gDW*h)]']);
    
    %plot(objFlux,flux_LDH_L,'-o','MarkerSize',2,'color','blue','LineWidth',3)
    xlabel('Basal insulin synthesis rate [mmol/(gDW*h)]');
    ylabel('Flux [mmol/(gDW*h)]');
    legend('PC','PDH');%,'LDH');
    set(gcf,'color','w');
    set(gca,'FontSize',11)
    hold off
    

    
    
    figure
    plot(objFlux,flux_ATP_TCA,'-o','MarkerSize',2,'Color',[98/255,197/255,218/255],'LineWidth',4)
    xlabel('Basal insulin synthesis rate [mmol/(gDW*h)]');
    axis tight;
    ylabel('TCA reactions which produce ATP [mmol/(gDW*h)]');
    set(gcf,'color','w');
    set(gca,'FontSize',11)
    
    figure
    plot(objFlux,flux_CSm,'-o','MarkerSize',2,'Color',[98/255,197/255,218/255],'LineWidth',4)
    xlabel('Basal insulin synthesis rate [mmol/(gDW*h)]');
    axis tight;
    ylabel('CSm flux [mmol/(gDW*h)]');
    set(gcf,'color','w');
    set(gca,'FontSize',11)
    
    figure
    plot(objFlux,flux_ASPTA,'-o','MarkerSize',2,'Color','[0.6350, 0.0780, 0.1840]','LineWidth',4)
    xlabel('Basal insulin synthesis rate [mmol/(gDW*h)]');
    axis tight;
    ylabel('ASPTA flux [mmol/(gDW*h)]');
    set(gcf,'color','w');
    set(gca,'FontSize',11)
    
    
    
    figure
    plot(objFlux, flux_MDHm,'-o','MarkerSize',2,'Color','[0.6350, 0.0780, 0.1840]','LineWidth',4)
    xlabel('Basal insulin synthesis rate [mmol/(gDW*h)]');
    axis tight;
    ylabel(['MDHm flux [mmol/(gDW*h)]']);
    set(gcf,'color','w');
    set(gca,'FontSize',11)
    
    
    figure
    plot(objFlux,flux_VALTA,'-o','MarkerSize',2,'Color','[0.6350, 0.0780, 0.1840]','LineWidth',4)
    xlabel('Basal insulin synthesis rate [mmol/(gDW*h)]');
    axis tight;
    ylabel(['VALTA flux [mmol/(gDW*h)]']);
    set(gcf,'color','w');
    set(gca,'FontSize',11)
    
    
    
    figure
    plot(objFlux,flux_LEUTA,'-o','MarkerSize',2,'Color','[0.6350, 0.0780, 0.1840]','LineWidth',4)
    xlabel('Basal insulin synthesis rate [mmol/(gDW*h)]');
    axis tight;
    ylabel(['LEUTA flux [mmol/(gDW*h)]']);
    set(gcf,'color','w');
    set(gca,'FontSize',11)
    
    figure
    plot(objFlux,flux_co2,'-o','MarkerSize',2,'Color',[98/255,197/255,218/255],'LineWidth',4)
    xlabel('Basal insulin synthesis rate [mmol/(gDW*h)]');
    axis tight;
    ylabel(['Mitochondrial CO2 production[mmol/(gDW*h)]']);
    set(gcf,'color','w');
    set(gca,'FontSize',11)
%     figure
%     plot(objFlux,flux_Glc,'c','LineWidth',3)
%     xlabel('Insulin synthesis[mmol/(gDW*h)]');
%     axis tight;
%     ylabel(['Glucose uptake flux[mmol/(gDW*h)]']);
%     
 
   
       
end