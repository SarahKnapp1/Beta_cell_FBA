%% Sarah Knapp

function Max_Fraction_EAA_ins=fraction_of_carbon_EAA_to_ins_syn(model) 
%   input: model - GEM model
%   output: Max_Fraction_EAA_ins - Number of carbon flux from (Essential AA
%   / Insulin synthesis)
%
%   The function checks if carbon come from glucose or EAA in order to
%   synthezie insulin 
%   Gives a plot of Number of carbon flux from (Essential AA / Insulin
%   synthesis) as a function of Insulin synthesis 
%
%   Note: need to take into acount that not all what is uptaken is relly been
%   carried by the reaction.

    EAA ={'EX_his_L_e','EX_ile_L_e','EX_leu_L_e','EX_lys_L_e','EX_met_L_e','EX_phe_L_e',...
        'EX_thr_L_e','EX_trp_L_e','EX_val_L_e'};
    carbon_EAA = [6,6,6,6,5,9,4,11,5];
    %nitrogen_EAA = [3,1,1,2,1,1,1,2,1];
    Glucose_carbon=6;
    proInsulin_carbon=410;
    %proInsulin_nitrogen=114;
    
    Max_Fraction_EAA_ins=[];
    
    INSsyn= zeros(81,1);
    
    for i=0:80
        N=-i*0.05;;
        %essential amino acid change bounds
        model = changeRxnBounds_essential_AA(model,N)


        %change objective function
        model=changeObjective(model,{'INSsyn','ATPS4m'},[0.999,0.001]);

        %optimize model
        FBAsolution=optimizeCbModel(model,'max',0);


        EAA_uptake__flux=[];

        for k=1:length(EAA)
            flux_rec=FBAsolution.x(find(ismember(model.rxns,EAA(k))));
            EAA_uptake__flux=[EAA_uptake__flux;flux_rec];        
        end

        %change from coulnm to rows
        EAA_uptake__flux=EAA_uptake__flux';
        %multiply the number of influx of each EAA in there number of carbons then sum it up
        carbon_EAA_flux=sum(abs(EAA_uptake__flux.*carbon_EAA));%nitrogen_EAA

        INSsyn(i+1)=FBAsolution.x(find(ismember(model.rxns,'INSsyn')));

        carbon_ins_syn=(INSsyn(i+1))*proInsulin_carbon; %proInsulin_nitrogen

        Max_Fraction_EAA_ins= [Max_Fraction_EAA_ins; carbon_EAA_flux/carbon_ins_syn];

    end
    
    
    figure
    %controlFlux=[0:0.01:2];

    plot(INSsyn,Max_Fraction_EAA_ins,'-o','MarkerSize',2,'color',[0 0.5 0],'LineWidth',3)
    
    axis tight;
    xlabel(['Insulin synthesis [mmol/(gDW*h)]']);
    ylabel(['Number of carbon flux from (Essential AA / Insulin synthesis)']);
    
     
end
