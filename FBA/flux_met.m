%%Sarah Knapp

%input: model
%output: matrix with the fluxs of interested reactions

%The function finds the flux for intersted reactions when optimizing the objective
%function in two states (norminsulinmia and hyperinsulinemia)

function matrix=flux_met(model)


    %change objective function
    model=changeObjective(model,{'INSsyn','ATPS4m'},[0.999,0.001]);

    %optimize model
    FBAsolution=optimizeCbModel(model,'max',0);
    
    % DM_atp_c_, Demand for ATP, Cytosolic, atp[c] + h2o[c] -> adp[c] + h[c] + pi[c]
    % ATPS4m, ATP synthase (four protons for one ATP), adp[m] + 4.0 h[i] + pi[m] -> atp[m] + 3.0 h[m] + h2o[m]
    % ATPtm ADP/ATP Transporter, Mitochondrial adp[c] + atp[m] -> adp[m] + atp[c]
    % ASPTA,Aspartate Transaminase, akg[c] + asp_L[c] <=> glu_L[c] + oaa[c]
    % PCm,Pyruvate Carboxylase,atp[m] + hco3[m] + pyr[m] -> adp[m] + h[m] + oaa[m] + pi[m]
    % PDHm Pyruvate Dehydrogenase, coa[m] + nad[m] + pyr[m] -> accoa[m] + co2[m] + nadh[m]
    % GLUDxm, Glutamate Dehydrogenase (NAD), Mitochondrial, glu_L[m] + h2o[m] + nad[m] <=> akg[m] + h[m] + nadh[m] + nh4[m]
    % OAADC, Oxaloacetate Decarboxylase, h[c] + oaa[c] -> co2[c] + pyr[c]
    % MDH, Malate Dehydrogenase, mal_L[c] + nad[c] <=> h[c] + nadh[c] + oaa[c]
    % LDH_L, L-Lactate Dehydrogenase, lac_L[c] + nad[c] <=> h[c] + nadh[c] + pyr[c]
    % GLNS, Glutamine Synthetase, atp[c] + glu_L[c] + nh4[c] -> adp[c] + gln_L[c] + h[c] + pi[c]
    % ASNS1, Asparagine Synthase (Glutamine-Hydrolysing), asp_L[c] + atp[c] + gln_L[c] + h2o[c] -> amp[c] + asn_L[c] + glu_L[c] + h[c] + ppi[c]
    %r0165, UTP:Pyruvate O2-Phosphotransferase, h[c] + pep[c] + udp[c] -> pyr[c] + utp[c]
    
    %Each time change to the reaction of your intrest.
    rec_of_interest={'LEUTA'};%'PCm','PDHm','LDH_L','CSm','ACONTm','ICDHxm','ICDHyrm','AKGDm','SUCOAS1m','SUCOASm','SUCD1m',...
       %'FUMm','MDHm','ASPTA','ASNS1','ALATA_L','GLNS','GLUDxm','GLUDym','r0081','r0399','VALTA','VALTAm','r0173','MDH','LEUTAm','LEUTA'};%,...
       % 'CITL','GLUt2m','MDH','FUMSO3tm','FUMSO4tm','FUMtm','r0822','MALSO3tm','MALSO4tm','MALtm'};
 %'EX_his_L_e','EX_ile_L_e','EX_leu_L_e','EX_lys_L_e','EX_met_L_e','EX_phe_L_e','EX_thr_L_e','EX_trp_L_e','EX_val_L_e'};
        %'r0355','HEX1','r0354','ASPTA','INSsyn','PYK','r0165','PCm','PDHm','CSm','GLUDxm', 'GLUDym','OAADC','LDH_L','L_LACDcm','GLNS','ASNS1','EX_glc_e','GLCt1r','GLCt2_2'};%,...
        %'GLCt2_2','GLCt1r','HEX1','PGI','PFK','FBA','GAPD','PGK','PGM','ENO','PYK','PYRt2m','CSm','ACONTm','ICDHxm','AKGDm','SUCOASm','SUCD1m','FUMm','MDHm'};
   
   
   %Reverse flux to reactions that register that enter the TCA circle but we want to see their exit.
   list_rec={'SUCOAS1m','SUCOASm','GLUDxm','GLUDym','GLUt2m','ALATA_L','ASPTAm','ASPTA','FUMSO3tm','FUMSO4tm','FUMtm'...
        'r0822','MALSO3tm','MALSO4tm','MALtm'};
    
   
    len_rec_of_interest=length(rec_of_interest);
    flux_rec_of_interest = zeros(len_rec_of_interest,1);
    
    for i=1:len_rec_of_interest
        if sum(strcmp(rec_of_interest(i),list_rec))
            flux_rec_of_interest(i)=-1*FBAsolution.x(find(ismember(model.rxns,rec_of_interest(i))));
        else
            flux_rec_of_interest(i)=FBAsolution.x(find(ismember(model.rxns,rec_of_interest(i))));
        end
    end
    matrix=[rec_of_interest;num2cell(((flux_rec_of_interest')))];
end