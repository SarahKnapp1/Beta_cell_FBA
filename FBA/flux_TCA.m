%% Sarah Knapp

function Flux_in_out_TCA=flux_TCA(model)

%   input: model- GEM model
%   output: Flux_in_out_TCA- an array with the amount of flux whitch enters and leaves the TCA
%   The function checks the flux which enters the TCA cycle and leaves 

    Met_TCA={'cit_m','icit_m','akg_m','succoa_m','succ_m','fum_m','mal_L_m','oaa_m'};
    Rec_TCA={'ACONTm','ICDHxm','ICDHyrm','AKGDm','SUCOAS1m','SUCOASm','SUCD1m','FUMm','MDHm'};
    
    %change objective function
    model = changeObjective(model,{'INSsyn','ATPS4m'},[0.999,0.001]);

    %optimize model
    FBAsolution=optimizeCbModel(model,'max',0);
     
    Flux_into_TCA=0;
    Flux_out_TCA=0;
    Flux_out_TCA_CO2=0;
    met_out_rec_neg=0;
    for i=1:length(FBAsolution.x)
        %exclude the rections in the TCA cycle itself
        %if ~sum(find(ismember(Rec_TCA,model.rxns(i))))
            if FBAsolution.x(i)>0 
                %find the indexes of the stoichiometric coefficient values of
                %each reaction, then find the metabolites which are produced as prodacts in the reaction 
                %by postive stoichiometric coefficient values
                %Enters the TCA
                met_in_rec=model.mets(find(model.S(:,i)>0));
                met_out_rec_neg=model.mets(find(model.S(:,i)<0));
               
                for j=1:length(met_in_rec)
                    %if the met' in the recation is in the TCA cycle add 
                    if sum(find(ismember(Met_TCA,met_in_rec(j))))>0 
                        met_in_rec
                         if sum(find(ismember(met_in_rec,'co2_m')))>0 
                             Flux_out_TCA_CO2=Flux_out_TCA_CO2+FBAsolution.x(i);
                         else
                            Flux_into_TCA=Flux_into_TCA+FBAsolution.x(i);
                         end
                    end
                end   
                for w=1:length(met_out_rec_neg)
                    %if the met' in the recation is in the TCA cycle add 
                    if sum(find(ismember(Met_TCA,met_out_rec_neg(w))))>0 
                        Flux_out_TCA=Flux_out_TCA+FBAsolution.x(i);
                        
                    end
                end 
            elseif FBAsolution.x(i)<0
                %find the indexes of the stoichiometric coefficient values of
                %each reaction, then find the metabolites which are precursor is a react in the reaction 
                %by negtive stoichiometric coefficient values
                %exitets the TCA 
                met_in_rec_neg=model.mets(find(model.S(:,i)<0));
                met_out_rec_pos=model.mets(find(model.S(:,i)>0));
                %CO2=model.mets(find(model.S(:,i)>0));
                for j=1:length(met_in_rec_neg)
                    if sum(find(ismember(Met_TCA,met_in_rec_neg(j))))>0 
                        
                        Flux_into_TCA=Flux_into_TCA+FBAsolution.x(i);
                        
                    
                    
                    end
                end 
                for k=1:length(met_out_rec_pos)
                    %if the met' in the recation is in the TCA cycle add 
                    if sum(find(ismember(Met_TCA,met_out_rec_pos(k))))>0  
                        
                        Flux_out_TCA=Flux_out_TCA+FBAsolution.x(i);
                    end
                end
                
            end
      %  end
    end
    Flux_in_out_TCA={Flux_out_TCA,Flux_into_TCA,Flux_out_TCA_CO2};
end
