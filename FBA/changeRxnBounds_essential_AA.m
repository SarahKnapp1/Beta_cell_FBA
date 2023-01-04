%% Sarah Knapp

function model = changeRxnBounds_essential_AA(model,N)

%   input: model - GEM model, N = the amount of uptake flux 
%   output: model- updated model
%   The function changes the bounds of the uptake reactions of the essential amino acids 

%   essential amino acids    
         
% model= changeRxnBounds(model,'EX_his_L_e',-1.2+N,'l');
% model= changeRxnBounds(model,'EX_lys_L_e',-2.9+N,'l');
% model= changeRxnBounds(model,'EX_met_L_e',-0.30+N,'l');
% model= changeRxnBounds(model,'EX_phe_L_e',-0.86+N,'l');
% model= changeRxnBounds(model,'EX_thr_L_e',-1.5+N,'l');
% model= changeRxnBounds(model,'EX_trp_L_e',-0.5+N,'l');
% 
% model= changeRxnBounds(model,'EX_val_L_e',-3.5+N,'l');
% model= changeRxnBounds(model,'EX_ile_L_e',-1.4+N,'l');
% model= changeRxnBounds(model,'EX_leu_L_e',-1.7+N,'l');


        model= changeRxnBounds(model,'EX_his_L_e',-0.7+N,'l');
        model= changeRxnBounds(model,'EX_lys_L_e',-0.7+N,'l');
        model= changeRxnBounds(model,'EX_met_L_e',-0.7+N,'l');
        
        model= changeRxnBounds(model,'EX_phe_L_e',-0.7+N,'l');
        model= changeRxnBounds(model,'EX_thr_L_e',-0.7+N,'l');
        model= changeRxnBounds(model,'EX_trp_L_e',-0.7+N,'l');
        
        model= changeRxnBounds(model,'EX_val_L_e',-0.7+N,'l');
        model= changeRxnBounds(model,'EX_ile_L_e',-0.7+N,'l');
        model= changeRxnBounds(model,'EX_leu_L_e',-0.7+N,'l');



end