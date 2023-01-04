%%Sarah Knapp


function changed_reactions=flux_change_in_hyperinsulinmia(model)

%   input: model- GEM model. 
%   output: changed_reactions - The reaction with a change in there flun 
%   The function checks in which reaction the flux changes in hyperinsulinmia by 2 units of flux.
%   change objective function norminsulinemia

    model = changeRxnBounds_essential_AA(model,-0.3)
    model=changeObjective(model,{'INSsyn','ATPS4m'},[0.999,0.001]);

    %optimize model
    FBAsolution=optimizeCbModel(model,'max',0);

    %change objective function hyperinsulinemia
    model = changeRxnBounds_essential_AA(model,-2)
    model=changeObjective(model,{'INSsyn','ATPS4m'},[0.999,0.001]);

    %optimize model
    FBAsolution2=optimizeCbModel(model,'max',0);

    changed_reactions=model.rxns((FBAsolution2.x-FBAsolution.x)>=2)
    
end