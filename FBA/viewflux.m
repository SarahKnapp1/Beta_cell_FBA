%%Sarah Knapp

function viewflux(model)

%   input: model- GEM model
%   output: a map with the flux in the TCA cycle when performing FBA.
%
%   The function visualizes the flux in the TCA cycle by the CellDesigner


    %change objective function
    model = changeObjective(model,{'INSsyn','ATPS4m'},[0.999,0.0001]);%'INSsyn' )%

    %model=changeObjective(model,'INSsyn',1);

    %optimize model
    FBAsolution=optimizeCbModel(model,'max',0);
    
    %make the reactions of intrest go out from the TCA cycle postive (clockwise).
    
   rec_of_interest={'SUCOAS1m','SUCOASm','GLUDxm','GLUDym','GLUt2m','ALATA_L','ASPTAm','ASPTA','FUMSO3tm','FUMSO4tm','FUMtm'...
        'r0822','MALSO3tm','MALSO4tm','MALtm'};
    
    FBAsolution.x(find(ismember(model.rxns,rec_of_interest)))=(-1*FBAsolution.x(find(ismember(model.rxns,rec_of_interest))));
    
    FBAsolution.v(find(ismember(model.rxns,rec_of_interest)))=(-1*FBAsolution.v(find(ismember(model.rxns,rec_of_interest))));

    

    
    
    %read xml file to map
    [xmlGly, mapGly] = transformXML2Map('TCA.xml');

    % [diffReactions, diffMetabolites, diffReversibility, diffFormula] = ...
    %     checkCDerrors(mapGly, model);
    % 
    % correctMets = diffMetabolites.extraMetsModel;
    % mapGly = correctMetNameCD(mapGly, ...
    %     diffMetabolites, correctMets);

    mapSpecificATP = addFluxFBAdirectionAndColor(mapGly, ...
        model, FBAsolution);
    transformMap2XML(xmlGly, ...
        mapSpecificATP, 'flux_in_TCA.xml');
end