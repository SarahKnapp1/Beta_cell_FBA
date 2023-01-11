%% Sarah Knapp 

function matrix_INSsyn=Glucose_ins_syn_change(model)
%   input: model - GEM model 
%   output: matrix_INSsyn - the values of flux of insulin synthesis when
%   changing glucose and EAA
%
%   The function checked where glucose causes insulin synthesis when there is a change in EAA
%   plot an heatmaap of where the y axis-glucose, x axis - EAA, cell- the flux of
%   insulin synthesis.
%


    matrix_INSsyn=[];
    for i=0:30
        %defining the allowed carbon sources Glucose.
        model= changeRxnBounds(model,'EX_glc_e',-i,'l')
        ins_syn_glc_row=[];
        for j=0:50
            N=-j*0.1;
            %essential amino acid change bounds
            model = changeRxnBounds_essential_AA(model,N)


            %change objective function
            model=changeObjective(model,{'INSsyn','ATPS4m'},[0.999,0.001]);

            %optimize model
            FBAsolution=optimizeCbModel(model,'max',0);
            
            INSsyn=FBAsolution.x(find(ismember(model.rxns,'INSsyn')))
            
            ins_syn_glc_row=[ins_syn_glc_row,INSsyn]
           
        end
        matrix_INSsyn=[matrix_INSsyn;ins_syn_glc_row];
        
    end


    xvalues = [0:0.1:5];
    yvalues = [0:30];
    h = heatmap(xvalues,yvalues,matrix_INSsyn);

    h.MissingDataLabel = 'No Flux';
    h.MissingDataColor='white'
%     
    mycolormap = customcolormap([0 0.8 1],[0 1 0.8;0 0.6 0.8;0 0 1;0 0 0.8]);
    colorbar;
    %colormap(mycolormap);
    colormap(summer);
   
    h.XLabel = 'Essential AA uptake[mmol/(gDW*h)]';
    h.YLabel = 'Glucose uptake';

end