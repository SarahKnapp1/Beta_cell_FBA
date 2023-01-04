%% Sarah Knapp

function matric_AA_reactions=amino_acid_from_TCA(model)

%   input: model -GEM model.
%   output: matric_AA_reactions- rows - EAA, column- reaction, each cell each 
%   cell contains the flux of the reaction that produces the amino acid.
%   
%
%   The function checks from where the nonessential amino acids are produced
%   in the mitochondria or the cytoplasm, what reaction produce them.
%   Then the function calculate how much flux there is through the reactions
%   that produce the non-essential amino acids. The function also generates
%   a heatmap of the fluxes of the reactions the produce the nonessential
%   amino acids from the TCA cycle.


    rec_of_interest={'AGTim','ALATA_L','ASNS1','ASNNm','ASPTA','ASPTAm','GHMT2r',...
        'GHMT2rm','GLNS','GLUDxm','GLUDym','GLUNm','PSERT','PSP_L','TYRTA','r0081','r0127','r0156',...
        'r0399','CYSO','P5CRxm','P5CRm','P5CR','ARGSL','ARGDCm','CYSTA'};
    %nonessential amino acids from the mitochondria or the cytoplasm
    nonessential_AA_met= {'ala_L','arg_L','asp_L','asn_L','cys_L','glu_L','gln_L','gly','pro_L','ser_L','tyr_L'};

    %change objective function
    model = changeObjective(model,{'INSsyn','ATPS4m'},[0.999,0.001]);

    %optimize model
    FBAsolution=optimizeCbModel(model,'max',0);
    len_rec_of_interest=length(rec_of_interest);
    len_nonessential_AA_met=length(nonessential_AA_met);
    matric_AA_reactions= zeros(len_nonessential_AA_met,len_rec_of_interest);
    flux_rec_of_interest = zeros(len_rec_of_interest,1);
    
    
    for i=1:len_rec_of_interest

        %find the indexes of the stoichiometric coefficient values of each reaction
        vec_met_stoichiometric_coefficient=find(model.S(:,find(ismember(model.rxns,rec_of_interest(i)))));
        for j=1:length(vec_met_stoichiometric_coefficient)
            for w=1:len_nonessential_AA_met
                if contains(char(model.mets(vec_met_stoichiometric_coefficient(j))),char(nonessential_AA_met(w)))
                    %flux of reaction of interest for the metabolite in the
                    %reaction
                    %if flux is - this means the precursor is a react
                    con1=model.S(vec_met_stoichiometric_coefficient(j),find(ismember(model.rxns,rec_of_interest(i))));
                    if (con1<0) &&((FBAsolution.x(find(ismember(model.rxns,rec_of_interest(i)))))<0)
                        matric_AA_reactions(w,i)=FBAsolution.x(find(ismember(model.rxns,rec_of_interest(i))));
                    end
                    %if it is postive it means product
                    if (con1>0) && (FBAsolution.x(find(ismember(model.rxns,rec_of_interest(i)))))>0
                        matric_AA_reactions(w,i)=FBAsolution.x(find(ismember(model.rxns,rec_of_interest(i))));
                    end
                end
            end
        end
    end
    
   
   figure
   mat__AA_reactions=log10(abs(matric_AA_reactions));
   
   mat__AA_reactions(mat__AA_reactions == -inf) = NaN;
%    mat__AA_reactions=abs(matric_AA_reactions)
   %heatmap of the nonessential amino acids from the TCA cycle
    xvalues = rec_of_interest;
    yvalues = nonessential_AA_met;
    h = heatmap(xvalues,yvalues,mat__AA_reactions);

    h.MissingDataLabel = 'No Flux';
    h.MissingDataColor='white'
    
    mycolormap = customcolormap([0 0.8 1],[0 1 0.8;0 0.6 0.8;0 0 0.8]);
    colorbar;
    colormap(mycolormap);
    %colormap(summer);
    h.Title = 'Nonessential AA from the TCA cycle';
    h.XLabel = 'Reactions';
    h.YLabel = 'nonessential AA';
    
end