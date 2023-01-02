%% Sarah Knapp 

%input: model, objective_function and the reactionof intrest(glucose
%exchange reaction).
%output: 3 figures of: 1. Glucose uptake rate as a function of glucose input flux.
% 2. Hexokinase as a function of glucose uptake.
% 3. Basal insulin synthesis rate as a function of glucose uptake rate. 

%The function checks how much flux does the exchange glucose reaction relly
%carries. Also, shows if the glucose which enters the cell continues into the 
%glycolysis pathway through the hexokinase reaction (HEX1). Addtionly,it determines 
%the effect of varying glucose uptake on insulin synthesis


function plot_objective_as_function_glu(model,objective_function,reaction)

    model = changeObjective(model,{'INSsyn','ATPS4m'},[0.999,0.001]);
    FBAsolution = optimizeCbModel(model,'max',0);
    if FBAsolution.stat == 1
    disp('Relaxed model is feasible');
    end
    % set both the upper and lower bounds of the glucose exchange reaction to values between 0 and -40 mmol gDW-1hr-1
    % use optimizeCbModel to perform FBA with each uptake rate.
    objFlux = zeros(401,1);
    flux_HEX=zeros(401,1);
    flux_EX_glc_e=zeros(401,1);
    for i=0:400
        N=i*0.05
        model = changeRxnBounds(model,reaction,-N,'l');% 
        FBAsolution = optimizeCbModel(model,'max',0);%optimize model
        
        flux_EX_glc_e(i+1)=-FBAsolution.x(find(ismember(model.rxns,'EX_glc_e')));
        
        %if i<51
            objFlux(i+1) = FBAsolution.x(find(ismember(model.rxns,'INSsyn')));%FBAsolution.f;
           
            flux_r0355=FBAsolution.x(find(ismember(model.rxns,'r0355')));
            flux_HEX1=FBAsolution.x(find(ismember(model.rxns,'HEX1')));
            flux_r0354=FBAsolution.x(find(ismember(model.rxns,'r0354')));
            flux_HEX(i+1)=flux_r0355+flux_HEX1+flux_r0354;
       % end
    end
   
    figure
    controlFlux=[0:0.1:40];
    objFlux = objFlux';
    length(objFlux);
    plot(flux_EX_glc_e,objFlux,'-o','MarkerSize',2,'color',[128 128 128]/255,'LineWidth',4)
    xlabel('Glucose Uptake rate [mmol/(gDW*h)]');
    axis tight;
    if strcmp('INSsyn',objective_function)
        ylabel(['Basal insulin synthesis rate [mmol/(gDW*h)]']);
    
    else 
        ylabel(['ATP production rate [mmol/(gDW*h)]']);
    end
    
    set(gcf,'color','w');
    set(gca,'FontSize',11)
    
    figure
    plot(flux_EX_glc_e,flux_HEX,'-o','MarkerSize',2,'color','c','LineWidth',4)
    xlabel('Glucose uptake rate [mmol/(gDW*h)]');
    axis tight;
    ylabel(['Hexokinase flux [mmol/(gDW*h)]']);
    set(gcf,'color','w');
    set(gca,'FontSize',11)
    
    figure
    plot(controlFlux,flux_EX_glc_e,'-o','MarkerSize',2,'color','c','LineWidth',4)
    xlabel('Glucose input flux [mmol/(gDW*h)]');
    axis tight;
    ylabel(['Glucose uptake rate [mmol/(gDW*h)]']);
    
    set(gcf,'color','w');
    set(gca,'FontSize',11)
           
     
       
end