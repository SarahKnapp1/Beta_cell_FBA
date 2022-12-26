%Sarah Knapp
%Yael's experiment with 13C-glucose.
%The deta is from islets
%Their values are from LC-MC analysis - Intensity. is). 
%Each metabolite has a value of intensity.

%Purpose: to show that labeled glucose has gone to non-essential amino acids.

close all

%Read the excel file of the islets 
[value_Medium_experiment,text_Medium_experiment,Medium_experiment_data]=xlsread('Medium_experiment.xlsx');

[value_cells_experiment,text_cells_experiment,cells_experiment_data]=xlsread('Experiment cells.xlsx');

X = categorical(Medium_experiment_data(1,2:end));
X = reordercats(X,Medium_experiment_data(1,2:end));
Y_Medium = cell2mat(Medium_experiment_data(2,2:end));
Y_cells = cell2mat(cells_experiment_data(2,2:end));

%bar plot, Y- axis:Peak area,X-axis:Essential amino acids
figure
bh=bar(X,[Y_Medium;Y_cells],'FaceColor',[0 .5 .5],'EdgeColor',[0 .5 .5])

bh(1).FaceColor=[ 0, 0.9,0.9]

legend('Beta cells' ,'Islets')
xlabel('Non-essential amino acids')
ylabel('Peak area')
set(gca,'FontSize',14)
set(gcf,'color','w');


 
