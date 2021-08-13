% Code for making Figure 4
%%
your_save_path = ['Z:\home\jen\Analysis\HVA Recording\Interneuron properties\'];
your_read_path = ['Z:\home\jen\Output\Manuscripts\2020 HVA Inhibition\Text\'];
addpath(your_read_path);

%% 
PV = load([your_save_path,'PV_grouped.mat']);
SOM = load([your_save_path,'SOM_grouped.mat']);
Pyr = load([your_save_path,'Pyr_grouped.mat']);

cell_names = {'PV','SOM','Pyr'};
area_names = {'LM','AL','PM','AM'};
int_color = [0.3 0.5 0.8;
    0.3 .8 0.4;
    0.7 0.7 0.7];
cmap = brewermap(4,'Set2');
cmap_cell = mat2cell(cmap,[1 1 1 1],3);
disp('analyzed and loaded');
%% Figure 4b-e
if isfield(PV,'your_read_path')
PV = rmfield(PV,'your_read_path');
end
plotStruct = [PV; SOM; Pyr];
font_size = 10;
figure;
posBins = [0,100,200,300,400,500];

for type_i = 1:3
    figure;
    subplot(1,4,1); hold on;
    plotSpread(plotStruct(type_i).AP_adapt,'distributionMarkers','o','binWidth',1,'spreadWidth',0.75,'distributionColors',int_color(type_i,:));
    fast_errbar(1:4,plotStruct(type_i).AP_adapt,2,'continuous',false);
    title(['AP adapt  p = ',num2str(plotStruct(type_i).p_APadapt)]);
    xlim([0 5]);xticks(1:4);xticklabels(area_names);ylim([0 1]);
    fix_axes(gcf,font_size,'Area','Adaptation');
    axis square

    subplot(1,4,2); hold on;
    plotSpread(plotStruct(type_i).Rin_grouped,'distributionMarkers','o','binWidth',1,'spreadWidth',0.75,'distributionColors',int_color(type_i,:));
    fast_errbar(1:4,plotStruct(type_i).Rin_grouped,2,'continuous',false);
    title(['Rin p = ',num2str(plotStruct(type_i).p_Rin)]);
    xlim([0 5]);xticks(1:4);xticklabels(area_names);ylim([0 200]);
    fix_axes(gcf,font_size,'Area','Rin (MOhms)');
    axis square


    subplot(1,4,3); hold on;
    plotY = plotStruct(type_i).Ih_slope;
    plotSpread(plotY,'distributionMarkers','o','binWidth',1,'spreadWidth',0.75,'distributionColors',int_color(type_i,:));
    fast_errbar(1:4,plotY,2,'continuous',false);
    title(['Ih  p = ',num2str(plotStruct(type_i).p_Ih)]);
    xlim([0 5]);xticks(1:4);xticklabels(area_names);
    ylim([0 20]);
    fix_axes(gcf,font_size,'Area','Ih (mV/pA)');
    axis square; 

    subplot(1,4,4); hold on;
    for area_i = 1:4
        fast_errbar(posBins,plotStruct(type_i).FI{area_i},1,'color',cmap(area_i,:),'shaded',true);
    end
%     legend(area_names,'Location','NorthEastOutside');
    fix_axes(gcf,font_size,'I inj','FR'); 
    if type_i == 3
        ylim([0 20]); 
        yticks(0:10:20);
    else
        ylim([0 200]);
    end
    title(['FI']);
    axis square;
    suptitle(cell_names{type_i});
end
