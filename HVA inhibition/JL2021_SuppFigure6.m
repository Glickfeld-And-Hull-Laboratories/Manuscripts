% Supplementary figure 6

ccc;
load('DEC_analysis_grouped.mat')
%%

cell_names = {'PV';'SOM';'Pyr'};
for cell_type = 1:3
subplot(1,3,cell_type); hold on;
plot(timeVector_trace,Template{cell_type});
axis square; fix_axes(gcf,10,'time (s)','pA');
title(cell_names{cell_type});
end
%% 
timeVector = make_time(PV_DEC_binned{1},4,1);
figure; 
plotStruct_DEC = {PV_DEC_norm SOM_DEC_norm Pyr_DEC_norm};
plotStruct_trace = {PV_TraceByArea SOM_TraceByArea Pyr_TraceByArea};

for cell_type = 1:3

subplot(3,3,cell_type);hold on;
fast_errbar(timeVector_trace,[plotStruct_trace{cell_type}{1}./max(abs(plotStruct_trace{cell_type}{1}(30:120,:))) plotStruct_trace{cell_type}{2}./max(abs(plotStruct_trace{cell_type}{2}(30:120,:)))],2,'shaded',true)
fast_errbar(timeVector_trace,[plotStruct_trace{cell_type}{3}./max(abs(plotStruct_trace{cell_type}{3}(30:120,:))) plotStruct_trace{cell_type}{4}./max(abs(plotStruct_trace{cell_type}{4}(30:120,:)))],2,'shaded',true,'color',[0.4 0.4 0.4])
fix_axes(gcf,8,'time (s)','norm pA');
axis square;
subplot(3,3,cell_type+3);hold on;
fast_errbar(timeVector,[plotStruct_DEC{cell_type}{1} plotStruct_DEC{cell_type}{2}],2,'shaded',true)
fast_errbar(timeVector,[plotStruct_DEC{cell_type}{3} plotStruct_DEC{cell_type}{4}],2,'shaded',true,'color',[0.4 0.4 0.4])
fix_axes(gcf,8,'time (s)','norm. deconv.');
ylim([0 1]);
axis square;
subplot(3,3,cell_type+6);hold on;
plotSpread((cellfun(@(x) max(x(35:55,:)),plotStruct_DEC{cell_type},'un',0)),'distributionMarker','o')

fast_errbar(1:4,(cellfun(@(x) max(x(35:55,:)),plotStruct_DEC{cell_type},'un',0)),2,'cells_as_x',true,'stats',true);
fix_axes(gcf,8,'area','late max'); ylim([0 1]);
xticks(1:4);xticklabels({'LM';'AL';'PM';'AM'}); xlim([0.5 4.5]);
axis square;
end

%%
area_i = 2; cell_i = 1;
figure;
subplot(2,1,1);
plot(PV_DEC_norm{area_i}(:,cell_i));
subplot(2,1,2);
plot(PV_TraceByArea{area_i}(:,cell_i))