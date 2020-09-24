% ccc;
clearvars -except PV_* SOM_*
[~,sheet_names] = xlsfinfo('HVA_IN_pyr_pair_connectivity.xlsx');
disp(sheet_names);
u_input = input('Choose data set type: ');

summary_table = readtable('HVA_IN_pyr_pair_connectivity.xlsx','Sheet',u_input);
area_names = [{'LM'},{'AL'},{'PM'},{'AM'}];
try 
    do_analysis = logical(cell2mat(cellfun(@isempty,summary_table.Analyzed,'un',0)));
catch
    do_analysis = isnan(summary_table.Analyzed);
end

if any(~do_analysis)
    temp_disp = summary_table.Date(~do_analysis);
    disp(temp_disp(:));
    disp('found analysis files')
end

arrayfun(@(x) JL_2020_pairStim_connectivityParams(x,u_input),find(do_analysis),'un',0);

%%
summary_table = readtable('HVA_IN_pyr_pair_connectivity.xlsx','Sheet',u_input);

for area_i = 1:numel(area_names)
    areaIdx = find(cell2mat(cellfun(@(x) strcmp(area_names{area_i},x),summary_table.Area,'un',0)));
    nCells(area_i) = numel(areaIdx);
    if nCells(area_i) > 0
        for cell = 1:nCells(area_i)
            group_data{area_i}(cell) = load([summary_table.Analyzed{areaIdx(cell)}, '.mat']);
        end
        pyrIN_connect{area_i} = cell2mat(cellfun(@(x) ~isempty(x),summary_table.Pyr__IN(areaIdx),'un',0));
        INpyr_connect{area_i} = cell2mat(cellfun(@(x) ~isempty(x),summary_table.IN__Pyr(areaIdx),'un',0));
        recip_connect{area_i} = cell2mat(cellfun(@(x) ~isempty(x),summary_table.Reciprocal(areaIdx),'un',0));
    end
    tempCat = vertcat(group_data{area_i}.data);
    IN_amp_group{area_i} = vertcat(tempCat.IN_EPSC_amp);
    pyr_amp_group{area_i} = vertcat(tempCat.pyr_IPSC_amp);
        
    IN_amp_connected{area_i} = IN_amp_group{area_i}(pyrIN_connect{area_i},:);
    pyr_amp_connected{area_i} = pyr_amp_group{area_i}(INpyr_connect{area_i},:);
    
    pyr_amp_norm{area_i} = pyr_amp_connected{area_i}./pyr_amp_connected{area_i}(:,1);
    IN_amp_norm{area_i} = IN_amp_connected{area_i}./IN_amp_connected{area_i}(:,1);
    
    recip_connect_IN_amp{area_i} = IN_amp_group{area_i}(recip_connect{area_i},:);
    recip_connect_pyr_amp{area_i} = pyr_amp_group{area_i}(recip_connect{area_i},:);

end

pyrIN_count_byArea = cell2mat(cellfun(@sum,pyrIN_connect,'un',0));
pyrIN_byArea = pyrIN_count_byArea./nCells;
pyrIN_lat = sum(pyrIN_count_byArea(1:2))/sum(nCells(1:2));
pyrIN_med = sum(pyrIN_count_byArea(3:4))/sum(nCells(3:4));
pyrIN_pos = sum(pyrIN_count_byArea([1,3]))/sum(nCells([1,3]));
pyrIN_ant = sum(pyrIN_count_byArea([2,4]))/sum(nCells([2,4]));

INpyr_count_byArea = cell2mat(cellfun(@sum,INpyr_connect,'un',0));
INpyr_byArea = INpyr_count_byArea./nCells;
INpyr_lat = sum(INpyr_count_byArea(1:2))/sum(nCells(1:2));
INpyr_med = sum(INpyr_count_byArea(3:4))/sum(nCells(3:4));
INpyr_pos = sum(INpyr_count_byArea([1,3]))/sum(nCells([1,3]));
INpyr_ant = sum(INpyr_count_byArea([2,4]))/sum(nCells([2,4]));

recip_count_byArea = cell2mat(cellfun(@sum,recip_connect,'un',0));
recip_byArea = recip_count_byArea./nCells;
recip_lat = sum(recip_count_byArea(1:2))/sum(nCells(1:2));
recip_med = sum(recip_count_byArea(3:4))/sum(nCells(3:4));
recip_pos = sum(recip_count_byArea([1,3]))/sum(nCells([1,3]));
recip_ant = sum(recip_count_byArea([2,4]))/sum(nCells([2,4]));

%%

IN_PPR_group = cellfun(@(x) x(:,2)./x(:,1),IN_amp_connected,'un',0);
pyr_PPR_group = cellfun(@(x) x(:,2)./x(:,1),pyr_amp_connected,'un',0);
IN_P1P10_group = cellfun(@(x) x(:,10)./x(:,1),IN_amp_connected,'un',0);
pyr_P1P10_group = cellfun(@(x) x(:,10)./x(:,1),pyr_amp_connected,'un',0);

figure;
subplot(1,2,1);
fast_errbar(1:4,IN_PPR_group,1,'continuous',false,'stats',true);
xticks(1:4);xticklabels(area_names);
fix_axes(gcf,10,'Area','PPR'); ylim([0 2]);
title('IN PPR');
subplot(1,2,2);
fast_errbar(1:4,pyr_PPR_group,1,'continuous',false,'stats',true);
fix_axes(gcf,10,'Area','PPR');
xticks(1:4);xticklabels(area_names);
title('pyr PPR'); ylim([0 2]);

figure;
subplot(1,2,1);
fast_errbar(1:4,IN_P1P10_group,1,'continuous',false,'stats',true);
xticks(1:4);xticklabels(area_names);
fix_axes(gcf,10,'Area','PPR'); ylim([0 2]);
title('IN PPR');
subplot(1,2,2);
fast_errbar(1:4,pyr_P1P10_group,1,'continuous',false,'stats',true);
fix_axes(gcf,10,'Area','PPR');
xticks(1:4);xticklabels(area_names);
title('pyr PPR'); ylim([0 2]);
%%
SOM_IN_pyr_PPR = pyr_PPR_group;
SOM_pyr_IN_PPR = IN_PPR_group;
SOM_IN_pyr_P1P10 = pyr_P1P10_group;
SOM_pyr_IN_P1P10 = IN_P1P10_group;
% 
% PV_IN_pyr_PPR = pyr_PPR_group;
% PV_pyr_IN_PPR = IN_PPR_group;
% PV_IN_pyr_P1P10 = pyr_P1P10_group;
% PV_pyr_IN_P1P10 = IN_P1P10_group;

figure;
subplot(2,2,1)
fast_errbar(1:2,[{cell2mat(PV_pyr_IN_PPR')} {cell2mat(SOM_pyr_IN_PPR')}],1,'continuous',false,'stats',true)
subplot(2,2,2)
fast_errbar(1:2,[{cell2mat(PV_IN_pyr_PPR')} {cell2mat(SOM_IN_pyr_PPR')}],1,'continuous',false,'stats',true)
subplot(2,2,3)
fast_errbar(1:2,[{cell2mat(PV_pyr_IN_P1P10')} {cell2mat(SOM_pyr_IN_P1P10')}],1,'continuous',false,'stats',true)
subplot(2,2,4)
fast_errbar(1:2,[{cell2mat(PV_IN_pyr_P1P10')} {cell2mat(SOM_IN_pyr_P1P10')}],1,'continuous',false,'stats',true)
