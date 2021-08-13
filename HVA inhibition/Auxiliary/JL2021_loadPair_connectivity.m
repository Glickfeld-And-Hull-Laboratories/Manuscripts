% analyze all files in the excel sheet

ccc;
[~,sheet_names] = xlsfinfo([your_read_path,'HVA_IN_pyr_pair_connectivity.xlsx']);
disp(sheet_names);

for u_input = 1:2
summary_table = readtable([your_read_path,'HVA_IN_pyr_pair_connectivity.xlsx'],'Sheet',u_input);
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

arrayfun(@(x) JL2021_pairStim_connectivityParams(x,u_input,your_read_path,your_save_path),find(do_analysis),'un',0);
end
%%
ccc;
cell_names = {'PV';'SOM'};
area_names = {'LM';'AL';'PM';'AM'};
for u_input = 1:2
summary_table = readtable([your_read_path,'HVA_IN_pyr_pair_connectivity.xlsx'],'Sheet',u_input);

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

IN_P1P2_group = cellfun(@(x) x(:,2)./x(:,1),IN_amp_connected,'un',0);
pyr_P1P2_group = cellfun(@(x) x(:,2)./x(:,1),pyr_amp_connected,'un',0);
IN_P1P10_group = cellfun(@(x) x(:,10)./x(:,1),IN_amp_connected,'un',0);
pyr_P1P10_group = cellfun(@(x) x(:,10)./x(:,1),pyr_amp_connected,'un',0);
pyr_amp_group = cellfun(@(x) x(:,1),pyr_amp_connected,'un',0);
IN_amp_group = cellfun(@(x) x(:,1),IN_amp_connected,'un',0);
clear group_data summary_table tempCat
save([your_save_path,cell_names{u_input},'_grouped.mat']);
end