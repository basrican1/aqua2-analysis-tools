function mOut = bwaBuildCellMetricOutFromTable(T, metricCol, epochTable)
% Convert a cell-level metric table into canonical-like summary/event structs.

req = {'epoch_idx','epoch_name'};
if ~all(ismember(req, T.Properties.VariableNames))
    error('Cell metric table missing epoch columns.');
end
if ~ismember('primary_cell_uid', T.Properties.VariableNames)
    if ismember('cell_uid', T.Properties.VariableNames)
        T.primary_cell_uid = string(T.cell_uid);
    else
        error('Cell metric table missing primary_cell_uid/cell_uid.');
    end
end
if ~ismember('cell_uid', T.Properties.VariableNames)
    T.cell_uid = string(T.primary_cell_uid);
end
mOut = bwaBuildMetricOutFromTable(T, metricCol, epochTable, 'entityKeys', {'primary_cell_uid'});
end
