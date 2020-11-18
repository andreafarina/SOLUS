function DispQuantifyTable(Q,REC,FieldsToDisplay)
StructFieldsToDisplay = cellfun(@(Str) strsplit(Str,'.'),FieldsToDisplay,'UniformOutput',false);
FieldNames = fieldnames(Q);
MergedTable = [];
for ifn = 1:numel(FieldNames)
    if contains(FieldNames{ifn},'mu')
        FieldVarNames = REC.radiometry.lambda';
        FieldVarNames = arrayfun(@(il) strcat(FieldNames{ifn},num2str(FieldVarNames(il))),1:REC.radiometry.nL,'UniformOutput',false);
        nFields = REC.radiometry.nL;
    end
    if strcmpi(FieldNames{ifn},'cromo')
        FieldVarNames = REC.spe.cromo_label;
        nFields = REC.spe.nCromo;
    end
    if strcmpi(FieldNames{ifn},'a')||strcmpi(FieldNames{ifn},'b')
        FieldVarNames = FieldNames(ifn);
        nFields = 1;
    end
    Data = zeros(numel(FieldsToDisplay),nFields);
    for ifd = 1:numel(FieldsToDisplay)
        ActFields = StructFieldsToDisplay{ifd};
        if numel(StructFieldsToDisplay{ifd}) == 1
            Data(ifd,:) = arrayfun(@(il) Q.(FieldNames{ifn})(il).(ActFields{1})',1:nFields);
        elseif numel(StructFieldsToDisplay{ifd}) == 2
            Data(ifd,:) = arrayfun(@(il) Q.(FieldNames{ifn})(il).(ActFields{1}).(ActFields{2})',1:nFields);
        elseif numel(StructFieldsToDisplay{ifd}) == 3
            Data(ifd,:) = arrayfun(@(il) Q.(FieldNames{ifn})(il).(ActFields{1}).(ActFields{2}).(ActFields{3})',1:nFields);
        end
    end
    Table = table;
    for ic = 1:nFields
        Table(:,ic) = table(Data(:,ic));
    end
    Table.Properties.RowNames = FieldsToDisplay;
    Table.Properties.VariableNames = FieldVarNames;
    if any(strcmpi(FieldNames,'a'))&&any(strcmpi(FieldNames,'cromo'))&&(strcmpi(FieldNames{ifn},'a')||strcmpi(FieldNames{ifn},'b')||strcmpi(FieldNames{ifn},'cromo'))
        if ~isempty(MergedTable)
            MergedTable = join(MergedTable,Table,'Keys','Row');
        else
            MergedTable = Table;
        end
        if ifn == numel(FieldNames)
            disp(MergedTable);
        end
    else
        disp(Table);
    end
end
end