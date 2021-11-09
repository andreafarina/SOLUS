function [babs,bsca,habs,hsca] = getTruthfromLog(MeasName,logbook,spectrabook)
% [babs,bsca,habs,hsca] = getTruthfromLog(MeasName,logbook,spectrabook)
% This function extracts the values of optical properties of a given
% phantom as specified from the file logbook and spectrabook
% MeasName: is a string name of the file containing the data 
% logbook: is a cell of 2 strings, the fist one is the name of the .xlsx
% file, the second is the sheet. The file logbook needs to contain the
% fields MeasureName, LesionType, BreastType
% spectrabook: is a cell of 2 strings, the fist one is the name of the .xlsx
% file containing the truth values for the media, the second is the sheet. The file logbook needs to contain the
% fields "Name" and two fields for each wavelegth written as A+wavlength
% for absorption and S+wavelegth for scattering

[~,MeasName] = fileparts(strrep(strrep(strrep(MeasName,'_Reconstruction',''),'_FirstGate',''),'.mat',''));

tablog = readtable(logbook{1},'Sheet',logbook{2});
idx_row = find(strcmpi(tablog{:,'MeasureName'},MeasName));
incl_name = tablog{idx_row,'LesionType'}{1,1};
bulk_name = tablog{idx_row,'BreastType'}{1,1};

tabspe = readtable(spectrabook{1},'Sheet',spectrabook{2});

lam = [640,670,830,905,930,970,1020,1050];
prop = {'A','S'};
k = 1;
for i = 1:2
    k = 1;
    for l=lam
        if i == 1
            Astrlam{k} = [prop{i},num2str(l)];
        elseif i ==2
            Sstrlam{k} = [prop{i},num2str(l)];
        end
        k = k + 1;
    end
end
%% incl
idxspe = find(strcmpi(tabspe{:,'Name'},incl_name));
habs = tabspe{idxspe,Astrlam}/10;
hsca = tabspe{idxspe, Sstrlam}/10;
%% bulk
idxspe = find(strcmpi(tabspe{:,'Name'},bulk_name));
babs = tabspe{idxspe,Astrlam}/10;
bsca = tabspe{idxspe, Sstrlam}/10;

end