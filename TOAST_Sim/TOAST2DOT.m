%%
clearvars;
close all
%%
prefix = 'mus1_05-1_56_';
d = dir([prefix,'*.mat']);
Nm = 8;
Nq = 8;
dmask = logical(1 - eye(Nm,Nq));
% dmask = false(Nm,Nq); 
% dmask(1,6) = true;
for file = {d.name}
    strfile = char(file);
    load(strfile);
    RefTD = zeros(test.prm.n_step,Nm*Nq);
    DataTD = zeros(test.prm.n_step,Nm*Nq);
    l = 1;
    for fn = fieldnames(test.results.gammahom)'
        RefTD(:,l:l+Nm - 1) = full(test.results.gammahom.(fn{1}));
        DataTD(:,l:l+Nm - 1) = full(test.results.gamma.(fn{1}));
        l = l + Nq;
    end
    RefCW = squeeze(sum(RefTD));
    DataCW = squeeze(sum(DataTD));
    figure,imagesc(reshape((DataCW-RefCW),[Nm,Nq]));
    pause(0.1)
    RefTD = RefTD(:,dmask(:));
    DataTD = DataTD(:,dmask(:));
    RefCW = RefCW(dmask(:));
    DataCW = DataCW(dmask(:));
    
    strfile = char(file);
    save([strfile(1:end-4),'_singl_FwdTeo'],'RefTD','DataTD','RefCW','DataCW','test');
end