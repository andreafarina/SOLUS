d = dir('*.mat');
prefix = 'mus1_05-1_56';
for file = {d.name}
    strfile = char(file);
    load(strfile);
    if isfield(test.prm.mesh,'toastmesh')
        test.prm.mesh = rmfield(test.prm.mesh,'toastmesh');
    end
    disp(['mua = ',num2str(struct2array(test.prm.mua))]);
    disp(['mus = ',num2str(struct2array(test.prm.mus))]);
    
    newstr = [prefix,'_',strfile(23:end)];
    save(newstr,'test');
end