function [time_domain, meshdata, DataTD, RefTD, DataCW, RefCW, nominalcoefficients ,...
    refindex, Muatoast, Mustoast, prior] = toast2dot(toastfilename)
%-
%-
%-
%-
%-
% [time_domain, meshdata, DataTD, RefTD, nominalcoefficients , Muatoast, Mustoast, prior] = toast2dot(toastfilename)
% meshdata: meshdata.vtx, meshdata.idx, meshdata.eltp
% nominalcoefficients =  [mua_i, mus_i, mua_o, mus_o]
% time_domain = [dt, steps]

%load([cd ,'/', toastfilename]);
load(toastfilename);


Nm = 8;
Nq = 8;
dmask = logical(1 - eye(Nm,Nq));
% dmask = false(Nm,Nq); 
% dmask(1,6) = true;


RefTD = zeros(test.prm.time.n_step,Nm*Nq);
DataTD = zeros(test.prm.time.n_step,Nm*Nq);
        
l = 1;
for fn = fieldnames(test.results.gamma)'        
    RefTD(:,l:l+Nm - 1) = full(test.results.gammahom.(fn{1}));
    DataTD(:,l:l+Nm - 1) = full(test.results.gamma.(fn{1}));
    l = l + Nq;
end


RefCW = squeeze(sum(RefTD));
DataCW = squeeze(sum(DataTD));
%       figure,imagesc(reshape((DataCW-RefCW),[Nm,Nq]));
RefTD = RefTD(:,dmask(:));
DataTD = DataTD(:,dmask(:));
RefCW = RefCW(dmask(:));
DataCW = DataCW(dmask(:));
%sdCW = sqrt(DataCW);
%sdTD = sqrt(DataTD);

dt = test.prm.time.dt;
nsteps = test.prm.time.n_step;
prior = test.prm.mesh.toastmeshdata.mask_toast;
meshdata.vtx = test.prm.mesh.toastmeshdata.vtx;
meshdata.idx =  test.prm.mesh.toastmeshdata.idx;
meshdata.eltp = test.prm.mesh.toastmeshdata.eltp;
Muatoast = test.prm.mesh.toastmeshdata.Mua;
Mustoast = test.prm.mesh.toastmeshdata.Mus;
refindex = test.prm.ref.out;
nominalcoefficients = [test.prm.mua.in, test.prm.mus.in, test.prm.mua.out, test.prm.mua.out];
time_domain = [dt, nsteps];

end