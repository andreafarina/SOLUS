file = 'mus1_05-1_56_benign_1.mat';

%% load
load(file);
%% create TOAST mesh
hmesh = toastMesh('slab_test.msh');
hmesh.Make(test.prm.mesh.toastmeshdata.vtx,...
    test.prm.mesh.toastmeshdata.idx,...
    test.prm.mesh.toastmeshdata.eltp);
%% create Q/M
ds = test.prm.mqwidth;
ds = 1
hmesh.SetQM(test.prm.mesh.QM.Qpos,test.prm.mesh.QM.Mpos);
Q = hmesh.Qvec('Neumann','Gaussian',ds);
M = hmesh.Mvec('Gaussian',ds,0);
%% plot
figure,hmesh.Display(sum(Q,2));
hmesh.delete;
