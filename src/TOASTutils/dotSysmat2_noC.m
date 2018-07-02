function S = dotSysmat2_noC(hmesh, mua, musp, n)
% compute S matrix for CW with kap = 1/(3mus) and 'Contini' Boundary term

kap = 1./(3*musp);
alpha = toastDotBndterm(n,'Contini');
S1 = hmesh.SysmatComponent('PFF',mua);  % absorption term
S2 = hmesh.SysmatComponent('PDD',kap);  % diffusion component
S3 = hmesh.SysmatComponent('BndPFF',ones(size(mua))./(2*alpha)); % boundary term
S = S1 + S2 + S3;