function [J_semi_PCBC] = J_Semi_Inf_PCBC_CW_3D(mua,mus,A,ri,rj,XX,YY,ZZ,vi)
% function [J_semi_PCBC] = J_Semi_Infinite_PCBC_TR_3D(mua,mus,A,ri,rj,XX,YY,ZZ,vi)
% CORRETTO DA ANDREA FARINA
% Jacobian for the Semi-infinite homogeneous space geometry
% The Jacobian is calculated for the 'perturbation Contrast, i.e. dR/R'
% with R DE green's function in the Time Domain for the reflectance
% The Jacobian is calculated for an inclusion placed in the medium where is
% assumed an absorption variation dmua with respect to the background
% If C=dR/R is the contrast, the Jacobian J is the coefficient J=dC/dmua 
% for the flux (Reflectance) calculated with the PCBC
% mua absorption coefficient (mm^-1)
% mus reduced scattering coefficient (mm^-1)
% cs speed of light (mm/ps)
% A factor that accounts Fresnel reflections
% kap diffusion coefficient
% vi volume of the inclusion (mm^3)
% t time (ps)
% ri position vector of the source (mm)
% rj position vector of the detector (mm)
% rk position vector of the inclusion (mm)
% rhoik distance between source and inclusion
% rhojk distance between detector and inclusion
% rhoji distance between detector and source
% It is used the the Born approximation applied to
% the DE solved with the Extrapolated Boundary Condition
% G0_semi_PCBC: Unperturbed Fluence, Reflectance=(1/2A)*Fluence
% dG_semi_PCBC: Perturbation for the Fluence, Reflectance=(1/2A)*Fluence
% G_semi_PCBC: Contrast or relative perturbation for fluence and
% refletcance
%---------------------------------------------------------------

kap = 1/(3*(mus));
mueff=sqrt(3*mua*mus);
ze=2*A*kap;
n=numel(XX);

 if ri(3) > 0,
    z0=ri(3);     
 elseif ri(3)==0,
    z0=1/mus;
    ri(3)=z0;
 end
      

z12plus=z0;
z12minus=-2*ze-z0;
%z23plus=rk(3);
z23plus=ZZ;
%z23minus=-2*ze-rk(3);
z23minus=-2*ze-ZZ;

rhoij_plus=sqrt((ri-rj)*(ri-rj)');
ri_min=ri;ri_min(3)=z12minus;
rhoij_minus=sqrt((ri_min-rj)*(ri_min-rj)');

r12plus=sqrt((XX-ri(1)).^2+(YY-ri(2)).^2+(ZZ-z12plus).^2);
r12minus=sqrt((XX-ri(1)).^2+(YY-ri(2)).^2+(ZZ-z12minus).^2);
r23plus=sqrt((rj(1)-XX).^2+(rj(2)-YY).^2+(rj(3)-z23plus).^2);
r23minus=sqrt((rj(1)-XX).^2+(rj(2)-YY).^2+(rj(3)-z23minus).^2);

rho12plus=reshape(r12plus,[1 n]);
rho12minus=reshape(r12minus,[1 n]);
rho23plus=reshape(r23plus,[1 n]);
rho23minus=reshape(r23minus,[1 n]);

%fact=exp(rhoij.*mueff);

J_semi_PCBC=-1/(4*pi*kap)^2.*(...
    +exp(-mueff*(rho12plus+rho23plus))./(rho12plus.*rho23plus)...
    -exp(-mueff*(rho12plus+rho23minus))./(rho12plus.*rho23minus)...
    -exp(-mueff*(rho12minus+rho23plus))./(rho12minus.*rho23plus)...
    +exp(-mueff*(rho12minus+rho23minus))./(rho12minus.*rho23minus));

%G0_semi_PCBC=exp(-mueff*rhoij_plus)./rhoij_plus-exp(-mueff*rhoij_minus)./rhoij_minus;

%J_semi_PCBC=dG_semi_PCBC./(2*A)%G0_semi_PCBC;
%G0_semi_PCBC/(2*A)
J_semi_PCBC=vi*J_semi_PCBC./(2*A);%./G0_semi_PCBC;
J_semi_PCBC(isnan(J_semi_PCBC)) = 0;
%J_semi_PCBC(J_semi_PCBC<0.01)=0;

end



