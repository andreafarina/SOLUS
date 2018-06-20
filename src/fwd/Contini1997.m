function [Rrhot,Trhot,Rrho,Trho,Rt,Tt,lrhoR,lrhoT,R,T,A,Z] = Contini1997(rho,t,s,mua,musp,n1,n2,phantom,DD,m)
%
% [Rrhot,Trhot,Rrho,Trho,Rt,Tt,lrhoR,lrhoT,R,T,A,Z] = Contini1997(rho,t,s,mua,musp,n1,n2,phantom,DD,m)
%
% From:
% Contini D, Martelli F, Zaccanti G.
% Photon migration through a turbid slab described by a model
% based diffusion approximation. I. Theory
% Applied Optics Vol 36, No 19, 1997, pp 4587-4599
%
% phantom: if phantom='semiinf' then semi-infinite medium (please set s=inf)
%          if phantom='slab' then slab
% rho: radial position of the detector at distance s
% t: time (ns)
% s: slab thickness (mm)
% m: maximum number of positive or negative sources
%    m automatically set to 200 if m=[] for a slab
%    and automatically set to m=0 for a semi-infinite medium
% mua: absorption coefficient mm^(-1)
% musp: reduced scattering coefficient
% n1: external medium
% n2: diffusing medium
%
% The diffusion coefficient utilized for the calculation is
% DD: if DD = 'Dmuas' then mua dependent, D = 1/(3*(musp+mua))  (equation (19))
%     if DD = 'Dmus' then                 D = 1/(3*musp)         (page 4588)
% Rrhot: time resolved reflectance mm^(-2) ps^(-1) (equation (36))
%        N x M matrix (N rho values and M t values)
% Trhot: time resolved transmittance mm^(-2) ps^(-1) (equation (39))
%        N x M matrix (N rho values and M t values)
% Rrho: reflectance mm^(-2) (equation (45))
%       N rho values
% Trho: transmittance mm^(-2) (equation (46))
%       N rho values
% Rt: equation (40) ps^(-1)
% Tt: equation (41) ps^(-1)
% lrhoR: equation (47) mm
% lrhoT: equation (48) mm
% A : equations (27) and (29)
% Z : equations (37) for -m:m (first line -m, last line m)
%
% 22/3/2013
% Tiziano BINZONI (University of Geneva)
% Fabrizio MARTELLI (University of Firenze)
% Alessandro TORRICELLI (Politecnico di Milano)
%
% 18/4/2013
% 1) Now it is possible to use m as input of Contini1997
% 2) Contini1997 generates now a warning message if any of the generated
%    variables values changes more than 1e-6 percent when going from m to m+1


if nargin < 10,
    error('number of input arguments for Contini1997 must be 10');
end

if strcmp(phantom,'slab'),
    if isempty(m),
        m=200;
    end
elseif strcmp(phantom,'semiinf'),
    m=0;
    s=1;
else
    error('define phantom type as: ''slab'' or ''semiinf''');
end

t=t*1e-9;
rho=rho*1e-3;
s=s*1e-3;
mua=mua*1e3;
musp=musp*1e3;

% max accepted error on compuzed data
err=1e-6;

% Generates equations (36) and (39)
[Rrhot,Trhot] = RTrhotfunction(rho,t,s,m,mua,musp,n1,n2,DD);
if ~strcmp(phantom,'semiinf'),
    [Rrhot1,Trhot1] = RTrhotfunction(rho,t,s,m+1,mua,musp,n1,n2,DD);
    WarningPrcfunction(Rrhot,Rrhot1,'Rrhot',err);
    WarningPrcfunction(Trhot,Trhot1,'Trhot',err);
end

% Generates equations (45) and (46)
[Rrho,Trho] = RTrhofunction(rho,s,m,mua,musp,n1,n2,DD);
if ~strcmp(phantom,'semiinf'),
    [Rrho1,Trho1] = RTrhofunction(rho,s,m+1,mua,musp,n1,n2,DD);
    WarningPrcfunction(Rrho,Rrho1,'Rrho',err);
    WarningPrcfunction(Trho,Trho1,'Trho',err);
end

% Generates equations (40) and (41)
[Rt,Tt] = RTtfunction(t,s,m,mua,musp,n1,n2,DD);
if ~strcmp(phantom,'semiinf'),
    [Rt1,Tt1] = RTtfunction(t,s,m+1,mua,musp,n1,n2,DD);
    WarningPrcfunction(Rt,Rt1,'Rt',err);
    WarningPrcfunction(Tt,Tt1,'Tt',err);
end

% Generates equations (47) and (48)
[lrhoR,lrhoT] = lrhoTRfunction(rho,s,m,mua,musp,n1,n2,DD);
if ~strcmp(phantom,'semiinf'),
    [lrhoR1,lrhoT1] = lrhoTRfunction(rho,s,m+1,mua,musp,n1,n2,DD);
    WarningPrcfunction(lrhoR,lrhoR1,'lrhoR',err);
    WarningPrcfunction(lrhoT,lrhoT1,'lrhoT',err);
end

% Generates equations (49) and (50)
[R,T] = TRfunction(s,m,mua,musp,n1,n2,DD);
if ~strcmp(phantom,'semiinf'),
    [R1,T1] = TRfunction(s,m+1,mua,musp,n1,n2,DD);
    WarningPrcfunction(R,R1,'R',err);
    WarningPrcfunction(T,T1,'T',err);
end

% Generates equations (27) and (29)
A = Afunction(n1,n2);

% Generates equations (37) for the utilized m values
Z=[];
for i=-m:m,
    [z1,z2,z3,z4]=zfunction(s,i,mua,musp,n1,n2,DD);
    Z=[Z;z1,z2,z3,z4];
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A=Afunction(n1,n2)
% Computes parameter A; equation (27)
%
% n1: external medium
% n2: diffusing medium

% page 4590
n=n2/n1;

if n>1,
    
    % equations (30)
    t1 = 4*(-1-n^2+6*n^3-10*n^4-3*n^5+2*n^6+6*n^7-3*n^8-(6*n^2+9*n^6)*(n^2-1)^(1/2))/...
        (3*n*(n^2-1)^2*(n^2+1)^3);
    t2 = (-8+28*n^2+35*n^3-140*n^4+98*n^5-13*n^7+13*n*(n^2-1)^3*(1-(1/n^2))^(1/2))/...
        (105*n^3*(n^2-1)^2);
    t3 = 2*n^3*(3+2*n^4)*log(...
        ((n-(1+n^2)^(1/2))*(2+n^2+2*(1+n^2)^(1/2))*(n^2+(n^4-1)^(1/2)))/...
        (n^2*(n+(1+n^2)^(1/2))*(-n^2+(n^4-1)^(1/2)))...
        )/...
        ((n^2-1)^2*(n^2+1)^(7/2));
    t4 = ( (1+6*n^4+n^8)*log((-1+n)/(1+n))+4*(n^2+n^6)*log((n^2*(1+n))/(n-1)) )/...
        ((n^2-1)^2*(n^2+1)^3);
    
    % equation (29)
    B = 1+(3/2)*( 2*(1-1/n^2)^(3/2)/3+t1+t2+( (1+6*n^4+n^8)*(1-(n^2-1)^(3/2)/n^3) )/( 3*(n^4-1)^2) +t3 );
    C = 1-( (2+2*n-3*n^2+7*n^3-15*n^4-19*n^5-7*n^6+3*n^7+3*n^8+3*n^9)/(3*n^2*(n-1)*(n+1)^2*(n^2+1)^2) )-t4;
    A = B/C;
    
elseif n==1,
    
    %page 4591
    A=1;
    
else
    
    % equations (28)
    r1 = (-4+n-4*n^2+25*n^3-40*n^4-6*n^5+8*n^6+30*n^7-12*n^8+n^9+n^11)/...
        (3*n*(n^2-1)^2*(n^2+1)^3);
    r2 = (2*n^3*(3+2*n^4))/((n^2-1)^2*(n^2+1)^(7/2))*...
        log( (n^2*(n-(1+n^2)^(1/2)))*(2+n^2+2*(1+n^2)^(1/2))/...
        (n+(1+n^2)^(1/2))/(-2+n^4-2*(1-n^4)^(1/2)) );
    r3 = (4*(1-n^2)^(1/2)*(1+12*n^4+n^8))/...
        (3*n*(n^2-1)^2*(n^2+1)^3);
    r4 = ( (1+6*n^4+n^8)*log((1-n)/(1+n))+4*(n^2+n^6)*log((1+n)/(n^2*(1-n)))  )/...
        ((n^2-1)^2*(n^2+1)^3);
    
    % equation (27)
    A = (1+(3/2)*(8*(1-n^2)^(3/2)/(105*n^3))-(((n-1)^2*(8+32*n+52*n^2+13*n^3))/(105*n^3*(1+n)^2)+r1+r2+r3) )/...
        (1-(-3+7*n+13*n^2+9*n^3-7*n^4+3*n^5+n^6+n^7)/(3*(n-1)*(n+1)^2*(n^2+1)^2)-r4);
    
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Rrhot,Trhot] = RTrhotfunction(rho,t,s,m,mua,musp,n1,n2,DD)
% Computes equations (36) and (39)
%
% rho: radial position of the detector at distance s
% t: time
% s: slab thickness
% m: maximum number of positive or negative sources
% mua: absorption coefficient
% musp: reduced scattering coefficient
% n1: external medium
% n2: diffusing medium

% Generates equations (36) and (39)
Rrhot=[];
Trhot=[];
for i=1:length(rho);
    
    [Rrhottmp,Trhottmp] = RTrhotfunctionPartial(rho(i),t,s,m,mua,musp,n1,n2,DD);
    Rrhot=[Rrhot;Rrhottmp];
    Trhot=[Trhot;Trhottmp];
    
end

if m==0, Trhot=[];end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Rrhot,Trhot] = RTrhotfunctionPartial(rho,t,s,m,mua,musp,n1,n2,DD)

% speed of light
c=299792458; % m s^(-1)

v=c/n2;

% equation (19)
D=Dfunction(DD,mua,musp);

% equations (36) and (39)
Rrhottmp = 0;
Trhottmp = 0;
for i=-m:m,
    [z1,z2,z3,z4]=zfunction(s,i,mua,musp,n1,n2,DD);
    Rrhottmp = Rrhottmp + ...
        z3*exp(-z3^2./(4*D*v*t)) - z4*exp(-z4^2./(4*D*v*t));
    Trhottmp = Trhottmp + ...
        z1*exp(-z1^2./(4*D*v*t)) - z2*exp(-z2^2./(4*D*v*t));
end

Rrhot = -exp(-mua*v*t-rho^2./(4*D*v*t))./(2*(4*pi*D*v)^(3/2)*t.^(5/2)) .* Rrhottmp;
Trhot =  exp(-mua*v*t-rho^2./(4*D*v*t))./(2*(4*pi*D*v)^(3/2)*t.^(5/2)) .* Trhottmp;

Rrhot=Rrhot*1e-6*1e-12;
Trhot=Trhot*1e-6*1e-12;

Rrhot(t<=0)=0;
Trhot(t<=0)=0;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Rrho,Trho] = RTrhofunction(rho,s,m,mua,musp,n1,n2,DD)
% Computes equations (45) and (46)
%
% rho: radial position of the detector at distance s
% s: slab thickness
% m: maximum number of positive or negative sources
% mua: absorption coefficient
% musp: reduced scattering coefficient
% n1: external medium
% n2: diffusing medium

% equation (19)
D=Dfunction(DD,mua,musp);

% equation (45) and (46)
Rrhotmp=0;
Trhotmp=0;
for i=-m:m,
    [z1,z2,z3,z4]=zfunction(s,i,mua,musp,n1,n2,DD);
    Rrhotmp=Rrhotmp+...
        z3*( D^(-1/2)*mua^(1/2)*(rho.^2+z3^2).^(-1) + (rho.^2+z3^2).^(-3/2) ).*exp( -(mua*(rho.^2+z3^2)/D).^(1/2) )-...
        z4*( D^(-1/2)*mua^(1/2)*(rho.^2+z4^2).^(-1) + (rho.^2+z4^2).^(-3/2) ).*exp( -(mua*(rho.^2+z4^2)/D).^(1/2) );
    Trhotmp=Trhotmp+...
        z1*( D^(-1/2)*mua^(1/2)*(rho.^2+z1^2).^(-1) + (rho.^2+z1^2).^(-3/2) ).*exp( -(mua*(rho.^2+z1^2)/D).^(1/2) )-...
        z2*( D^(-1/2)*mua^(1/2)*(rho.^2+z2^2).^(-1) + (rho.^2+z2^2).^(-3/2) ).*exp( -(mua*(rho.^2+z2^2)/D).^(1/2) );
end
Rrho=-1/(4*pi)*Rrhotmp;
Trho= 1/(4*pi)*Trhotmp;

Rrho=Rrho*1e-6;
Trho=Trho*1e-6;

if m==0, Trho=[];end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Rt,Tt] = RTtfunction(t,s,m,mua,musp,n1,n2,DD)
% Computes equations (40) and (41)
%
% t: time
% s: slab thickness
% m: maximum number of positive or negative sources
% mua: absorption coefficient
% musp: reduced scattering coefficient
% n1: external medium
% n2: diffusing medium

% speed of light
c=299792458; % m s^(-1)

v=c/n2;

% equation (19)
D=Dfunction(DD,mua,musp);

% equations (40) and (41)
Rttmp=0;
Tttmp=0;
for i=-m:m,
    [z1,z2,z3,z4]=zfunction(s,i,mua,musp,n1,n2,DD);
    Rttmp=Rttmp+z3*exp(-z3^2./(4*D*v*t))-z4*exp(-z4^2./(4*D*v*t));
    Tttmp=Tttmp+z1*exp(-z1^2./(4*D*v*t))-z2*exp(-z2^2./(4*D*v*t));
end
Rt=-exp(-mua*v*t)./(2*(4*pi*D*v)^(1/2)*t.^(3/2)).*Rttmp;
Tt= exp(-mua*v*t)./(2*(4*pi*D*v)^(1/2)*t.^(3/2)).*Tttmp;

Rt=Rt*1e-12;
Tt=Tt*1e-12;

if m==0, Tt=[];end

Tt(t<=0)=0;
Rt(t<0)=0;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lrhoR,lrhoT] = lrhoTRfunction(rho,s,m,mua,musp,n1,n2,DD)
% Computes equations (47) and (48)
% These equations are valid only for mua>0
%
% rho: radial position of the detector at distance s
% s: slab thickness
% m: maximum number of positive or negative sources
% mua: absorption coefficient
% musp: reduced scattering coefficient
% n1: external medium
% n2: diffusing medium

% equation (19)
D=Dfunction(DD,mua,musp);

% equations (45) and (46)
[Rrho,Trho] = RTrhofunction(rho,s,m,mua,musp,n1,n2,DD);
Rrho=Rrho*1e6;
Trho=Trho*1e6;

% equations (47) and (48)
if mua<=0,
    
    lrhoR=[];
    lrhoT=[];
    
else
    
    lrhoRtmp=0;
    lrhoTtmp=0;
    for i=-m:m,
        [z1,z2,z3,z4]=zfunction(s,i,mua,musp,n1,n2,DD);
        lrhoRtmp=lrhoRtmp+...
            z3*(rho.^2+z3^2).^(-1/2).*exp(-(mua*(rho.^2+z3^2)/D).^(1/2))-...
            z4*(rho.^2+z4^2).^(-1/2).*exp(-(mua*(rho.^2+z4^2)/D).^(1/2));
        lrhoTtmp=lrhoTtmp+...
            z1*(rho.^2+z1^2).^(-1/2).*exp(-(mua*(rho.^2+z1^2)/D).^(1/2))-...
            z2*(rho.^2+z2^2).^(-1/2).*exp(-(mua*(rho.^2+z2^2)/D).^(1/2));
    end
    lrhoR=(-1./(8*pi*D*Rrho)).*lrhoRtmp;
    if ~isempty(Trho),
        lrhoT= (1./(8*pi*D*Trho)).*lrhoTtmp;
    else
        lrhoT=[];
    end
    
    lrhoR=lrhoR*1e3;
    lrhoT=lrhoT*1e3;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R,T] = TRfunction(s,m,mua,musp,n1,n2,DD)
% Computes equations (49) and (50)
% These equations are valid only for mua>0
%
% s: slab thickness
% m: maximum number of positive or negative sources
% mua: absorption coefficient
% musp: reduced scattering coefficient
% n1: external medium
% n2: diffusing medium

% equation (19)
D=Dfunction(DD,mua,musp);

% equations (49) and (50)
if mua>0,
    
    Rtmp=0;
    Ttmp=0;
    for i=-m:m,
        [z1,z2,z3,z4]=zfunction(s,i,mua,musp,n1,n2,DD);
        Rtmp=Rtmp+...
            (sign(z3)*exp(-(mua/D)^(1/2)*abs(z3))-sign(z4)*exp(-(mua/D)^(1/2)*abs(z4)));
        Ttmp=Ttmp+...
            (sign(z1)*exp(-(mua/D)^(1/2)*abs(z1))-sign(z2)*exp(-(mua/D)^(1/2)*abs(z2)));
    end
    R=-(1/2)*Rtmp;
    T= (1/2)*Ttmp;
    
else
    
    R=[];
    T=[];
end

if m==0, T=[];end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z1,z2,z3,z4]=zfunction(s,m,mua,musp,n1,n2,DD)
% Computes equation (37)
%
% s: slab thickness
% m: maximum number of positive or negative sources
% mua: absorption coefficient
% musp: reduced scattering coefficient
% n1: external medium
% n2: diffusing medium

% page 4592
z0 = 1/musp;
%z0 = 1/(musp+mua);

% equation (27)
A = Afunction(n1,n2);
% A=504.332889-2641.0021*(n2/n1)+5923.699064*(n2/n1)^2-...
%     7376.355814*(n2/n1)^3+5507.53041*(n2/n1)^4-...
%     2463.357945*(n2/n1)^5+610.956547*(n2/n1)^6-...
%     64.8047*(n2/n1)^7;

% equation (19)
D=Dfunction(DD,mua,musp);

%page 4592
ze = 2*A*D;

z1 = s*(1-2*m) - 4*m*ze - z0;
z2 = s*(1-2*m) - (4*m-2)*ze + z0;
z3 = -2*m*s -4*m*ze - z0;
z4 = -2*m*s -(4*m-2)*ze + z0;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D=Dfunction(DD,mua,musp)
% equation (19)

if strcmp(DD,'Dmuas'),
    D = 1/(3*(musp+mua));
elseif strcmp(DD,'Dmus'),
    D = 1/(3*(musp));
else
    error('DD must be equal to ''Dmuas'' or ''Dmus''');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function WarningPrcfunction(x,y,label,err)
%size(x)
%size(y)
if ~isempty(x) && ~isempty(y),
    prcMax=max(max(abs( (x-y)./x*100)));
    if prcMax>err,
        warning(['increase m (error larger than ',num2str(err),'% for ',label,')']);
    end
end

end
