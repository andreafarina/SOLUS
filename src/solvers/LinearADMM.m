function [fout,cout,objFunVal] = LinearADMM(ops,g,alpha,K,tau,rho,LSQRtol,LSQRit,FLAGS)
% LinearADMM minimizes the objective function 
% 1/2 * || A f - g ||_2^2 + alpha || W f ||_1
% where W is the Haar wavelet transform by appling ISTA to
% 1/2 * || A W^{-1} c - g ||_2^2 + alpha || c ||_1
%
% falpha = LinearADMM(ops,g,alpha,K,tau,rho,tol,maxit,FLAGS)
%
%  INPUT:
%   ops : Structure for operators
%   ops.A : function A(x) for forward operator
%   ops.ATr : function ATr(x) for adjoint forward operator
%   ops.W : function W(x) for transform operator
%   ops.WTr : function WTr(x) for adjoint transform operator
%   ops.Winv : function Winc(x) for inverse transform operator
%   ops.B : function B(x)mapping from image to solution basis.
%   ops.ndat : length of data
%   ops.nsol : length of solution
%   ops.Dims : dimensions of image basis
%   ops.ntc  : number of coefficients after transform
%
%   g - data, as 1D vector
%   alpha - the regularization parameter
%   K - number of iterations
%   tau - the stepsize parameter
%   rho - the ADMM equality penalty 
%   tol - tolerance for primal solver
%   maxit - maximum iterations for primal solver
%  FLAGS :
%    pos - impose positivity (image space) on each iteration
%    FISTA - use FISTA acceleration
%    Wavelets - use Wavelets. If not, use pixels.
%    Iterates - output all iterations, if not, just the last.
%  OUTPUTS:
%   fout - the solution in 'sol' basis 
%   cout - the Haar wavelet coefficients of the image
%   objFunVAl - value of the objective function at the solution falpha
%
% Author: Felix Lucka
% Extended to DOT by Simon Arridge

% determine if the additional work of computing the objective function is
% necessary 
computeObjVal = false; 
if(nargout > 2)
   computeObjVal = true; 
   objFunVal = zeros(K,1);
end

if(FLAGS.pos)
    disp('Positivity Constrained');
end
if (FLAGS.Wavelets)
    disp('Using Wavelet Transform');
else
    disp('Using Pixel Basis');
end

ndat = ops.ndat;
nsol = ops.nsol;
% get the size of the image and check if it is square
bd = ops.Dims;
Nx = bd(1); Ny= bd(2);
disp(['Image basis ',num2str(Nx),'X',num2str(Ny)])

if(FLAGS.Wavelets & (Nx ~= Ny))
   error('wavelet debluring requires images with Nx=Ny');
end
if(FLAGS.Iterates)
    disp('Outputing all iterations');
    cout = zeros(ops.ntc,K);
    fout = zeros(nsol,K);
end
% define some function handles

    W    = @(x) ops.W(x);
    Winv = @(c) ops.Winv(c);
    WTr  = @(c) ops.WTr(c);
    A    = @(x) ops.A(x);
    ATr  = @(y) ops.ATr(y);
    B    = @(x) ops.B(x);

S   = @(x,alpha) max(abs(x)-alpha,0).*sign(x);


% start the iteration
sqrho = sqrt(rho);
v = zeros(ops.ntc,1);
w = v;
%gd = ATr(g);  % gradient descent direction
xk = ops.Binv(zeros(nsol,1));
tic;
for i = 1:K
    % primal step
    
    xk = lsqr(@(x,tflag)lsqr_wafun(ops,sqrho,x,tflag),[g ; sqrho*(v-w) ],LSQRtol,LSQRit,[],[],xk);
    xkb = B(xk); % this cuts exterior pixels. A bit clunky...
    
    % positivity
        
    xkb(find(xkb < 0)) = 0;
    
    % dual step
    ww = W(xk);
%    v = SoftThresh(ww+w,alpha/rho);
    v = S(ww+w,alpha/rho);
    % update step
    w = w + ww - v;
    if(computeObjVal)
        % compute objective function
        objFunVal(i) = 1/2 * sum(sum((A(xk) - g).^2)) + alpha * sum(abs(ww));
    end
    if(FLAGS.Iterates)
        %        fout(:,i) = ops.hBasis.Map('B->S',Winv(c));
        fout(:,i) = xkb;
        cout(:,i) = ww;
    end
end

% retransform the optimal coefficients to get the image back
if(~FLAGS.Iterates)
    cout = ww;
%    fout = ops.hBasis.Map('B->S',Winv(cout));
    fout = B(Winv(cout));
end
end