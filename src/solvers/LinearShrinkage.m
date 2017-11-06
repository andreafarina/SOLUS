function [fout,cout,objFunVal] = LinearShrinkage(ops,g,alpha,K,tau0,FLAGS)
% DOTWAVDECONVBYISTA minimizes the objective function 
% 1/2 * || A f - g ||_2^2 + alpha || W f ||_1
% where W is the Haar wavelet transform by appling ISTA to
% 1/2 * || A W^{-1} c - g ||_2^2 + alpha || c ||_1
%
% falpha = LinearShrinkage(ops,g,alpha,K,tau,FLAGS)
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


if FLAGS.FISTA
    disp('FISTA Shrinkage');
else
    disp('ISTA Shrinkage');
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
if (FLAGS.Rayleighstep)
    disp('step length uses Rayleigh criteria');
end
% define some function handles

    W    = @(x) ops.W(x);
    Winv = @(c) ops.Winv(c);
    WTr  = @(c) ops.WTr(c);
    A    = @(x) ops.A(x);
    ATr  = @(y) ops.ATr(y);
    B    = @(x) ops.B(x);

S   = @(x,alpha) max(abs(x)-alpha,0).*sign(x);

% initialize inner variables
c = zeros(ops.ntc,1);
if FLAGS.FISTA
    c0  = c;
    t1 = 1;
end
% start the iteration
tau = tau0;
for i=1:K
    % compute the gradient
    p = W(ATr(A(Winv(c)) - g));
    % apply the soft thresholding
    if FLAGS.Rayleighstep;
        tau = (norm(Winv(p))/norm(A(Winv(p) )))^2;
    end
    if FLAGS.FISTA
        c1 = S(c - tau * p,alpha * tau);
        % update t
        t2 = (1+sqrt(1+4*t1^2))/2;
        % update y
        c = c1 + (t1-1)/t2*(c1-c0);
        c0 = c1;
        t1 = t2;
    else
        c = S(c - tau * p,alpha * tau);
    end
    if(FLAGS.pos)  % impose positivity at each iteration.
        ff = Winv(c); ff(find(ff < 0)) = 0;
        c = W(ff);
    end
    if(computeObjVal)
        % compute objective function
        objFunVal(i) = 1/2 * sum(sum((A(Winv(c)) - g).^2)) + alpha * sum(abs(c));
    end
    if(FLAGS.Iterates)
%        fout(:,i) = ops.hBasis.Map('B->S',Winv(c));
        fout(:,i) = B(Winv(c));
        cout(:,i) = c;
    end
end

% retransform the optimal coefficients to get the image back
if(~FLAGS.Iterates)
    cout = c;
%    fout = ops.hBasis.Map('B->S',Winv(cout));
    fout = B(Winv(cout));
end
end