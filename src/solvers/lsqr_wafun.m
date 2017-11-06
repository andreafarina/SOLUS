function y = lsqr_wafun(ops,sqrho,x,transp_flag)
% forward and adjoint function for LSQR for the wavelet ADMM problem

if strcmp(transp_flag,'transp')
    xdat = x(1:ops.ndat);
    xsol = x(ops.ndat+1:end);
    y = ops.ATr(xdat) + sqrho*ops.WTr(xsol);
    y = reshape(y,[],1);
elseif strcmp(transp_flag,'notransp')
    ydat = ops.A(x);
    ysol = sqrho*ops.W(x);
    y = [ydat;ysol];
end