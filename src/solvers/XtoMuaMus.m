function [bmua,bmus] = XtoMuaMus(x,mua0,mus0,type)

switch lower(type)
    case 'mua'
        bmua = x;
        bmus = ones(size(bmua)) * mus0;
    case 'd'
        bmus = 1./(3*x);
        bmua = ones(size(bmus)) * mua0;
    case 'muad'
        bmua = x(1:numel(x)/2);
        bmus = 1./(3*x(numel(x)/2+1:end));
end
