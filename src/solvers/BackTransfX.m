function y = BackTransfX(x,x0,xtransf)

switch lower(xtransf)
    case 'log(x/x0)'
        y = x0.*exp(x);
    case 'log(x)'
        y = exp(x);
    case '(x/x0)'
        y = x.*x0;
    case 'x'
        y = x;
    otherwise
        disp('error in BackTransfX');
end
end