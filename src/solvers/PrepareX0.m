function x0 = PrepareX0(mu,n,type)

switch lower(type)
    case 'mua'
        x0 = mu(1) * ones(n,1);
    case 'd'
        x0 = mu(2) * ones(n,1);
    case 'muad'
        x0 = [mu(1)*ones(n,1);mu(2)*ones(n,1)];
end

