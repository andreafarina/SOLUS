function [nvar,str_var] = ExtractVariables_spectral(x,spe)
nvar = numel(x);
str_var = '';
for i = 1:nvar
    dummy = x{i};
    if strcmpi(dummy,'mus')
        dummy = 'D';
    end
    str_var = strcat(str_var,dummy);
end
if nvar == 2
    nvar = spe.nCromo+2;
elseif strcmpi(str_var,'mua')
    nvar = spe.nCromo;
elseif strcmpi(str_var,'D')
    nvar = 2;
end
