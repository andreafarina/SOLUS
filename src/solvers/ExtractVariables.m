function [nvar,str_var] = ExtractVariables(x)
nvar = numel(x);
str_var = '';
for i = 1:nvar
    dummy = x{i};
    if strcmpi(dummy,'mus')
        dummy = 'D';
    end
    str_var = strcat(str_var,dummy);
end
