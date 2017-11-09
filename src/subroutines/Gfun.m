function y = Gfun(p,space,size_space)
mu = [p(3),p(4),p(5)];
    cov_matrix = [p(6),p(7),p(8);
    p(7),p(9),p(10);
    p(8),p(10),p(11)];
%% with automatic implicit expansion, compatible with R0216a and later
% G = @(A,b,mu,cov_matrix)(b+(exp(-1/2*sum((space-mu)/(cov_matrix).*(space-mu),2))*...
%     1./(sqrt((2*pi)^3*det(cov_matrix)))*A));
%% compatible with earlier
space_off = space - repmat(mu,[size(space,1),1]);
G = @(A,b,mu,cov_matrix)(b+(exp(-1/2*sum((space_off)/(cov_matrix).*(space_off),2))*...
    1./(sqrt((2*pi)^3*det(cov_matrix)))*A));
%%
%p(2) = 0;
y = G(p(1),p(2),mu,cov_matrix);
y = reshape(y,size_space);
end

