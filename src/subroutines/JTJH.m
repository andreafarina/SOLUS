function y = JTJH(x,J,H,lambda)
z = J*x;
y = J' * z + lambda*H*x;

