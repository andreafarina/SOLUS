function P = RasterScan(x1,x2,y1,y2,Nx,Ny,z)

dx = (x2-x1)/Nx;
dy = (y2-y1)/Ny;
[X,Y] = ndgrid(x1:dx:(x2-dx),y1:dy:(y2-dy),2);
P = zeros(Nx*Ny,3);
P(:,1) = X(:);
P(:,2) = Y(:);
P(:,3) = z;
