function B = GenerateMontageMatrix(A,res_fact)
%% Generate a block matrix of the volume A 
% based the resize vector in the 3 directions
% A. Farina - CNR - 02/11/2020

if nargin<2
    res_fact = 1;
end

[Nx,Ny,Nz] = size(A);
Ar = imresize3(A,res_fact.*[Nx,Ny,Nz]);
%% resize
[Nx,Ny,Nz] = size(Ar);

Nrow = ceil(sqrt(Nz));
Ncol = round(sqrt(Nz));

il = 1;
delta = 2;

B = ones(Nrow * Ny + (Nrow-1)*delta,Ncol * Nx + (Ncol-1)*delta);
B = (max(A(:))+min(A(:)))/2*B;
%B = mean(A(:))*B;
for i = 1:Nrow
    for j = 1:Ncol
        offsetx = (j-1)*(Nx+delta);
        offsety = (i-1)*(Ny+delta);
        
        if il<=Nz
            B(offsety + (1:Ny),offsetx + (1:Nx)) = flipud(Ar(:,:,il)');
            il = il + 1;
            %figure(100),imagesc(B),pause(0.1)
        end
    end
end

%figure,imagesc(B,[min(A(:)),max(A(:))]),axis image colorbar;
        
        
    