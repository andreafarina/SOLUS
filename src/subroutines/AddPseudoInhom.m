function [mat, noisy] = AddPseudoInhom( hom , dim,  corr, perc)

% function [mat, noisy] = AddPseudoInhom(hom, dim,  corr, perc)
% given as input:
% - hom: an input matrix of values or a single value (if so the hom will be interpreted
%           as the matrix hom = hom * ones(dim) )
% - size of the input matrix
% - corr : Characteristic correlation distance of the noise expressed in number of pixels 
% - perc : Maximum variation expressed in percentage
% Returns:
% - mat: output matrix of dimensions dim, obtained adding correlated noise, defined by corr and perc, to hom
% - noisy: output matrix of dimensions dim displaying the percentage
%       variation of noise applied to hom so that, mat = noisy * ones(dim)


%%
if numel(hom) == 1
    hom = ones(dim);
end
x = (-0.5*(dim(1)-1)):1:(0.5*(dim(1)-1));
y = (-0.5*(dim(2)-1)):1:(0.5*(dim(2)-1));
z = (-0.5*(dim(3)-1)):1:(0.5*(dim(3)-1));


lcar = corr;
lcar = 2 * sqrt((3/(8*pi^3))*(lcar^2));
sigsq = (lcar*lcar);
factexp =1/(sqrt(2*pi*sigsq^3));

for i = 1:1:length(x)
    for j = 1:1:length(y)
        for k = 1:1:length(z)
            
        expin(i,j,k) = (x(i).^2) + (y(j).^2) + (z(k).^2);
        expin(i,j,k) = expin(i,j,k) / (2 * sigsq);
        
        end
    end
end



noise = rand(dim(1), dim(2), dim(3));
FFnoise = fftn(noise);
smoother(:,:,:) =  fftn((factexp) *exp( - expin));
FFnoise = ifftn(smoother .* FFnoise);
FFnoise = abs(FFnoise);
meanFFnoise = mean(FFnoise(:));


FFnoise = FFnoise - meanFFnoise;

maxFFnoise = max(FFnoise(:));
minFFnoise = min(FFnoise(:));
perc = perc/( maxFFnoise - minFFnoise);

mat = hom + perc * hom .* (FFnoise);
%mat = smooth3(mat, 'box',5);
noisy = ones(dim(1), dim(2), dim(3)) + perc *(FFnoise);
%noisy = smoother;
%noisy = smooth3(noisy, 'box', 5);
end

