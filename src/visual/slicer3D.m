function h =  slicer3D(im, REC, position,str_plot)

% h =  slicer3D(im, REC, position,str_plot)
% plots the image im in 3D, requires the toolbox MatImage 
% im : 3D image to be shown
% REC : REC structure as outputted from the reconstruction
% position : position where to centre the visualisation of im. If empty display will be centre in the middle of im
% str_plot: string governing the output.
%           if '3D' it outputs a 3D visualisation 
%           if 'slice' it outputs a 3D orthogonal view of im
%           if 'both' it outputs a combined plot



% load /scratch0/NOT_BACKED_UP/gdisciac/SOLUS_spectralfitUCL_bis/SOLUS/example/example_200202_series/US_cylprior_exp_t0.01/ex_bulk10_0.1_incl10_0.05/Test_Standard_REC.mat
h = figure();
% a= reshape(REC.opt.bmua(:,3), REC.grid.dim);
% close all
a = im;

x = REC.grid.x; 
y = REC.grid.y;
z = REC.grid.z;

if nargin < 3 || isempty(position)  
    position = round([numel(x)/2,numel(y)/2,numel(z)/2]);
end
if nargin<4
    str_plot = 'both';
end
%position = permute(position, [2,1,3]);
position(1) = numel(x) - position(1);

spacing = [REC.grid.dx,REC.grid.dy, REC.grid.dz];
origin = [REC.grid.x(1),REC.grid.y(1), REC.grid.z(1)];


if nargin < 4 || ( exist('str_plot','var') && strcmpi(str_plot, 'both')==1)
    sb1 = subplot(1,2,1);
    lims = [prctile(a(:),0.5), prctile(a(:),99.5)];
    to_plot = permute(flip(a), [2,1,3]);
    c = OrthoSlicer3d(to_plot,'spacing', spacing,...
        'origin',[0,0,0], 'origin', origin,'position', position,...
        'colormap', 'parula', 'displayRange', lims);
    title('Reconstruction', 'FontSize', 20)
    caxis(lims);
    xlabel('X (mm)');
    ylabel('Y (mm)');
    zlabel('Z (mm)');
    set(gca, 'Zdir', 'reverse');
    cbh = colorbar();
    pause(3)


    sb2 =subplot(1,2,2);
    %out = identify_inclusion(a);
    idxb = identify_inclusion(a);
    b = zeros(size(a));
    b(idxb) = 1;
    cmap = jet(256);
    p = patch(isosurface(REC.grid.x,REC.grid.y,REC.grid.z,permute(flip(b),[2,1,3]), 0.1, 'noshare'));
    pause(1)
    xlim([x(1), x(end)])
    ylim([y(1), y(end)])
    zlim([z(1), z(end)])
    drawnow
    pause(2)
    view(25,30)
    %colormap([0.3,0.3,0.3])
    my_cmap = [0.3,0.3,0.3];
    grid on
    set(p, 'FaceColor', my_cmap,'EdgeColor', 'none'); %'FaceColor', cmap(end,:)
    extent = stackExtent(size(b), spacing, origin);
    % setup display
    axis image;
    axis(extent);
    view(3);
    axis('vis3d');
    title('Region of Inclusion', 'FontSize', 20)
    light;
    camlight(30,-75)
    set(gca, 'zdir', 'reverse');
    xlabel('X (mm)');
    ylabel('Y (mm)');
    zlabel('Z (mm)');
    view(-25,30)
    % 
else
    if strcmpi(str_plot, '3D')
        
        idxb = identify_inclusion(a);
        if ~isempty(idxb)
            b = zeros(size(a));
            b(idxb) = 1;
        else
            b=a;
        end
        cmap = jet(256);
        p = patch(isosurface(REC.grid.x,REC.grid.y,REC.grid.z,permute(flip(b),[2,1,3]), mean(b(:)), 'noshare'));
        pause(1)
        xlim([x(1), x(end)])
        ylim([y(1), y(end)])
        zlim([z(1), z(end)])
        
        drawnow
        pause(2)
        view(25,30)
        colormap([0.3,0.3,0.3])
        my_cmap = [0.3,0.3,0.3];
        grid on
        set(p, 'FaceColor', my_cmap, 'EdgeColor', 'none');
        extent = stackExtent(size(b), spacing, origin);
        % setup display
        axis image;
        axis(extent);
        view(3);
        axis('vis3d');
        %title('Region of Inclusion', 'FontSize', 20)
        light;
        camlight(30,-75)
        set(gca, 'zdir', 'reverse');
        xlabel('X (mm)');
        ylabel('Y (mm)');
        zlabel('Z (mm)');
        view(-25,30)
        
    elseif strcmpi(str_plot, 'slice')    
            a = im;
        x = REC.grid.x; 
        y = REC.grid.y;
        z = REC.grid.z;
        
        if nargin < 3 || isempty(position)  
            position = round([numel(x)/2,numel(y)/2,numel(z)/2]);
        end
        
        %position = permute(position, [2,1,3]);
        position(1) = numel(x) - position(1);
        spacing = [REC.grid.dx,REC.grid.dy, REC.grid.dz];
        origin = [REC.grid.x(1),REC.grid.y(1), REC.grid.z(1)];
        %sb1 = subplot(1,2,1);
        lims = [prctile(a(:),0.5), prctile(a(:),99.5)];
        to_plot = permute(flip(a), [2,1,3]);
        c = OrthoSlicer3d(to_plot,'spacing', spacing,...
            'origin',[0,0,0], 'origin', origin,'position', position,...
            'colormap', 'parula', 'displayRange', lims);
        %title('Reconstruction', 'FontSize', 20)
        if diff(lims) ==0
            lims(1) = 0.99*lims(1);
            lims(2) = 1.01*lims(2);
        end
        caxis(lims);
        xlabel('X (mm)');
        ylabel('Y (mm)');
        zlabel('Z (mm)');
        set(gca, 'Zdir', 'reverse');
        cbh = colorbar();

    end
 %%%%%   %set(gcf, 'Position', get(0, 'Screensize'));

    % pos1 = get(sb1, 'Position') % gives the position of current sub-plot
    % new_pos1 = pos1 +[0 0 +0.05 +0.05];
    % set(sb1, 'Position', new_pos1)
    % %cbh.Position = [new_pos1(1) + , ];

%     pos2 = get(h, 'Position') % gives the position of current sub-plot
%     new_pos2 = pos2 +[0.025 0.025 -0.025 -0.025];
%     set(h, 'Position', new_pos2)

    
end
