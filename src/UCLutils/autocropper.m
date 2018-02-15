function cropped = autocropper(im_in)
sens = 0.1;
tresh = 1;
im(:,:) = sum(im_in(:,:,1) + im_in(:,:,2) + im_in(:,:,3),3); 

[siz_x, siz_y] = size(im);

x_max = 1;
y_max = 1;
x_min = siz_x;
y_min = siz_y;

j = 1;
while( j <= ceil((1 - sens)*siz_y))
    i = 1;
    while( i <= ceil((1 - sens)*siz_x)) 
  
        
        tocheck = im( i: (i - 1 + ceil((sens * siz_x))), j: (j - 1 + ceil((sens * siz_y))));
        tocheck = (tocheck ~= 0);
        par = mean(mean(tocheck));
        
        
        if( par >= tresh)
            
            if( i < x_min)
                x_min = i;
            end
            if( (i + ceil((sens * siz_x))) > x_max ) 
                x_max = i + ceil((sens * siz_x));
            end
            if( j < y_min)
                y_min = j;
            end
            if( (j + ceil((sens * siz_y))) > y_max) 
               y_max = j + ceil((sens * siz_y));
            end
            
                       
        end
        
        i = i + 1; 
    end
    j = j + 1;
end


% % if ( x_max > x_min)
% %     if(y_max > y_min)
        cropped = im_in(x_min:x_max, y_min:y_max, :);
%     else
%         cropped = im_in(x_min:x_max, y_max:y_min, :);
%     end
% else
%     if(y_max > y_min)
%         cropped = im_in(x_max:x_min, y_min:y_max, :);
%     else
%         cropped = im_in(x_max:x_min, y_max:y_min, :);
%     end
    
end