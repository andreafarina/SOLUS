function mask=roispline(I,kind,tension)

% function [mask,perimeter,area]=roispline(I,str,tension)
%
% INPUT :
% - I       : Grayscale or color image
% - kind    : String specifying the kind of spline ('natural' for natural 
%             cubic spline, 'cardinal' for cardinal cubic spline). 
%             DEFAULT = 'natural'
% - tension : Tension parameter between 0 and 1 for cardinal spline.
%             DEFAULT = 0
%
% OUTPUT :
% - mask    : Binary mask of segmented ROI


warning off

if nargin<2
    kind='natural';
elseif nargin<3
    tension = 0;
end

siz=size(I);
imshow(I, [0, 255]);
title('Click right button or close curve clicking near first point');
hold on;


% initializes number of points for each step
npoints=sqrt(siz(1)*siz(2))*.1;

dim=10/835*sqrt(siz(1)*siz(2));
if dim<8
    dim=8;
end

bottone=1;
i=0;
flag=1;
points=[];
switch kind
    case 'natural'
        while flag                   
            i=i+1;
            [cor(1,i),cor(2,i),bottone]=ginput(1); %[x, y, button of mouse (1,2,3)]
            if bottone ~= 1 % if bigger than 1 (i.e. left) then the last value should be equal to the first one
                if i>1 % if bigger than 1 thenan arry of points is present and the last value should be equal to the first one
                    cor(:,i)=cor(:,1);
                    flag=0;
                else
                    break;
                end
            end
            plot(cor(1,i),cor(2,i),'y.'); %draw the selected point with a yellow dot
    
            if i>1      
                if (norm(cor(:,1)-cor(:,i))<dim & i>2) & flag % if you click close enough to the initial point then you identify that as your last point
                    cor(:,i)=cor(:,1);
                    flag=0;
                end
                imshow(I,[0, 255])
                hold on;
                %spcv = cscvn(cor);              
                %points=myfnplt(spcv,npoints*i);
                points = Lee_spline(cor, npoints*i);
                plot(points(1,:),points(2,:),'b','LineWidth',2);
                plot(cor(1,:),cor(2,:),'y.');
                drawnow;
            end
        end
    case 'cardinal'                  
        while flag                         
            i=i+1;
            [x(i),y(i),bottone]=ginput(1);
            if bottone ~= 1
                if i>1
                    x(i)=[];
                    y(i)=[];
                    Px=[x(end) x x(1) x(2)];
                    Py=[y(end) y y(1) y(2)];
                    i=i-1;
                    flag=0;
                else
                    break;
                end
            end            
            plot(x(i),y(i),'y.');  
            
            if i>1 
                if flag
                    if (norm([x(1) y(1)]-[x(i) y(i)])<dim & i>2)
                        x(i)=[];
                        y(i)=[];
                        Px=[x(end) x x(1) x(2)];
                        Py=[y(end) y  y(1) y(2)];
                        flag=0;
                    else
                        Px=[x(1) x x(end)];
                        Py=[y(1) y y(end)];
                    end
                end
                pointsx=[];
                pointsy=[];
                for k=1:length(Px)-3
                    [xvec,yvec]=EvaluateCardinal2DAtNplusOneValues([Px(k),Py(k)],[Px(k+1),Py(k+1)],[Px(k+2),Py(k+2)],[Px(k+3),Py(k+3)],tension,npoints);
                    pointsx=[pointsx, xvec];
                    pointsy=[pointsy, yvec];    
                end
                imshow(I,[0, 255]);
                hold on;
                plot(pointsx,pointsy,'LineWidth',2);
                plot(x,y,'y.');
                drawnow;
            end
        end
        points=[pointsx; pointsy];
end

% mask calculation
points=points';
if i>2
    mask=roipoly(I,points(:,1),points(:,2));
else
    mask=logical(zeros(siz(1),siz(2)));
end



end

