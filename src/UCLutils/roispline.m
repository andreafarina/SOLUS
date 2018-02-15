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
imshow(I);
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
            [cor(1,i),cor(2,i),bottone]=ginput(1);
            if bottone ~= 1
                if i>1
                    cor(:,i)=cor(:,1);
                    flag=0;
                else
                    break;
                end
            end
            plot(cor(1,i),cor(2,i),'y.');
    
            if i>1      
                if (norm(cor(:,1)-cor(:,i))<dim & i>2) & flag
                    cor(:,i)=cor(:,1);
                    flag=0;
                end
                imshow(I)
                hold on;
                spcv = cscvn(cor);              
                points=myfnplt(spcv,npoints*i);
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
                imshow(I);
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




% internal functions

function [points,t] = myfnplt(f,npoints,varargin)

% interpret the input:
symbol=''; interv=[]; linewidth=[]; jumps=0;
for j=3:nargin
   arg = varargin{j-1};
   if ~isempty(arg)
      if ischar(arg)
         if arg(1)=='j', jumps = 1;
         else, symbol = arg;
         end
      else
         [ignore,d] = size(arg);
         if ignore~=1
	    error('SPLINES:FNPLT:wrongarg',['arg',num2str(j),' is incorrect.']), end
         if d==1
            linewidth = arg;
         else
            interv = arg;
         end
      end
   end
end

% generate the plotting info:
if ~isstruct(f), f = fn2fm(f); end

% convert ND-valued to equivalent vector-valued:
d = fnbrk(f,'dz'); if length(d)>1, f = fnchg(f,'dim',prod(d)); end

switch f.form(1:2)
case 'st'
   if ~isempty(interv), f = stbrk(f,interv);
   else
      interv = stbrk(f,'interv');
   end
%    npoints = 150; 
   d = stbrk(f,'dim');
   switch fnbrk(f,'var')
   case 1
      x = linspace(interv{1}(1),interv{1}(2),npoints);
      v = stval(f,x);
   case 2
      x = {linspace(interv{1}(1),interv{1}(2),npoints), ...
                    linspace(interv{2}(1),interv{2}(2),npoints)};
      [xx,yy] = ndgrid(x{1},x{2});
      v = reshape(stval(f,[xx(:),yy(:)].'),[d,size(xx)]);
   otherwise
      error('SPLINES:FNPLT:atmostbivar', ...
            'Cannot handle st functions with more than 2 variables.')
   end
otherwise
   if ~strcmp(f.form([1 2]),'pp')
      givenform = f.form; f = fn2fm(f,'pp'); basicint = ppbrk(f,'interval');
   end
   
   if ~isempty(interv), f = ppbrk(f,interv); end
      
   [breaks,coefs,l,k,d] = ppbrk(f);
   if iscell(breaks)
      m = length(breaks);
      for i=m:-1:3
         x{i} = (breaks{i}(1)+breaks{i}(end))/2;
      end
%       npoints = 150;
      ii = [1]; if m>1, ii = [2 1]; end
      for i=ii
         x{i}= linspace(breaks{i}(1),breaks{i}(end),npoints);
      end
      v = ppual(f,x);
      if exist('basicint','var') 
                         % we converted from B-form to ppform, hence must now
                         % enforce the basic interval for the underlying spline.
         for i=ii
            temp = find(x{i}<basicint{i}(1)|x{i}>basicint{i}(2));
            if d==1
               if ~isempty(temp), v(:,temp,:) = 0; end
               v = permute(v,[2,1]);
            else
               if ~isempty(temp), v(:,:,temp,:) = 0; end
               v = permute(v,[1,3,2]);
            end
         end
      end
   else     % we are dealing with a univariate spline
%       npoints = 500;
      x = [breaks(2:l) linspace(breaks(1),breaks(l+1),npoints)];
      v = ppual(f,x); 
      if l>1 % make sure of proper treatment at jumps if so required
         if jumps
            tx = breaks(2:l); temp = repmat(NaN, d,l-1);
         else
            tx = []; temp = zeros(d,0);
         end
         x = [breaks(2:l) tx x]; 
         v = [ppual(f,breaks(2:l),'left') temp v];
      end
      [x,inx] = sort(x); v = v(:,inx);
   
      if exist('basicint','var') 
                         % we converted from B-form to ppform, hence must now
                         % enforce the basic interval for the underlying spline.
                         % Note that only the first d components are set to zero
                         % outside the basic interval, i.e., the (d+1)st 
                         % component of a rational spline is left unaltered :-)
         if jumps, extrap = repmat(NaN,d,1); else, extrap = zeros(d,1); end
         temp = find(x<basicint(1)); ltp = length(temp);
         if ltp
            x = [x(temp),basicint([1 1]), x(ltp+1:end)];
            v = [zeros(d,ltp+1),extrap,v(:,ltp+1:end)];
         end
         temp = find(x>basicint(2)); ltp = length(temp);
         if ltp
            x = [x(1:temp(1)-1),basicint([2 2]),x(temp)];
            v = [v(:,1:temp(1)-1),extrap,zeros(d,ltp+1)];
         end
         %   temp = find(x<basicint(1)|x>basicint(2));
         %   if ~isempty(temp), v(temp) = zeros(d,length(temp)); end
      end
   end
   
   if exist('givenform','var')&&givenform(1)=='r' 
                                           % we are dealing with a rational fn:
                                           % need to divide by last component
      d = d-1;
      sizev = size(v); sizev(1) = d;
      % since fnval will replace any zero value of the denominator by 1,
      % so must we here, for consistency:
      v(d+1,find(v(d+1,:)==0)) = 1;
      v = reshape(v(1:d,:)./repmat(v(d+1,:),d,1),sizev);
   end 
end

%  use the plotting info, to plot or else to output:
if nargout==0
   if iscell(x)
      switch d
      case 1
         [yy,xx] = meshgrid(x{2},x{1});
         surf(xx,yy,reshape(v,length(x{1}),length(x{2})))
      case 2
         v = squeeze(v); roughp = 1+(npoints-1)/5;
         vv = reshape(cat(1,...
              permute(v(:,1:5:npoints,:),[3,2,1]),...
              repmat(NaN,[1,roughp,2]),...
              permute(v(:,:,1:5:npoints),[2,3,1]),...
              repmat(NaN,[1,roughp,2])), ...
                    [2*roughp*(npoints+1),2]);
         plot(vv(:,1),vv(:,2))
      case 3
         v = permute(reshape(v,[3,length(x{1}),length(x{2})]),[2 3 1]);
         surf(v(:,:,1),v(:,:,2),v(:,:,3))
      otherwise
      end
   else
      if isempty(symbol), symbol = '-'; end
      if isempty(linewidth), linewidth = 2; end
      switch d
      case 1, plot(x,v,symbol,'linew',linewidth)
      case 2, plot(v(1,:),v(2,:),symbol,'linew',linewidth)
      otherwise
         plot3(v(1,:),v(2,:),v(3,:),symbol,'linew',linewidth)
      end
   end
else
   if iscell(x)
      switch d
      case 1
         [yy,xx] = meshgrid(x{2},x{1});
         points = {xx,yy,reshape(v,length(x{1}),length(x{2}))};
      case 2
         [yy,xx] = meshgrid(x{2},x{1});
         points = {xx,yy,reshape(v,[2,length(x{1}),length(x{2})])};
      case 3
         points = {squeeze(v(1,:)),squeeze(v(2,:)),squeeze(v(3,:))};
         t = {x{1:2}}
      otherwise
      end
   else
      if d==1, points = [x;v];
      else, t = x; points = v([1:min([d,3])],:); end
   end
end


function [xvec,yvec]=EvaluateCardinal2DAtNplusOneValues(P0,P1,P2,P3,T,N)
% Evaluate cardinal spline at N+1 values for given four points and tesion.
% Uniform parameterization is used.

% P0,P1,P2 and P3 are given four points.
% T is tension.
% N is number of intervals (spline is evaluted at N+1 values).

xvec=[]; yvec=[];

% u vareis b/w 0 and 1.
% at u=0 cardinal spline reduces to P1.
% at u=1 cardinal spline reduces to P2.

u=0;
[xvec(1),yvec(1)] =EvaluateCardinal2D(P0,P1,P2,P3,T,u);
du=1/N;
for k=1:N
    u=k*du;
    [xvec(k+1),yvec(k+1)] =EvaluateCardinal2D(P0,P1,P2,P3,T,u);
end


function [xt,yt] =EvaluateCardinal2D(P0,P1,P2,P3,T,u)
% Evaluates 2D Cardinal Spline at parameter value u

% INPUT
% P0,P1,P2,P3  are given four points. Each have x and y values.
% P1 and P2 are endpoints of curve.
% P0 and P3 are used to calculate the slope of the endpoints (i.e slope of P1 and
% P2).
% T is tension (T=0 for Catmull-Rom type)
% u is parameter at which spline is evaluated

% OUTPUT
% cardinal spline evaluated values xt,yt,zt at parameter value u

s= (1-T)./2;

% MC is cardinal matrix
MC=[-s     2-s   s-2        s;
    2.*s   s-3   3-(2.*s)   -s;
    -s     0     s          0;
    0      1     0          0];

GHx=[P0(1);   P1(1);   P2(1);   P3(1)];
GHy=[P0(2);   P1(2);   P2(2);   P3(2)];

U=[u.^3    u.^2    u    1];


xt=U*MC*GHx;
yt=U*MC*GHy;



