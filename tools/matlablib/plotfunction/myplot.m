function myplot(varargin)

if nargin==1&&isstr(varargin{1})
    varargin{1}=str2double(varargin{1});
end

if nargin==3
    nrows=varargin{1};
    ncols=varargin{2};
    thisplot=varargin{3};
elseif nargin==1
    nrows=floor(varargin{1}/100);
    ncols=floor(rem(varargin{1},100)/10);
    thisplot=rem(rem(varargin{1},100),10);
else
    error('wrong input argument')
end

f=gca;
p=get(f,'position');
h=(0.9250-0.1100)/nrows;
p(2)=0.9250-(floor((thisplot-1)/ncols)+1)*h;
p(4)=h-0.03;
set(f,'position',p);
%y=get(f,'yticklabel');
%y(end,:)=char(32);
%set(f,'yticklabel',y);
if thisplot~=nrows*ncols
    set(f,'xticklabel',[]);
end