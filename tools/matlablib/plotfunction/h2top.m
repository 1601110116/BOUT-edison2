function h2top(varargin)

h=findobj(gcf,'type','axes');
h=sort(h);
if nargin>1
    error('wrong input argument!')
elseif nargin==0
    varargin={1};
end

if cell2mat(varargin)>length(h)
    error('such handle not exist!')
else
    set(gcf,'currentaxes',h(cell2mat(varargin)));
end