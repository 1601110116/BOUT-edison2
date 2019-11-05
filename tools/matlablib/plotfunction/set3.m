function set3(varargin)

h=findobj(gcf,'type','axes');

for i=1:length(h)
        set(h(i),varargin{:});
end