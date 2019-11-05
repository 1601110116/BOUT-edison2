function set2(varargin)

h_axes=findobj(gcf,'type','axes');

if strcmp(varargin{1},'color')
     for i=1:length(h_axes)
        h=findobj(h_axes(i),'type','line');   %handles of plots under one axe
        for j=1:length(h)
            x=get(h(j),'xdata');
            y=get(h(j),'ydata');
            if sum(abs(diff(x)))~=0&&sum(abs(diff(y)))~=0
                set(h(j),varargin{:});
            end
        end
     end
else
    for i=1:length(h_axes)
        h=findobj(h_axes(i),'type','line');   %handles of plots under one axe
        for j=1:length(h)
            set(h(j),varargin{:});
        end
    end
end