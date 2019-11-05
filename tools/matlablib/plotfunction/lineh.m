function lineh(ypoint,c)
%a horizontal line

if nargin==0
    ypoint=0;
end
if nargin<=1
    c='k';
end

h=findobj(gcf,'type','axes');
len=length(h);

for i=1:len
    if isempty(get(h(i),'tag'))
        xrange=get(h(i),'xlim');
        set(gcf,'currentaxes',h(i));
        for j=1:length(ypoint)
            Y=[ypoint(j) ypoint(j)];
            X=[xrange(1) xrange(2)];
            line(X,Y,'linestyle','--','linewidth',2,'color',c);
        end
        temp=get(gca,'children');
        new=[temp(2:end) ;temp(1)];
%         set(temp(1),'edgecolor','white')
        set(gca,'children',new);
        set(gca,'layer','top')
    end
end