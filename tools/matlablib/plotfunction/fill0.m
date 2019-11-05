function fill0(xrange,c)

if nargin==1
    c='r';
end

h=findobj(gcf,'type','axes','-not','tag','legend');
len=length(h);

for i=1:len
        yrange=get(h(i),'ylim');
        set(gcf,'currentaxes',h(i));
        X=[xrange(1) xrange(1) xrange(2) xrange(2)];
        Y=[yrange(1) yrange(2) yrange(2) yrange(1)];
        fill2(X,Y,c);
end