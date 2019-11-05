function fill2(X,Y,C)
hold on;
axis_h=gca;

fill(X,Y,C)
temp=get(axis_h,'children');
new=[temp(2:end) ;temp(1)];
set(temp(1),'edgecolor','white')
set(axis_h,'children',new);
set(axis_h,'layer','top')
