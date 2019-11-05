function plotmir(h)

if nargin==0
    h=gca;
end

h_plot=findobj(h,'type','line');

for i=1:length(h_plot)
    temp=h_plot(i);
    xdata=get(temp,'xdata');
    ydata=get(temp,'ydata');
    set(h,'nextplot','add');
    plot(-xdata,ydata,'color','r','linestyle','--');
end