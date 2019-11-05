function plotmis(bis,bis_phase,f,r)
%for mislin
%r: cm

if nargin==3
    r=1;
end

len=length(bis);

num=max(get(0,'children'));
if isempty(num)
    num=0;
end
for i=1:len
    eval(['figure(' num2str(num+i) ')'])
    subplot(2,1,1)
    plot(f{i},bis{i});
    
    myplot(2,1,1)
    ylabel2('\Gamma^{2}(f_{1},f_{2})')
    
    subplot(2,1,2)
    plot(f{i},bis_phase{i}/r,'linestyle','none','marker','^','markersize',4)
    myplot(2,1,2)
    ylabel2('\Delta k(cm^{-1})')
    xlabel2('f_{1}(kHz)')
    xtick2
    set2('linewidth',2.5)
end