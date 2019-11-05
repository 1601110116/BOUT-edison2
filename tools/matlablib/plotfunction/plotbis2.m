function plotbis2(bis,bis_phase,f)
%for bispeclin

num=max(get(0,'children'));
if isempty(num)
    num=0;
end

if ~iscell(bis);
    eval(['figure(' num2str(num+1) ')'])
    subplot(2,1,1)
    plot(f,bis);

    myplot(2,1,1)
    ylabel2('bi^{2}(f_{1},f_{3}-f_{1})')

    subplot(2,1,2)
    plot(f,bis_phase/pi,'linestyle','none','marker','^','markersize',4)
    myplot(2,1,2)
    ylabel2('\Delta\theta(f_{1},f_{3}-f_{1})/\pi')
    xlabel2('f_{1}')
    set2('linewidth',2.5)
else 
    len=length(bis);
    for i=1:len
        eval(['figure(' num2str(num+i) ')'])
        subplot(2,1,1)
        plot(f{i},bis{i});

        myplot(2,1,1)
        ylabel2('bi^{2}(f_{1},f_{3}-f_{1})')

        subplot(2,1,2)
        plot(f{i},bis_phase{i}/pi,'linestyle','none','marker','^','markersize',4)
        myplot(2,1,2)
        ylabel2('\Delta\theta(f_{1},f_{3}-f_{1})/\pi')
        xlabel2('f_{1}')
        set2('linewidth',2.5)
    end
end