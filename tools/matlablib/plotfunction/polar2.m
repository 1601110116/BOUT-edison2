function hpol = polar2(thr,varargin)
% thr is a threshold value for plot: lower limiter- thr(1)
%                                              upper limiter- thr(2)

% i.e.    polar2(0.042, bis, bis_pha, '*r')
%          polar2([0.042 0.2], bis, bis_pha, 'Ok' )

%   Copyright 1984-2006 The MathWorks, Inc.
%   $Revision: 5.22.4.7 $  $Date: 2006/05/27 18:07:42 $

% Parse possible Axes input

% edit by LAD 2008/04/20

[cax,args,nargs] = axescheck(varargin{:});
error(nargchk(1,3,nargs));

figure

if nargs < 1 || nargs > 3
    %error('Requires 2 or 3 data arguments.')
    error('MATLAB:polar:InvalidInput', 'Requires 2 or 3 or 4 data arguments.')
elseif nargs == 2 
    theta = args{1};
    rho = args{2};
    if ischar(rho)
        line_style = rho;
        rho = theta;
        [mr,nr] = size(rho);
        if mr == 1
            theta = 1:nr;
        else
            th = (1:mr)';
            theta = th(:,ones(1,nr));
        end
    else
        line_style = 'auto';
    end
elseif nargs == 1
    theta = args{1};
    line_style = 'auto';
    rho = theta;
    [mr,nr] = size(rho);
    if mr == 1
        theta = 1:nr;
    else
        th = (1:mr)';
        theta = th(:,ones(1,nr));
    end
else % nargs == 3
    [theta,rho,line_style] = deal(args{1:3});
end
if ischar(theta) || ischar(rho)
    %error('Input arguments must be numeric.');
    error('MATLAB:polar:InvalidInputType', 'Input arguments must be numeric.');
end
if ~isequal(size(theta),size(rho))
    %error('THETA and RHO must be the same size.');
    error('MATLAB:polar:InvalidInput', 'THETA and RHO must be the same size.');
end

    index=find(rho<=thr(1));
    rho(index)=[];
    theta(index)=[];
if length(thr)==2
    index=find(rho>=thr(2));
    rho(index)=[];
    theta(index)=[];
end

% get hold state
cax = newplot(cax);

next = lower(get(cax,'NextPlot'));
hold_state = ishold(cax);

% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');
ls = get(cax,'gridlinestyle');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
fUnits  = get(cax, 'DefaultTextUnits');
set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
    'DefaultTextFontName',   get(cax, 'FontName'), ...
    'DefaultTextFontSize',   get(cax, 'FontSize'), ...
    'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
    'DefaultTextUnits','data')

% only do grids if hold is off
if ~hold_state

% make a radial grid
    hold(cax,'on');
    maxrho = max(abs(rho(:)));
    hhh=line([-maxrho -maxrho maxrho maxrho],[-maxrho maxrho maxrho -maxrho],'parent',cax);
    set(cax,'dataaspectratio',[1 1 1],'plotboxaspectratiomode','auto')
    v = [get(cax,'xlim') get(cax,'ylim')];
    ticks = sum(get(cax,'ytick')>=0);
    delete(hhh);
% check radial limits and ticks
    rmin = 0; rmax = v(4); rticks = max(ticks-1,2);
    if rticks > 5   % see if we can reduce the number
        if rem(rticks,2) == 0
            rticks = rticks/2;
        elseif rem(rticks,3) == 0
            rticks = rticks/3;
        end
    end

% define a circle
    th = 0:pi/50:2*pi;
    xunit = cos(th);
    yunit = sin(th);
% now really force points on x/y axes to lie on them exactly
    inds = 1:(length(th)-1)/4:length(th);
    xunit(inds(2:2:4)) = zeros(2,1);
    yunit(inds(1:2:5)) = zeros(3,1);
% plot background if necessary
    if ~ischar(get(cax,'color')),
       patch('xdata',xunit*rmax,'ydata',yunit*rmax, ...
             'edgecolor',tc,'facecolor',get(cax,'color'),...
             'handlevisibility','off','parent',cax);
    end

   err=thr(1);
% draw radial circles
    c82 = cos(60*pi/180);
    s82 = sin(60*pi/180);
     rinc = (rmax-rmin)/rticks;
    text((rmax+rinc/20)*c82,(rmax+rinc/20)*s82, ...
            ['max=' num2str(rmax)],'verticalalignment','bottom',...
            'interpreter','latex','handlevisibility','off','parent',cax)
    line(xunit*err,yunit*err,'linestyle',ls,'color',tc,'linewidth',1,...
                  'handlevisibility','off','parent',cax);
%     text(-err*4/5,-err*4/5, ...
%             'noise','verticalalignment','bottom',...
%             'handlevisibility','off','parent',cax)
%     

    rt = 1.1*rmax;
    th=((1:12)*30-15)*pi/180;
    cst = cos(th); snt = sin(th);
    nu=histc(phaseturn(theta),(0:30:360)*pi/180);
    num=100*nu(1:12)/length(theta);
    for i = 1:length(th)
        if num(i)>1
            text(rt*cst(i),rt*snt(i),[int2str(num(i)) '%'],...
                 'color','b','horizontalalignment','center',...
                 'handlevisibility','off','parent',cax);     
        elseif num(i)<1&&num(i)~=0
             text(rt*cst(i),rt*snt(i),'<1%',...
                 'color','b','horizontalalignment','center',...
                 'handlevisibility','off','parent',cax);     
        end
    end
    
%  plot spokes
    th = (1:6)*2*pi/12;
    cst = cos(th); snt = sin(th);
    cs = [-cst; cst];
    sn = [-snt; snt];
    line(rmax*cs,rmax*sn,'linestyle',ls,'color',tc,'linewidth',1,...
         'handlevisibility','off','parent',cax)

   
    
% % % plot spokes
%     th = (1:2)*2*pi/4;
%     cst = cos(th); snt = sin(th);
%     cs = [-cst; cst];
%     sn = [-snt; snt];
%     line(rmax*cs,rmax*sn,'linestyle',ls,'color',tc,'linewidth',1,...
%          'handlevisibility','off','parent',cax)
%  

% annotate spokes in degrees
         th = (1:2)*2*pi/4;
        cst = cos(th); snt = sin(th);
        text(rt*cst(1),rt*snt(1),'$$\pi/2$$',...
             'horizontalalignment','center',...
             'interpreter','latex','handlevisibility','off','parent',cax);
         text(-rt*cst(1),-rt*snt(1),'$$-\pi/2$$','horizontalalignment','center',...
             'interpreter','latex','handlevisibility','off','parent',cax);
          text(rt*cst(2),rt*snt(2),'$$\pm\pi$$','horizontalalignment','center',...
              'interpreter','latex','handlevisibility','off','parent',cax);
         text(-rt*cst(2),-rt*snt(2),'0','horizontalalignment','center',...
             'handlevisibility','off','parent',cax);


% set view to 2-D
    view(cax,2);
% set axis limits
    axis(cax,rmax*[-1 1 -1.15 1.15]);
end

% Reset defaults.
set(cax, 'DefaultTextFontAngle', fAngle , ...
    'DefaultTextFontName',   fName , ...
    'DefaultTextFontSize',   fSize, ...
    'DefaultTextFontWeight', fWeight, ...
    'DefaultTextUnits',fUnits );

% transform data to Cartesian coordinates.
xx = rho.*cos(theta);
yy = rho.*sin(theta);

% plot data on top of grid
if strcmp(line_style,'auto')
    q = plot(xx,yy,'parent',cax);
else
    q = plot(xx,yy,line_style,'parent',cax);
end

if nargout == 1
    hpol = q;
end

if ~hold_state
    set(cax,'dataaspectratio',[1 1 1]), axis(cax,'off'); set(cax,'NextPlot',next);
end
set(get(cax,'xlabel'),'visible','on')
set(get(cax,'ylabel'),'visible','on')

if ~isempty(q) && ~isdeployed
    makemcode('RegisterHandle',cax,'IgnoreHandle',q,'FunctionName','polar');
end
