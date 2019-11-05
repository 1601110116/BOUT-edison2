function [M]=data2movie(var,scale,delay,filename)
% DATA2MOVIE show the plot of var and shore it into the filename,if there
% were no filename in input, this function is just like showdata coded by Minwoo Kim
% SYNTAX
% [output] = DATA2MOVIE(var,scale,delay,filename)
%
% INPUT PARAMETERS
% A: mesh in patch format with the following additional subfields
% - A.FA:
% - A.FN:
%
% OUTPUT PARAMETERS
% output: output description
%
% DESCRIPTION
% Exhaustive and long description of functionality of rms2movie.
%
% Examples:
% description of example for rms2movie
% >> [output] = data2movie(A);
%
% See also:
% showdata, getframe,movie,movie2avi
%
% Copyright (c) ASIPP 1985-2015 Defeng Kong
% All Rights Reserved.
% email: dfkong@ipp.ac.cn
% $Revision: 1.0$ Created on: 05-Jun-2015 17:04:21

% write down your codes from here.
delay_def = 0.05;

% Check input arguments
if ( nargin < 1 )
    fprintf('\tIt needs input arguments.\n');    
    return
end
if nargin<4
    filename='filename.avi';
end

if ( ischar(var) )
    fprintf('\tFirst input must be numeric format.\n');    
end
var = squeeze(var);
var_size = size(var);
if ( length(var_size) == 4 ) % 4 dimension
    fprintf('4-D sturcture visulization is not implemented yet.');
    return
elseif ( length(var_size) == 3 ) % 3 dimension
    option = 3;
elseif ( length(var_size) == 2 && sum(var_size) > 2 ) % 2 dimension
    option = 2;
else % Constant number
    fprintf('Input variable is just constant value.');
    return
end

if ( nargin < 2 )
    scale = 'off';
    delay =delay_def;
end
if ( nargin == 2 )
    if ( ischar(scale) )
        delay = delay_def;
    else
        delay = scale;
        scale = 'off';
    end
end
if ( ischar(delay) )
    delay = delay_def;
end

if ( strcmp(scale, 'on') )
    switch option
        case 2;
            axis_limit_max = max(max(var));
            axis_limit_min = min(min(var));
        case 3;
            axis_limit_max = max(max(max(var)));
            axis_limit_min = min(min(min(var)));
    end
    option = -option;
end

figure;

for i = 1:var_size(end)
   
    switch option
        
        case 2;
            ax1=subplot(2,1,1);
            set(ax1,'Position',[0.1300 0.6838 0.7750 0.2412]);
            plot(log(squeeze(var(324,:))));
            hold on;
            ylim0=get(gca,'Ylim');
            ylabel('log(P)','FontSize',20)
            plot([i,i],ylim0,'r--');
            set(gca, 'FontSize', 16);
            hold off;
            
            ax2=subplot(2,1,2);
            set(ax2,'Position',[0.1300 0.1300 0.7750 0.4212]);
            
            plot( var( :, i ), 'LineWidth', 2.5);
            %plot( var( :, i ), 'LineWidth', 2.5, 'Marker', 'x', 'MarkerSize', 10 );
            set(gca, 'FontSize', 16);
            text_time = ['t = ' num2str(i)];
            xlabel(text_time, 'FontSize', 20);
            M(i)=getframe(gcf);
            pause(delay);
         
        case -2;
            ax1=subplot(2,1,1);
            set(ax1,'Position',[0.1300 0.6838 0.7750 0.2412]);
            plot(log(squeeze(var(324,:))));
            hold on;
            ylim0=get(gca,'Ylim');
            ylabel('log(P)','FontSize',20)
            plot([i,i],ylim0,'r--');
            set(gca, 'FontSize', 16);
            set(gca,'XMinorTick','on','YMinorTick','on')
            hold off;
            
            ax2=subplot(2,1,2);
            set(ax2,'Position',[0.1300 0.1100 0.7750 0.4412]);
            
            plot( var( :, i ), 'LineWidth', 2.5);
            %plot( var( :, i ), 'LineWidth', 2.5, 'Marker', 'x', 'MarkerSize', 10 );
            set(gca, 'FontSize', 16);
            text_time = ['t = ' num2str(i)];
            xlabel(text_time, 'FontSize', 20);
            ylim( [axis_limit_min axis_limit_max] );
            set(gca,'XMinorTick','on','YMinorTick','on')
            M(i)=getframe(gcf);
            pause(delay);
            
        case 3;
            ax1=subplot(2,1,1);
            set(ax1,'Position',[0.1300 0.6838 0.7750 0.2412]);
            plot(log(squeeze(var(174,34,:))));
            hold on;
            ylim0=get(gca,'Ylim');
            ylabel('log(P)','FontSize',20)
            plot([i,i],ylim0,'r--');
            set(gca,'XMinorTick','on','YMinorTick','on')
            hold off;
            
            ax2=subplot(2,1,2);
            set(ax2,'Position',[0.1300 0.1100 0.7750 0.4412]);
            surf(var( :, :, i) );shading interp;
            ylim([0,var_size(1)]);
            xlim([0,var_size(2)]);
            set(gca, 'FontSize', 18);
            text_time = ['t = ' num2str(i)];
            xlabel(text_time, 'FontSize', 25);
            ylabel('r position', 'FontSize', 18);
            set(gca,'XMinorTick','on','YMinorTick','on')
            M(i)=getframe(gcf); %get the frame of current figure;
            %zlim( [axis_limit_min axis_limit_max] );
            pause(delay);
            
        case -3;
            ax1=subplot(2,1,1);
            set(ax1,'Position',[0.1300 0.6838 0.7750 0.2412]);
            plot(log(squeeze(var(174,34,:))));
            hold on;
            ylim0=get(gca,'Ylim');
            ylabel('log(P)','FontSize',20)
            plot([i,i],ylim0,'r--');
            set(gca,'XMinorTick','on','YMinorTick','on')
            hold off;
            
            ax2=subplot(2,1,2);
            set(ax2,'Position',[0.1300 0.1100 0.7750 0.4412]);
            surf(var( :, :, i) );shading interp;
            ylim([0,var_size(1)]);
            xlim([0,var_size(2)]);
            set(gca, 'FontSize', 18);
            text_time = ['t = ' num2str(i)];
            xlabel(text_time, 'FontSize', 25);
            ylabel('r position', 'FontSize', 18);
            set(gca,'XMinorTick','on','YMinorTick','on')
            M(i)=getframe(gcf); %get the frame of current figure;
            zlim( [axis_limit_min axis_limit_max] );
            pause(delay);
    
    end
end
if nargin==4
     movie2avi(M,filename);
end
end