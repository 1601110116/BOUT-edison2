function [var_name,varval] = readgridnc(filename,varargin)
% READGRIDNC read the BOUT++ input grid nv file
% SYNTAX
% [output] = READGRIDNC(filename,varargin)
%
% INPUT PARAMETERS
% -filename: the full path of the grid nc file
% -varargin: variables in the grid nv file 
%
% OUTPUT PARAMETERS
% output: output description
%
% DESCRIPTION
% Exhaustive and long description of functionality of readgridnc.
%
% Examples:
% description of example for readgridnc
% >> [output] = readgridnc('examples/xu2/data_t3020_zs5-p/cbm18_dens6_ne9_0.1_1.4.grid.nc','Ti0');
% >> [output] = readgridnc('examples/xu2/data_t3020_zs5-p/cbm18_dens6_ne9_0.1_1.4.grid.nc','Ti0','Te0');
%
% See also:
% info_file
%
% References:
% [1] info_file:Coded by Minwoo Kim(Mar. 2012) 

% Copyright (c) ASIPP 1985-2015 Defeng Kong
% All Rights Reserved.
% email: dfkong@ipp.ac.cn
% $Revision: 1.0$ Created on: 18-May-2015 15:32:56

% write down your codes from here.

j=length(varargin);
if ( nargin < 1 )
    fprintf('\tIt needs input arguments. Please input the full path of grid nc file.\n');
    return
elseif j==0
    if ( ~ischar(filename) )    
        fprintf('\tThe path of grid nc file must be string format.\n');    
    return
    end
elseif j>=1
    if ( ~ischar(filename) )    
        fprintf('\tThe filename must be string format.\n');    
    return
    end
    if ( ~ischar(varargin{1}) )    
        fprintf('\tThe varname must be string format.\n');    
    return
    end
end

% Check netCDF file existence
try
    netcdf.open(filename, 'nc_nowrite');
catch exception
    if ( strcmp(exception.identifier, 'MATLAB:netcdf:open:noSuchFile') )
        fprintf('\tCan not open the file. Please, check file name and path.\n');
        return
    end
end

% Open netCDF file
fid = netcdf.open(filename, 'nc_nowrite');
[~, num_var, ~, ~ ] = netcdf.inq(fid);
% [num_dim, num_var, num_global_atts, unlim_dim_ID ] = netcdf.inq(fid);

% fprintf( '\n' );
% fprintf( 'Variable list of netCDF file\n\n' );
% fprintf( '\t Var_ID\t\t Var_Name\t\tDimension\n' );

var_name=cell(num_var,1);
varval=cell(num_var,1);
for i=1 : num_var

    % Inquire saved variable information
    [var_name{i}, ~, ~, ~] = netcdf.inqVar(fid, i-1);
%     [var_name, xtype, dim_ids, num_atts] = netcdf.inqVar(fid, i-1);
%     if ( length(var_name) < 7 )
%         fprintf( '\t %d\t\t %s\t\t', i-1, var_name );        
%     else
%         fprintf( '\t %d\t\t %s\t', i-1, var_name );
%     end
        
    var_tmp = netcdf.getVar( fid, i-1 );
    if ( length(size(var_tmp)) == 4 ) % 4 dimension
        var_tmp = permute(var_tmp, [3 2 1 4]); % [Z, Y, X, T] -> [X, Y, Z, T]
    elseif ( length(size(var_tmp)) == 3 ) % 3 dimension
        var_tmp = permute(var_tmp, [3 2 1]); % [Z, Y, X] -> [X, Y, Z]
    elseif ( length(size(var_tmp)) == 2 && sum(size(var_tmp)) > 2 ) % 2 dimension
        var_tmp = permute(var_tmp, [2 1]); % [Y, X] -> [X, Y]   
    end
    varval{i}=var_tmp;
%     var_dim = size(var_tmp);
%     fprintf( '\t%d', var_dim(1) );
%     for j= 2:(length(var_dim))
%         fprintf( ' X %d', var_dim(j) );
%     end
%     fprintf( '\n' );  
end

netcdf.close(fid);

%output
if j==0&&nargout==0
    assignin('base','var_name',var_name);
    assignin('base','varval',varval);
elseif j>=1&&nargout==0
    for ii=1:num_var
        for jj=1:j
            if strcmp(varargin{jj},var_name{ii})
                assignin('base',regexprep(var_name{ii},'-','_'),varval{ii});
            end
        end
    end
end
