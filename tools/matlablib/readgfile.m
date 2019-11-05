function  readgfile(filename,varargin)
% READGFILE short description
% SYNTAX
% [output] = READGFILE(A)
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
% Exhaustive and long description of functionality of readgfile.
%
% Examples:
% description of example for readgfile
% >> [output] = readgfile(A);
%
% See also:
% ALSO1, ALSO2
%
% References:
% [1] A. Einstein, Die Grundlage der allgemeinen Relativitaetstheorie,
%     Annalen der Physik 49, 769-822 (1916).

% Copyright (c) ASIPP 1985-2015 dfkong
% All Rights Reserved.
% email: dfkong@ipp.ac.cn
% $Revision: 1.0$ Created on: 21-May-2015 03:01:05

% write down your codes from here.
j=length(varargin);
if ( nargin < 1 )
    fprintf('\tIt needs input arguments. Please input the full path of pfile.\n');
    return
elseif j==0
    if ( ~ischar(filename) )    
        fprintf('\tThe path of pfile must be string format.\n');    
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

%readdata
fid=fopen(filename,'r');
firstline=fgetl(fid);
expression=' ';
s=regexp(firstline,expression,'split');
ns=length(s);
if ns<3
    fprintf('\tError: Expecting at least 3 numbers on first line');
    return
end

idum=str2num(s{ns-2});
nxefit=str2num(s{ns-1});
nyefit=str2num(s{ns});
fprintf('\tidum=%d\tnxefit=%d\tnyefit=%d\n',idum,nxefit,nyefit);
fpol   = read_1d(nxefit,fid);
nihao=1;
