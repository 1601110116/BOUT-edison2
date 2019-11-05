function d = read_1d(n,fid)
% READ_1D read 1-D array data of gfile or pfile
% SYNTAX
% d = READ_1D(n)
%
% INPUT PARAMETERS
% d: mesh in patch format with the following additional subfields
% - n: length of the a
% - A.FN:
%
% OUTPUT PARAMETERS
% output: output description
%
% DESCRIPTION
% Exhaustive and long description of functionality of read_1d.
%
% Examples:
% description of example for read_1d
% >> [output] = read_1d(A);
%
% See also:
% ALSO1, ALSO2
%
% References:
% [1] A. Einstein, Die Grundlage der allgemeinen Relativitaetstheorie,
%     Annalen der Physik 49, 769-822 (1916).

% Copyright (c) ASIPP 1958-2015 Defeng Kong
% All Rights Reserved.
% email: dfkong@ipp.ac.cn
% $Revision: 1.0$ Created on: 26-May-2015 17:26:21

% write down your codes from here.
% global data
% global nleft
d=zeros(n,1);
data=0;
nleft=0;
pattern = '[+-]?([0-9])+($| |(\.[0-9]+([eE][+-]?[0-9]*)?))';
for ii=1:n
    if nleft==0
        line=fgetl(fid);
        data=regex_extract(line,pattern);
        nleft=length(data);
    end
    d(ii)=str2num(data(1));
    
    nleft=nleft-1;
    if nleft>0
        data=data(2:end);
    end
end