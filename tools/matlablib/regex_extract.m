function result=regex_extract(line,pattern)
% RESULT=REGEX_EXTRACT(LINE,PATTERN,NMATCH): get the data from string
% varable "line".
% SYNTAX
% [output] = RESULT
%
% INPUT PARAMETERS
% A: mesh in patch format with the following additional subfields
% - line: a string
% - pattern: for gfile and pfile, 
%           pattern = '[+-]?([0-9])+($| |(\.[0-9]+([eE][+-]?[0-9]*)?))'
%
% Examples:
% description of example for result=regex_extract(line,pattern,nmatch)
% >> result=regex_extract(line,pattern,nmatch);
%
% See also:
% read_1d, read_2d, next_double
%
% References:
% [1] IDL function: regex_extract.pro 
%               by Ben Dudson, University of York, Feb 2010

% Copyright (c) ASIPP 1958-2015 Defeng Kong
% All Rights Reserved.
% email: dfkong@ipp.ac.cn
% $Revision: 1.0$ Created on: 26-May-2015 14:28:19

% write down your codes from here
result=0;
nmatch=0;
str=line;
ind=regexp(str,pattern,'start');
%s=regexp(str,pattern,'match');
if length(ind)==0
    result='';
end
while(ind(1)<=2)
    if length(ind)>1
        s=regexp(str,pattern,'match');
        str=strcat(s{2:end});      
        if nmatch==0
            result=s{1};
        else
            result=[result,s{1}];
            nmatch=nmatch+1;
        end
        break
    elseif length(ind)==1
        if nmatch==0
            result=s{1};
        else
            result=[result,s{1}];
            nmatch=nmatch+1;
        end
        str='';
        break
    end
    ind=regexp(str,pattern,'start');
    if length(ind)==0
        break;
    end
end
