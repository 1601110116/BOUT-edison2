function arr4d = arrto4d(arr,varargin)
% ARRTO4D change an array to a 4D array, e.g. Ti0(m,n) to Ti0(m,n,k,l) 
%
% INPUT PARAMETERS
% arr: input array, the array should be an array no more than 4D;
% varargin: dimension contral, for 1d input, varargin should be 'n,k,l'; for
% 2d input,varargin should be 'k,l',for 3d input, varargin should be 'l' 
%
%
% DESCRIPTION
% Exhaustive and long description of functionality of arrto4d.
%
% Examples:
% description of example for arrto4d
% arr4d = ARRTO4D(num,m,n,k,l);
% arr4d = ARRTO4D(arr1d,n,k,l);
% arr4d = ARRTO4D(arr2d,k,l);
% arr4d = ARRTO4D(arr3d,l);
%
% See also:
% ones

% Copyright (c) ASIPP 1985-2015 Defeng Kong
% All Rights Reserved.
% email: dfkong@ipp.ac.cn
% $Revision: 1.0$ Created on: 04-Jun-2015 11:32:51

if length(size(arr))>4
    error('dimensions of arr should be no more than 4');
end

len=length(varargin);
% if len==1&&iscell(varargin{1})
%     len=length(varargin{1});
%     varargin=varargin{1};
% end

switch len
    case 1;
        arr4d=arr(:,:,:,ones(1,1,1,varargin{1}));
    case 2;
        arr4d=arr(:,:,ones(1,1,varargin{1},1),ones(1,1,1,varargin{2}));
    case 3;
        arr4d=arr(:,ones(1,varargin{1},1,1),ones(1,1,varargin{2},1),ones(1,1,1,varargin{3}));
    case 4
        arr4d=arr(ones(varargin{1},1,1,1),ones(1,varargin{2},1,1),ones(1,1,varargin{3},1),ones(1,1,1,varargin{4}));
end
    

