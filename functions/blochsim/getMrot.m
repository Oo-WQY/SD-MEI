function Mrot = getMrot(M,theta,varargin)
%getMrot transforms the magnetization from the laboratory frame of
%reference into the rotating frame of reference
%
% Syntax:
%       getMrot(M,theta,varargin)
%
% Inputs:
%       M - magnetization in the lab-frame
%       theta - instantaneous phase omega*t [rad]
%       phi - [optional] additional phase [rad]
%
% Outputs:
%       Mrot - magnetization in the rot-frame
%
% Example:
%       getMrot(M,theta)
%
% Other m-files required:
%       none
%
% Subfunctions:
%       none
%
% MAT-files required:
%       none
%
% See also: BLOCHUS
% Author: Thomas Hiller
% email: thomas.hiller[at]leibniz-liag.de
% License: GNU GPLv3 (at end)

%------------- BEGIN CODE --------------

phi = 0;
if nargin > 2
    phi = varargin{1};
end
Mrot = zeros(size(M));

Mrot(:,1) =  M(:,1).*cos(theta + phi) + M(:,2).*sin(theta + phi);
Mrot(:,2) = -M(:,1).*sin(theta + phi) + M(:,2).*cos(theta + phi);
Mrot(:,3) =  M(:,3);

return

%------------- END OF CODE --------------

%% License:
% GNU GPLv3
%
% BLOCHUS
% Copyright (C) 2019 Thomas Hiller
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
