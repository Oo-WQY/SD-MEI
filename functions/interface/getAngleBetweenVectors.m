function [theta,sgn] = getAngleBetweenVectors(x,y)
%getAngleBetweenVectors calculates the angle theta between two vectors 'x' and 'y'
%
% Syntax:
%       getAngleBetweenVectors(x,y)
%
% Inputs:
%       x - vector
%       y - vector
%
% Outputs:
%       theta - angle between x and y [rad]
%       sgn - sign of theta
%
% Example:
%       getAngleBetweenVectors([1 0 0],[0 0 1])
%
% Other m-files required:
%       none;
%
% Subfunctions:
%       none
%
% MAT-files required:
%       none
%
% See also BLOCHUS
% Author: Thomas Hiller
% email: thomas.hiller[at]leibniz-liag.de
% License: GNU GPLv3 (at end)

%------------- BEGIN CODE --------------

if numel(x)<=3 % vector treatment
    % if x is a vector make x and y column vectors
    x = x(:);
    y = y(:);
    % angle [rad]
    theta = acos(dot(x,y)./(norm(x).*norm(y)));    
    % sign
    sgn = sign(cross(x,y));
    sgn = sgn(3);
else % matrix treatment
    % angle [rad]
    normx = sqrt(x(:,1).^2+x(:,2).^2+x(:,3).^2);
    normy = sqrt(y(:,1).^2+y(:,2).^2+y(:,3).^2);
    theta = acos(dot(x,y,2)./(normx.*normy));    
    % sign
    sgn = sign(cross(x,y));
    sgn = sgn(:,3);
end

if ~isreal(theta)
    theta = real(theta);
end

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
