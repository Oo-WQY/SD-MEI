function omega0 = getOmega0(gamma,B)
%getOmega0 calculates the angular frequency from a given B-field and
%gyromagnetic ratio
%
% Syntax:
%       getOmega0(gamma,B)
%
% Inputs:
%       gamma - gyromagnetic ration [rad/s/T]
%       B - magnetic field [T]
%
% Outputs:
%       omega0 - angular frequency [rad/s]
%
% Example:
%       getOmega0(getGyroRatio('1H'),50e-6)
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

% this formula is sign sensitive:
% omega0 is negative for positive gyromagnetic ratio
omega0 = -gamma.*B;

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
