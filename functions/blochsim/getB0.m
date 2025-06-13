function B = getB0(gamma,omega0)
%getB0 calculates the magnetic field from a given angular frequency and
%gyromagnetic ratio
%
% Syntax:
%       getB0(gamma,omega0)
%
% Inputs:
%       gamma - gyromagnetic ration [rad/s/T]
%       omega0 - angular frequency [rad/s]
%
% Outputs:
%       B - magnetic field [T]
%
% Example:
%       getB0(getGyroRatio('1H'),2000*2*pi)
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

% even though the original formula is sign-sensitive, as a precaution
% I make the output B-field always positive
B = abs(omega0/-gamma);

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
