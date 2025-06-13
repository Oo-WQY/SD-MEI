function ax_phi = getPulseAxisPhase(ax_dir)
%getPulseAxisPhase provides the phase angle of the given pulse axis
%
% Syntax:
%       getPulseAxisPhase(ax_dir)
%
% Inputs:
%       ax_dir - string holding the pulse axis direction
%
% Outputs:
%       ax_phi - phase angle of the pulse axis [rad]
%
% Example:
%       getPulseAxisPhase('+x')
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

switch ax_dir
    case '+x'
        ax_phi = 0;
    case '+y'
        ax_phi = pi/2;
    case '-x'
        ax_phi = pi;
    case '-y'
        ax_phi = 3*pi/2;
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
