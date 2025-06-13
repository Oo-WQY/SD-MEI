function phi_ref = getReferencePhase(gamma)
%getReferencePhase provides the reference phase after Levitt depending on
%the sign of the gyromagnetic ratio
%
% Syntax:
%       getReferencePhase(gamma)
%
% Inputs:
%       gamma - gyromagnetic ratio [rad/s/T]
%
% Outputs:
%       phi_ref - reference phase [rad]
%
% Example:
%       getReferencePhase(getGyroRatio('1H'))
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

% reference phase depending on the sign of gamma
% after M. H. Levitt, Spin Dynamics - Basics of Nuclear Magnetic Resonance
% (John Wiley & Sons, LTD, 2002) page 244 eq. 10.17 
if gamma > 0
    phi_ref = pi;
elseif gamma <= 0
    phi_ref = 0;
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
