function gamma = getGyroRatio(nuc)
%getGyroRatio provides the gyromagnetic ratio of different nuclei
%all values are from https://en.wikipedia.org/wiki/Gyromagnetic_ratio
%
% Syntax:
%       getGyroRatio(nuc)
%
% Inputs:
%       nuc - nucleus type [string]
%
% Outputs:
%       gamma - gyromagnetic ration [rad/s/T]
%
% Example:
%       getGyroRatio('1H')
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

%% all freq values are in [MHz/T] and /2pi
switch nuc
    case '1H'
        freq = 42.57747892;
    case '2H'
        freq = 6.536;
    case '3He'
        freq = -32.434;
    case '7Li'
        freq = 16.546;
    case '13C'
        freq = 10.705;
    case '14N'
        freq = 3.077;
    case '15N'
        freq = -4.316;
    case '17O'
        freq = -5.772;
    case '19F'
        freq = 40.053;
    case '23Na'
        freq = 11.262;
    case '27Al'
        freq = 11.103;
    case '29Si'
        freq = -8.465;    
    case '31P'
        freq = 17.235;
    case '57Fe'
        freq = 1.382;
    case '63Cu'
        freq = 11.319;
    case '67Zn'
        freq = 2.669;    
    case '129Xe'
        freq = -11.777;
    otherwise
        freq = NaN;
end

% transform to [rad/s/T]
gamma = 2*pi*freq*1e6;

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
