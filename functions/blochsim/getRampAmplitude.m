function Bp = getRampAmplitude(t,param)
%getRampAmplitude provides pre-polarization switch-off B-field amplitude
%
% Syntax:
%       getRampAmplitude(t,param)
%
% Inputs:
%       t - time [s]
%       param - struct containing extra settings:
%               ramp  - struct containing ramp shape: 'exp', 'linexp',
%                       'halfcos', 'lin'
%               gamma  - gyromagnetic ratio [rad/s/T]
%               B0     - Earth magnetic field amplitude [T]
%               Bmax   - maximum pre-polarization amplitude [T]
%               Bstar  - switch magnetic field between linear and
%                        exponential part of the 'linexp' ramp
%               Tramp  - switch-off ramp time [s]
%               Tslope - switch time between linear and
%                        exponential part of the 'linexp' ramp [s]
%
% Outputs:
%       Bp - pre-polarization B-field amplitude
%
% Example:
%       getRampAmplitude(t,param)
%
% Other m-files required:
%       none
%
% Subfunctions:
%       getLinExpAmp
%
% MAT-files required:
%       none
%
% See also: BLOCHUS
% Author: Thomas Hiller
% email: thomas.hiller[at]leibniz-liag.de
% License: GNU GPLv3 (at end)

%------------- BEGIN CODE --------------

% get the different parameter
ramp = param.ramp;
gamma = param.gamma;
B0 = param.B0;
Bmax = param.Bmax;
Bstar = param.Bstar;
Tramp = param.Tramp;
Tslope = param.Tslope;

switch ramp
    case 'exp' % exponential
        Bp = Bmax .* exp(-t/Tslope);
    case 'linexp' % linear + exponential
        Bp = getLinExpAmp(Bmax,Bstar,Tslope,t);
    case 'halfcos' % half cosine
        Bp = Bmax .* (0.5+(cos(pi*t./Tramp)./2));
    case 'lin' % linear
        Bp = Bmax.*(1-t./Tramp);
end

return
% ---

% ---
function Bp = getLinExpAmp(Bmax,Bstar,T,t)
% linear + exponential ramp after:
% Conradi et al., 2017, Journal of Magnetic Resonance 281, p.241-245
% https://doi.org/10.1016/j.jmr.2017.06.001

if numel(t)>1
    % linear part
    Bplin = (-Bmax/T)*t + Bmax;
    % exponential part
    Bpexp = exp(-t /(Bstar*T/Bmax));
    % find change
    index = find(abs(Bplin-Bstar)==min(abs(Bplin-Bstar)),1,'first');
    % merge the lin- and exp-part and scale the amplitude of the exp-part
    % to that of the lin-part at the switch-over time t(index)
    scale_point = Bplin(index)/Bpexp(index);
    % in case something goes south due to very small numbers set the
    % amplitude to 0
    if isinf(scale_point) || isnan(scale_point)
        scale_point = 0;
    end
    % the final amplitude vector
    Bp = [Bplin(1:index-1); scale_point * Bpexp(index:end)];
else
    % linear part
    Bplin = (-Bmax/T)*t + Bmax;
    % exponential part
    Bpexp = exp(-t /(Bstar*T/Bmax));
    % Bstar time tstar
    tstar = (Bstar-Bmax)/(-Bmax/T);
    % amplitude at tstar for scaling
    Btstar = exp(-tstar /(Bstar*T/Bmax));
    % apply
    if t<tstar
        Bp = Bplin;
    else
        Bp = (Bstar/Btstar) * Bpexp;
    end
end

% if due to division by "0" the value is NaN ... set it to 0
Bp(isnan(Bp)) = 0;


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
