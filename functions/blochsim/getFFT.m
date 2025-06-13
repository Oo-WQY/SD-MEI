function [Y,f] = getFFT(t,s)
%getFFT calculates the FFT for a given time and signal pair
%
% Syntax:
%       getFFT(t,m)
%
% Inputs:
%       t - time [s]
%       s - signal
%
% Outputs:
%       Y - amplitudes of the FFT
%       f - frequency [Hz]
%
% Example:
%       getFFT(t,s)
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

% length of the time series [s]
Tmax = t(end);
% NOTE: because the ode-solver uses adaptive time-stepping, the time vector
% may be irregularly spaced, so we use the shortest time step
dt = diff(t);
dt = max([1e-6 abs(min(dt))]); % [s] max. 1MHz
% dt = 1e-7;
% sampling frequency [Hz]
FS = 1/dt;
% new time vector with regular time stepping [s]
tt=0:1/FS:Tmax;
% number of time steps [-]
Nt = numel(tt);
% if the signal is complex we need to interpolate both parts
if size(s,2)== 2
    % interpolate signal to new time steps
    SxI = interp1(t,s(:,1),tt);
    SyI = interp1(t,s(:,2),tt);
    % new complex signal
    y = complex(SxI,SyI);
else
    % interpolate signal to new time steps
    y = interp1(t,s,tt);
end
% complex FFT
YY = fft(y);
% apply fftshift to center the spectrum correctly
Y = fftshift(YY)./Nt;
% frequency sampling [Hz]
f = FS/2*linspace(-1,1,Nt);

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
