function getRampParameters(fig)
%getRampParameters updates all relevant switch-off ramp settings
%
% Syntax:
%       getRampParameters(fig)
%
% Inputs:
%       fig - figure handle
%
% Outputs:
%       none
%
% Example:
%       getRampParameters(gui.figh)
%
% Other m-files required:
%       getAngleBetweenVectors
%       getRampAmplitude
%       getRotationMatrixFromAngleandAxis
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

% get GUI data
gui  = getappdata(fig,'gui');
data = getappdata(fig,'data');

% z-axis unit vector
zunit = [0 0 1]';

% --- switch-off ramp settings ---
% switch-off ramp type [string]
rampparam.ramp = data.prepol.Ramp;
% gyromagnetic ratio [rad/s/T]
rampparam.gamma = data.basic.gamma;
% primary (Earth's) magnetic field amplitude B0 [T]
rampparam.B0 = data.basic.B0;
% maximum pre-polarization amplitude (factor*B0) [T]
rampparam.Bmax = data.basic.B0*data.prepol.Factor;
% switch-over amplitude for the "linexp" ramp (factor*B0) [T]
rampparam.Bstar = data.basic.B0*data.prepol.SwitchFactor;
% switch-off ramp time [s]
rampparam.Tramp = data.prepol.Tramp/1e3;
% switch-over time for the "linexp" ramp [s]
rampparam.Tslope = data.prepol.Tslope/1e3;
% switch-off ramp time vector [s] discretized with 500kHz
t = (0:1/500000:rampparam.Tramp)';
% get the amplitude of the pre-polarization field [T]
% this amplitude decreases over time due to the particular switch-off ramp
Bp = getRampAmplitude(t,rampparam);

% now adjust the direction of the pre-polarization field
% 1.) rotation by angle theta [deg] around y-axis for the z-unit vector
% this means an angle of 90° around the y-axis will turn the z-unit into
% the x-unit vector
RM = getRotationMatrixFromAngleandAxis(deg2rad(data.prepol.Theta),[0 1 0]);
orient = RM * zunit;
% now this new orientation vector gets rotated by  angle phi [deg] around
% the z-axis for the x-unit vector
% this means an angle of 90° around the z-axis will turn the x-unit into
% the y-unit vector
RM = getRotationMatrixFromAngleandAxis(deg2rad(data.prepol.Phi),[0 0 1]);
orient = RM * orient;

% primary (Earth's) magnetic field vector [T]
Be = rampparam.B0.*zunit;
% pre-polarization field vector oriented into the correct direction [T]
Bpre = Bp.*repmat(orient',[numel(Bp) 1]);
% the effective B-field vector as a combination of primary and
% pre-polarization field [T]
Beff = Be' + Bpre;

% angle between primary and pre-polarization field [rad]
alpha = getAngleBetweenVectors(repmat(Be',[size(Beff,1) 1]),Beff);
% amplitude of the effective B-field [T]
Beffn = sqrt(Beff(:,1).^2+Beff(:,2).^2+Beff(:,3).^2);
% angular frequency of the effective B-field [rad/s]
omega = rampparam.gamma.*Beffn;
% rate of change of the angle alpha [rad/s]
dt = t(2)-t(1);
dadt = abs(diff(alpha)./dt);

% save data
data.results.prepol.orient = orient;
data.results.prepol.t = t;
data.results.prepol.Bmax = rampparam.Bmax;
data.results.prepol.Bstar = rampparam.Bstar;
data.results.prepol.Bp = Bp;
data.results.prepol.Beff = Beff;
data.results.prepol.alpha = alpha;
data.results.prepol.omega = omega;
data.results.prepol.dadt = dadt;

% because the ramp data changed, deactivate the "Animate" button
set(gui.push_handles.Animate,'Enable','off');

% update GUI data
setappdata(fig,'data',data);
end

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
