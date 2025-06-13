function dM = fcn_BLOCHUS_ode_y(t,m,param)
%fcn_BLOCHUS_ode is the objective function for the ode-solver which solves
%the Bloch equation in the laboratory frame of reference
%
% Syntax:
%       fcn_BLOCHUS_ode(t,m,param)
%
% Inputs:
%       t - instantaneous time
%       m - instantaneous magnetization
%       param - struct that holds additional settings:
%                   type       : flag between 'std', 'prepol',
%                               'pulse', 'prepolpulse'
%                   M0         : equilibrium magnetization [A/m]
%                   B0         : Earth magnetic field [T]
%                   T1         : relaxation times [s]
%                   T2         : relaxation times [s]
%                   gamma      : gyromagnetic ratio [rad/s/T]
%                   rampparam  : struct holding switch-off settings
%                   pulseparam : struct holding pulse settings
%
% Outputs:
%       dM - time derivative of m
%
% Example:
%       fcn_BLOCHUS_ode(t,m,param)
%
% Other m-files required:
%       getRampAmplitude
%       getPulseTimeSeries
%       getRotationMatrixFromVectors
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

% unit vectors
xunit = [1 0 0]';
yunit = [0 1 0]';
zunit = [0 0 1]';

% basic parameters
M0 = param.M0;
B0 = param.B0;
T1 = param.T1;
T2 = param.T2;
gamma = param.gamma;

switch param.type
    case 'std'
        B = B0*zunit;        
        % dM/dt
        dM = gamma*cross(m,B) - ( (m(1)*xunit+m(2)*yunit) / T2 ) - ( ((m(3)-M0(3))*zunit) / T1 );
        
    case 'prepol'
        % switch-off ramp parameters
        rampparam = param.rampparam;
        rampparam.B0 = B0;
        rampparam.gamma = gamma;
        % decreasing Bp amplitude over time
        Bp = getRampAmplitude(t,rampparam);
        % B0 - Earth field
        Be = B0*zunit;
        % check if we are during the switch-off or not
        if t <= rampparam.Tramp
            % during the switch-off we need to account for the decreasing
            % Bp-field amplitude
            B1 = Bp*param.orient;
            B = Be + B1;            
            % if RDS = 0 (switched-off) ignore T1 and T2 during switch-off
            % ramp
            if param.RDS == 0
                T1 = T1*1e6;
                T2 = T2*1e6;
            end            
            % if the B-field vector is not parallel to B0
            % rotate B and m into z-axis to apply relaxation
            if any(B./norm(B) ~= zunit)
                % get rotation matrix from Beff to z-axis
                R = getRotationMatrixFromVectors(B,zunit);
                B = R*B; m = R*m;
                % dM/dt
                dM = gamma*cross(m,B) - ( (m(1)*xunit+m(2)*yunit) / T2 ) - ( ((m(3)-M0(3))*zunit) / T1 );
                % rotate dM back
                dM = R'*dM;
            else
                % dM/dt
                dM = gamma*cross(m,B) - ( (m(1)*xunit+m(2)*yunit) / T2 ) - ( ((m(3)-M0(3))*zunit) / T1 );
            end
        else
            % pure B0 evolution after the switch-off
            B = Be;
            % dM/dt
            dM = gamma*cross(m,B) - ( (m(1)*xunit+m(2)*yunit) / T2 ) - ( ((m(3)-M0(3))*zunit) / T1 );
        end
        
    case 'pulse'
        pulseparam = param.pulseparam;
        pulseparam.t = t;
        % check if we are during the pulse or not
        if t <= param.Ttau       %更改为加载仿真时间!!!
            % get pulse amplitude over time
            B1 = getPulseTimeSeries(pulseparam);   
            Bx = B1(2);
            By = B1(1);
            B_df = param.pulseparam.B_df;
            B = [gamma*Bx,gamma*By,B_df]';
            % if RDP = 0 (switched-off) ignore T1 and T2 during pulse
            if param.RDP == 0
                T1 = T1*1e6;
                T2 = T2*1e6;
            end
        else
            % pure B0 evolution after the pulse
            B = B0*zunit;
        end        
        % dM/dt
        dM = cross(m,B) - ( (m(1)*xunit+m(2)*yunit) / T2 ) - ( ((m(3)-M0(3))*zunit) / T1 );
end

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