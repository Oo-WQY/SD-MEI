function theta = getPulsePhase(t,df,param,flag)
%getPulsePhase provides the instantaneous phase of a pulse; this is
%needed for frequency modulated pulses (e.g. AHP) because the frequency is
%actually modulated via the phase
%NOTE: because f = dphi/dt*2pi the time domain phase is the integral of the
%frequency:
%phi(t) = phi0 + 2pi*int_0^t f(tau) dtau
%
%this means e.g. for a linear frequency chirp from f0 to f1 like:
% f(t) = k*t + f0, with slope k = (f1-f0)/(t1-t0)
%
% the instantaneous phase is given as
% phi(t) = phi0 + 2pi*(k/2*t^2 + f0*t)
%
% Syntax:
%       getPulsePhase(t,df,param,flag)
%
% Inputs:
%       t - time [s]
%       df - off-resonance [Hz]
%       param - struct containing extra settings:
%               fmod   : struct containing frequency modulation settings
%               gamma  : gyromagnetic ratio [rad/s/T]
%               omega0 : angular frequency [rad/s]
%       flag - 0 -> phase during pulse
%              1 -> phase after pulse (only relaxation)
%
% Outputs:
%       theta - instantaneous phase angle of the pulse
%
% Example:
%       getPulsePhase(t,0,param,0)
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

% modulation parameter
fmod = param.fmod;
% Larmor freq. [Hz]
f0 = param.omega0/2/pi;

% manipulating the phase angle theta modulates the frequency of the pulse
switch fmod.shape
    case 'const'
        switch flag
            case 0 % during pulse
                theta = (f0-df)*2*pi.*t;
            case 1 % after pulse
                theta = -fmod.v1*2*pi*t;
        end
        
    case 'lin'
        k = (fmod.v0-fmod.v1)./(fmod.t1-fmod.t0);
        switch flag
            case 0 % during pulse
                theta = 2*pi*((f0-fmod.v0).*t + k./2.*t.^2);
            case 1 % after pulse
                theta = -2*pi*(fmod.v1*t + k/2*t^2);
        end
        
    case 'tanhMIDI'
        delta_f = fmod.v0-fmod.v1;
        delta_t = fmod.t1-fmod.t0;
        A = 2*pi*fmod.A;
        B = fmod.B;
        C = tanh( (A./fmod.t1).*(fmod.t0-B.*(fmod.t1/2)));
        D = tanh( (A./fmod.t1).*(fmod.t1-B.*(fmod.t1/2)));
        E = fmod.v1;
        F = -delta_f;
        T = delta_t;
        % sign switch
        delta_val = -delta_f;
        switch flag
            case 0 % during pulse
                theta = 2.*pi.*(f0.*t + (( A.*B.*C.*F.*T + 2.*F.*T.*log(cosh((A.*(t-(B.*T./2)))./T)) + 2.*A.*t.*(C.*E-D.*(E+F)) ) ./ (2.*A.*(C-D))) );
            case 1 % after pulse
                theta = 2.*pi.*(fmod.v1*t + (( A.*B.*C.*F.*T + 2.*F.*T.*log(cosh((A.*(t-(B.*T./2)))./T)) + 2.*A.*t.*(C.*E-D.*(E+F)) ) ./ (2.*A.*(C-D))) );
        end
        
    case 'tanhGMR'
        delta_f = fmod.v0-fmod.v1;
        delta_t = fmod.t1-fmod.t0;
        switch flag
            case 0 % during pulse
                theta = 2*pi*( (f0-fmod.v0).*t + (delta_t/3).*delta_f.*log(cosh(3.*t./delta_t)) );
            case 1 % after pulse
                theta = -2*pi*(fmod.v0*t - (delta_t/3)*delta_f*log(cosh(3*t/delta_t)) );
        end
        
    case 'exp'
        delta_f = fmod.v0-fmod.v1;
        delta_t = fmod.t1-fmod.t0;
        switch flag
            case 0 % during pulse
                theta = 2*pi*(f0.*t + (delta_t*delta_f./fmod.A).*exp(fmod.A.*(-t/delta_t)));
            case 1 % after pulse
                theta = 2*pi*(fmod.v1*t + (delta_t*delta_f./fmod.A)*exp(fmod.A*(-t/delta_t)));
        end
        
    otherwise
        switch flag
            case 0 % during pulse
                theta = (f0-df)*2*pi.*t;
            otherwise
                % nothing to do
        end
        
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
