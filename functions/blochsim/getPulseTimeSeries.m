function [Bout,df,I,theta] = getPulseTimeSeries(param)
%getPulseTimeSeries returns the B-field amplitudes of the pulse either for a
%single point or a complete time series
%
% Syntax:
%       getPulseTimeSeries(param)
%
% Inputs:
%       param - struct containing the pulse settings
%               PulseType         : pulse type [string]
%               gamma             : gyromagnetic ratio [rad/s/T]
%               omega0            : angular frequency [rad/s]
%               t                 : time [s]
%               Amp               : pulse amplitude
%               phi               : optional phase [rad]
%               PulseAxis         : pulse axis direction [string]
%               PulsePolarization : pulse polarization [string]
%               fmod              : struct containing the frequency modulation
%                                   settings
%               Imod              : struct containing the current modulation
%                                   settings
%               MIDI              : [optional] struct containing the predetermined
%                                   discrete MIDI pulse
%
% Outputs:
%       Bout - B-field amplitudes
%       df - frequency modulation [Hz]
%       I - current modulation [A]
%       theta - instantaneous phase [rad]
%
% Example:
%       getPulseTimeSeries(param)
%
% Other m-files required:
%       modulatePulse
%       getPulseAxisPhase
%       getReferencePhase
%       getPulsePhase
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

PulseType = param.PulseType;
gamma = param.gamma;
t = param.t;
Amp = param.Amp;
phi = param.phi;
fmod = param.fmod;
Imod = param.Imod;
Trep = param.basic.Trep/1e3;
Tp = param.basic.Tp/1e3;
[T_single,N_now,Flag_on] = getsingletime(t,Trep,Tp);
% if the pulse axis is not given, it is set to '+x'
if isfield(param,'PulseAxis')
    PulseAxis = param.PulseAxis;
else
    PulseAxis = '+x';
end
% if the pulse polarization is not given, it is set to 'circular'
if isfield(param,'PulsePolarization')
    PulsePolarization = param.PulsePolarization;
else
    PulsePolarization = 'circular';
end

switch PulseType
    % handle the 'discrete' MIDI pulses
    case {'MIDI_OR','MIDI_AP'}
        t = t(:);
        fmod.t = t;
        Imod.t = t;
        df = modulatePulse(fmod,'df');
        I = modulatePulse(Imod,'I');        
        
        % interpolate amplitudes at current time t from the pulse stored in
        % param.MIDI
        xMIDI = interp1(param.MIDI.t,param.MIDI.y(:,1),t,'linear');
        yMIDI = interp1(param.MIDI.t,param.MIDI.y(:,2),t,'linear');
        
        % get instantaneous phase due to frequency modulation
        % here the phase angle theta is simply the integral of the frequency
        % modulation; this obviously only works for a vector (so this is
        % not done during the ode calculation)
        if numel(t)>1
            ft = -param.MIDI.mparam.fL-df;
            theta = cumtrapz(t,ft.*2*pi);
        end        
        
        % I is already incorporated in the MIDI fit
        switch PulsePolarization
            case 'circular'
                Bout = Amp.*[xMIDI yMIDI];
            case 'linear'
                Bout = 2.*Amp.*[xMIDI 0.*t];
        end
        
    otherwise % all other (continuous) pulse types
        
        % direction of pulse axis
        phi_ax = getPulseAxisPhase(PulseAxis);
        % reference phase due to gyromagnetic ratio
        % NOTE: due to the reference phase, M flips according to the
        % right-hand-rule in the rotating frame of reference
        % -> also protons with negative gyromagnetic ratio!
        % refer to Levitt, 2002
        phi_ref = getReferencePhase(gamma);

        % pulse modulation
        t = t(:);
        fmod.t = t;
        Imod.t = t;
        df = modulatePulse(fmod,'df');
        I = modulatePulse(Imod,'I');
        
        % check for quality factor tuning 用于绝热脉冲，调节动态变化
%         if Imod.useQ && Imod.Q > 0
%             % get line (band) width -> f_L/Q (simple bandwidth for bandpass)
%             Lwidth = abs(param.omega0/2/pi) / Imod.Q;
%             % apply "Breit-Wigner" formula (here already normalized to 1 by
%             % multiplying with pi*Lwidth)
%             L = 1 ./ ( ((df+Imod.Qdf).^2 ./ Lwidth.^2) + 1 );
%             I = I.*L;
%         end
           
        % get instantaneous phase due to frequency modulation
        theta = getPulsePhase(T_single,df,param,0);
        
        switch PulsePolarization
            case 'circular'
                if  Flag_on==1
                if  T_single>0.01 & T_single<0.01
                Bout = 0.*[cos(theta + phi + phi_ax + phi_ref) ...
                    sin(theta + phi + phi_ax + phi_ref)];
                else
                Bout = I.*Amp.*[cos(theta + phi + phi_ax + phi_ref) ...
                    sin(theta + phi + phi_ax + phi_ref)];
                end
                else
                Bout = 0.*[cos(theta + phi + phi_ax + phi_ref) ...
                    sin(theta + phi + phi_ax + phi_ref)];
                end
            case 'linear'
                if  Flag_on==1
                Bout = 2.*I.*Amp.*[cos(theta + phi + phi_ax + phi_ref) ...
                    0.*t];
                else
                Bout = 0.*[cos(theta + phi + phi_ax + phi_ref) ...
                    0.*t];
                end
            case 'normal'
                if  Flag_on==1
                    Bout = I.*Amp.*[cos(theta + phi + phi_ax + phi_ref) ...
                    0.*t];
                else
                    Bout = 0.*[cos(theta + phi + phi_ax + phi_ref) ...
                    0.*t];
                end
            case 'rot'
                if  Flag_on==1
                    Bout = I.*Amp.*[1.*ones(size(t,1),size(t,2)) ...
                    0.*t];
                else
                    Bout = 0.*[1.*t ...
                    0.*t];
                end
        end
end
% if B-field values are NaN, set them to zero
% this can happen due to the interpolation at the end of the MIDI-pulses
Bout(isnan(Bout)) = 0;

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
