function getPulseParameters(fig)
%getPulseParameters updates all relevant pulse settings
%
% Syntax:
%       getPulseParameters(fig)
%
% Inputs:
%       fig - figure handle
%
% Outputs:
%       none
%
% Example:
%       getPulseParameters(gui.figh)
%
% Other m-files required:
%       none
%
% Subfunctions:
%       none
%
% MAT-files required:
%       getMIDI_Tx
%       getOmega0
%       getPulseTimeSeries
%
% See also: BLOCHUS
% Author: Thomas Hiller
% email: thomas.hiller[at]leibniz-liag.de
% License: GNU GPLv3 (at end)

%------------- BEGIN CODE --------------

% get GUI data
gui = getappdata(fig,'gui');
data = getappdata(fig,'data');

% dummy time vector with sampling freq 50 kHz used mainly for plotting
pulse_t = 0:1/50000:data.pulse.Ttau/1e3; % [s]

% --- frequency modulation struct ---
% this is also used when no modulation is done (then all settings are
% basically constant)
% pulse type [string]
fmod.PulseType = data.pulse.Type;
% shape of the modulation function [string]
fmod.shape = data.pulse.DFmode;
% time vector [s]
fmod.t = pulse_t;
% pulse start time [s]
fmod.t0 = pulse_t(1);
% pulse end time [s]
fmod.t1 = pulse_t(end);
% start frequency [Hz] 
fmod.v0 = data.pulse.DFstart;
% end frequency [Hz]
fmod.v1 = data.pulse.DFend;
% modulation parameter (MIDI)
fmod.A = data.pulse.DFA;
% modulation parameter (MIDI)
fmod.B = data.pulse.DFB;

% --- current modulation struct ---
% this is also used when no modulation is done (then all settings are
% basically constant)
% pulse type [string]
Imod.PulseType = data.pulse.Type;
% shape of the modulation function [string]
Imod.shape = data.pulse.Imode;
% check if quality factor tuning is activated [0/1]
Imod.useQ = get(gui.check_handles.PulseQ,'Value');
% quality factor
Imod.Q = data.pulse.Q;
% quality factor off-resonance frequency [Hz]
Imod.Qdf = data.pulse.Qdf;
% time vector [s]
Imod.t = pulse_t;
% pulse start time [s]
Imod.t0 = pulse_t(1);
% pulse end time [s]
Imod.t1 = pulse_t(end);
% start current [A] 
Imod.v0 = data.pulse.Istart;
% end current [A] 
Imod.v1 = data.pulse.Iend;
% modulation parameter (MIDI)
Imod.A = data.pulse.IA;
% modulation parameter (MIDI)
Imod.B = data.pulse.IB;

% adjust some pulse settings depending on the pulse type
switch data.pulse.Type    
    case {'pi_half','pi','free'}
        % set constant frequency modulation shape 
        fmod.shape = 'const';
        % equal start and end frequency
        fmod.v1 = fmod.v0;
        % set constant current modulation shape 
        Imod.shape = 'const';
        % equal start and end current
        Imod.v1 = Imod.v0;
        
    case 'MIDI_OR'
        % set constant frequency modulation shape
        fmod.shape = 'const';
        % equal start and end frequency
        fmod.v1 = fmod.v0;
        % set constant current modulation shape 
        Imod.shape = 'const';
        % equal start and end current
        Imod.v1 = Imod.v0;

        % MIDI parameter settings for discrete pulses
        % pulse type [string]
        mparam.Tx = data.pulse.Type;
        % Larmor frequency [Hz] WIHTOUT SIGN
        mparam.fL = abs(data.basic.Omega0);
        % gyromagnetic ratio [rad/s/T]
        mparam.gamma = data.basic.gamma;
        % frequency offset [Hz]
        mparam.df = data.pulse.DFstart;
        % MIDI sampling frequency [Hz] (default is 50kHz)
        mparam.sf = data.pulse.MIDIsf;
        % number of periods
        mparam.P = data.pulse.MIDINP;
        % Tx-current [A] (within BLOCHUS always "1" because the pulse
        % amplitude is scaled from the GUI)
        mparam.I = 1;
        % duty cycle width - start
        mparam.DCmin = data.pulse.Istart;
        % duty cycle width - end set equal to start
        mparam.DCmax = mparam.DCmin;
        % pulse axis [string]
        mparam.PulseAxis = data.pulse.Axis;
        % assemble the discrete pulse
        [t,y,AP] = getMIDI_Tx(mparam);
        % because the ODE-solver does not like long rows of "zeros" due to
        % the adaptive time stepping, the "zeros" in the pulse are set to a
        % very small number
        ind = find(y==0);
        y(ind) = (2.*rand(numel(ind),1)-1)*1e-7;
        % output struct containing the MIDI pulses
        param.MIDI.mparam = mparam;
        param.MIDI.t = t';
        param.MIDI.y = y';
        param.MIDI.AP = AP;
        % update data
        data.pulse_MIDI = param.MIDI;
        % update dummy time vector
        pulse_t = t;
        
    case 'MIDI_AP'
        % MIDI parameter settings for discrete pulses
        % pulse type [string]
        mparam.Tx = data.pulse.Type;
        % Larmor frequency [Hz] WIHTOUT SIGN
        mparam.fL = abs(data.basic.Omega0);
        % gyromagnetic ratio [rad/s/T]
        mparam.gamma = data.basic.gamma;
        % start frequency [Hz]
        mparam.df = data.pulse.DFstart;
        % MIDI sampling frequency [Hz] (default is 50kHz)
        mparam.sf = data.pulse.MIDIsf;
        % number of periods
        mparam.P = data.pulse.MIDINP;
        % Tx-current [A] (within BLOCHUS always "1" because the pulse
        % amplitude is scaled from the GUI)
        mparam.I = 1;
        % duty cycle width - start
        mparam.DCmin = data.pulse.Istart;
        % duty cycle width - end
        mparam.DCmax = data.pulse.Iend;
        % pulse axis [string]
        mparam.PulseAxis = data.pulse.Axis;
        % frequency modulation struct
        mparam.fmod = fmod;
        % current modulation struct
        mparam.Imod = Imod;
        % assemble the discrete pulse
        [t,y,AP] = getMIDI_Tx(mparam);
        % because the ODE-solver does not like long rows of "zeros" due to
        % the adaptive time stepping, the "zeros" in the pulse are set to a
        % very small number
        ind = find(y==0);
        y(ind) = (2.*rand(numel(ind),1)-1)*1e-7;
        % output struct containing the MIDI pulses
        param.MIDI.mparam = mparam;
        param.MIDI.t = t';
        param.MIDI.y = y';
        param.MIDI.AP = AP;
        % update data
        data.pulse_MIDI = param.MIDI;
        % update dummy time vector
        pulse_t = t;
        
        % update the time vector within the modulation functions
        fmod.t = t;
        fmod.t1 = t(end);
        Imod.t = t;
        Imod.t1 = t(end);
        
        % adjust the pulse length to the actual value
        data.pulse.Ttau = data.pulse_MIDI.t(end)*1e3;
        % in case of pure "Pulse" update the total simulation time
        switch data.basic.type
            case 'pulse'
                data.basic.Tsim = data.pulse.Ttau;
                set(gui.edit_handles.Tsim,'String',num2str(data.basic.Tsim));
            otherwise
                % nothing to do
        end
    otherwise
        % nothing to do        
end

% --- pulse settings ---
% pulse type [string]
param.PulseType = data.pulse.Type;
% gyromagnetic ratio [rad/s/T]
param.gamma = data.basic.gamma;
% Larmor frequency [rad/s]
param.omega0 = getOmega0(data.basic.gamma,data.basic.B0);
% pulse amplitude [B0]
param.Amp = data.basic.B0*data.pulse.B1Factor;
% pulse frequency modulation [struct]
param.fmod = fmod;
% pulse current modulation [struct]
param.Imod = Imod;
% auxiliary pulse phase [rad]
param.phi = 0;
% pulse axis [string]
param.PulseAxis = data.pulse.Axis;
% pulse polarization [string]
param.PulsePolarization = data.pulse.Polarization;
% time vector [s]
param.t = pulse_t;
% get the actual pulse time series (including frequency and current modulation)
[pulse_Bxy,df,I] = getPulseTimeSeries(param);
% dummy time vector [ms]
data.results.pulse.t = pulse_t'.*1e3;

% get FFT of the pulse
[Xb,fbx] = getFFT(pulse_t(:),pulse_Bxy(:,1:2));
                
% update data
data.results.pulse.fmod = fmod;
data.results.pulse.Imod = Imod;
data.results.pulse.df = df;
data.results.pulse.I = I;
data.results.pulse.Bxy = pulse_Bxy;
data.results.pulse.Bspec.fx = fbx;
data.results.pulse.Bspec.X = Xb;

% because the pulse data changed, deactivate the "Animate" button
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
