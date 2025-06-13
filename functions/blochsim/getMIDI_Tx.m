function [t,Bxy,AP] = getMIDI_Tx(param)
%getMIDI_Tx creates discrete on-resonant or adiabatic pulses
%
% Syntax:
%       getMIDI_Tx(param)
%
% Inputs:
%       param - struct containing the pulse settings
%
% Outputs:
%       t - time vector [s]
%       Bxy - pulse amplitudes (x,y)
%       AP - struct containing adiabatic pulse settings
%
% Example:
%       [t,y,~] = getMIDI_Tx(param)
%
% Other m-files required:
%       getPulseTimeSeries
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

% input data
% Larmor freq. [Hz]
fL = param.fL;
% gyromagnetic ratio [rad/s/T]
gamma = param.gamma;
% sampling freg. [Hz]
sf = param.sf;
% number of periods
P = param.P;
% Tx current
I = param.I;
% duty cycle min
DCmin = param.DCmin;
% duty cycle max
DCmax = param.DCmax; 

% time sampling due to sampling frequency [s]
dt = 1/sf;

switch param.Tx
    case 'MIDI_OR'
        % creating on-resonant pulses is straight forward
        
        % Larmor sampling + offset [s]
        % this determines the time of a full period with the given
        % frequency
        dtL = 1/(fL+param.df); 
        
        % time vector [s]
        % total length of pulse depending on the number of periods and the
        % time length of one single period
        t = 0:dt:P*dtL; 
        
        % positive (up) and negative (down) "Larmor" peaks
        % we always start with a negative (down) one
        tLup = 3*dtL/4:dtL:P*dtL;
        tLdown = dtL/4:dtL:P*dtL;
        
        % peak width of positive and negative peaks depending on the duty
        % cycle
        peakwidthup = DCmax*dtL/2.*ones(size(tLup));
        peakwidthdown = DCmax*dtL/2.*ones(size(tLup));
        
        % no adiabatic pulse data
        AP = 0;
        
    case 'MIDI_AP'
        % creating the adiabatic pulses is more tricky because due to the
        % frequency modulation it is not clear at the beginning how long
        % the actual pulse is given the number of periods
        
        % Larmor sampling + offset [s]
        % this determines the time of a full period with the given
        % frequency
        dtL = 1/(fL); 
        
        % time vector [s]
        % total length of pulse depending on the number of periods and the
        % time length of one single period
        t = 0:dt:P*dtL; 

        % the adiabatic pulse is assembled in an iterative approach
        % maybe there is a more neat way, but I could not think of one
        search = true;
        count = 0;
        while search
            % first create a continuous adiabatic pulse with the modulation
            % settings from the GUI
            
            % temporary time vector [s]
            tt = t(:);
            
            % get the modulation functions
            fmod = param.fmod;
            Imod = param.Imod;
            
            % update the time vector and start and end point
            fmod.t = tt;
            fmod.t0 = tt(1);
            fmod.t1 = t(end);
            
            Imod.t = tt;
            Imod.t0 = tt(1);
            Imod.t1 = tt(end);
            
            % standard adiabatic pulse settings
            pparam.PulseType = 'AHP';
             % gyromagnetic ratio [rad/s/T]
            pparam.gamma = gamma;
            % angular frequency [rad/s]
            if gamma > 0
                pparam.omega0 = -fL*2*pi;
            else
                pparam.omega0 = fL*2*pi;
            end
            % normalized pulse amplitude
            pparam.Amp = 1;
            % pulse frequency modulation [struct]
            pparam.fmod = fmod;
            % pulse current modulation [struct]
            pparam.Imod = Imod;
            % auxiliary pulse phase [rad]
            pparam.phi = 0;
            % pulse axis [string]
            if isfield(param,'PulseAxis')
                pparam.PulseAxis = param.PulseAxis;
            else
                pparam.PulseAxis = '+y';
            end
            % pulse polarization [string]
            pparam.PulsePolarization = 'linear';            
            % temporary time vector [s]
            pparam.t = tt;
            % create the continuous adiabatic pulse
            [AP.Bout,AP.df,AP.I,AP.theta] = getPulseTimeSeries(pparam);
            
            % find the indices of all peaks within the pulse
            [~,locs] = findpeaks(abs(AP.Bout(:,1)));
            % take only the positive one
            indup = AP.Bout(locs,1)>0;
            % if there are more peaks than number of periods we are done
            if sum(indup) >= P
                search = false;
            else
                % other wise extend the time vector and start again
                count = count + 2;
                t = 0:dt:(P+count)*dtL;
            end
        end
        
        % find the last positive peak that is equal to the number of
        % periods
        ind = find(cumsum(indup)==P,1,'first');
        % keep all peaks (up & down) within this range
        locs = locs(1:ind);
        % set the pulse length to the corresponding time and add a quarter
        % of a period
        tmax = tt(locs(end)) + dtL/4;
        
        % trim the temporary time vector
        tt = tt(tt<=tmax);
        % trim the output data
        AP.Bout = AP.Bout(1:numel(tt),:);
        AP.df = AP.df(1:numel(tt),:);
        AP.I = AP.I(1:numel(tt),:);
        AP.theta = AP.theta(1:numel(tt),:);
        
        % update the original time vector
        t = tt';
        
        % find all positive peaks
        indup = AP.Bout(locs,1)>0;
        % find all negative peaks
        inddown = AP.Bout(locs,1)<0;
        % find the corresponding time samples
        tLup = tt(locs(indup));
        tLdown = tt(locs(inddown));
        
        % the "width" of the peak in [s]
        dtLup = diff([0;tLup]);
        dtLdown = diff([0;tLdown]);
        
        % the actual peak width is the width in [s] scaled by the
        % corresponding current amplitude of the continuous adiabatic pulse
        % therewith the duty cycle increases as given by the current
        % modulation function
        peakwidthup = AP.I(locs(indup)).*(dtLup./2);
        peakwidthdown = AP.I(locs(inddown)).*(dtLdown/2);
        
        % save all output data
        AP.locs = locs;
        AP.indup = indup;
        AP.inddown = inddown;
        AP.tLup = tLup;
        AP.tLdown = tLdown;
        AP.peakwidthup = peakwidthup;
        AP.peakwidthdown = peakwidthdown;
end

% loop over all negative peaks
for i = 1:numel(tLdown)
    % find the discrete time steps that are within the range of peak
    % width
    ind = find(t>=tLdown(i)-peakwidthdown(i)/2 & t<tLdown(i)+peakwidthdown(i)/2);
    % store the indices within the time vector
    if i == 1
        locd = ind;
    else
        locd = [locd ind]; %#ok<AGROW>
    end
end

% loop over all negative peaks
for i = 1:numel(tLup)
    % find the discrete time steps that are within the range of peak
    % width
    ind = find(t>=tLup(i)-peakwidthup(i)/2 & t<tLup(i)+peakwidthup(i)/2);
    % store the indices within the time vector
    if i == 1
        locu = ind;
    else
        locu = [locu ind]; %#ok<AGROW>
    end
end

% create the amplitude vector
x = zeros(size(t));
% all negative peak positions get the corresponding Tx amplitude (scaled
% by the Tx current I)
x(locd) = -I*ones(size(locd));
% all positive peak positions get the corresponding Tx amplitude (scaled
% by the Tx current I)
x(locu) = I*ones(size(locu));

% now create an artificial y-component by shifting the original signal
% (x-component) 90° (quarter of a full period)
% get the length of a quarter period
t90 = dtL/4;
% cut out these first samples and put them to the end of the original signal
y = [x(t>=t90) x(t<t90)];
% merge the x- and y-components
Bxy = [x;y];

% because the created pulse is a +y-pulse adjust its pulse axis here if
% necessary (only the on-resonant case here because for the adiabatic
% pulse this is already taken care of during assembling)
if strcmp(param.Tx,'MIDI_OR') && isfield(param,'PulseAxis')
    % depending on the frequency, there is a different amount of time
    % steps for the "90° phase shift"
    shiftind = sum(t<t90);
    % now shift both components by the necessary amount of samples
    switch param.PulseAxis
        case '+x'
            tmp = circshift(Bxy',-shiftind);
        case '+y'
            tmp = Bxy'; % nothing to do
        case '-x'
            tmp = circshift(Bxy',+shiftind);
        case '-y'
            tmp = circshift(Bxy',+2*shiftind);
    end
    % update the pulse amplitudes
    Bxy = tmp';
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
