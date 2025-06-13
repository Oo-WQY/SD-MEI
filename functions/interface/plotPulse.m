function plotPulse(fig)
%plotPulse plots different pulse parameter
%
% Syntax:
%       plotPulse(fig)
%
% Inputs:
%       fig - figure handle
%
% Outputs:
%       none
%
% Example:
%       plotPulse(gui.figh)
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

% get GUI data
data = getappdata(fig,'data');
gui = getappdata(fig,'gui');
myui = gui.myui;

% time vector in [ms]
t = data.results.pulse.t;

%% pulse frequency and current modulation
value = data.results.pulse.df;
% min-max-spread
d = max(value)-min(value);
if d == 0 % check if min=max
    d = 20;
end
ax = gui.axes_handles.PulseSetupF;
cla(ax);
plot(t,value,'LineWidth',gui.myui.linewidth,'Color',myui.color.pulse,...
    'Parent',ax);
set(ax,'Xlim',[min(t) max(t)],'Ylim',[min(value)-d/20 max(value)+d/20]);
grid(ax,'on');
set(get(ax,'XLabel'),'String','t [ms]');
set(get(ax,'YLabel'),'String','df [Hz]');
legend(ax,'df mod.','Location','SouthEast');
set(ax,'FontSize',myui.axfontsize);

ax = gui.axes_handles.PulseSetupI;
cla(ax);
plot(t,data.results.pulse.I,'LineWidth',gui.myui.linewidth,'Color',myui.color.pulse,...
    'Parent',ax);
set(ax,'Xlim',[min(t) max(t)],'Ylim',[-0.05 1.05]);
grid(ax,'on');
set(get(ax,'XLabel'),'String','t [ms]');
set(get(ax,'YLabel'),'String','I [A]');
legend(ax,'I mod.','Location','SouthEast');
set(ax,'FontSize',myui.axfontsize);
% in case of discrete pulses the modulation is done via the duty cycle
if strcmp(data.pulse.Type,'MIDI_OR') || strcmp(data.pulse.Type,'MIDI_AP')
    set(ax,'Ylim',[-0.03 0.63]);
    set(get(ax,'YLabel'),'String','duty cycle');
end


%% pulse
value = data.results.pulse.Bxy./data.basic.B0;
% min-max-spread
d = max(value(:))-min(value(:));
ax = gui.axes_handles.PulseB;
cla(ax);
hold(ax,'on');
plot(t,value(:,1),'r','LineWidth',gui.myui.linewidth,'Parent',ax);
plot(t,value(:,2),'g','LineWidth',gui.myui.linewidth,'Parent',ax);
hold(ax,'off');
set(ax,'Xlim',[min(t) max(t)]);
set(ax,'Ylim',[min(value(:))-d/20 max(value(:))+d/20]);
grid(ax,'on');
set(get(ax,'XLabel'),'String','t [ms]');
set(get(ax,'YLabel'),'String','B_{1} [B_0]');
legend(ax,'x','y','Location','SouthWest');
set(ax,'FontSize',myui.axfontsize);

%% FFT
% Larmor freq. [Hz]
fL = getOmega0(data.basic.gamma,data.basic.B0)/2/pi;

ax = gui.axes_handles.PulseFFT;
f = data.results.pulse.Bspec.fx;
X = data.results.pulse.Bspec.X;
% plot data
clearSingleAxis(ax);
hold(ax,'on');
plot(f,abs(X),'r','Parent',ax);
% vertical line indicating Larmor frequency
line([fL fL],[0 max(abs(X))],'Color','k','LineStyle','--',...
    'LineWidth',0.75,'Parent',ax);
hold(ax,'off');
% axis settings
set(ax,'XLim',[-abs(2*fL) abs(2*fL)],'YLim',[0 max(abs(X))].*1.1);
grid(ax,'on');
set(get(ax,'XLabel'),'String','F [Hz]');
set(get(ax,'YLabel'),'String','amplitude');
% legend
legend(ax,'B','\omega_0','Location','NorthEast');
% font size
set(ax,'FontSize',myui.axfontsize);

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
