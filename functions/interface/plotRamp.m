function plotRamp(fig)
%plotRamp plots different pre-polarization switch-off ramp parameter
%
% Syntax:
%       plotRamp(fig)
%
% Inputs:
%       fig - figure handle
%
% Outputs:
%       none
%
% Example:
%       plotRamp(fig)
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
t = data.results.prepol.t.*1e3;

%% pre-polarization switch-off ramp amplitude
% Bp amplitude in units of the Earth's magnetic field B0
value = data.results.prepol.Bp./data.basic.B0;
% min-max-spread
d = max(value)-min(value);
ax = gui.axes_handles.Bpre;
cla(ax);
hold(ax,'on');
plot(t,value,'LineWidth',gui.myui.linewidth,'Color',myui.color.prepolB,...
    'Parent',ax);
line([min(t) max(t)],[1 1],'Linestyle','--','Color','k','LineWidth',1,...
    'Parent',ax);
hold(ax,'off');
set(ax,'XLim',[min(t) max(t)],'YLim',[min(value)-d/20 max(value)+d/20]);
grid(ax,'on');
set(get(ax,'XLabel'),'String','t [ms]');
set(get(ax,'YLabel'),'String','Bp [B0]');
legend(ax,'Bp','Bp=B0','Location','NorthEast');
set(ax,'FontSize',myui.axfontsize);

%% angle alpha between B and B0
value = rad2deg(data.results.prepol.alpha);
% min-max-spread
d = max(value)-min(value);
ax = gui.axes_handles.alpha;
cla(ax);
plot(t,value,'LineWidth',gui.myui.linewidth,'Color',myui.color.prepolB,...
    'Parent',ax);
set(ax,'XLim',[min(t) max(t)],'YLim',[min(value)-d/20 max(value)+d/20]);
grid(ax,'on');
set(get(ax,'XLabel'),'String','t [ms]');
set(get(ax,'YLabel'),'String','\alpha [deg]');
legend(ax,'\alpha=\angle B0 B','Location','best');
set(ax,'FontSize',myui.axfontsize);

%% time derivative of alpha
value = data.results.prepol.dadt;
% min-max-spread
d = max(value)-min(value);
ax = gui.axes_handles.dadt;
cla(ax);
semilogy(t(1:length(data.results.prepol.dadt)),value,'LineWidth',gui.myui.linewidth,...
    'Color',myui.color.prepolB,'Parent',ax);
set(ax,'XLim',[min(t) max(t)],'YLim',[min(value)-d/20 max(value)+d/20]);
grid(ax,'on');
set(get(ax,'XLabel'),'String','t [ms]');
set(get(ax,'YLabel'),'String','d\alpha / dt');
set(ax,'FontSize',myui.axfontsize);

%% ratio of da/dt and omega=gamma*B
% ratio = dadt./omega
% ratio << 1 -> adiabatic condition (see Melton 1995, eq. 1)
value = data.results.prepol.dadt./data.results.prepol.omega(1:numel(data.results.prepol.dadt));
% min-max-spread
d = max(value)-min(value);
ax = gui.axes_handles.wda;
cla(ax);
hold(ax,'on');
plot(t(1:numel(data.results.prepol.dadt)),value,'LineWidth',gui.myui.linewidth,...
    'Color',myui.color.prepolB,'Parent',ax);
line([min(t) max(t)],[1 1],'Linestyle','--','Color','k','LineWidth',1,...
    'HandleVisibility','off','Tag','MarkerLines','Parent',ax);
hold(ax,'off');
set(ax,'XLim',[min(t) max(t)],'YLim',[min(value)-d/20 max(value)+d/20]);
grid(ax,'on');
set(get(ax,'XLabel'),'String','t [ms]');
set(get(ax,'YLabel'),'String','(d\alpha/dt) / \gammaB');
lh = legend(ax,'$$\ll$$ 1 $$\rightarrow$$ adiabatic cond.','Location','best');
set(lh,'Interpreter','Latex');
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
