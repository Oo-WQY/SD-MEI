function plotResults(fig)
%plotResults plots results depending on the chosen settings
%
% Syntax:
%       plotResults(fig)
%
% Inputs:
%       fig - figure handle
%
% Outputs:
%       none
%
% Example:
%       plotResults(fig)
%
% Other m-files required:
%       plotBSphere
%
% Subfunctions:
%       plotMag
%       plotSphere
%       plotFFT
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

% update the different plots depending on the simulation type
switch data.basic.type
    
    case 'std' % only relaxation
        plotMag(data,gui,'lab');
        plotMag(data,gui,'rot');
        plotFFT(data,gui,'M');
        plotSphere(data,gui,'lab');
        plotSphere(data,gui,'rot');
        
    case 'prepol' % pre-polarization switch-off + (relaxation)
        plotMag(data,gui,'lab');
        plotSphere(data,gui,'lab');
        
        % clear rot-frame axis M
        clearSingleAxis(gui.axes_handles.MagR);
        % clear rot-frame axis sphere
        clearSingleAxis(gui.axes_handles.SphereR);
        
    case 'pulse' % pulse + (relaxation)
        plotMag(data,gui,'lab');
        plotMag(data,gui,'rot');
        plotFFT(data,gui,'M');
        plotSphere(data,gui,'lab');
        plotSphere(data,gui,'rot');
        plotFFT(data,gui,'B');
        
    case 'prepolpulse' % pre-polarization switch-off + pulse + (relaxation)
        plotMag(data,gui,'lab');
        plotMag(data,gui,'rot');
        plotFFT(data,gui,'M');
        plotSphere(data,gui,'lab');
        plotSphere(data,gui,'rot');
        plotFFT(data,gui,'B');
end

end

%% magnetization components
function plotMag(data,gui,frame)
myui = gui.myui;

% for plotting everything is in [ms]
T = data.results.basic.T.*1e3;

% all relevant time marker
Tsim = data.basic.Tsim;
Ttau = data.pulse.Ttau;
Tramp = data.prepol.Tramp;
Twait = data.pulse.Twait;

% choose reference frame to plot
switch frame
    case 'lab'
        M = data.results.basic.M;
        ax = gui.axes_handles.MagL;
    case 'rot'
        M = data.results.basic.Mrot;
        ax = gui.axes_handles.MagR;
end
% xy-component of magnetization vector
Mxy = sqrt(M(:,1).^2+M(:,2).^2);
% norm of magnetization vector
Mamp = sqrt(M(:,1).^2+M(:,2).^2+M(:,3).^2);

% normalize M in case of pre-polarization
% and adjust axis limits
switch data.basic.type
    case {'prepol','prepolpulse'}
        M = M./data.basic.B0;
        Mxy = Mxy./data.basic.B0;
        Mamp = Mamp./data.basic.B0;
        d = max([M(:);Mxy(:)])-min([M(:);Mxy(:)]);
        ymin = min([M(:);Mxy(:)])-d/20;
        ymax = max([M(:);Mxy(:)])+d/20;
    otherwise
        ymin = -1.05;
        ymax = 1.05;
end

% plot the data
clearSingleAxis(ax);
hold(ax,'on');
plot(T,M(:,1),'LineWidth',myui.linewidth,'Color','r','Parent',ax);
plot(T,M(:,2),'LineWidth',myui.linewidth,'Color','g','Parent',ax);
plot(T,M(:,3),'LineWidth',myui.linewidth,'Color','b','Parent',ax);
plot(T,Mxy,'LineWidth',myui.linewidth,'Color','m','Parent',ax);
plot(T,Mamp,'LineWidth',myui.linewidth,'Color','k',...
    'LineStyle','--','Parent',ax);
set(ax,'XLim',[min(T) max(T)],'YLim',[ymin ymax]);

% -------------------------------------------------------------------------
% plot vertical line at end of switch-off ramp
if (strcmp(data.basic.type,'prepol') || strcmp(data.basic.type,'prepolpulse'))...
        && Tramp < Tsim
    plot([Tramp Tramp],get(ax,'YLim'),'LineStyle','--',...
        'Color',myui.color.prepol,'LineWidth',1,'Tag','MarkerLines',...
        'HandleVisibility','off','Parent',ax);
end
% plot vertical line at end of pulse
if strcmp(data.basic.type,'pulse') && Ttau < Tsim
    plot([Ttau Ttau],get(ax,'YLim'),'LineStyle','--',...
        'Color',myui.color.pulse,'LineWidth',1,'Tag','MarkerLines',...
        'HandleVisibility','off','Parent',ax);
end
% plot vertical line at end of pulse
if strcmp(data.basic.type,'prepolpulse')
    plot([Ttau+Tramp+Twait Ttau+Tramp+Twait],get(ax,'YLim'),'LineStyle','--',...
        'Color',myui.color.pulse,'LineWidth',1,'Tag','MarkerLines',...
        'HandleVisibility','off','Parent',ax);
    if Twait > 0
         plot([Tramp+Twait Tramp+Twait],get(ax,'YLim'),'LineStyle','--',...
             'Color',myui.color.wait,'LineWidth',1,'Tag','MarkerLines',...
             'HandleVisibility','off','Parent',ax);
    end
end
% -------------------------------------------------------------------------
hold(ax,'off');

% axis settings
grid(ax,'on');
set(get(ax,'XLabel'),'String','t [ms]');
set(get(ax,'YLabel'),'String','magnetization M / M0');
% legend
lh = legend(ax,'x','y','z','|xy|','|M|','Location','SouthWest');

% adiabatic quality of switch-off ramp
if strcmp(data.basic.type,'prepol') || strcmp(data.basic.type,'prepolpulse')
%     set(get(lh,'Title'),'String',{'adiab. qual.',['p = ',sprintf('%4.3f',data.results.prepol.p)]})
    set(get(ax,'Title'),'String',['adiabatic quality p = ',sprintf('%4.3f',data.results.prepol.p)])
else
    set(get(ax,'Title'),'String','');
end
% font size
set(ax,'FontSize',myui.axfontsize);

end

%% magnetization components on Bloch sphere
function plotSphere(data,gui,frame)
myui = gui.myui;

% for plotting everything is in [ms]
T = data.results.basic.T.*1e3;

% all relevant time marker
Tsim = data.basic.Tsim;
Ttau = data.pulse.Ttau;
Tramp = data.prepol.Tramp;
Twait = data.pulse.Twait;

% choose reference frame to plot
switch frame
    case 'lab'
        M = data.results.basic.M;
        ax = gui.axes_handles.SphereL;
    case 'rot'
        M = data.results.basic.Mrot;
        ax = gui.axes_handles.SphereR;
end
% norm of magnetization vector
Mamp = sqrt(M(:,1).^2+M(:,2).^2+M(:,3).^2);

% for visualization purposes M always gets normalized on Bloch sphere
M = M./Mamp(1);

% chhose what to plot depending on simulation type
switch data.basic.type
    case 'std'
        indS1 = 1:numel(T); % relaxation
        indS2 = []; % relaxation at end of "prepolpulse"-type
        indR = []; % switch-off ramp
        indB = []; % switch-off B-field
        indP = []; % pulse
        indE = numel(T); % end point
        Beffn = [0 0 0];
        
    case 'prepol'
        indS1 = T>Tramp;
        indS2 = [];
        indR = T<=Tramp;
        indB = [];
        indP = [];
        indE = numel(T);
        % switch-off B-field trajectory
        Bamp = sqrt(data.results.prepol.Beff(:,1).^2+...
            data.results.prepol.Beff(:,2).^2+...
            data.results.prepol.Beff(:,3).^2);
        Beffn = data.results.prepol.Beff./Bamp;
        
    case 'pulse'
        indS1 = T>Ttau;
        indS2 = [];
        indR = [];
        indB = [];
        indP = T<=Ttau;
        indE = numel(T);
        Beffn = [0 0 0];
        
    case 'prepolpulse'
        indS2 = T>Tramp & T<=Tramp+Twait;
        indS1 = T>Tramp+Twait+Ttau;
        indR = T<=Tramp;
        indB = [];
        indP = T>Tramp+Twait & T<=Tramp+Twait+Ttau;
        indE = numel(T);
        % switch-off B-field trajectory
        Bamp = sqrt(data.results.prepol.Beff(:,1).^2+...
            data.results.prepol.Beff(:,2).^2+...
            data.results.prepol.Beff(:,3).^2);
        Beffn = data.results.prepol.Beff./Bamp;
end

% plot data
clearSingleAxis(ax);
hold(ax,'on');
plot3(M(indS1,1),M(indS1,2),M(indS1,3),'LineWidth',myui.linewidth,...
    'Color',myui.color.basic,'Parent',ax);
plot3(M(indS2,1),M(indS2,2),M(indS2,3),'LineWidth',myui.linewidth,...
    'Color',myui.color.wait,'Parent',ax);
plot3(M(indP,1),M(indP,2),M(indP,3),'LineWidth',myui.linewidth,...
    'Color',myui.color.pulse,'Parent',ax);
plot3(M(indR,1),M(indR,2),M(indR,3),'LineWidth',myui.linewidth,...
    'Color',myui.color.prepol,'Parent',ax);
plot3(Beffn(:,1),Beffn(:,2),Beffn(:,3),'LineWidth',myui.linewidth,...
    'Color',myui.color.prepolB,'Parent',ax)
plot3(M(indE,1),M(indE,2),M(indE,3),'ko','MarkerSize',8,'Parent',ax);
% plot actual Bloch sphere
plotBSphere(18,18,ax);
% axis settings
view(ax,[-35 30]);
hold(ax,'off');
set(ax,'XLim',[-1.05 1.05],'YLim',[-1.05 1.05],'ZLim',[-1.05 1.05])
set(ax,'Color','w','XColor','none','YColor','none','ZColor','none');
axis(ax,'equal');
axis(ax,'tight');

end

%% FFT of magnetization and pulse
function plotFFT(data,gui,type)

% Larmor freq. [Hz]
fL = getOmega0(data.basic.gamma,data.basic.B0)/2/pi;

% choose what to plot
switch type
    case 'M'
        ax = gui.axes_handles.MagFFT;
        f = data.results.basic.Mspec.fx;
        X = data.results.basic.Mspec.X;
        lgdstr = {'Mxy','\omega_0'};
    case 'B'
        ax = gui.axes_handles.PulseFFT;
        f = data.results.pulse.Bspec.fx;
        X = data.results.pulse.Bspec.X;
        lgdstr = {'B','\omega_0'};
end

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
legend(ax,lgdstr,'Location','NorthEast');
% font size
set(ax,'FontSize',gui.myui.axfontsize);

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
