function bsh = plotBSphere(dlat,dlong,ax,varargin)
%plotBSphere plots a (Bloch)-Sphere (basically a unit sphere) with
%increments 'dlat' [deg] and 'dlong' [deg] for latitude and longitude,
%respectively
%
% Syntax:
%       plotBSphere
%
% Inputs:
%       dlat - latitude increment [deg]
%       dlong - longitude increment [deg]
%       ax - axes handle
%		varargin - optional radius R
%
% Outputs:
%       none
%
% Example:
%       plotBSphere(30,30)
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

% radius of the sphere
R = 1;
if nargin > 3
    R = varargin{1};
end
% gray scale factor
sf = 0.85;

% lines along longitude
[lon1,lat1] = meshgrid(-180:dlong:180,linspace(-90,90,181));
% lines along latitude:
[lat2,lon2] = meshgrid(-90:dlat:90,linspace(-180,180,361));

% spherical to Cartesian coordinate transform
[x1,y1,z1] = sph2cart(deg2rad(lon1),deg2rad(lat1),R);
[x2,y2,z2] = sph2cart(deg2rad(lon2),deg2rad(lat2),R);

% plotting the lines in the current axes
bsh1 = plot3(x1,y1,z1,'-','Color',sf*[1 1 1],'LineWidth',1,'Parent',ax);
bsh2 = plot3(x2,y2,z2,'-','Color',sf*[1 1 1],'LineWidth',1,'Parent',ax);

bsh3 = line([-R R],[0 0],[0 0],'Color',sf*[1 1 1],'LineWidth',1,'Parent',ax);
bsh4 = line([0 0],[-R R],[0 0],'Color',sf*[1 1 1],'LineWidth',1,'Parent',ax);
bsh5 = line([0 0],[0 0],[-R R],'Color',sf*[1 1 1],'LineWidth',1,'Parent',ax);

bsh6 = line([0 R*1.2],[0 0],[0 0],'Color','r','LineWidth',1.3,'Parent',ax);
bsh7 = line([0 0],[0 R*1.2],[0 0],'Color','g','LineWidth',1.3,'Parent',ax);
bsh8 = line([0 0],[0 0],[0 R*1.2],'Color','b','LineWidth',1.3,'Parent',ax);

t1 = text(R*1.4,0,0,'X','HorizontalAlignment','center','Color','r','Parent',ax);
t2 = text(0,R*1.4,0,'Y','HorizontalAlignment','center','Color','g','Parent',ax);
t3 = text(0,0,R*1.4,'Z','HorizontalAlignment','center','Color','b','Parent',ax);

% output handles of all lines and text
bsh.grid = [bsh1; bsh2];
bsh.axes1 = [bsh3; bsh4; bsh5];
bsh.axes2 = [bsh6; bsh7; bsh8];
bsh.label = [t1; t2; t3];

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
