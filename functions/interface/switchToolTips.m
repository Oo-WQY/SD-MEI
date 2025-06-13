function switchToolTips(gui,onoff)
%switchToolTips switches GUI tool tips either "on" or "off"
%
% Syntax:
%       switchToolTips(gui,onoff)
%
% Inputs:
%       gui - figure gui elements structure
%       onoff - 'on' or 'off'
%
% Outputs:
%       none
%
% Example:
%       switchToolTips(gui,'on')
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

%% list of handles that have switchable tool tips
h = {'check_handles','edit_handles','popup_handles',...
    'push_handles'};

%% process all handles
switch lower(onoff)
    case 'on' % switch tool tips on
        for i = 1:numel(h)
            if isfield(gui,h{i})
                eval(['fnames = fieldnames(gui.',h{i},');']);
                for j = 1:numel(fnames)
                    eval(['ud = get(gui.',h{i},'.',fnames{j},...
                        ',''UserData'');']);
                    if isfield(ud,'Tooltipstr')
                        tstr = ud.Tooltipstr;
                        eval(['set(gui.',h{i},'.',fnames{j},...
                            ',''ToolTipString'',tstr);']);
                    end
                end
            end
        end
        
    case 'off' % switch tool tips off
        for i = 1:numel(h)
            if isfield(gui,h{i})
                eval(['fnames = fieldnames(gui.',h{i},');']);
                for j = 1:numel(fnames)
                    eval(['ud = get(gui.',h{i},'.',fnames{j},...
                        ',''UserData'');']);
                    if isfield(ud,'Tooltipstr')
                        tstr = ud.Tooltipstr;
                        eval(['set(gui.',h{i},'.',fnames{j},...
                            ',''ToolTipString'','''');']);
                    end
                end
            end
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
