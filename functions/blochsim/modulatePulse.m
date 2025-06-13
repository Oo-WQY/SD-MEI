function out = modulatePulse(mod,type)
%modulatePulse generates the frequency and current modulation functions
%
% Syntax:
%       modulatePulse(mod,type)
%
% Inputs:
%       mod - struct containing modulation settings
%             shape : shape of the modulation function [string]
%             t     : time (scalar / vector) [s]
%             t0    : pulse start time [s]
%             t1    : pulse end time [s]
%             v0    : start value 
%             v1    : end value
%             A     : modulation parameter (MIDI)
%             B     : modulation parameter (MIDI)
%       type - switch for frequency 'df' or current 'I' modulation [string]
%
% Outputs:
%       out - modulated signal (either 'df' or 'I' as scalar / vector)
%
% Example:
%       modulatePulse(fmod,'df')
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

% general parameters
t = mod.t(:);
% t0 = mod.t0;
% t1 = mod.t1;
% v0 = mod.v0;
v1 = mod.v1;
% A = mod.A; 
% B = mod.B;
% 
% % pulse duration [s]
% tau = t1-t0;
% % modulation range
% delta = v1-v0;

switch mod.shape
    
    case 'const' % no modulation
        out = v1.*ones(size(t));
        
    case 'lin' % linear modulation
        out = v0 + (t.*delta./tau);
        
    case 'tanhMIDI' % tanh-modulation with two parameters A & B
        % slope parameter
        N = tanh( ((2.*pi.*A)./tau).*(t - B.*(tau/2)) );
        N0 = tanh( ((2.*pi.*A)./tau).*(t0 - B.*(tau/2)) );
        N1 = tanh( ((2.*pi.*A)./tau).*(tau - B.*(tau/2)) );
        
        switch type
            case 'df'
                % sign switch (MMP ;-))
                delta = -delta;
                out = v1 + delta*(1-((N-N0)/(N1-N0)));
            case 'I'
                out = v0 + delta*(N-N0)/(N1-N0);
        end
        
    case 'tanhGMR' % tanh-modulation GMR style
        switch type
            case 'df'
                tau = 3*t./tau; % (RD: pers. comm. Grunewald 13.10.2016)
            case 'I'
                tau = pi*t./tau; % (RD: pi is arbitrary.)
        end
        out = v0 + (delta * tanh(tau));
        
    case 'exp' % exponential modulation
        out = v1 - delta .* exp(A.*(-t./tau));
        
    case 'custom' % custom modulation (interpolation)
        switch type
            case 'df'
                out = interp1(mod.custom_t,mod.custom_df,t,'linear');
            case 'I'
                out = interp1(mod.custom_t,mod.custom_I,t,'linear');
        end
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
