function P = parabola_TMRM_parameters(data, t, model)
%PARABOLA_TMRM_PARAMETERS defined parameter properties for the parabola TMRM model
%
% Copyright © 2019 Daniel Woschée <daniel.woschee@physik.lmu.de>
% Faculty of Physics / Ludwig-Maximilians-Universität München
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License version 2
% as published by the Free Software Foundation.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.

% Get data properties
t_min = min(t);
t_max = max(t);
data_min = min(data);
data_max = max(data);
data_amp = data_max - data_min;

% Define parameter minima
par_min = model.par_min;
par_min(2) = t_min - 10;
par_min([3, 6]) = data_min - .5 * data_amp;
par_min(4) = t_min;
par_min(7) = t_min;
par_min(8) = t_max;

% Define parameter maxima
par_max = model.par_max;
par_max(2) = t_max + 10;
par_max([3, 6]) = data_max + .5 * data_amp;
par_max(4) = t_max;
par_max(7) = t_min;
par_max(8) = t_max;

% Return parameter limits
P = struct( ...
	'min', par_min, ...
	'max', par_max ...
	);
end
