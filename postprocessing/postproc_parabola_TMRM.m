function varargout = postproc_parabola_TMRM(~, ~, R)
% This function contains the postprocessing routine for TMRM_interactive.
%
% Required arguments:
% ===================
% %D					structure, containing these fields:
% %D.data					array of measured data
% %D.index					index of current trace in datafile
% %F					structure, containing these fields:
% %F.t_sim					time vector of simulated data (in hours)
% %F.modelname				name of the model used
% %F.simulate				simulation function of the model
% R					structure, containing these fields:
% R.data_sim				array of simulated data
% R.max_val					maximum value
% R.max_ind					index of maximum
% R.min_val					minimum value
% R.min_ind					index of minimum
% R.amplitude				amplitude (R.max_val - R.min_val)
% R.params					vector of fit parameters
%
% Return value:
% =============
% varargout{1}			breakdown time of the signal
% varargout{2}			structure with model-specific information
% varargout{3}			derivative at breakdown
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

%% Define thresholds
min_decay_steepness = 0.1;
min_limit_dist = .25;
min_breakdown_amp = .05 * R.amplitude;
min_deriv_change = 3;

%% Test breakdown time
a_parab = R.params(1);
x_parab = R.params(2);
t_breakdown = R.params(4);
a_break = R.params(5);
y_end = R.params(6);
t_start = R.params(7);
t_end = R.params(8);
y_break = parabola_TMRM_simulate(t_breakdown, R.params);


breakdown_deriv = -a_break * (y_break - y_end);
parab_deriv = 2 * a_parab * (t_breakdown - x_parab);
deriv_change = parab_deriv - breakdown_deriv;

if a_break < min_decay_steepness
	% Decay must be steep enough
	t_breakdown = NaN;
	
elseif t_breakdown < t_start + min_limit_dist || t_breakdown > t_end - min_limit_dist
	% Breakdown must be within limit
	t_breakdown = NaN;

elseif y_break - y_end < min_breakdown_amp
	% Breakdown must be high enough
	t_breakdown = NaN;
	
elseif deriv_change < min_deriv_change
	% Too small change in derivative at breakdown
	t_breakdown = Inf;
end

%% Calculate breakdown derivative
if ~isfinite(t_breakdown)
	breakdown_deriv = NaN;
end

%% Return values
% Breakdown time
varargout{1} = t_breakdown;

% Fit type
if nargout > 1
	varargout{2} = int8(5); % type: parabola_TMRM
end

% Linear parameter of parabola in normal form
if nargout > 2
	varargout{3} = breakdown_deriv;
end

% Custom data
if nargout > 3
	varargout{4} = 1;
	varargout{5} = deriv_change;
end
