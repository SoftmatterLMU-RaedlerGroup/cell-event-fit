function varargout = postproc_mRNA_trivial(~, F, R)
% This function contains the postprocessing routine for the
% trivial model of mRNA expression.
%
%
% Required arguments:
% ===================
% D					structure, containing these fields:
% D.data					array of measured data
% D.index					index of current trace in datafile
% R					structure, containing these fields:
% R.data_sim				array of simulated data
% R.max_val					maximum value
% R.max_ind					index of maximum
% R.min_val					minimum value
% R.min_ind					index of minimum
% R.amplitude				amplitude (R.max_val - R.min_val)
% R.amp_max					file-wide maximum value of amplitude
% R.params					vector of fit parameters
% F					structure, containing these fields:
% F.t_sim					time vector of simulated data (in hours)
% F.fullpath				absolute path to the initial file
%
% Return value:
% =============
% varargout{1}			onset time of the signal
% varargout{2}			structure with model-specific information
% varargout{3}			derivative of the rising edge
%
% Copyright © 2018 Daniel Woschée <daniel.woschee@physik.lmu.de>
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

%% Definitions
% Define threshold values
assume_const_delta = 1E-2;

% Allocate memory
middle_amplitude = R.min_val + R.amplitude / 2;
mid_amp_ind = NaN;

%% Get 1st derivative of flourescence signal w.r.t. time
deriv = NaN(1,size(R.data_sim,1));

for i = 2:length(R.data_sim)
	deriv(i) = ( R.data_sim(i) - R.data_sim(i-1) ) ...
		/ ( F.t_sim(i) - F.t_sim(i-1) );
end
deriv(1) = deriv(2);

% Get maximum of derivative of rising edge
rising_deriv = deriv(2);
rising_ind = 2;

for i = 3:R.max_ind
	if deriv(i) > rising_deriv
		rising_deriv = deriv(i);
		rising_ind = i;
	elseif deriv(i) < 0
		break
	end
end

% Define rising edge coefficients
rising_edge_xoff = R.data_sim(rising_ind);
rising_edge_yoff = F.t_sim(rising_ind);
rising_edge_scale = 1 / deriv(rising_ind);

% Find index next to middle amplitude in rising edge
min_diff = Inf;
for i = R.max_ind:-1:1
	temp_diff = abs(R.data_sim(i) - middle_amplitude);

	if temp_diff < min_diff
		mid_amp_ind = i;
		min_diff = temp_diff;
	elseif isfinite(min_diff) && temp_diff > min_diff
		break
	end
end

%% Determine onset time
if max(deriv) <= assume_const_delta
	t_onset = NaN;

elseif deriv(1) < -assume_const_delta || R.max_ind < R.min_ind
	t_onset = -Inf;

else
	% Define coefficients for rising edge
	mid_amp_deriv = ...
		( F.t_sim(mid_amp_ind + 1) - F.t_sim(mid_amp_ind - 1) ) ...
		/ ( R.data_sim( mid_amp_ind + 1 ) - R.data_sim( mid_amp_ind - 1 ) );

	rising_edge_xoff = middle_amplitude;
	rising_edge_yoff = F.t_sim(mid_amp_ind);
	rising_edge_scale = mid_amp_deriv;

	ascent_time = @(Y) ( Y - rising_edge_xoff ) * rising_edge_scale  ...
		+ rising_edge_yoff;

	% Calculate onset time
	t_onset = ascent_time(R.data_sim(1));
end

%% Create output
varargout{1} = t_onset;

varargout{2}.type = int8(3);	% mRNA transfection
varargout{2}.rising_edge_xoff = rising_edge_xoff;
varargout{2}.rising_edge_yoff = rising_edge_yoff;
varargout{2}.rising_edge_scale = rising_edge_scale;

% Main derivative
if nargout > 2
	if isfinite(t_onset)
		varargout{3} = 1 / mid_amp_deriv;
	else
		varargout{3} = 1 / deriv(rising_ind);
	end
end