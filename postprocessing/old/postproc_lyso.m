function stat_mat = postproc_lyso(D, F, output_info)
% This function contains the postprocessing routine for LysoTracker.
%
% Required arguments:
% ===================
% D					structure, containing these fields:
% D.data					array of measured data
% D.data_sim				array of simulated data
% D.max_val					maximum value
% D.max_ind					index of maximum
% D.min_val					minimum value
% D.min_ind					index of minimum
% D.amplitude				amplitude (D.max_val - D.min_val)
% D.amp_max					file-wide maximum value of amplitude
% D.amp_max_restricted		D.amp_max, cleaned for runaway values
% D.					vector of fit parameters
% D.logPost					log-likelihood for parameter set D.params
% D.index					index of current trace in datafile
% F					structure, containing these fields:
% F.t						time vector of measured data (in hours)
% F.t_sim					time vector of simulated data (in hours)
% F.t_sim_res				resolution (stepsize) of F.t_sim
% F.t_unit					time unit (for plot)
% F.target_dir				directory in which to save data
% F.name					name of the initial file
% F.format_cell_number		formatted string of cell number
% F.modelname				name of the model used
% output_info		structure, containing these fields:
% output_info.debug			boolean whether debug figure should be created
% output_info.marker		marker model name
% output_info.time_now		formatted string of current time
%
% Return value:
% =============
% stat_mat			status (line) vector; contains onset times and amplitudes
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

%% Define thresholds
noise = D.amp_max * .1;
assume_const_delta = 1E-2;
parab_threshold = 0.01;

%% Set default values
parab_threshold = parab_threshold * D.amplitude;
deriv_min_ind = NaN;
t_breakdown = NaN;
switch_to_extrapolate = 0;
parabola = zeros(length(F.t_sim), 1);
parabola(:) = NaN;
fit_type = 0;
start_deriv = NaN;

%% Evaluate differences to parabola

% Get parabola coefficients (index depends on model)
switch F.modelname
	case 'Lyso1'
		vertex_pos = D.params(3);
		coeff_const = D.params(4);
		coeff_lin = D.params(5);
		coeff_quad = - D.params(6);

	case 'Lyso2'
		vertex_pos = D.params(1);
		coeff_const = D.params(2);
		coeff_lin = D.params(3);
		coeff_quad = - D.params(4);

	case 'Lyso3'
		vertex_pos = D.params(1);
		coeff_const = D.params(2);
		coeff_lin = 0;
		coeff_quad = - D.params(3);

	case 'CellROX'
		vertex_pos = D.params(6);
		coeff_const = D.params(1) + D.params(2) * D.params(7);
		coeff_lin = 0;
		coeff_quad = - D.params(2) * D.params(5);

	otherwise
		error(['Unknown model for ' output_info.marker]);
end

% Calculate parabola
for i = 1:length(F.t_sim)
	x = F.t_sim(i) - vertex_pos;
	parabola(i) = coeff_const + x * ( coeff_lin + x * coeff_quad );

	if parabola(i) < D.min_val - parab_threshold
		parabola(i) = NaN;
	end
end

%% Sort out badly-behaved cells; calculate breakdown time
if D.amplitude < noise
	% Cell does not die at all
	t_breakdown = NaN;

else
	% Get 1st derivative of flourescence signal w.r.t. time
	deriv = zeros(1,size(D.data_sim,1));
	deriv(1) = NaN;

	for i = 2:size(D.data_sim,1)
		deriv(i) = ( D.data_sim(i) - D.data_sim(i-1) ) ...
			/ ( F.t_sim(i) - F.t_sim(i-1) );
	end

	% Find slope minimum
	[~,deriv_min_ind] = min(deriv(D.max_ind:end));
	deriv_min_ind = deriv_min_ind + D.max_ind - 1;
		
	start_deriv = deriv(2);

	% Sort out cells dying too early or too late
	if start_deriv <= assume_const_delta
		% Cell died before measurement
		t_breakdown = -Inf;

	elseif D.data_sim(end) == D.max_val || deriv_min_ind == length(D.data_sim)
		% Cell died to late
		%disp(['Cell ' num2str(D.index) ': end_deriv=' num2str(end_deriv) ', end_val=' num2str(D.data_sim(end))]);
		t_breakdown = Inf;

	else
		% Cell dies within observed time

		% Locate kink by deviation from parabola
		for i = 1:length(F.t_sim)
			parab_devi = abs(D.data_sim(i) - parabola(i));

			if isnan(parab_devi) || parab_devi <= parab_threshold
				% Nothing interesting happens
				continue;

			elseif deriv(i) >= 0 && i < D.max_ind
				% Too flat breakdown, switch to extrapolation mode
				switch_to_extrapolate = 1;
				break;

			else
				if i > D.max_ind && D.data_sim(i-1) <= D.min_val + noise
					% Signal is too low, still no breakdown, switch
					% to extrapolation mode
					switch_to_extrapolate = 1;
					break;
				elseif i > 1
					% Breakdown time found
					t_breakdown = F.t_sim(i);
					fit_type = 1;
					break;
				else
					continue;
				end

			end
		end

		if switch_to_extrapolate
			% Get normalized 1st derivative

			% Set up extrapolation function
			breakdowntime = @(y) (y - D.data_sim(deriv_min_ind)) / deriv(deriv_min_ind) ...
				+ F.t_sim(deriv_min_ind);

			% Calculate onset time by extrapolating from slope minimum
			t_breakdown = breakdowntime(D.max_val);
			fit_type = -1;
		end
	end
end

%% Write status matrix row
%	1st column: cell number
%	2nd column: breakdown time
%	3rd column: absolute amplitude
%	4th column: relative amplitude
%	5th column: log-likelihood
%	6th column: fit_type (1: kink, 0: no fit, -1: extrapolation)
%	7th column: initial slope
stat_mat = [ D.index t_breakdown D.amplitude D.amplitude/D.amp_max D.logPost fit_type start_deriv ];

%% Plot debug figures
if output_info.debug
	fh = figure;
	plot(F.t_sim,D.data_sim,'r-', F.t,D.data,'b-');

	% Minimum amplitude level:
	line(xlim, [D.min_val+noise D.min_val+noise], 'LineStyle','-.', 'Color',[1,.5,0]);

	if isfinite(t_breakdown)
		hold on;

		if switch_to_extrapolate
			% Extrapolation:
			line(breakdowntime(ylim), ylim, 'LineStyle','-', 'Color',[1,.5,0]);
			% Minimum of slope:
			line([F.t_sim(deriv_min_ind), F.t_sim(deriv_min_ind)], ylim, 'LineStyle',':', 'Color',[1,.5,0]);
		else
			plot(F.t_sim,parabola,'m-');
		end

		% Onset time:
		line([t_breakdown t_breakdown], ylim, 'LineStyle','-', 'Color',[1,.5,0]);
	end

	title([output_info.marker ' Fluorescence in Cell ' num2str(D.index)]);
	xlabel(['Time ' F.t_unit]); xlim([0,max(F.t)]);
	ylabel('Fluorescence Intensity [a.u.]');
	set(fh, 'PaperUnits','centimeters', 'PaperSize',[10 7], 'PaperPosition',[0 0 10 7])
	print(fh, '-dpdf', fullfile(F.target_dir, [F.name '_DEBUG_' num2str(D.index, F.format_cell_number) '_' output_info.time_now '.pdf']));
	close(fh);
end
