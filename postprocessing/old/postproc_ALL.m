function stat_mat = postproc_ALL(D, F, output_info)
% This function contains the postprocessing routine for PI (depracated).
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
% D.params					vector of fit parameters
% D.logPost					log-likelihood for parameter set D.params
% D.index					index of current trace in datafile
% F					structure, containing these fields:
% F.t						time vector of measured data (in hours)
% F.t_sim					time vector of simulated data (in hours)
% F.t_unit					time unit (for plot)
% F.target_dir				directory in which to save data
% F.name					name of the initial file
% F.format_cell_number		formatted string of cell number
% output_info		structure, containing these fields:
% output_info.debug			boolean whether debug figure should be created
% output_info.marker		marker model name
% output_info.time_now		formatted string of current time
%
% Return value:
% =============
% stat_mat		status (line) vector; contains onset times and amplitudes
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
noise = D.amp_max_restricted * .05;
onset_level = .05;
assume_const_delta = 1E-2;

%% Populate status matrix
% Calculate onset time
if D.amplitude < noise || D.amplitude > D.amp_max_restricted
	% Cell does not die at all / signal too high
	die_level = NaN;
	t_onset = NaN;
else
	% Calculate temporal derivative
	start_deriv = D.data_sim(2) - D.data_sim(1);
	end_deriv = D.data_sim(length(F.t_sim)) - D.data_sim(length(F.t_sim)-1);

	if D.max_ind <= D.min_ind || start_deriv < 0 || start_deriv > assume_const_delta
		% Cell died before measurement
		die_level = NaN;
		t_onset = -Inf;
	elseif end_deriv > assume_const_delta
		% Cell died to late; reaches maximum after measurement
		die_level = NaN;
		t_onset = Inf;
	else
		% Cell dies within observed time
		die_level = D.min_val + onset_level * D.amplitude;
		die_delta = D.amplitude;

		for x = 1:D.max_ind
			die_delta_temp = abs( D.data_sim(x) - die_level );

			if die_delta_temp < die_delta && x < D.max_ind
				die_delta = die_delta_temp;
				t_onset = F.t_sim(x);
			end
		end
	end
end

% Write status matrix row:
%	1st column: cell number
%	2nd column: onset time
%	3rd column: absolute amplitude
%	4th column: relative amplitude
%	5th column: log-likelihood
stat_mat = [ D.index t_onset D.amplitude D.amplitude/D.amp_max_restricted D.logPost ];

%% Plot debug figures
if output_info.debug
	fh = figure;
	p = plot(F.t_sim,D.data_sim,'r-', F.t,D.data,'b-');
		% Minimum amplitude level:
		line(xlim, [D.min_val+noise D.min_val+noise], 'LineStyle', '-.', 'Color', [1,.5,0]);
		% Onset level:
		line(xlim, [die_level die_level], 'LineStyle', '-', 'Color', [1,.5,0]);
		% Onset time:
		line([t_onset t_onset], ylim, 'LineStyle', '-', 'Color', [1,.5,0]);
	uistack(p, 'top');
	title([output_info.marker ' Fluorescence in Cell ' num2str(D.index)]);
	xlabel(['Time ' F.t_unit]); xlim([0,max(F.t)]);
	ylabel('Fluorescence Intensity [a.u.]');
	set(fh, 'PaperUnits','centimeters', 'PaperSize',[10 7], 'PaperPosition',[0 0 10 7])
	print(fh, '-dpdf', fullfile(F.target_dir, [F.name '_DEBUG_' num2str(D.index, F.format_cell_number) '_' output_info.time_now '.pdf']));
	close(fh);
end

