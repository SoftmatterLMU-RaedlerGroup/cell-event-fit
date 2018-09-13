%% postproc_psiva.m
% This script contains the postprocessing routine for pSIVA signal.
%
% Required variables:
% ===================
% data					array of measured data
% t						time vector of measured data (in hours)
% data_sim				array of simulated data
% t_sim					time vector of simulated data (in hours)
% t_sim_res				time resolution of time vector for simulated data
% max_val				vector of maximum values of all traces
% max_ind				vector of indices of maxima in all traces
% min_val				vector of minimum values of all traces
% min_ind				vector of indices of minima in all traces
% amplitude				vector of amplitudes of all traces
% model.simulate		simulation function (fitting function)
% outmode.debug			boolean whether debug figure should be created
% outmode.states		boolean whether state values should be created
% model.marker			marker model name
% t_unit				time unit (for plot)
% target_dir			directory in which to save data
% filename				name of the initial file
% format_cell_number	formatted string of cell number
% time_now				formatted string of current time
% 
% Provided variables:
% ===================
% stat_mat				status matrix; contains onset times and amplitudes
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

%% Allocate memory
stat_mat = zeros(size(data,2), 4);
die_level = NaN;

%% Define thresholds
% Find global maximum (exclude misfits by comparison with median)
[amps_sort,amps_sort_ind] = sort(amplitude, 'descend');
amp_median = amps_sort(floor( numel(amps_sort) / 2 ));

for i = 1:numel(amps_sort)
	amp_max = amps_sort(i);

	if amps_sort(i) < 10 * amp_median
		break
	else
		% Probably misfit; prevent it from being evaluated
		amplitude(amps_sort_ind(i)) = -Inf;
	end
end

% Define threshold values
noise = amp_max * .05;
onset_level = .05;
assume_const_delta = 1E-2;

%% Populate status matrix
for c = 1:size(data,2)

	% Calculate onset time
	if amplitude(c) < noise
		% Cell does not die at all
		die_level = NaN;
		t_onset = NaN;
	else
		% Calculate temporal derivative
		start_deriv = data_sim(2,c) - data_sim(1,c);
		end_deriv = data_sim(length(t_sim),c) - data_sim(length(t_sim)-1,c);

		if start_deriv > assume_const_delta
			% Cell died before measurement
			die_level = NaN;
			t_onset = -Inf;
		elseif end_deriv > assume_const_delta
			% Cell died to late; reaches maximum after measurement
			die_level = NaN;
			t_onset = Inf;
		else
			% Cell dies within observed time
			if min_ind(c) < max_ind(c)
				% Minimum lies before maximum
				die_level = min_val(c) + onset_level * amplitude(c);
				die_delta = amplitude(c);
			else
				% Minimum comes after maximum; find peak
				[peak_height,peak_loc,~,peak_prom] = findpeaks(data_sim(:,c), 'SortStr', 'descend');

				if numel(peak_prom) == 0 || peak_prom(1) < noise
					% No significant peak found
					die_level = NaN;
					if start_deriv < 0
						t_onset = -Inf;
					else
						t_onset = NaN;
					end
				else
					% We seem to have a significant peak
					die_level = peak_height(1) - (1 - onset_level) * peak_prom(1);
					die_delta = amplitude(c) + 1;
					max_ind(c) = peak_loc(1);
					max_val(c) = peak_height(1);
				end
			end

			% Find the time where signal is nearest to die_level
			if isfinite(die_level)
				% Initialize variables
				late_die_value = max_val(c);
				early_die_index = NaN;
				late_die_index = NaN;

				% Search value nearest to die_level, starting at maximum
				for x = max_ind(c):-1:1
					die_delta_temp = abs( data_sim(x,c) - die_level );

					if die_delta_temp < die_delta
						% Got nearer to die_level than before
						die_delta = die_delta_temp;
						early_die_index = x;

						if data_sim(x,c) <= late_die_value
							% Signal decreases monotonically
							late_die_value = data_sim(x,c);
							late_die_index = x;
						end
					elseif die_delta_temp > die_delta
						% Going away from die_level
						break
					end
				end

				if early_die_index == late_die_index
					% Had only monotonic descent
					t_onset = t_sim(late_die_index);
				elseif isfinite(late_die_index)
					% Signal increased again
					time_step = t_sim_res / 2;
					t_temp = t_sim(late_die_index) + time_step;
					die_delta_temp = abs(model.simulate(t_temp, params(c,:)) - die_level);

					while die_delta_temp >= die_delta
						time_step = time_step / 2;

						t_temp_plus = t_temp + time_step;
						t_temp_minus = t_temp - time_step;

						die_delta_temp_plus = abs(model.simulate(t_temp_plus, params(c,:)) - die_level);
						die_delta_temp_minus = abs(model.simulate(t_temp_minus, params(c,:)) - die_level);

						if die_delta_temp_plus < die_delta_temp_minus
							die_delta_temp = die_delta_temp_plus;
							t_temp = t_temp_plus;
						else
							die_delta_temp = die_delta_temp_minus;
							t_temp = t_temp_minus;
						end
					end

					t_onset = t_temp;
				end
			end
		end
	end

	% Write status matrix row:
	%	 1st column: cell number
	%	 2nd column: onset time
	%	 3rd column: absolute amplitude
	%	 4th column: relative amplitude
	stat_mat(c,:) = [ c amplitude(c) amplitude(c)/amp_max t_onset ];

	%% Plot debug figures
	if outmode.debug
		fh = figure;
		p = plot(t_sim,data_sim(:,c),'r-', t,data(:,c),'b-');
			% Minimum amplitude level:
			line(xlim, [min_val(c)+noise min_val(c)+noise], 'LineStyle', '-.', 'Color', [1,.5,0]);
			% Onset level:
			line(xlim, [die_level die_level], 'LineStyle', '-', 'Color', [1,.5,0]);
			% Onset time:
			line([t_onset t_onset], ylim, 'LineStyle', '-', 'Color', [1,.5,0]);
		uistack(p, 'top');
		title([model.marker ' Fluorescence in Cell ' num2str(c)]);
		xlabel(['Time ' t_unit]); xlim([0,max(t)]);
		ylabel('Fluorescence Intensity [a.u.]');
		set(fh, 'PaperUnits','centimeters', 'PaperSize',[10 7], 'PaperPosition',[0 0 10 7])
		print(fh, '-dpdf', fullfile(target_dir, [filename '_DEBUG_' num2str(c, format_cell_number) '_' time_now '.pdf']));
		close(fh);
	end
end

%% Export status matrix
if outmode.states
	csvwrite(fullfile(target_dir, [filename '_ALL_STATE_' time_now '.csv']), stat_mat);
end