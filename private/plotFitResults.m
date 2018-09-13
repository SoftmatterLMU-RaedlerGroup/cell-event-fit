function [] = plotFitResults(Ffile, loop_start)
%plotFitResults plots the results found in the mat-file.
%
% Input parameters
% ================
% Ffile			String containing the path of the results mat-file
% loop_start	(optional) Vector with:
%					loop_start(1)	first trace index to be used
%					loop_start(2)	first file index to be used
%					Set to 1 to loop over all traces/files.
%					Set to NaN to skip corresponding loop.
%					`loop_start(2)` may be omitted to export all traces.
%
% If `loop_start` is not given, a plot log file `plot.log` in the temporary
% directroy defined in `Ffile` is queried. If the log file is found,
% exporting is continued where the log file indicates. If no log file is
% found and the output directories are empty, all traces and files are
% exported.
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

%% Open matfile
mf = matfile(Ffile);
ntraces = mf.ntraces;
ndatafiles = mf.ndatafiles;
outmode = mf.outmode;
plot_log = fullfile(mf.temp_dir, 'plot.log');

%% Parse input
% Check if `loop_start` is defined and valid
is_start_defined = false;
is_start_valid = false;

if nargin > 1 && ~isempty(loop_start)
	is_start_defined = true;
	if isnumeric(loop_start)
		for i_ls = 1:min(2, length(loop_start))
			lsi = loop_start(i_ls);
			if isnan(lsi)
				% NaN is allowed; ignores traces/files
				is_start_valid = true;
			elseif lsi > 0 && floor(lsi) == lsi
				% Check for positive integer (or Inf)
				is_start_valid = true;
			else
				% Other values are not allowed
				is_start_valid = false;
				break
			end
		end
	end
end

% Get starting indices according to `loop_start`:
% If `loop_start` is defined and valid, use these values.
% If `loop_start` is defined but not valid, throw an error.
% If `loop_start` is not defined, check if `plot.log` exists.
% If `plot.log` exists, use these values or throw an error.
% If `plot.log` does not exist, check if there is already content there.
% If content is there, throw an error.
% If there is no content yet, export all traces and files.
first_trace_index = 1;
first_file_index = 1;

if is_start_valid
	if length(loop_start) >= 1
		if loop_start(1) > 0 && loop_start(1) <= ntraces
			% Loop over traces, beginning at loop_start(1)
			first_trace_index = loop_start(1);

		elseif ~isfinite(loop_start(1))
			% Skip loop over traces
			first_trace_index = inf;
		end
	end

	if length(loop_start) >= 2
		if loop_start(2) > 0 && loop_start(2) <= ndatafiles
			% Loop over datafiles, beginning at loop_start(2)
			first_file_index = loop_start(2);

		elseif ~isfinite(loop_start(2))
			% Skip loop over datafiles
			first_file_index = inf;
		end
	end

elseif ~is_start_defined
	% Try to read start index from plot.log
	if exist(plot_log, 'file')
		is_valid_log = false;
		log = regexpi(fileread(plot_log), ...
			'^(:?plotting)?\s*(?<kind>\w+):?\s*(?<index>\d+)', ...
			'names', 'once');

		if ~isempty(log) && isfield(log, 'kind') && isfield(log, 'index')
			index = str2double(log.index);

			if strcmpi(log.kind, 'trace')
				% `plot.log` specifies trace
				first_trace_index = index;
				is_valid_log = true;

			elseif strcmpi(log.kind, 'file')
				% `plot.log` specifies file
				first_file_index = index;
				is_valid_log = true;
			end
		end

		if ~is_valid_log
			error(['`plot.log` has invalid content. ' ...
				'Please specify `loop_start` to define which data ' ...
				'should be plotted; see help.']);
		end

	else
		% Check if files exist in directories
		is_content_found = false;
		for d = dir(mf.temp_dir)'
			if ~d.isdir
				continue
			elseif strcmp(d.name(1), '.')
				continue
			end
			for d2 = dir(fullfile(d.folder, d.name))'
				if ~strcmp(d2.name(1), '.')
					is_content_found = true;
					break
				end
			end
			if is_content_found
				break
			end
		end

		% Do not overwrite existing content by default
		if is_content_found
			error(['There was content found in "%s". ' ...
				'Please specify `loop_start` to define which data ' ...
				'should be plotted; see help.'], mf.temp_dir);
		end
	end
else
	error('Invalid argument encountered. `loop_start` must be numeric; see help.');
end

% Print info if not plotting all files and traces
if first_trace_index ~= 1 || first_file_index ~= 1
	fprintf('%sExport start values: for trace %.0f, for file %.0f\n', ...
		get_time, first_trace_index, first_file_index);
end

%% Plot the single cell data
for d = first_trace_index:ntraces

	% Write log file
	write_plot_log(sprintf('Plotting trace: %d\n', d));

	% Get data from matfile
	f = mf.file_ind(d,1);
	model = mf.models(mf.model_ind(f,1),1);

	% Ensure that target directory exists
	target_dir = mf.target_dir(mf.target_dir_ind(f,1),1:mf.target_dir_len(f,1));
	if ~exist(target_dir, 'dir')
		mkdir(target_dir);
	end

	%% Print postprocessing and prepare state matrix
	if outmode.debug

		% Make options structure
		info = struct( ...
			'name',		mf.name_F(mf.name_F_ind(f,1),1:mf.name_F_len(f,1)), ...
			'target_dir',	target_dir, ...
			't',		mf.t(1:mf.data_len(f,1), mf.t_ind(f,1)), ...
			't_sim',	mf.t_sim(1:mf.data_sim_len(f,1), mf.t_sim_ind(f,1)), ...
			't_unit',	mf.t_unit, ...
			'data',		mf.data(1:mf.data_len(f,1),d), ...
			'data_sim',	mf.data_sim(1:mf.data_sim_len(f,1), d), ...
			'noise',	mf.min_val(d,1) + mf.noise(f,1), ...
			't_event',	mf.t_event(d,1), ...
			'fit_type',	mf.fit_type(d,1), ...
			'rising_edge_xoff',		mf.rising_edge_xoff(d,1), ...
			'rising_edge_yoff',		mf.rising_edge_yoff(d,1), ...
			'rising_edge_scale',	mf.rising_edge_scale(d,1), ...
			'falling_edge_xoff',	mf.falling_edge_xoff(d,1), ...
			'falling_edge_yoff',	mf.falling_edge_yoff(d,1), ...
			'falling_edge_scale',	mf.falling_edge_scale(d,1), ...
			'index',				mf.index_F(d,1), ...
			'marker',				model.marker, ...
			'format_cell_number',	mf.format_cell_number(mf.format_cell_number_ind(f,1),:), ...
			'time_now',				mf.time_now ...
			);

		if mf.fit_parabola_ind(d,1)
			info.fit_parabola = mf.fit_parabola(1:mf.data_sim_len(f,1), mf.fit_parabola_ind(d,1));
		else
			info.fit_parabola = NaN(mf.data_sim_len(f,1), 1);
		end

		% Print postprocessing results
		plot_postproc(info);
	end

	% Plot the multistart results
	if outmode.ms && mf.fh_MS_ind(d,1) && isgraphics(mf.fh_MS(mf.fh_MS_ind(d,1),1))
		fh_MS = mf.fh_MS(mf.fh_MS_ind(d,1),1);
		set(fh_MS, 'PaperUnits','centimeters', 'PaperSize',[25 15], 'PaperPosition',[0 0 25 15])
		print(fh_MS, '-dpdf', fullfile(target_dir, ...
			strcat(mf.name_F(mf.name_F_ind(f,1),1:mf.name_F_len(f,1)), '_MS_', num2str(mf.index_F(d,1)), ...
			mf.format_cell_number(mf.format_cell_number_ind(f,1),:), '_', mf.time_now, '.pdf')));
		delete(fh_MS);
	end

	% Plot the fitting result
	if outmode.single
		fh = figure('Visible','off');
		plot( mf.t( 1:mf.data_len(f,1), mf.t_ind(f,1) ), ...
				mf.data( 1:mf.data_len(f,1), d ), 'b-', ...
			mf.t_sim( 1:mf.data_sim_len(f,1), mf.t_sim_ind(f,1) ), ...
				mf.data_sim( 1:mf.data_sim_len(f,1), d ), 'r-' ...
			);
		title([model.marker ' Fluorescence in Cell ' num2str(mf.index_F(d,1))]);
		xlabel(['Time ' mf.t_unit]); ...
			xlim([0,max( mf.t( 1:mf.data_len(f,1), mf.t_ind(f,1) ) )]);
		ylabel('Fluorescence Intensity [a.u.]');
		set(fh, 'PaperUnits','centimeters', 'PaperSize',[10 7], 'PaperPosition',[0 0 10 7])
		print(fh, '-dpdf', fullfile(target_dir, ...
			strcat( mf.name_F(mf.name_F_ind(f,1),1:mf.name_F_len(f,1)), '_FIT_', num2str( ...
			mf.index_F(d,1), mf.format_cell_number(mf.format_cell_number_ind(f,1), :) ), ...
			'_', mf.time_now, '.pdf' ) ));
		delete(fh);
	end
end

%% Plot the whole dataset
for f = first_file_index:ndatafiles

	% Write log file
	write_plot_log(sprintf('Plotting file: %d\n', f));

	% Get data from matfile
	model = mf.models(mf.model_ind(f,1),1);

	% Ensure that target directory exists
	target_dir = mf.target_dir(mf.target_dir_ind(f,1),1:mf.target_dir_len(f,1));
	if ~exist(target_dir, 'dir')
		mkdir(target_dir);
	end

	% Write status matrix to file
	if outmode.states
		stat_mat = zeros(mf.size_F(f,1), 7);

		for i = mf.dataindices(f, 1:mf.size_F(f,1))
			j = mf.index_F(i,1);
			%% Write status matrix row:
			%	1st column: cell number
			%	2nd column: event time
			%	3rd column: absolute amplitude
			%	4th column: relative amplitude
			%	5th column: log-likelihood
			%	6th column: fit_type (one of the following)
			%				1:	extrapol_rising
			%				2:	kink_parabola
			%				-2:	nokink_parabola
			%				3: mRNA transfection
			%				4: min between peaks
			%				other: currently not defined
			%	7th column: event slope
			stat_mat(j,1) = mf.index_F(i,1);
			stat_mat(j,2) = mf.t_event(i,1);
			stat_mat(j,3) = mf.amplitude(i,1);
			stat_mat(j,4) = mf.amplitude(i,1) / mf.amp_max(f,1);
			stat_mat(j,5) = mf.logPost(i,1);
			stat_mat(j,6) = mf.fit_type(i,1);
			stat_mat(j,7) = mf.event_deriv(i,1);
		end

		csvwrite( fullfile(target_dir, ...
			strcat(mf.name_F( mf.name_F_ind(f,1),1:mf.name_F_len(f,1)), '_ALL_STATE_', mf.time_now, '.csv' ) ), ...
			stat_mat );
	end

	% Calculate parameters and simulated values
	if outmode.params

		% Allocate parameter matrix
		par_mat = NaN(mf.size_F(f,1), model.par_num, 'single');

		% Populate parameter matrix from matfile
		for i = 1:mf.size_F(f,1)
			par_mat(i,:) = mf.params(mf.dataindices(f,i), 1:model.par_num);
		end

		% Export the parameter matrix
		csvwrite(fullfile(target_dir, ...
			strcat(mf.name_F(mf.name_F_ind(f,1),1:mf.name_F_len(f,1)), '_ALL_PARAMS_', mf.time_now, '.csv')), par_mat);
		clear par_mat;
	end

	% Write simulated values
	if outmode.simulated
		csvwrite(fullfile( target_dir, ...
				strcat( mf.name_F(mf.name_F_ind(f,1),1:mf.name_F_len(f,1)), '_ALL_SIMULATED_', mf.time_now, '.csv' ) ...
			), ...
			[	mf.t_sim( 1:mf.data_sim_len(f,1), mf.t_sim_ind(f,1) ), ...
				mf.data_sim( 1:mf.data_sim_len(f,1), mf.dataindices(f,1:mf.size_F(f,1)) ) ...
			]);
	end

	% Plot whole datafile
	if outmode.total
		fh = figure('Visible','off');
		col_map = colormap(fh, jet(double(mf.size_F(f,1))));

		for i = 1:mf.size_F(f,1)
			plot( mf.t(1:mf.data_len(f,1), mf.t_ind(f,1)), ...
				mf.data(1:mf.data_len(f,1), mf.dataindices(f,i) ), ...
				'-', 'linewidth',1, 'color',col_map(i,:) );
			hold on;
			plot( mf.t_sim(1:mf.data_sim_len(f,1), mf.t_sim_ind(f,1)), ...
				mf.data_sim(1:mf.data_sim_len(f,1), mf.dataindices(f,i)), ...
				'-', 'linewidth',2, 'color',col_map(i,:));
		end

		title([model.marker ' Fluorescence']);
		xlim([ 0, max( mf.t(1:mf.data_len(f,1), mf.t_ind(f,1)) ) ]);
		xlabel([ 'Time ' mf.t_unit ]);
		ylabel('Fluorescence Intensity [a.u.]');
		set(fh, 'PaperUnits','centimeters', 'PaperSize',[25 14], ...
			'PaperPosition',[0 0 25 14], 'PaperPositionMode','manual');
		print(fh, '-dpdf', fullfile( target_dir, ...
			strcat(mf.name_F(mf.name_F_ind(f,1),1:mf.name_F_len(f,1)), '_ALL_FIT_', mf.time_now, '.pdf')));
		delete(fh);
	end
end

%% Delete logfile
delete(plot_log);

%% Auxiliary log function
function write_plot_log(msg)
	% write_plot_log writes the message msg to a log file
	fh_log = fopen(plot_log, 'w');

	if fh_log == -1
		warning('Cannot open plot logfile “%s”\n', plot_log);
	else
		fprintf(fh_log, msg);
		fclose(fh_log);
	end
end
end