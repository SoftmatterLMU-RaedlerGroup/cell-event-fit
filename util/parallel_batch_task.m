function parallel_batch_task(mf, indices)
%PARALLEL_BATCH_TASK performs non-interactive fitting in parallel
%
% Input arguments:
% ================
%	mf			handle to matfile holding the data
%	indices		(optional) indices of traces to fit; else, fit all traces
%
% The Parallel Processing Toolbox is required for standing to benefit from
% the parallelization; else, the fitting is performed serially.
%
%
% Copyright © 2018-2019 Daniel Woschée <daniel.woschee@physik.lmu.de>
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

if ~exist('indices', 'var') || isempty(indices)
	indices = 1:mf.ntraces;
end

parfor i = indices
%for i = 1:ntraces
	% Read in temporary files
	file_ind = mf.file_ind(i,1);
	model = mf.models(mf.model_ind(file_ind,1), 1);

	% Write log message
	w_id = 0;
	if exist('getCurrentTask', 'file') == 2
		worker = getCurrentTask();
		if isprop(worker, 'ID')
			w_id = worker.ID;
		end
	end

	disp([get_time 'Worker ' num2str(w_id) ' processing trace ' ...
		num2str(mf.index_F(i,1)) ' from file ' num2str(file_ind)])

	% Allocate memory
	data_indices = [];

	%% Do marker-specific pre-processing
	if isa(model.preproc, 'function_handle')
		data_indices = model.preproc( ...
			mf.data(1:mf.data_len(file_ind,1),i), ...
			mf, i);
	end

	% Define data range to be fitted (Default: fit whole data)
	if isempty(data_indices)
		data_indices = 1:mf.data_len(file_ind,1);
	end

	%% Perform the actual fit
	% Get data points for fitting
	mf_data = mf.data(1:mf.data_len(file_ind, 1), i);
	mf_data = mf_data(data_indices);
	mf_t = mf.t(1:mf.data_len(file_ind, 1), mf.t_ind(file_ind, 1));
	mf_t = mf_t(data_indices);

	% Get parameters
	[params,logPost] = do_fitting(model, mf_t, mf_data);

	% Calculate real parameters (not logarithmic ones)
	if any(model.par_log)
		if length(model.par_log) == 1
			params = 10 .^ params;
		else
			params(model.par_log) = 10 .^ params(model.par_log);
		end
	end

	% Ensure that params vector has right dimensions
	if iscolumn(params)
		params = params';
	end

	% Evaluate fitted function
	data_sim = model.simulate(mf.t_sim(1:mf.data_sim_len(file_ind,1),mf.t_sim_ind(file_ind,1)), params);

	%% Postprocessing
	% Predefine default result values
	max_val = NaN;
	max_ind = NaN;
	min_val = NaN;
	min_ind = NaN;
	amplitude = NaN;
	data_amp = NaN;
	t_event = NaN;
	event_deriv = NaN;
	fit_parabola = [];
	fit_type = NaN;
	rising_edge_xoff = NaN;
	rising_edge_yoff = NaN;
	rising_edge_scale = NaN;
	falling_edge_xoff = NaN;
	falling_edge_yoff = NaN;
	falling_edge_scale = NaN;

	% Execute postprocessing function
	if isfield(model, 'postproc') && isa(model.postproc, 'function_handle')
		% Determine extrema and amplitude of fit
		[max_val, max_ind] = max(data_sim);
		[min_val, min_ind] = min(data_sim);
		amplitude = max_val - min_val;

		% Determine amplitude of experimental data
		data = mf.data(1:mf.data_len(mf.file_ind(i,1),1),i);
		data_amp = max(data) - min(data);

		% Write structure as input arguments for postprocessing
		D = struct( 'data',mf.data(:,i), 'index',mf.index_F(i,1) );
		F = struct( 't_sim', mf.t_sim( 1:mf.data_sim_len(file_ind,1), mf.t_sim_ind(file_ind,1) ), ...
			'modelname', mf.model_name(mf.model_name_ind(file_ind,1),1:mf.model_name_len(file_ind,1)), ...
			'simulate', model.simulate );
		R = struct( 'data_sim', data_sim, 'max_val', max_val, 'params', params, ...
				'max_ind', max_ind, 'min_val', min_val, 'min_ind', min_ind, ...
				'amplitude', amplitude, 'data_amp', data_amp );

		% Run model specific postprocessing routine
		if isfield(model, 'postproc') && isa(model.postproc, 'function_handle')
			[t_event, fit_info, event_deriv] = model.postproc(D, F, R);
		end

		fit_type = fit_info.type;
		if isfield(fit_info, 'parabola')
			fit_parabola = fit_info.parabola;
		end

		if isfield(fit_info, 'rising_edge_xoff')
			rising_edge_xoff = fit_info.rising_edge_xoff;
		end
		if isfield(fit_info, 'rising_edge_yoff')
			rising_edge_yoff = fit_info.rising_edge_yoff;
		end
		if isfield(fit_info, 'rising_edge_scale')
			rising_edge_scale = fit_info.rising_edge_scale;
		end

		if isfield(fit_info, 'falling_edge_xoff')
			falling_edge_xoff = fit_info.falling_edge_xoff;
		end
		if isfield(fit_info, 'falling_edge_yoff')
			falling_edge_yoff = fit_info.falling_edge_yoff;
		end
		if isfield(fit_info, 'falling_edge_scale')
			falling_edge_scale = fit_info.falling_edge_scale;
		end

	else
		disp(['No postprocessing routine defined for model ' model.name]);
	end

	%% Write results to temporary matfile
	append_worker_file(mf.temp_dir, mf.time_now, w_id, i, ...
		max_val, ...
		max_ind, ...
		min_val, ...
		min_ind, ...
		amplitude, ...
		data_amp, ...
		logPost, ...
		t_event, ...
		event_deriv, ...
		fit_type, ...
		rising_edge_xoff, ...
		rising_edge_yoff, ...
		rising_edge_scale, ...
		falling_edge_xoff, ...
		falling_edge_yoff, ...
		falling_edge_scale, ...
		data_sim, ...
		params, ...
		fit_parabola)
end
