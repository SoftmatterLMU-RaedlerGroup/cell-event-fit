function parallel_interactive_task(mf, indices)
%PARALLEL_INTERACTIVE_TASK performs interactive fitting in parallel
%
% Input arguments:
% ================
%	mf			handle to matfile holding the data
%	indices		(optional) indices of traces to fit; else, fit all traces
%
% The Parallel Processing Toolbox is required for standing to benefit from
% the parallelization; else, the fitting is performed serially.
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

if ~exist('indices', 'var') || isempty(indices)
	indices = 1:mf.ntraces;
end

if isempty(gcp)
	is_parallel = false;
else
	is_parallel = true;
end
futures = parallel.FevalFuture.empty;
future_map = containers.Map('KeyType','uint8', 'ValueType','uint32');
tempinter_map = containers.Map('KeyType','uint32', 'ValueType','any');
MAX_FUTURES = 100;

while ~isempty(indices) || ~isempty(futures)
	if ~isempty(futures)
		[idxFut, S] = fetchNext(futures);
		f, idx = remove_future(idxFut);
		% TODO
	end

	if ~isempty(indices)
		% Pop index of a trace to fit
		idx = indices(1);
		indices = indices(2:end);

		% Get trace data
		file_ind = mf.file_ind(idx, 1);
		index_F = mf.index_F(idx, 1);
		model = mf.models(mf.model_ind(file_ind, 1), 1);
		datalen = mf.data_len(file_ind, 1);

		mf_data = mf.data(1:datalen, idx);
		mf_t = mf.t(1:datalen, mf.t_ind(file_ind, 1));
		
		% User interaction
		private = inter_private(idx);
		[is_to_fit, params, data_indices, par_fun, private] = ...
			model.interactive(mf_t, mf_data, model, params, mf, private);
		if is_to_fit
			inter_private(idx, private);
		else
			inter_private(idx, []);
		end

		% Get parameters
		if is_to_fit

			% Define data range to be fitted (Default: fit whole data)
			if isempty(data_indices) && isa(model.preproc, 'function_handle')
				data_indices = model.preproc( ...
					mf.data(1:datalen, idx), ...
					mf, idx);
			end
			if isempty(data_indices)
				data_indices = 1:datalen;
			end
			mf_t = mf_t(data_indices);
			mf_data = mf_data(data_indices);
			mf_tsim = mf.t_sim(1:mf.data_sim_len(file_ind,1), mf.t_sim_ind(file_ind,1));

			% Perform optimization
			add_future(idx_trace, index_F, model, mf_t, mf_data, mf_tsim, par_fun)
		end
		
		
	end
	
end

	function [F, idx_trace] = remove_future(idx)
		%remove_future removes a future from the future array
		if isempty(idx)
			F = [];
			idx_trace = [];
			return
		end

		F = futures(idx);
		idx_trace = future_map(idx);

		futures = futures([1:idx-1, idx+1:end]);
		remove(future_map, idx);
	end

	function add_future(idx_trace, index_F, model, t, data, tsim, par_fun)
		%add_future adds a new future to the future array
		futures(end+1) = parfeval(parevalfun, 1, ...
			index_F, model, t, data, tsim, par_fun);
		future_map(length(futures)) = idx_trace;
	end

	function p = inter_private(idx, p_new)
		%inter_private writes, reads or removes data in `tempinter_map`
		p = [];
		if nargin == 1 && tempinter_map.iskey(idx)
			p = tempinter_map(idx);
		elseif isempty(p_new)
			tempinter_map.remove(idx);
		else
			tempinter_map(idx) = p_new;
		end
	end

end


function S = parevalfun(index_F, model, t, data, t_sim, par_fun)
	% Initialize struct of results
	S = struct( ...
		'params', single.empty(0, model.par_num), ...
		'logPost', NaN('single'), ...
		'data_sim', NaN(size(t_sim), 'single'), ...
		'max_val', NaN('single'), ...
		'max_ind', uint32(0), ...
		'min_val', NaN('single'), ...
		'min_ind', uint32(0), ...
		'amplitude', NaN('single'), ...
		'data_amp', NaN('single'), ...
		't_event', NaN('single'), ...
		'event_deriv', NaN('single'), ...
		'fit_type', int8(0), ...
		'custom_data_labels', int32.empty(0,1), ...
		'custom_data_values', single.empty(0,1) ...
		);

	% Perform fitting
	[S.params(:), S.logPost(1)] = do_fitting(model, t, data, par_fun);
	if any(model.par_log)
		if length(model.par_log) == 1
			S.params = 10 .^ S.params;
		else
			S.params(model.par_log) = 10 .^ S.params(model.par_log);
		end
	end

	% Calculate fitted function
	S.data_sim(:) = model.simulate(t_sim, S.params);

	%% Postprocessing
	% Determine extrema and amplitude of fit
	[S.max_val(1), S.max_ind(1)] = max(S.data_sim);
	[S.min_val(1), S.min_ind(1)] = min(S.data_sim);
	S.amplitude(1) = S.max_val - S.min_val;

	% Determine amplitude of experimental data
	S.data_amp(1) = max(data) - min(data);

	if isfield(model, 'postproc') && isa(model.postproc, 'function_handle')
		% Write structure as input arguments for postprocessing
		D = struct('data', data, 'index', index_F);
		F = struct('t_sim', t_sim, 'modelname', model.name, ...
			'simulate', model.simulate);
		R = struct('data_sim', S.data_sim, 'max_val', S.max_val, 'params', S.params, ...
				'max_ind', S.max_ind, 'min_val', S.min_val, 'min_ind', S.min_ind, ...
				'amplitude', S.amplitude, 'data_amp', S.data_amp);

		% Run model specific postprocessing routine
		[t_event, event_deriv, fit_type, custom_data_labels, custom_data_values] = ...
			model.postproc(D, F, R);

		% Write postprocessing results to results structure
		if ~isempty(t_event)
			S.t_event(1) = t_event(1);
		end
		if ~isempty(event_deriv)
			S.event_deriv(1) = event_deriv(1);
		end
		if ~isempty(fit_type)
			S.fit_type(1) = fit_type(1);
		end

		lbl_len = length(custom_data_labels);
		val_len = length(custom_data_values);
		custom_len = max(lbl_len, val_len);
		if custom_len
			if lbl_len < custom_len
				S.custom_data_labels = NaN(custom_len, 1, 'single');
				S.custom_data_labels(1:lbl_len) = custom_data_labels;
				S.custom_data_values = custom_data_values;
			elseif val_len < custom_len
				S.custom_data_labels = custom_data_labels;
				S.custom_data_values = NaN(custom_len, 1, 'single');
				S.custom_data_values(1:val_len) = custom_data_values;
			else
				S.custom_data_labels = custom_data_labels;
				S.custom_data_values = custom_data_values;
			end
		end

	else
		disp(['No postprocessing routine defined for model ' model.name]);
	end

end