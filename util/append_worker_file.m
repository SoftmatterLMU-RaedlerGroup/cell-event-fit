function append_worker_file(temp_dir, time_now, worker_id, i, ...
	params, ...
	data_sim, ...
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
	custom_data_labels, ...
	custom_data_values)
%APPEND_WORKER_FILE appends new data to the worker-specific matfile
%
% Input parameters
% ================
% temp_dir    Char vector, path of the directory of the worker matfile
% time_now    Char vector, timestamp of session start (only for filenames)
% worker_id   Non-negative integer, ID of worker (only for filenames)
% i           Index of current trace in global matfile
% [others]    parameters returned by postprocessing
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

% Open matfile related to this worker and save computation results
wfname = fullfile(temp_dir, [ 'TEMP_' time_now '_worker' num2str(worker_id) '.mat' ]);

% Copy data to temporary matfile
if exist(wfname, 'file')
	% Rename old temporary matfile
	[path,name,ext] = fileparts(wfname);
	wfoname = fullfile(path, [name '_old' ext]);
	movefile(wfname, wfoname);

	% Open old and new temporary matfile
	wfo = matfile(wfoname);
	wf = matfile(wfname);

	% Write old and new data to new temporary matfile
	wf.version = wfo.version;
	append_to_mf_vec('map', i)
	append_to_mf_vec('max_val', max_val);
	append_to_mf_vec('max_ind', max_ind);
	append_to_mf_vec('min_val', min_val);
	append_to_mf_vec('min_ind', min_ind);
	append_to_mf_vec('amplitude', amplitude);
	append_to_mf_vec('data_amp', data_amp);
	append_to_mf_vec('logPost', logPost);
	append_to_mf_vec('t_event', t_event);
	append_to_mf_vec('event_deriv', event_deriv);
	append_to_mf_vec('fit_type', fit_type);

	old_size = size(wfo, 'params');
	new_length = length(params);
	total_size = [(old_size(1) + 1), max(old_size(2), new_length)];
	wf.params = NaN(total_size, 'single');
	wf.params(1:old_size(1), 1:old_size(2)) = wfo.params;
	wf.params(total_size(1), 1:new_length) = params;

	old_size = size(wfo, 'data_sim');
	new_length = length(data_sim);
	total_size = [max(old_size(1), new_length), (old_size(2) + 1)];
	wf.data_sim = NaN(total_size, 'single');
	wf.data_sim(1:old_size(1), 1:old_size(2)) = wfo.data_sim;
	wf.data_sim(1:new_length, total_size(2)) = data_sim;

	% DEPRECATED: copy rising/falling edge and fit parabola data
	if ~isempty(who(wfo, 'rising_edge_xoff'))
		wf.rising_edge_xoff = wfo.rising_edge_xoff;
	end
	if ~isempty(who(wfo, 'rising_edge_yoff'))
		wf.rising_edge_yoff = wfo.rising_edge_yoff;
	end
	if ~isempty(who(wfo, 'rising_edge_scale'))
		wf.rising_edge_scale = wfo.rising_edge_scale;
	end
	if ~isempty(who(wfo, 'falling_edge_xoff'))
		wf.falling_edge_xoff = wfo.falling_edge_xoff;
	end
	if ~isempty(who(wfo, 'falling_edge_yoff'))
		wf.falling_edge_yoff = wfo.falling_edge_yoff;
	end
	if ~isempty(who(wfo, 'falling_edge_scale'))
		wf.falling_edge_scale = wfo.falling_edge_scale;
	end
	if ~isempty(who(wfo, 'fit_parabola'))
		wf.fit_parabola = wfo.fit_parabola;
	end
	if ~isempty(who(wfo, 'fit_parabola_ind'))
		wf.fit_parabola_ind = wfo.fit_parabola_ind;
	end

	% Delete old matfile
	delete(wfoname);

else
	% Create new matfile and write data to it
	wf = matfile(wfname);
	wf.version = [2,0];
	wf.map = uint32(i);

	wf.max_val = single(max_val);
	wf.max_ind = uint32(max_ind);
	wf.min_val = single(min_val);
	wf.min_ind = uint32(min_ind);
	wf.amplitude = single(amplitude);
	wf.data_amp = single(data_amp);
	wf.logPost = single(logPost);
	wf.t_event = single(t_event);
	wf.event_deriv = single(event_deriv);
	wf.fit_type = int8(fit_type);

	wf.data_sim(1:length(data_sim),1) = single(data_sim(:));
	wf.params(1,1:length(params)) = single(params);

	% Custom data
	if isempty(custom_data_labels) && isempty(custom_data_values)
		wf.custom_data_map = zeros(1, 2, 'uint32');
		wf.custom_data_labels = int32.empty(0,1);
		wf.custom_data_values = single.empty(0,1);
	else
		len_lbl = length(custom_data_labels);
		len_val = length(custom_data_values);
		len_custom = max(len_lbl, len_val);
		wf.custom_data_map = uint32([1, len_custom]);
		wf.custom_data_labels = zeros(len_custom, 1, 'int32');
		wf.custom_data_labels(1:len_lbl, 1) = custom_data_labels;
		wf.custom_data_values = NaN(len_custom, 1, 'single');
		wf.custom_data_values(1:len_val, 1) = custom_data_values;
	end
end


	function append_to_mf_vec(name, value, dtype)
		if ~exist('value', 'var') || isempty(value)
			wf.(name) = wfo.(name);
		end
		s = whos(wfo, name);
		if ~exist('dtype', 'var')
			dtype = s.class;
		end
		if ~strcmp(class(value), dtype)
			value = cast(value, dtype);
		end
		value_len = length(value);
		new_len = max(s.size) + value_len;
		wf.(name) = zeros(new_len, 1, dtype);
		wf.(name)(1:new_len-value_len, 1) = wfo.(name);
		wf.(name)(new_len-(value_len-1:-1:0), 1) = value;
	end

end