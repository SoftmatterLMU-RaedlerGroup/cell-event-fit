function append_worker_file(temp_dir, time_now, worker_id, i, ...
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
%COMBINE_WORKER_FILES appends new data to the worker-specific matfile
%   Detailed explanation goes here
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
	wf.map(1:length(wfo.map)+1,1) = [wfo.map; i];

	wf.max_val = [wfo.max_val; max_val];
	wf.max_ind = [wfo.max_ind; max_ind];
	wf.min_val = [wfo.min_val; min_val];
	wf.min_ind = [wfo.min_ind; min_ind];
	wf.amplitude = [wfo.amplitude; amplitude];
	wf.data_amp = [wfo.data_amp; data_amp];
	wf.logPost = [wfo.logPost; logPost];
	wf.t_event = [wfo.t_event; t_event];
	wf.event_deriv = [wfo.event_deriv; event_deriv];
	wf.fit_type = [wfo.fit_type; fit_type];
	wf.rising_edge_xoff = [wfo.rising_edge_xoff; rising_edge_xoff];
	wf.rising_edge_yoff = [wfo.rising_edge_yoff; rising_edge_yoff];
	wf.rising_edge_scale = [wfo.rising_edge_scale; rising_edge_scale];
	wf.falling_edge_xoff = [wfo.falling_edge_xoff; falling_edge_xoff];
	wf.falling_edge_yoff = [wfo.falling_edge_yoff; falling_edge_yoff];
	wf.falling_edge_scale = [wfo.falling_edge_scale; falling_edge_scale];

	d_len = length(data_sim) - size(wfo.data_sim, 1);
	if d_len > 0
		wf.data_sim = [ [wfo.data_sim; NaN(d_len, size(wfo.data_sim, 2))], ...
			data_sim ];
	elseif d_len < 0
		wf.data_sim = [ wfo.data_sim, [data_sim; NaN(-d_len, 1)] ];
	else
		wf.data_sim = [ wfo.data_sim, data_sim ];
	end

	d_len = length(params) - size(wfo.params, 2);
	if d_len > 0
		wf.params = [ [wfo.params, NaN(size(wfo.params, 1), d_len)]; params ];
	elseif d_len < 0
		wf.params = [ wfo.params; [params, NaN(1, -d_len)] ];
	else
		wf.params = [ wfo.params; params ];
	end

	fit_parabola_new = wfo.fit_parabola;
	fit_parabola_ind = uint32(0);
	if ~isempty(fit_parabola)
		fit_parabola_ind(1) = size(fit_parabola_new, 2) + 1;
		fit_parabola_new( 1:length(fit_parabola), fit_parabola_ind ) = fit_parabola;
	end
	wf.fit_parabola = fit_parabola_new;
	wf.fit_parabola_ind = [ wfo.fit_parabola_ind; fit_parabola_ind ];

	% Delete old matfile
	delete(wfoname);

else
	% Create new matfile and write data to it
	wf = matfile(wfname);
	wf.map(1,1) = uint32(i);

	wf.max_val(1,1) = single(max_val);
	wf.max_ind(1,1) = uint32(max_ind);
	wf.min_val(1,1) = single(min_val);
	wf.min_ind(1,1) = uint32(min_ind);
	wf.amplitude(1,1) = single(amplitude);
	wf.data_amp(1,1) = single(data_amp);
	wf.logPost(1,1) = single(logPost);
	wf.t_event(1,1) = single(t_event);
	wf.event_deriv(1,1) = single(event_deriv);
	wf.fit_type(1,1) = int8(fit_type);
	wf.rising_edge_xoff(1,1) = single(rising_edge_xoff);
	wf.rising_edge_yoff(1,1) = single(rising_edge_yoff);
	wf.rising_edge_scale(1,1) = single(rising_edge_scale);
	wf.falling_edge_xoff(1,1) = single(falling_edge_xoff);
	wf.falling_edge_yoff(1,1) = single(falling_edge_yoff);
	wf.falling_edge_scale(1,1) = single(falling_edge_scale);

	wf.data_sim(1:length(data_sim),1) = single(data_sim(:));
	wf.params(1,1:length(params)) = single(params);

	if isempty(fit_parabola)
		wf.fit_parabola = single([]);
		wf.fit_parabola_ind(1,1) = uint32(0);
	else
		wf.fit_parabola = single(fit_parabola);
		wf.fit_parabola_ind(1,1) = uint32(1);
	end
end
