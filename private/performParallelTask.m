function [ Ffile ] = performParallelTask( Ffile, ntraces )
%performParallelTask performs the parallel computing task.
%
% Input parameters
% ================
% Ffile		String containing the path of the input mat-file
% ntraces	The number of traces to be processed
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

%% Execute parallel job
% Suppress warning on parallel pool startup
warning('off','MATLAB:datetime:NonstandardSystemTimeZone')

parfor i = 1:ntraces
%for i = 1:ntraces
	% Read in temporary files
	mf = matfile(Ffile);
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
	if strcmp(model.name, 'pSIVA_old') || strcmp(model.name, 'pSIVA')
		% For pSIVA, only fit to the minimum after the first peak if
		% there are two large peaks
		data_indices = preproc_psiva(mf.data(1:mf.data_len(file_ind,1),i));
	end

	% Define data range to be fitted (Default: fit whole data)
	if isempty(data_indices)
		data_indices = 1:mf.data_len(file_ind,1);
	end

	%% Perform the actual fit
	% Get parameters
	[params,logPost] = do_fitting(model, mf.t(data_indices,mf.t_ind(file_ind,1)), mf.data(data_indices,i) );

	% Calculate real parameters (not logarithmic ones)
	params = 10 .^ params;

	% Ensure that params vector has right dimensions
	if iscolumn(params)
		params = params';
	end

	% Evaluate fitted function
	data_sim = model.simulate(mf.t_sim(1:mf.data_sim_len(file_ind,1),mf.t_sim_ind(file_ind,1)), params);

	if isfield(model, 'postproc')
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
		else
			t_event = NaN;
			fit_info.type = NaN;
			event_deriv = NaN;
		end

		fit_type = fit_info.type;
		if isfield(fit_info, 'parabola')
			fit_parabola = fit_info.parabola;
		else
			fit_parabola = [];
		end

		if isfield(fit_info, 'rising_edge_xoff')
			rising_edge_xoff = fit_info.rising_edge_xoff;
		else
			rising_edge_xoff = NaN;
		end
		if isfield(fit_info, 'rising_edge_yoff')
			rising_edge_yoff = fit_info.rising_edge_yoff;
		else
			rising_edge_yoff = NaN;
		end
		if isfield(fit_info, 'rising_edge_scale')
			rising_edge_scale = fit_info.rising_edge_scale;
		else
			rising_edge_scale = NaN;
		end

		if isfield(fit_info, 'falling_edge_xoff')
			falling_edge_xoff = fit_info.falling_edge_xoff;
		else
			falling_edge_xoff = NaN;
		end
		if isfield(fit_info, 'falling_edge_yoff')
			falling_edge_yoff = fit_info.falling_edge_yoff;
		else
			falling_edge_yoff = NaN;
		end
		if isfield(fit_info, 'falling_edge_scale')
			falling_edge_scale = fit_info.falling_edge_scale;
		else
			falling_edge_scale = NaN;
		end

	else
		disp(['No postprocessing routine defined for model ' model.name]);
		max_val = NaN;
		max_ind = NaN;
		min_val = NaN;
		min_ind = NaN;
		amplitude = NaN;
		data_amp = NaN;
		t_event = NaN;
		event_deriv = NaN;
		fit_parabola = [];
		rising_edge_xoff = NaN;
		rising_edge_yoff = NaN;
		rising_edge_scale = NaN;
		falling_edge_xoff = NaN;
		falling_edge_yoff = NaN;
		falling_edge_scale = NaN;
	end

	%% Write results to temporary matfile

	% Open matfile related to this worker and save computation results
	wfname = fullfile(mf.temp_dir, [ 'TEMP_' mf.time_now '_worker' num2str(w_id) '.mat' ]);

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
end

%% Perform cleanup (combine single worker matfiles to global matfile)
disp([ get_time 'Cleaning up worker matfiles ...' ])

% Get number of workers
p = gcp('nocreate');
if isempty(p)
	nWorkers = 0;
else
	nWorkers = p.NumWorkers;
end

% Perform cleanup
reduceData(Ffile, nWorkers);

%% Continue post-processing
disp([ get_time 'Fitting finished; starting amplitude comparison ...' ])

% Open global matfile
mf = matfile(Ffile, 'Writable', true);
ndatafiles = mf.ndatafiles;

% Get file-wide maximum of amplitude
for f = 1:ndatafiles

	% Read data from matfile
	model = mf.models(mf.model_ind(f,1),1);

	% Test if postprocessing routine is defined
	if ~isempty( model.postproc )
		amp_max_restricted = -Inf;

		% Find global maximum (clean for runaways by comparison with median)
		[amps_sort,amps_sort_ind] = sort(mf.amplitude(mf.dataindices(f,1:mf.size_F(f,1)),1), 'descend');
		amp_median = amps_sort(ceil( length(amps_sort) / 2 ));

		for i = 1:length(amps_sort)
			j = mf.dataindices(f, amps_sort_ind(i));
			if ~isnan(mf.t_event(j,1)) && amps_sort(i) < 100 * amp_median
				amp_max_restricted = amps_sort(i);
				break
			end
		end

		% Write filewide amplitude data to global matfile
		mf.amp_max(f,1) = amp_max_restricted;
		mf.noise(f,1) = amp_max_restricted * .1;

		% Continue model-specific postprocessing
		for i = mf.dataindices(f, 1:mf.size_F(f,1))
			% Get event time
			if mf.amplitude(i,1) <= mf.noise(f,1) || ...
					mf.amplitude(i,1) > 1.25 * mf.data_amp(i,1) || ...
					mf.amplitude(i,1) < 0.4 * mf.data_amp(i,1) %|| ...
					%mf.amplitude(i,1) > mf.amp_max(f,1)
				mf.t_event(i,1) = NaN;
				% else leave mf.t_event(i,1) as is
			end
		end
	end
end

end