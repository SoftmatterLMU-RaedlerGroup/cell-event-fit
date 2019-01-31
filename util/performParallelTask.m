function [ Ffile ] = performParallelTask( Ffile )
%performParallelTask performs the parallel computing task.
%
% Input parameters
% ================
% Ffile		String containing the path of the input mat-file
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

%% Get session data
mf = matfile(Ffile);
has_interactivity = mf.has_interactivity;

%% Prepare parallel pool
% Suppress warning on parallel pool startup
warning('off','MATLAB:datetime:NonstandardSystemTimeZone')

% Get number of workers
p = gcp();
if isempty(p)
	nWorkers = 0;
else
	nWorkers = p.NumWorkers;
end

%% Check for traces with interactivity
if has_interactivity
	% Get interactive models
	n_models = size(mf, 'models');
	inter_models = false(1, n_models);
	for iM = 1:n_models
		m = mf.models(iM, 1);
		inter_models(iM) = m.interactive;
	end
	[~,~,inter_models] = find(inter_models);

	% Get files associated with interactive models
	[~,~,inter_files] = find(ismember(mf.model_ind, inter_models));

	% Get traces from files with and without interactivity
	inter_idx = ismember(file_ind, inter_files);

	[~,~,batch_idx] = find(~inter_idx);
	if ~isrow(batch_idx)
		batch_idx = batch_idx';
	end

	[~,~,inter_idx] = find(inter_idx);
	if ~isrow(inter_idx)
		inter_idx = inter_idx';
	end
else
	batch_idx = [];
end

%% Execute parallel job
if has_interactivity
	interactive_task(mf, inter_idx);
end
parallel_batch_task(mf, batch_idx);

%% Perform cleanup (combine single worker matfiles to global matfile)
disp([ get_time 'Cleaning up worker matfiles ...' ])

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