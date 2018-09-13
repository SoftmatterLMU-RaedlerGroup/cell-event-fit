function reduceData( Ffile, nWorkers, isDebug )
%REDUCEDATA Merges all worker mat-files into the general mat-file
%
% Input arguments:
%	Ffile			char vector of the general mat-file
%	nWorkers		number of workers that have created worker mat-files
%	isDebug			if true, enable verbose output
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

%% Get input parameters
if ~exist(Ffile, 'file')
	error(['File not found: ' Ffile])
end

if nargin < 2 || ~isnumeric(nWorkers) || nWorkers < 0
	nWorkers = 0;
end

if nargin < 3 || ~islogical(isDebug)
	isDebug = false;
end

%% Open global matfile
mf = matfile(Ffile, 'Writable', true);
ntraces = mf.ntraces;
ndatafiles = mf.ndatafiles;

if isDebug
	fprintf('Found %d traces in %d files.\n', ntraces, ndatafiles);
end

%% Allocate arrays for results in global matfile
% trace-wise indexing
mf.amplitude(1:ntraces,1) = single(NaN);
mf.data_amp(1:ntraces,1) = single(NaN);
mf.max_val(1:ntraces,1) = single(NaN);
mf.max_ind(1:ntraces,1) = uint32(NaN);
mf.min_val(1:ntraces,1) = single(NaN);
mf.min_ind(1:ntraces,1) = uint32(NaN);
mf.logPost(1:ntraces,1) = single(NaN);
mf.t_event(1:ntraces,1) = single(NaN);
mf.event_deriv(1:ntraces,1) = single(NaN);
mf.fit_type(1:ntraces,1) = int8(0);
mf.rising_edge_xoff(1:ntraces,1) = single(NaN);
mf.rising_edge_yoff(1:ntraces,1) = single(NaN);
mf.rising_edge_scale(1:ntraces,1) = single(NaN);
mf.falling_edge_xoff(1:ntraces,1) = single(NaN);
mf.falling_edge_yoff(1:ntraces,1) = single(NaN);
mf.falling_edge_scale(1:ntraces,1) = single(NaN);

mf.data_sim = NaN(max(mf.data_sim_len), ntraces, 'single');
mf.params = NaN(ntraces, mf.max_par_num, 'single');

mf.fit_parabola_ind(1:ntraces,1) = uint32(0);

% file-wise indexing
mf.noise(1:ndatafiles,1) = single(NaN);
mf.amp_max(1:ndatafiles,1) = single(-Inf);

% Special cases
fit_parabola = NaN(0, 'single');

%% Loop over worker matfiles
for f = 0:nWorkers
	% Open worker matfile
	wfname = fullfile(mf.temp_dir, [ 'TEMP_' mf.time_now '_worker' num2str(f) '.mat' ]);
	if isDebug
		fprintf('Opening worker matfile %s\n', wfname);
	end
	if ~exist(wfname, 'file')
		if isDebug
			fprintf('Matfile not found; continue.\n');
		end
		continue
	end
	wf = matfile(wfname);

	% Copy data from worker matfile to global matfile
	for i = 1:length(wf.map)
		
		if isDebug
			fprintf('Processing trace %d …\n', i);
		end

		% Get local copy of worker matfile array element
		j = wf.map(i,1);

		% Copy scalar data into global matfile
		mf.amplitude(j,1) = wf.amplitude(i,1);
		mf.data_amp(j,1) = wf.data_amp(i,1);
		mf.max_val(j,1) = wf.max_val(i,1);
		mf.max_ind(j,1) = wf.max_ind(i,1);
		mf.min_val(j,1) = wf.min_val(i,1);
		mf.min_ind(j,1) = wf.min_ind(i,1);
		mf.logPost(j,1) = wf.logPost(i,1);
		mf.t_event(j,1) = wf.t_event(i,1);
		mf.event_deriv(j,1) = wf.event_deriv(i,1);
		mf.fit_type(j,1) = wf.fit_type(i,1);
		mf.rising_edge_xoff(j,1) = wf.rising_edge_xoff(i,1);
		mf.rising_edge_yoff(j,1) = wf.rising_edge_yoff(i,1);
		mf.rising_edge_scale(j,1) = wf.rising_edge_scale(i,1);
		mf.falling_edge_xoff(j,1) = wf.falling_edge_xoff(i,1);
		mf.falling_edge_yoff(j,1) = wf.falling_edge_yoff(i,1);
		mf.falling_edge_scale(j,1) = wf.falling_edge_scale(i,1);

		% Copy non-scalar data into local struct array
		dsl = mf.data_sim_len(mf.file_ind(j, 1), 1);
		mf.data_sim(1:dsl, j) = wf.data_sim(1:dsl, i);
		mf.params(j, 1:size(wf.params, 2)) = wf.params(i,:);

		if wf.fit_parabola_ind(i,1)
			d_len = size(wf.fit_parabola, 1) - size(fit_parabola, 1);
			if d_len > 0
				fit_parabola = [ [fit_parabola; NaN(d_len, size(fit_parabola, 2))], ...
					wf.fit_parabola(:,wf.fit_parabola_ind(i,1)) ];
			elseif d_len < 0
				fit_parabola = [ fit_parabola, ...
					[wf.fit_parabola(:,wf.fit_parabola_ind(i,1)); NaN(-d_len, 1)] ];
			else
				fit_parabola = [ fit_parabola, ...
					wf.fit_parabola(:,wf.fit_parabola_ind(i,1)) ];
			end
			mf.fit_parabola_ind(j,1) = uint32( size(fit_parabola, 2) );
		end
	end

	% Delete worker matfile
	clear wf;
	
	if isDebug
		fprintf('Deleting worker matfile %s\n', wfname);
	end
	delete(wfname);
	clear wfname;
end

%% Write more complex arrays to matfile
mf.fit_parabola = fit_parabola;

if isDebug
	fprintf('Finished reducing.\n');
end
end
