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

if isempty(gcp)
	is_parallel = false;
else
	is_parallel = true;
end
futures = parallel.FevalFuture.empty;

while ~isempty(indices) || ~isempty(futures)
	% TODO
	
end