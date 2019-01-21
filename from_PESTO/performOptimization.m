function [logPost_best, par_best] = performOptimization(parameters, logPosterior, options)
%performOptimization performs the multistart optimization
%
% Copyright © 2018-2019 Daniel Woschée <daniel.woschee@physik.lmu.de>
% Faculty of Physics / Ludwig-Maximilians-Universität München
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License version 2
% as published by the Free Software Foundation.
% However, parts of this file may be subject to the license of PESTO. See
% the file "LICENSE" in this folder for details of the license of PESTO.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.

%% Define default options
defaultoptions.fmincon = optimset('algorithm', 'interior-point',...
							'GradObj', 'off',...
							'MaxIter', 300,...
							'MaxFunEvals', 300*4);
defaultoptions.n_starts = 100; % or 150

defaultoptions.logPost_options.sign = 'negative';
defaultoptions.logPost_options.grad_ind = (1:parameters.number)';
defaultoptions.logPost_options.logpar = true;

if nargin >= 3
	options = setdefault(options, defaultoptions);
else
	options = defaultoptions;
end

if ~isfield(parameters, 'guess')
	parameters.guess = [];
elseif size(parameters.guess, 2) > options.n_starts
	options.n_starts = size(parameters.guess, 2);
end

%% Sampling of starting points
% Sampling from latin hypercube
par0 = [parameters.guess; ...
	bsxfun(@plus, parameters.min, bsxfun(@times, parameters.max - parameters.min, ...
	lhsdesign(options.n_starts - size(parameters.guess, 2), parameters.number, 'smooth', 'off')))];

%% Prepare arrays of optimization results
logPost_opt = zeros(options.n_starts, 1);
par_opt = zeros(options.n_starts, parameters.number);

%% Loop: Multi-starts
for iMS = 1 : options.n_starts
	% Check if initial parameter values are valid
	if ~isfinite(logPosterior(par0(iMS,:), options.logPost_options))
		logPost_opt(iMS) = NaN;
		par_opt(iMS,:) = NaN(1, parameters.number);
		continue
	end

	% Do optimization using `fmincon`
	try
		% fmincon as local optimizer
		[logPost_opt(iMS), par_opt(iMS,:), ~, ~, ~] = performOptimizationFmincon( ...
				parameters, logPosterior, par0(iMS,:), options);

	catch ErrMsg
		warning(['Multi-start number ' num2str(iMS) ' failed. More details on the error:']);
		disp(['Last Error in function ' ErrMsg.stack(1).name ', line ' ...
			num2str(ErrMsg.stack(1).line) ', file ' ErrMsg.stack(1).file '.']);

		% Assign values to exitflag, gradient etc. to avoid
		% failure when assigning results
		logPost_opt(iMS) = NaN;
		par_opt(iMS,:) = NaN(1, parameters.number);
	end
end


%% Find best multistart result
[~,ind] = max(logPost_opt);
logPost_best = logPost_opt(ind);

if ~isfinite(logPost_best)
	logPost_best = NaN;
	par_best = NaN;
else
	par_best = par_opt(ind,:);
end

end % end of `performOptimization`


function [logPost_opt, par_opt, exitflag, n_objfun, n_iter] ...
	= performOptimizationFmincon(parameters, negLogPost, par0, options)

	[par_opt, negLogPost_opt, exitflag, results_fmincon] = fmincon( ...
			@(p) negLogPost(p, options.logPost_options), ...  % negative log-likelihood function
			par0, ...    % initial parameter
			[], [], ...  % linear inequality constraints
			[], [], ...  % linear equality constraints
			parameters.min, ...     % lower bound
			parameters.max, ...     % upper bound
			[], ...      % nonlinear constraints
			options.fmincon);   % options

	logPost_opt = -negLogPost_opt;

	% Assignment of disgnosis
	n_objfun = results_fmincon.funcCount;
	n_iter = results_fmincon.iterations;
end