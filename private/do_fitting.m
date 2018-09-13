function [varargout] = do_fitting(model, t, data, varargin)
%do_fitting configures the multistart optimization
%
% Input parameters
% ================
% model					model struct
%   model.name				name of the model
%   model.simulate			simulation function
%   model.par_num			number of fitting parameters
%   model.par_max			array of maximum values for fitting parameters
%   model.par_min			array of minimum values for fitting parameters
% t						time vector
% data					data vector
% varargin{1}			are parameters given logarithmically? (0: no, 1: yes)
%
% Output parameters
% =================
% P				fitted parameters
% logPost		logarithm of posterior
% toc			computation time needed for fitting
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

%% Define fitting options
% Options for multi-start optimization
options.fmincon = optimset('algorithm', 'interior-point', ...
                           'Display', 'off', ...
                           'GradObj', 'off', ...
                           'MaxIter', 300, ...
                           'MaxFunEvals', 300*4);
options.n_starts = 100; % or 150
options.plot = 'false';	% 'true','hidden'
options.mode = 'silent';

if numel(varargin) >= 1     % logarithmic parameters?
    options.logPost_options.logpar = varargin{1};
else
	options.logPost_options.logpar = 1;		% Per default: yes
end

% Profile likelihood options
%options_PL.fmincon = options.fmincon;

% Define parameter values
parameters.name = model.par_names;
parameters.min = model.par_min;
parameters.max = model.par_max;
parameters.number = model.par_num;

% Suppress annoying warnings
warning('off', 'MATLAB:nearlySingularMatrix'); %fmincon
warning('off', 'MATLAB:singularMatrix'); %fmincon

%% Perform multistart optimization
if nargout > 2
	tic
end

[logPost, par] = performOptimization(parameters, ...
            @(theta,options_lik) logLik( ...
                    model, theta, t, data, options_lik), ...
            options);

%% Return parameters
% Return parameters
varargout{1} = par;

% Return log-likelihood
varargout{2} = logPost;

% Return computing time
if nargout > 2
    varargout{3} = toc;
end
