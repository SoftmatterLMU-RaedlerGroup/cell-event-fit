function logL = logLik(model, theta, t, D, varargin)
%logLik calculates the logarithmic likelihood for given parameter estimates
%
% Parameters
% ==========
% model     Struct with the following fields:
%    model.name      Name of the model
%    model.simulate  Function handle of a model function with signature:
%                        Y = model.simulate(t, 10.^theta)
%    model.par_num   Number of parameters (the first model.par_num
%                    parameters will be passed to model.simulate)
% theta     Vector of parameters, where the last parameter is sigma
% t         Vector of timepoints
% D         Vector of measured data (single cell)
% options   Optional option struct
% max_ind	Maximum index to be evaluated
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

%% Initialization
theta = theta(:);
if length(theta) < model.par_num
    error('ERROR in calculation loglikelihood: too few parameters given')
end

options.sign = 'negative';
options.grad_ind = (1:length(theta))';
options.dist = 'normal'; % alternative: 'log-normal'; not implemented yet
options.logpar = 1; % 1: take 10^params, 0: take params

if nargin > 4
    options = setdefault(varargin{1},options);
end

% If options.logpar == true, consider theta as logarithmic values
if options.logpar
    theta = 10 .^ theta;
end

% Calculate values for given parameters
Y = model.simulate(t, theta);

% If a maximum index is given, evaluate logL only up to that index
if nargin > 5 && varargin{2} > 0
    max_ind = varargin{2};

    D = D(1:max_ind,:);
    Y = Y(1:max_ind,:);
end

%% Evaluation of likelihood function
sigma = theta(end); % non-system parameter (variance)

% Log-likelihood:
% log of gaussian distribution
logL = 1/2 * sum( log(2*pi*sigma^2) + ((D - Y) ./ sigma).^2 );

if isnan(logL)
    logL = inf;
end

%% Output
if strcmpi(options.sign, 'positive')
	logL = -logL;
end
