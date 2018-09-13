function Y = psiva_old_simulate(t, params)
% Provides a function with two peaks
%
% Arguments:
% t: time
% params:  coefficients of the sigmoid functions:
%		params(1):	offset of the function
%		params(2):	scaling factor of the function
%		params(3):	jump point of first peak
%		params(4):	steepness of step of first peak
%		params(5):	decay constant of first peak
%		params(6):	decay shift of first peak
%		params(7):	offset at beginning
%		params(8):	slope at beginning
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

Y = zeros(length(t), 1);

for i = 1:length(t)
    Y(i) = params(1) + params(2) * max([ ...
				sigmoid(t(i), params(3:6)) ...
				ramp(t(i), params(7:8)) ...
			]);
end

%% Sigmoid function
function S = sigmoid(t, params)

S = 1 / ( 1 + exp( (params(1)-t) * params(2) ) ) ...
    * exp(- params(3) * (t - params(4)) );


%% Ramp function
function R = ramp(t, params)
R = params(1) - t * params(2);