function Y = lyso_1_simulate(t, params)
% Provides a sigmoid function with parameters
%
% Arguments:
% t: time
% params:  coefficients of the sigmoid functions:
%	params(1): end time of linear, start time of non-linear behaviour
%	params(2): steepness of (initial) linear behaviour
%	params(3): temporal position of parabola vertex
%	params(4): parabola vertex height
%	params(5): linear coefficient of parabola
%	params(6): negative quadratic coefficient of parabola
%	params(7): position of sigmoid transition
%	params(8): steepness of sigmoid transition
%	params(9): constant end value
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

	%% Allocate memory
	Y = zeros(length(t), 1);

	%% Calculate value
	for i = 1:length(t)
		if t < params(1)
			Y(i) = curv(params(1)) + params(2) * (params(1) - t(i));
		else
			Y(i) = curv(t(i));
		end
	end

	%% Sigmoid function
	function Y = sigmoid(t)
		Y = 1 / (1 + exp( (t - params(7)) * params(8) ));
	end

	%% The whole curve
	function P = curv(t)
		x = t - params(3);
		P = ( params(4) + x * (params(5) - x * params(6)) ) * sigmoid(t) ...
			+ params(9) * ( 1 - sigmoid(t) );
	end

end