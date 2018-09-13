function Y = lyso_2_simulate(t, params)
% Provides a sigmoid function with parameters
%
% Arguments:
% t: time
% params:  coefficients of the sigmoid functions:
%	params(1): temporal position of parabola vertex
%	params(2): parabola vertex height
%	params(3): linear coefficient of parabola
%	params(4): negative quadratic coefficient of parabola
%	params(5): position of sigmoid transition
%	params(6): steepness of sigmoid transition
%	params(7): constant end value
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
		Y(i) = curv(t(i));
	end

	%% Sigmoid function
	function Y = sigmoid(t)
		Y = 1 / (1 + exp( (t - params(5)) * params(6) ));
	end

	%% The whole curve
	function P = curv(t)
		x = t - params(1);
		P = ( params(2) + x * (params(3) - x * params(4)) ) * sigmoid(t) ...
			+ params(7) * ( 1 - sigmoid(t) );
	end

end