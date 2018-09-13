function Y = polynomial_sigmoid_simulate(t, params)
% Provides a breakdown function
%
% Arguments:
% t: time
% params:  coefficients of the sigmoid functions:
%   params(1): height of parabola vertex
%   params(2): positive quadratic coefficient of parabola
%   params(3): temporal location of parabola vertex
%	params(4): exponent of nonlinear term
%   params(5): temporal location of breakdown
%   params(6): breakdown steepness
%   params(7): final constant value
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

	coeff_const = params(1);
	coeff_nonlin = params(2);
	T_vertex = params(3);
	expon = params(4);
	T_step = params(5);
	steepness = params(6);
	end_val = params(7);

	Y = zeros(length(t), 1);

	for i = 1:length(t)
		S = sigmoid(t(i));
		Y(i) = S * (coeff_const + coeff_nonlin * (t(i) - T_vertex)^expon ) ...
			+  (1 - S) * end_val;
	end

	function S = sigmoid(t)
		S = 1 / (1+exp( (t-T_step) * steepness ));
	end

end
