function Y = parabola_sigmoid_simulate(t, params)
% Provides a breakdown function
%
% Arguments:
% t: time
% params:  coefficients of the sigmoid functions:
%   params(1): height of parabola vertex
%   params(2): absolute value of quadratic coefficient of parabola
%	params(3): sign of quadratic coefficient of parabola
%   params(4): temporal location of parabola vertex
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
	coeff_quad = params(2);
	sgn = sign( log10(params(3)) );
	T_vertex = params(4);
	T_step = params(5);
	gamma = params(6);
	b = params(7);

	Y = zeros(length(t), 1);

	for i = 1:length(t)
		S = sigmoid(t(i));
		Y(i) = S * ( coeff_const + sgn * coeff_quad * (t(i) - T_vertex)^2 ) ...
			+  (1 - S) * b;
	end

	function S = sigmoid(t)
		S = 1 / (1+exp( (t-T_step) * gamma ));
	end

end
