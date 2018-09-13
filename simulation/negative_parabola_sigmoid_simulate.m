function Y = negative_parabola_sigmoid_simulate(t, params)
% Provides a breakdown function
%
% Arguments:
% t: time
% params:  coefficients of the sigmoid functions:
%   params(1): height of parabola vertex
%   params(2): negative quadratic coefficient of parabola
%   params(3): temporal location of parabola vertex
%   params(4): temporal location of breakdown
%   params(5): breakdown steepness
%   params(6): final constant value
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
	T_vertex = params(3);
	T_step = params(4);
	gamma = params(5);
	b = params(6);

	Y = zeros(length(t), 1);

	for i = 1:length(t)
		S = 1 / (1 + exp( (t(i)-T_step) * gamma ));
		Y(i) = S * ( coeff_const - coeff_quad * (t(i) - T_vertex)^2 ) ...
			+  (1 - S) * b;
	end
end
