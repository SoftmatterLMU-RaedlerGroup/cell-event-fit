function Y = lyso_exponential_sigmoid_simulate(t, params)
% Provides a sigmoid function with parameters
% 
% Arguments:
% t: time
% params:  coefficients of the sigmoid functions:
%    params(1): offset of the function
%    params(2): scaling factor of the function
%    params(3): breakdown point
%    params(4): breakdown steepness
%    params(5): size of gaussian
%    params(6): shift of gaussian
%    params(7): width of gaussian
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

% Y is approximately theta-shaped: Y = params(1) + params(2) * theta(x-params(3))

for i = 1:length(t)
    Y(i) = params(1) + params(2) * ( ...
        ( 1 / (1+exp( (t(i)-params(3)) * params(4) )) ) ...
        + params(5) * exp( (log(t(i)) + params(6))^2 * params(7) ) ...
        );
end

function Y = sig(t)
	Y = t * params(1);
end