function Y = pi_logistic_decay_simulate(t, params)
% Provides a sigmoid function with parameters
% 
% Arguments:
% t: time
% params:  coefficients of the sigmoid functions:
%    params(1): offset of the function
%    params(2): scaling factor of the function
%    params(3): jump point of the fuction (= inflection point)
%    params(4): steepness of step
%    params(5): decay constant
%    params(6): decay offset
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

Y = pdf('Logistic', t, params(3), params(4));

for i = 1:length(t)
    Y(i) = params(1) + params(2) * Y(i) ...
        * exp(- params(5) * (t(i) - params(6)) );
end
