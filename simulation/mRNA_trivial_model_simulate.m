function Y = mRNA_trivial_model_simulate(t, params)
% Provides a simulation function for the trivial model of mRNA expression,
% as shown in eq. (3) on p. 684 of:
%		C. Leonhardt, G. Schwake et al.:
%		"Single-cell mRNA transfection studies: Delivery, kinetics and
%		statistics by numbers",
%		Nanomedicine: NBM 2014;10:679-688,
%		http://dx.doi.org/10.1016/j.nano.2013.11.008
% 
% Arguments:
% ==========
% t: time
% params:  coefficients of the sigmoid functions:
%	params(1):	offset of the function
%	params(2):	scaling factor (product of illumination intensity,
%				translation rate (k_TL) and initial mRAN concentration (m_0)
%	params(3):	mRNA degradation rate (\delta)
%	params(4):	protein degradation rate (\beta)
%	params(5):	time of expression onset (t_0)
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
    Y(i) = params(1) + params(2) ...
        / ( params(3) - params(4) ) ...
		* ( 1 - exp( - (params(3) - params(4)) * (t(i) - params(5)) ) ) ...
		* exp( - params(4) * (t(i) - params(5)) );
		%* t(i)^20 / (params(5)^20 + t(i)^20) ...
	
	if Y(i) < params(1)
		Y(i) = params(1);
	end
end
