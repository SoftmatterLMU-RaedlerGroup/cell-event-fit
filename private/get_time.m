function timestring = get_time
%get_time builds a time string for log entries
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

timevec = clock;
hh = timevec(4);
mm = timevec(5);
ss = timevec(6);

timestring = [ '[' num2str(hh, '%02d') ':' num2str(mm, '%02d') ':' ...
	num2str(ss, '%02.0f') '] ' ];
