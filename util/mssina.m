function [arr, index] = mssina(arr, s)
%mssina "Make Sure that String `s` is contained IN Array of strings `arr`"
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
	index = issina(arr, s);

	% If string is not in array, append it
	if ~(index > 0)
		arr = astoa(arr, s);
		index = size(arr, 1);
	end
end

function index = issina(s_arr, s)
%issina tests if string is in array;
% if true return (vertical) index, else return 0
	eq = false;
	index = 0;

	% Iterate over array lines
	for i = 1:size(s_arr, 1)
		si = s_arr(i,:);

		% Determine string length differences
		d_len = length(si) - length(s);

		% Test strings for equality, ignoring trailing null bytes
		if d_len == 0
			eq = strcmp(si, s);

		elseif d_len > 0
			l_s = length(s);
			eq = ( strcmp( si(1:l_s), s ) && ~isempty(regexp( si(l_s+1:end), '^\0+$', 'once') ));

		elseif d_len < 0
			l_si = length(si);
			eq = ( strcmp( s(1:l_si), si ) && ~isempty(regexp( s(l_si+1:end), '^\0+$', 'once' ) ));
		end

		if eq
			% String is already contained in s_arr; return its index
			index = i;
			break
		end
	end
end

function arr_new = astoa(arr_old, s)
%astoa "Append String TO Array of strings"

	% Get dimensions of old array
	[size_v, size_h] = size(arr_old);

	% Determine dimension differences
	d_len = length(s) - size_h;

	% Calculate dimensions of fill arrays
	if d_len < 0
		fill_s = -d_len;
		fill_a = 0;
	else
		fill_s = 0;
		fill_a = d_len;
	end

	% Build new string array
	arr_new = [	arr_old, zeros(size_v, fill_a); ...
				s, zeros(1, fill_s) ];
end