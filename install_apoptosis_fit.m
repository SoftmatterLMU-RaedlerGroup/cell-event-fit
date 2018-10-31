%install_apoptosis_fit creates necessary entries in your Matlab path.
% To save the entries across sessions, save your Matlab path after running
% this script. Else, this script must be re-run in each session where you
% want to use this program.
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
base_dir = fileparts(mfilename('fullpath'));
sub_dirs = {'from_PESTO'; 'postprocessing'; 'preprocessing'; 'simulation'};
new_dirs = [ base_dir; fullfile(base_dir, sub_dirs) ];
addpath(new_dirs{:});
