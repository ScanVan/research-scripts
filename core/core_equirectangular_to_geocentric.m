
    %  research-scripts
    %
    %     Nils Hamel - nils.hamel@bluewin.ch
    %     Copyright (c) 2016-2019 EPFL, HES-SO Valais
    %
    %  This program is free software: you can redistribute it and/or modify
    %  it under the terms of the GNU General Public License as published by
    %  the Free Software Foundation, either version 3 of the License, or
    %  (at your option) any later version.
    %
    %  This program is distributed in the hope that it will be useful,
    %  but WITHOUT ANY WARRANTY; without even the implied warranty of
    %  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %  GNU General Public License for more details.
    %
    %  You should have received a copy of the GNU General Public License
    %  along with this program.  If not, see <http://www.gnu.org/licenses/>.

    function c_geocentric = core_equirectangular_to_geocentric( c_equirectangular, c_width, c_height )

        % initialise memory %
        c_geocentric = zeros( size( c_equirectangular ) );

        % coordinate conversion - equirectangular mapping to geocentric %
        c_geocentric(:,1) = ( c_equirectangular(:,1) / c_width ) * 2.0 * pi;

        % coordinate conversion - equirectangular mapping to geocentric %
        c_geocentric(:,2) = ( ( c_equirectangular(:,2) / ( c_height - 1 ) ) - 0.5 ) * pi;

    end
