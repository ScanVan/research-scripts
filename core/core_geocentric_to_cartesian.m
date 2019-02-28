
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

    function c_cartesian = core_geocentric_to_cartesian( c_geocentric )

        % initialise memory %
        c_cartesian = zeros( size( c_geocentric ) );

        % coordinates conversion - longitude/latitude to 3D unit sphere %
        c_cartesian(:,1) = cos( c_geocentric(:,2) ) .* cos( c_geocentric(:,1) );
        c_cartesian(:,2) = cos( c_geocentric(:,2) ) .* sin( c_geocentric(:,1) );
        c_cartesian(:,3) = sin( c_geocentric(:,2) );

    end
