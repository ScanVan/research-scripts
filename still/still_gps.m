
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

    function still_gps( s_path )

        % create listing %
        s_list = dir( [ s_path '/*.txt' ] );

        s_pos = zeros( size( s_list, 1 ), 3 );

        % parsing listing %
        for s_i = 1 : size( s_list, 1 )

            % avoid calibration file %
            if ( strcmp( s_list(s_i).name, 'calibration.txt' ) == 0 )

                % display information %
                fprintf( 2, 'Import %s ...\n', s_list(s_i).name );

                % import file content %
                s_data = textread( [ s_path '/' s_list(s_i).name ], '%s' );

                % extract position %
                s_pos(s_i,:) = still_gps_get( s_data );

            end

        end

        % convert to cartesian coordinates %
        s_cart = still_gps_convert( s_pos, 6378137.0, 298.257223563 )

        % mesaure step %
        s_step = 32;

        % parsing position %
        for s_i = s_step + 1 : size( s_cart, 1 )

            % compute components %
            s_comp_x = s_pos(s_i,1) - s_pos(s_i-s_step,1);
            s_comp_y = s_pos(s_i,2) - s_pos(s_i-s_step,2);
            s_comp_z = s_pos(s_i,3) - s_pos(s_i-s_step,3);

            % compute distance %
            s_dist(s_i - 1) = sqrt( s_comp_x .* s_comp_x + s_comp_y .* s_comp_y + s_comp_z .* s_comp_z );

        end

        % display distances %
        still_gps_show( s_dist, 'GPS position distances', 'export.png' );

    end

    function s_pos = still_gps_get( s_content )

        % parsing content rows %
        for s_r = 1 : length( s_content )

            % check tag %
            if ( strcmp( s_content( s_r ), 'latitude:' ) == 1 )

                % extract latitude %
                s_pos(2) = str2num( cell2mat( s_content( s_r + 1 ) ) ) * ( pi / 180.0 );

            elseif ( strcmp( s_content( s_r ), 'longitude:' ) == 1 )

                % extract latitude %
                s_pos(1) = str2num( cell2mat( s_content( s_r + 1 ) ) ) * ( pi / 180.0 );

            elseif ( strcmp( s_content( s_r ), 'altitude:' ) == 1 )

                % extract latitude %
                s_pos(3) =  str2num( cell2mat( s_content( s_r + 1 ) ) );

            end

        end

    end

    function s_cart = still_gps_convert( s_pos, s_rad, s_flat )

        % square eccentricity %
        s_e = 2 * ( 1.0 / s_flat ) - ( ( 1.0 / s_flat ) * ( 1.0 / s_flat ) );

        % trigonometric values %
        s_sin = sin( s_pos(:,2) );
        s_cos = cos( s_pos(:,2) );

        % normal curvature radii %
        s_r = s_rad ./ sqrt( 1.0 - s_e * s_sin .* s_sin );

        s_cart(:,1) = ( s_r + s_pos(:,3) ) .* s_cos .* cos( s_pos(:,1) );
        s_cart(:,2) = ( s_r + s_pos(:,3) ) .* s_cos .* sin( s_pos(:,1) );
        s_cart(:,3) = ( s_r * ( 1.0 - s_e ) + s_pos(:,3) ) .* s_sin;

    end

    function still_gps_show( s_dist, s_title, s_export )

        % create figure %
        figure;

        % figure configuration %
        hold on;
        grid on;
        box  on;

        % display criterion %
        plot( [2:size(s_dist,2)+1], s_dist, '-b' );

        % axis configuration %
        xlabel( 'Image index' );
        ylabel( 'Distance [Cartesian Metre]' );

        % ad-hoc axis scaling (avoid bad GPS signal) %
        ylim( [ 0, 2 ] );

        % figure title %
        title( s_title );

        % check exportation %
        if ( exist( 's_export', 'var' ) )

            % export figure %
            print( '-dpng', '-r300', s_export );

        end

    end

