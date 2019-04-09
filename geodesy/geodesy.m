

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

    function geodesy( g_root, g_path )

        % import source dataset path %
        g_source = [ g_root '/' strtrim( fileread( [ g_root '/' g_path '/input/input_dataset.txt' ] ) ) ];

        % create segment listing %
        g_list = dir( [ g_root '/' g_path '/output/8_models_derive/' ] );

        % parsing segment %
        for g_i = 1 : length( g_list )

            % avoid redirection %
            if ( strcmp( g_list(g_i).name, '.' ) == 1 )

                % continue parsing %
                continue

            end

            % avoid redirection %
            if ( strcmp( g_list(g_i).name, '..' ) == 1 )

                % continue parsing %
                continue

            end

            % process segment %
            geodesy_segment( [ g_root '/' g_path '/output/8_models_derive/' g_list(g_i).name ], g_source );

        end

    end

    function geodesy_segment( g_path, g_source )

        % create link listing %
        g_list = dir( [ g_path '/image/*' ] );

        % check segment size %
        if ( length( g_list ) < 16 )

            % abort process %
            return;

        end

        % import GPS information %
        g_gps = geodesy_segment_gps( g_path, g_source, g_list );

        % import segment track %
        dlmread( [ g_path '/path.xyz' ] );

        % convert GPS ellipsoidal to cartesian %
        g_gps = geodesy_cartesian( g_gps, 6378137.0, 298.257223563 );

        % extract reference %
        g_ref = g_gps(1,:);

        % express GPS to reference %
        g_gps(:,1) = g_gps(:,1) - g_ref(1);
        g_gps(:,2) = g_gps(:,2) - g_ref(2);
        g_gps(:,3) = g_gps(:,3) - g_ref(3);

    end

    function g_gps = geodesy_segment_gps( g_path, g_source, g_list )

        % initialise memory %
        g_gps = zeros( length( g_list ), 3 );

        % parsing image link %
        for g_i = 1 : length( g_list )

            % import gps information %
            g_raw = geodesy_segment_gps_get( [ g_source '/' g_list(g_i).name '.txt' ] );

            % compose wgs84 position %
            g_gps(g_i,1) = g_raw(1);
            g_gps(g_i,2) = g_raw(2);
            g_gps(g_i,3) = g_raw(3) + g_raw(4);

        end

    end

    function g_gps = geodesy_segment_gps_get( g_file )

        % import image information %
        g_content = textread( g_file, '%s' );

        % parsing content rows %
        for s_r = 1 : length( g_content )

            % check tag %
            if ( strcmp( g_content( s_r ), 'latitude:' ) == 1 )

                % extract latitude %
                g_gps(2) = str2num( cell2mat( g_content( s_r + 1 ) ) ) * ( pi / 180.0 );

            elseif ( strcmp( g_content( s_r ), 'longitude:' ) == 1 )

                % extract latitude %
                g_gps(1) = str2num( cell2mat( g_content( s_r + 1 ) ) ) * ( pi / 180.0 );

            elseif ( strcmp( g_content( s_r ), 'altitude:' ) == 1 )

                % extract latitude %
                g_gps(3) =  str2num( cell2mat( g_content( s_r + 1 ) ) );

            elseif ( strcmp( g_content( s_r ), 'geoid:' ) == 1 )

                % extract latitude %
                g_gps(4) =  str2num( cell2mat( g_content( s_r + 1 ) ) );

            end

        end

    end

    function g_cart = geodesy_cartesian( g_pos, g_rad, g_flat )

        % square eccentricity %
        g_e = 2 * ( 1.0 / g_flat ) - ( ( 1.0 / g_flat ) * ( 1.0 / g_flat ) );

        % trigonometric values %
        g_sin = sin( g_pos(:,2) );
        g_cos = cos( g_pos(:,2) );

        % normal curvature radii %
        g_r = g_rad ./ sqrt( 1.0 - g_e * g_sin .* g_sin );

        g_cart(:,1) = ( g_r + g_pos(:,3) ) .* g_cos .* cos( g_pos(:,1) );
        g_cart(:,2) = ( g_r + g_pos(:,3) ) .* g_cos .* sin( g_pos(:,1) );
        g_cart(:,3) = ( g_r * ( 1.0 - g_e ) + g_pos(:,3) ) .* g_sin;

    end

