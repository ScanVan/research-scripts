

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
            geodesy_segment( [ g_root '/' g_path ], [ g_root '/' g_path '/output/8_models_derive/' g_list(g_i).name ], g_list(g_i).name, g_source );

        end

    end

    function geodesy_segment( g_export, g_path, g_name, g_source )

        % create link listing %
        g_list = dir( [ g_path '/image/*' ] );

        % check segment size %
        if ( length( g_list ) < 16 )

            % abort process %
            return;

        end

        % compose exportation path %
        g_export = [ g_export '/output/9_geodesy_derive' ];

        % create directory %
        mkdir( g_export );

        % compose exportation path %
        g_export = [ g_export '/' g_name ];

        % create directory %
        mkdir( g_export );

        % import GPS information %
        g_gps = geodesy_segment_gps( g_path, g_source, g_list );

        % import segment track %
        g_track = dlmread( [ g_path '/path.xyz' ] );

        % convert GPS ellipsoidal to cartesian %
        g_gps = geodesy_cartesian( g_gps, 6378137.0, 298.257223563 );

        % temporary %
        g_gps(:,1) = g_gps(:,1) - 4.3925e+06;
        g_gps(:,2) = g_gps(:,2) - 5.5620e+05;
        g_gps(:,3) = g_gps(:,3) - 4.5763e+06;

        % compute alignement transformation %
        [ g_r, g_t, g_s ] = geodesy_align_detect( g_gps, g_track );

        % align visual odometry %
        g_track = geodesy_align( g_track, g_r, g_t, g_s );

        % export aligned visual odometry %
        dlmwrite( [ g_export '/path.xyz' ], g_track, ' ' );

        % import segment model %
        g_model = dlmread( [ g_path '/model.xyz' ] );

        % align model %
        g_model = geodesy_align( g_model, g_r, g_t, g_s );

        % export aligned model %
        dlmwrite( [ g_export '/model.xyz' ], g_model, ' ' );

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

    function [ g_r, g_t, g_s ] = geodesy_align_detect( g_ref, g_pts )

        % compute scale factor %
        g_s = geodesy_align_detect_scale( g_ref, g_pts );

        % compute alignment tranformation %
        [ g_r, g_t ] = geodesy_align_detect_rigid( g_ref, g_pts * g_s );

    end

    function g_s = geodesy_align_detect_scale( g_ref, g_pts )

        % compute distance %
        g_ref_dist = norm( g_ref(1,:) - g_ref(end,:) );

        % compute distance %
        g_pts_dist = norm( g_pts(1,:) - g_pts(end,:) );

        % compute scale factor %
        g_s = g_ref_dist / g_pts_dist;

    end

    function [ g_r, g_t ] = geodesy_align_detect_rigid( g_ref, g_pts )

        % memory management %
        g_m = zeros( 3, 3 );

        % points collection centroids computation %
        g_c = sum( g_ref ) / size( g_ref, 1 );
        g_d = sum( g_pts ) / size( g_pts, 1 );

        % compute estimator matrix %
        for i = 1 : size( g_ref, 1 ); g_m = g_m + ( g_ref(i,1:3) - g_c(1:3) )' * ( g_pts(i,1:3) - g_d(1:3) ); end

        % singular values decomposition of estimator matrix %
        [ g_U g_S g_V ] = svd( g_m );

        % compute rotation matrix %
        g_r = g_V * g_U';

        % check for reflections %
        if ( det( g_r ) < 0 ); g_V(:,3) = -g_V(:,3); g_r = g_V * g_U'; end

        % compute translation vector %
        g_t = g_d' - g_r * g_c';

    end

    function g_align = geodesy_align( g_pts, g_r, g_t, g_s )

        % apply scale factor %
        g_pts = g_pts * g_s;

        % apply translation %
        g_pts(:,1) = g_pts(:,1) - g_t(1);
        g_pts(:,2) = g_pts(:,2) - g_t(2);
        g_pts(:,3) = g_pts(:,3) - g_t(3);

        % apply transformation %
        g_align(:,1) = g_pts(:,1) * g_r'(1,1) + g_pts(:,2) * g_r'(1,2) + g_pts(:,3) * g_r'(1,3);
        g_align(:,2) = g_pts(:,1) * g_r'(2,1) + g_pts(:,2) * g_r'(2,2) + g_pts(:,3) * g_r'(2,3);
        g_align(:,3) = g_pts(:,1) * g_r'(3,1) + g_pts(:,2) * g_r'(3,2) + g_pts(:,3) * g_r'(3,3);

    end

