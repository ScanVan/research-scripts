

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

    function geodesy( g_path )

        % compose source dataset path %
        g_source = [ g_path '/../../../' strtrim( fileread( [ g_path '/input/input_dataset.txt' ] ) ) ];

        % create segment listing %
        g_list = dir( [ g_path '/output/8_models_derive/*' ] );

        % create directory %
        mkdir( [ g_path '/output/9_geodesy_derive' ] );

        % parsing segment %
        for g_i = 1 : length( g_list )

            % create segment directory %
            g_import = [ g_path '/output/8_models_derive/' g_list(g_i).name ];

            % create image listing %
            g_image = dir( [ g_import '/image/*' ] );

            % check image count %
            if ( length( g_image ) >= 16 )

                % display information %
                fprintf( 2, 'Processing segment %s ...\n', g_list(g_i).name );

                % create segment path %
                g_export = [ g_path '/output/9_geodesy_derive/' g_list(g_i).name ];

                % create segment directory %
                mkdir( g_export );

                % align segment %
                geodesy_segment( g_path, g_import, g_export, g_source, g_image );

            end

        end

    end

    function geodesy_segment( g_path, g_import, g_export, g_source, g_list )

        % import track : gps %
        g_gtrack = geodesy_segment_track_gps( g_source, g_list );

        % import track : odometry %
        g_otrack = dlmread( [ g_import '/path.xyz' ] );

        % translate track : gps % temporary
        g_gtrack = g_gtrack - [ 4.3925e+06, 5.5620e+05, 4.5763e+06 ];

        % compute transformation %
        [ g_r, g_t, g_s ] = geodesy_align( g_gtrack, g_otrack );

        % apply transformation %
        g_otrack = geodesy_align_apply( g_otrack, g_r, g_t, g_s );

        % export aligned track %
        dlmwrite( [ g_export '/path.xyz' ], g_otrack, ' ' );

        % read sparse model %
        g_model = dlmread( [ g_import '/model.xyz' ] );

        % apply transformation %
        g_model = geodesy_align_apply( g_model, g_r, g_t, g_s );

        % export sparse model %
        dlmwrite( [ g_export '/model.xyz' ], g_model, ' ' );

        % create image directory %
        mkdir( [ g_export '/image' ] );

        % parsing image list %
        for g_i = 1 : length( g_list )

            % import absolute transformation %
            g_trans = dlmread( [ g_import '/image/' g_list(g_i).name ] );

            % apply transformation %
            g_trans(1:3,1:3) = g_r' * g_trans(1:3,1:3);

            % apply transformation %
            g_trans(1:3,4) = g_r' * ( ( g_trans(1:3,4) * g_s ) - g_t );

            % apply transformation %
            g_trans(1:3,5) = g_trans(1:3,5) * g_s;

            % export corrected transformation %
            dlmwrite( [ g_export '/image/' g_list(g_i).name ], g_trans, ' ' );

        end

    end

    function g_track = geodesy_segment_track_gps( g_source, g_list )

        % initialise memory %
        g_track = zeros( length( g_list ), 3 );

        % parsing image list %
        for g_i = 1 : length( g_list )

            % import track coordinates %
            g_gps = geodesy_segment_track_gps_read( [ g_source '/' g_list(g_i).name '.txt' ] );

            % compose track position - wgs84-msl %
            g_track(g_i,1:3) = g_gps(1:3);

            % compose track position - wgs84-hae %
            g_track(g_i,3) = g_track(g_i,3) + g_gps(4);

        end

        % convert coordinates %
        g_track = geodesy_segment_track_gps_convert( g_track, 6378137.0, 298.257223563 );

    end

    function g_gps = geodesy_segment_track_gps_read( g_file )

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

    function g_cart = geodesy_segment_track_gps_convert( g_pos, g_rad, g_flat )

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

    function [ g_r, g_t, g_s ] = geodesy_align( g_ref, g_pts )

        % compute scale factor %
        g_s = geodesy_align_scale( g_ref, g_pts );

        % compute alignment tranformation %
        [ g_r, g_t ] = geodesy_align_rigid( g_ref, g_pts * g_s );

    end

    function g_s = geodesy_align_scale( g_ref, g_pts )

        % compute distance %
        g_ref_dist = norm( g_ref(1,:) - g_ref(end,:) );

        % compute distance %
        g_pts_dist = norm( g_pts(1,:) - g_pts(end,:) );

        % compute scale factor %
        g_s = g_ref_dist / g_pts_dist;

    end

    function [ g_r, g_t ] = geodesy_align_rigid( g_ref, g_pts )

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

    function g_pts = geodesy_align_apply( g_pts, g_r, g_t, g_s )

        % apply transformation %
        g_pts = ( g_r' * ( g_pts * g_s - g_t' )' )';

    end

