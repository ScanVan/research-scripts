
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

    % Notes : This densification process uses the result of the optical flows
    % performed on three successive images :
    %
    %   https://people.csail.mit.edu/celiu/OpticalFlow/
    %
    % The optical flows are provided as arrays (mapping the images) to the script
    % that uses them to compute dense matches set for the image triplet.
    %
    % The opical flows have to be computed following : image_2 -> image_1 and
    % image_2 -> image_3.
    %
    % Research : Point to keep in mind for further developments :
    %
    % - RGB image vs Grayscale image
    % - Blur or adaptative blur on image to improve optical flow
    % - Need a better filtering process (spherical epipolar + other)
    % - Expend histogram of image to the maximum range (+ blur ?)

    function densify( d_path, d_estimation, d_f21u_path, d_f21v_path, d_f23u_path, d_f23v_path, d_mask_path )

        % display information %
        fprintf( 2, 'Import pose estimation\n' );

        % import pose estimation %
        [ d_r12, d_t12, d_r23, d_t23 ] = densify_read_estimation( [ d_path '/' d_estimation ] );

        % display information %
        fprintf( 2, 'Import flow\n' );

        % import flow %
        d_f21u = dlmread( [ d_path '/' d_f21u_path ] );
        d_f21v = dlmread( [ d_path '/' d_f21v_path ] );

        % import flow %
        d_f23u = dlmread( [ d_path '/' d_f23u_path ] );
        d_f23v = dlmread( [ d_path '/' d_f23v_path ] );

        % display information %
        fprintf( 2, 'Import mask\n' );

        % import mask %
        d_mask = imread( [ d_path '/' d_mask_path ] );

        % display information %
        fprintf( 2, 'Extract matches\n' );

        % create matches set %
        d_match = densify_match( d_f21u, d_f21v, d_f23u, d_f23v, d_mask );

        % display information %
        fprintf( 2, 'Convert matches\n' );

        % convert matches %
        d_vector_1 = densify_cartesian( d_match(:,1:2), size( d_f21u, 2 ), size( d_f21u, 1 ) );
        d_vector_2 = densify_cartesian( d_match(:,3:4), size( d_f21u, 2 ), size( d_f21u, 1 ) );
        d_vector_3 = densify_cartesian( d_match(:,5:6), size( d_f21u, 2 ), size( d_f21u, 1 ) );

        % display information %
        fprintf( 2, 'Alignment to frame of sphere 1\n' );

        % expresse scene in sphere 1 frame %
        [ d_1_p, d_1_d, d_2_p, d_2_d, d_3_p, d_3_d ] = densify_frame_on_sphere_1( d_vector_1, d_vector_2, d_vector_3, d_r12, d_r23, d_t12, d_t23 );

        % display information %
        fprintf( 2, 'Compute scene and error\n' );

        % compute scene %
        [ d_1_r, d_2_r, d_3_r, d_1_e, d_2_e, d_3_e, d_scene ] = densify_compute( d_1_p, d_1_d, d_2_p, d_2_d, d_3_p, d_3_d );

        % display information %
        fprintf( 2, 'Filter scene\n' );

        % filter scene %
        d_scene = densify_filter( d_scene, d_1_e, d_2_e, d_3_e, d_t12, d_t23 );

        % export scene %
        dlmwrite( 'scene.xyz', d_scene, ' ' );

        %figure; hold on; plot( d_1_e, '-r' ); plot( d_2_e, '-g' ); plot( d_3_e, '-b' );
        %figure; hold on; plot3( d_scene(:,1), d_scene(:,2), d_scene(:,3), '.r' ); set( gca, 'ydir', 'reverse' );
        %figure; hold on; plot( d_match(:,1), d_match(:,2), '.r' ); set( gca, 'ydir', 'reverse' );

    end

    function [ d_rot12, d_tra12, d_rot23, d_tra23 ] = densify_read_estimation( d_estimation )

        % read pose estimation file %
        d_pose = dlmread( d_estimation );

        % extract rotation matrix %
        d_rot12 = d_pose(1:3,1:3);

        % extract translation vector %
        d_tra12 = d_pose(1:3,4);

        % extract rotation matrix %
        d_rot23 = d_pose(1:3,5:7);

        % extract translation vector %
        d_tra23 = d_pose(1:3,8);

    end

    function d_match = densify_match( d_f21u, d_f21v, d_f23u, d_f23v, d_mask )

        % initialise memory %
        d_match = zeros( prod( size( d_f21u ) ), 6 );

        % initialise index %
        d_k = 1;

        % parsing image %
        for d_y = 1 : size( d_f21u, 1 )

            % parsing image %
            for d_x = 1 : size( d_f21u, 2 )

                % check mask %
                if ( d_mask( d_y, d_x ) > 0 )

                    % compute match coordinates %
                    d_match(d_k,1) = d_x + d_f21u(d_y,d_x);
                    d_match(d_k,2) = d_y + d_f21v(d_y,d_x);

                    % compute match coordinates %
                    d_match(d_k,3) = d_x;
                    d_match(d_k,4) = d_y;

                    % compute match coordinates %
                    d_match(d_k,5) = d_x + d_f23u(d_y,d_x);
                    d_match(d_k,6) = d_y + d_f23v(d_y,d_x);

                    % update index %
                    d_k = d_k + 1;

                end

            end

        end

    end

    function d_point = densify_cartesian( d_match, d_width, d_height )

        % coordinates re-normalisation %
        d_match(:,1) = ( ( d_match(:,1) - 1 ) / d_width ) * 2.0 * pi;

        % coordinates re-normalisation %
        d_match(:,2) = ( ( d_match(:,2) / d_height ) - 0.5 ) * pi;

        % initialise memory %
        d_point = zeros( size( d_match, 1 ), 3 );

        % coordinates conversion %
        d_point( :, 1 ) = cos( d_match( :, 2 ) ) .* cos( d_match( :, 1 ) );
        d_point( :, 2 ) = cos( d_match( :, 2 ) ) .* sin( d_match( :, 1 ) );
        d_point( :, 3 ) = sin( d_match( :, 2 ) );

    end

    function d_rotate = densify_rotate( d_p, d_r )

        % initialise memory %
        d_rotate = zeros( size( d_p ) );

        % parsing points %
        for d_i = 1 : size( d_p, 1 )

            % apply rotation %
            d_rotate( d_i, : ) = ( d_r * d_p( d_i, : )' )';

        end

    end

    function [ d_1_p_, d_1_d_, d_2_p_, d_2_d_, d_3_p_, d_3_d_ ] = densify_frame_on_sphere_1( d_1_d, d_2_d, d_3_d, d_r12, d_r23, d_t12, d_t23 )

        % common frame - second camera - direction %
        d_1_d_ = d_1_d;

        % common frame - second camera - direction %
        d_2_d_ = densify_rotate( d_2_d, d_r12' );

        % common frame - second camera - direction %
        d_3_d_ = densify_rotate( densify_rotate( d_3_d, d_r23' ), d_r12' );

        % common frame - second camera - centers %
        d_1_p_ = + zeros( 3, 1 );

        % common frame - second camera - centers %
        d_2_p_ = - d_r12' * d_t12;

        % common frame - second camera - centers %
        d_3_p_ = + d_2_p_ - d_r12' * d_r23' * d_t23;

    end

    function d_inter = densify_intersect( d_1_p, d_1_d, d_2_p, d_2_d, d_3_p, d_3_d )

        % intermediate computation %
        d_w1 = eye(3,3) - d_1_d' * d_1_d;
        d_q1 = d_w1 * d_1_p;

        % intermediate computation %
        d_w2 = eye(3,3) - d_2_d' * d_2_d;
        d_q2 = d_w2 * d_2_p;

        % intermediate computation %
        d_w3 = eye(3,3) - d_3_d' * d_3_d;
        d_q3 = d_w3 * d_3_p;

        % compute best intersection point %
        d_inter = ( inv( d_w1 + d_w2 + d_w3 ) * ( d_q1 + d_q2 + d_q3 ) )';

    end

    function [ d_1_r, d_2_r, d_3_r, d_1_e, d_2_e, d_3_e, d_scene ] = densify_compute( d_1_p, d_1_d, d_2_p, d_2_d, d_3_p, d_3_d )

        % initialise memory %
        d_1_e = zeros( size( d_1_d, 1 ), 1 );
        d_2_e = zeros( size( d_2_d, 1 ), 1 );
        d_3_e = zeros( size( d_3_d, 1 ), 1 );

        % initialise memory %
        d_scene = zeros( size( d_1_d, 1 ), 3 );

        % parsing features %
        for d_i = 1 : size( d_1_d, 1 )

            % compute optimised intersection %
            d_inter = densify_intersect( d_1_p, d_1_d(d_i,:), d_2_p, d_2_d(d_i,:), d_3_p, d_3_d(d_i,:) );

            % push sence point %
            d_scene(d_i,1:3) = d_inter;

            % compute radius correction %
            d_1_r(d_i) = dot( d_1_d(d_i,:), d_inter - d_1_p' );
            d_2_r(d_i) = dot( d_2_d(d_i,:), d_inter - d_2_p' );
            d_3_r(d_i) = dot( d_3_d(d_i,:), d_inter - d_3_p' );

            % compute intersection disparity %
            d_1_e(d_i) = norm( d_1_p' + d_1_d(d_i,:) * d_1_r(d_i) - d_inter );
            d_2_e(d_i) = norm( d_2_p' + d_2_d(d_i,:) * d_2_r(d_i) - d_inter );
            d_3_e(d_i) = norm( d_3_p' + d_3_d(d_i,:) * d_3_r(d_i) - d_inter );

        end

    end

    function d_filter = densify_filter( d_scene, d_1_e, d_2_e, d_3_e, d_t12, d_t23 )

        % compute naive condition %
        d_norm = ( 0.5 * ( norm(d_t12) + norm(d_t23) ) ) * 0.05;

        % initialise index %
        d_k = 1;

        % parsing scene %
        for d_i = 1 : size( d_scene, 1 )

            % filtering condition %
            if ( d_1_e(d_i) < d_norm )
            if ( d_2_e(d_i) < d_norm )
            if ( d_3_e(d_i) < d_norm )

                % select scene point %
                d_filter(d_k,:) = d_scene(d_i,:);

                % update index %
                d_k = d_k + 1;

            end
            end
            end

        end

    end


