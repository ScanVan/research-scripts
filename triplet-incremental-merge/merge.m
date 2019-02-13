
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

    function merge( t_path, t_from, t_to )

        % create file list %
        t_pose = dir( [ t_path '/output/5_pose_3/*' ] );

        % create file list %
        t_sparse = dir( [ t_path '/output/6_sparse_3/*' ] );

        % initialise cumulative rotation matrix %
        t_rot = eye(3);

        % initialise cumulative position %
        t_pos = zeros(1,3);

        % initialise cumulative scale %
        t_scl = 1;

        % results - position array %
        t_vop = [];

        % results - model array %
        t_vom = [];

        % parsing files %
        for t_file = max( 1, t_from ) : min( size( t_pose, 1 ), t_to )

            % compute differential index - avoid ply file %
            t_index = ( ( t_file - 1 ) * 2 ) + 1;

            % check consistency %
            if ( t_sparse(t_index).bytes == 0 )

                % abort merging %
                break

            end

            % display information %
            fprintf( 2, 'merging %s ...\n', t_pose(t_file).name );

            % read estimated pose %
            t_data = dlmread( [ t_path '/output/5_pose_3/' t_pose(t_file).name ] );

            % extarct rotation 1-2 %
            t_r12 = t_data(1:3,1:3);

            % extract translation 1-2 %
            t_t12 = t_data(1:3,4)';

            % extarct rotation 1-2 %
            t_r23 = t_data(1:3,5:7);

            % extract translation 1-2 %
            t_t23 = t_data(1:3,8)';

            % extract sparse model - avoid ply files %
            t_model = dlmread( [ t_path '/output/6_sparse_3/' t_sparse(t_index).name ] );

            % compute scale factor %
            t_factor = t_scl / norm( t_t12 );

            % scale translation %
            t_t12 = t_t12 * t_factor;

            % scale translation %
            t_t23 = t_t23 * t_factor;

            % scale sparse model %
            t_model = t_model * t_factor;

            % compute position of sphere in triplet first sphere frame %
            [ t_p1, t_p2, t_p3 ] = merge_position( t_r12, t_t12, t_r23, t_t23 );

            % transfrom position %
            t_p1 = merge_rotation( t_p1, t_rot ) + t_pos;
            t_p2 = merge_rotation( t_p2, t_rot ) + t_pos;
            t_p3 = merge_rotation( t_p3, t_rot ) + t_pos;

            % transform model %
            t_model = merge_rotation( t_model, t_rot ) + t_pos;

            % store positions %
            t_vop = [ t_vop; t_p1; t_p2; t_p3 ];

            % store model %
            t_vom = [ t_vom; t_model ];

            % update cumulative matrix %
            t_rot = t_rot * ( t_r12' ); % to check %

            % update cumulative position %
            t_pos = t_p2;

            % update cumulative scale %
            t_scl = norm( t_t23 );

        end

        % create output stream %
        t_f = fopen( [ t_path '/output/7_odometry/nh_odometry.xyz' ], 'w' );

        % parsing model %
        for t_i = 1 : size( t_vom, 1 )

            % export scene point %
            fprintf( t_f, '%g %g %g 224 224 224\n', t_vom(t_i,:) );

        end

        % parsing position %
        for t_i = 1 : size( t_vop, 1 )

            % export position %
            fprintf( t_f, '%g %g %g 255 0 0\n', t_vop(t_i,:) );

        end

        % delete output stream %
        fclose( t_f );

        %figure;
        %hold on;
        %plot( t_vop(:,1), t_vop(:,2), '-ko', 'linewidth', 3 );
        %plot( t_vop(1:3,1), t_vop(1:3,2), '-xr' );
        %plot( t_vop(3:6,1), t_vop(3:6,2), '-ob' );

        %figure;
        %hold on;
        %plot3( t_vop(:,1), t_vop(:,2), t_vop(:,3), 'x-r' );
        %plot3( t_vom(:,1), t_vom(:,2), t_vom(:,3), '.k' );
        %axis( 'equal' );

    end

    function [ t_p1, t_p2, t_p3 ] = merge_position( t_r12, t_t12, t_r23, t_t23 )

        % compute triplet position of sphere one in frame of sphere one %
        t_p1 = zeros(1,3);

        % compute triplet position of sphere two in frame of sphere one %
        t_p2 = ( - t_r12' * t_t12' )';

        % compute triplet position of sphere three in frame of sphere one %
        t_p3 = ( t_p2' - t_r12' * t_r23' * t_t23' )';

    end

    function t_points = merge_rotation( t_points, t_rotation )

        % parsing point list %
        for t_i = 1 : size( t_points, 1 )

            % apply rotation %
            t_points(t_i,:) = ( t_rotation * t_points(t_i,:)' )';

        end

    end

