
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

    function merge( m_path, m_from, m_to )

        % create file list %
        m_pose = dir( [ m_path '/output/5_pose_3/*' ] );

        % create file list %
        m_sparse = dir( [ m_path '/output/6_sparse_3/*' ] );

        % initialise cumulative rotation matrix %
        m_rot = eye(3);

        % initialise cumulative position %
        m_pos = zeros(1,3);

        % initialise cumulative scale %
        m_scl = 1;

        % results - position array %
        m_vop = [];

        % results - model array %
        m_vom = [];

        % parsing files %
        for m_file = max( 1, m_from ) : min( size( m_pose, 1 ), m_to )

            % compute differential index - avoid ply file %
            m_index = ( ( m_file - 1 ) * 2 ) + 1;

            % check consistency %
            if ( m_sparse(m_index).bytes == 0 )

                % abort merging %
                break

            end

            % display information %
            fprintf( 2, 'merging %s ...\n', m_pose(m_file).name );

            % read estimated pose %
            m_data = dlmread( [ m_path '/output/5_pose_3/' m_pose(m_file).name ] );

            % extarct rotation 1-2 %
            m_r12 = m_data(1:3,1:3);

            % extract translation 1-2 %
            m_t12 = m_data(1:3,4)';

            % extarct rotation 1-2 %
            m_r23 = m_data(1:3,5:7);

            % extract translation 1-2 %
            m_t23 = m_data(1:3,8)';

            % extract sparse model - avoid ply files %
            m_model = dlmread( [ m_path '/output/6_sparse_3/' m_sparse(m_index).name ] );

            % compute scale factor %
            m_factor = m_scl / norm( m_t12 );

            % scale translation %
            m_t12 = m_t12 * m_factor;

            % scale translation %
            m_t23 = m_t23 * m_factor;

            % scale sparse model %
            m_model = m_model * m_factor;

            % compute position of sphere in triplet first sphere frame %
            [ m_p1, m_p2, m_p3 ] = merge_position( m_r12, m_t12, m_r23, m_t23 );

            % transfrom position %
            m_p1 = merge_rotation( m_p1, m_rot ) + m_pos;
            m_p2 = merge_rotation( m_p2, m_rot ) + m_pos;
            m_p3 = merge_rotation( m_p3, m_rot ) + m_pos;

            % transform model %
            m_model = merge_rotation( m_model, m_rot ) + m_pos;

            % store positions %
            m_vop = [ m_vop; m_p1; m_p2; m_p3 ];

            % store model %
            m_vom = [ m_vom; m_model ];

            % update cumulative matrix %
            m_rot = m_rot * ( m_r12' ); % to check %

            % update cumulative position %
            m_pos = m_p2;

            % update cumulative scale %
            m_scl = norm( m_t23 );

        end

        % create output stream %
        m_f = fopen( [ m_path '/output/7_odometry/nh_odometry.xyz' ], 'w' );

        % parsing model %
        for m_i = 1 : size( m_vom, 1 )

            % export scene point %
            fprintf( m_f, '%g %g %g 224 224 224\n', m_vom(m_i,:) );

        end

        % parsing position %
        for m_i = 1 : size( m_vop, 1 )

            % export position %
            fprintf( m_f, '%g %g %g 255 0 0\n', m_vop(m_i,:) );

        end

        % delete output stream %
        fclose( m_f );

        %figure;
        %hold on;
        %plot( m_vop(:,1), m_vop(:,2), '-ko', 'linewidth', 3 );
        %plot( m_vop(1:3,1), m_vop(1:3,2), '-xr' );
        %plot( m_vop(3:6,1), m_vop(3:6,2), '-ob' );

        %figure;
        %hold on;
        %plot3( m_vop(:,1), m_vop(:,2), m_vop(:,3), 'x-r' );
        %plot3( m_vom(:,1), m_vom(:,2), m_vom(:,3), '.k' );
        %axis( 'equal' );

    end

    function [ m_p1, m_p2, m_p3 ] = merge_position( m_r12, m_t12, m_r23, m_t23 )

        % compute triplet position of sphere one in frame of sphere one %
        m_p1 = zeros(1,3);

        % compute triplet position of sphere two in frame of sphere one %
        m_p2 = ( - m_r12' * m_t12' )';

        % compute triplet position of sphere three in frame of sphere one %
        m_p3 = ( m_p2' - m_r12' * m_r23' * m_t23' )';

    end

    function m_points = merge_rotation( m_points, m_rotation )

        % parsing point list %
        for m_i = 1 : size( m_points, 1 )

            % apply rotation %
            m_points(m_i,:) = ( m_rotation * m_points(m_i,:)' )';

        end

    end

