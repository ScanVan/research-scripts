
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

    function merge( m_path )

        % create image listing %
        m_list = dir( [ m_path '/output/1_features/*' ] );

        % create directory %
        mkdir( [ m_path '/output/7_odometry' ] );

        % initialise index %
        m_index = 1;

        % merging process %
        while ( m_index < ( length( m_list ) - 2 ) )

            % merge segment %
            [ m_index, m_count, m_vom, m_vop ] = merge_segment( m_path, m_list, m_index );

            % check merge count %
            if ( m_count > 1 )

                % export merged model %
                merge_export( m_path, m_list, m_index - m_count, m_index + 1, m_vom, m_vop );

            end

        end

    end

    function [ m_index, m_count, m_vom, m_vop ] = merge_segment( m_path, m_list, m_index )

        % push initial index %
        m_push = m_index;

        % initialise merge count %
        m_count = 0;

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

        % parsing listing %
        while ( m_index < ( length( m_list ) - 2 ) )

            % compose triplet name %
            m_name = [ m_list(m_index).name '_' m_list(m_index+1).name '_' m_list(m_index+2).name ];

            % check consistency %
            if ( exist( [ m_path '/output/5_pose_3/' m_name ], 'file' ) == 2 )

                % display information %
                fprintf( 2, 'Processing : %s\n', m_name );

            else

                % abort incremental merge %
                break;

            end

            % read estimated pose %
            m_data = dlmread( [ m_path '/output/5_pose_3/' m_name ] );

            % extarct rotation 1-2 %
            m_r12 = m_data(1:3,1:3);

            % extract translation 1-2 %
            m_t12 = m_data(1:3,4)';

            % extarct rotation 1-2 %
            m_r23 = m_data(1:3,5:7);

            % extract translation 1-2 %
            m_t23 = m_data(1:3,8)';

            % extract sparse model  %
            m_model = dlmread( [ m_path '/output/6_sparse_3/' m_name '.xyz' ] );

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

            % avoid initial triplet %
            if ( length( m_vop ) > 0 )

                % compute position agreement %
                m_check = norm( m_p2 - m_vop(end,:) ) / norm( m_t23 );

                % check position agreement %
                if ( m_check < 0.1 )

                    % update index %
                    m_index = m_index + 1;

                    % update merge count %
                    m_count = m_count + 1;

                else

                    % abort incremental merge %
                    break;

                end

            else

                % update index %
                m_index = m_index + 1;

                % update merge count %
                m_count = m_count + 1;

            end

            % store positions %
            m_vop = [ m_vop; m_p1; m_p2; m_p3 ];

            % store model %
            m_vom = [ m_vom; m_model ];

            % update cumulative matrix %
            m_rot = m_rot * ( m_r12' );

            % update cumulative position %
            m_pos = m_p2;

            % update cumulative scale %
            m_scl = norm( m_t23 );

        end

    end

    function merge_export( m_path, m_list, m_start, m_stop, m_vom, m_vop )

        % create exportation name %
        m_name = [ m_path '/output/7_odometry/' m_list(m_start).name '_' m_list(m_stop).name ];

        % export model data %
        dlmwrite( [ m_name '.xyz' ], m_vom, ' ' );

        % export visual odometry %
        dlmwrite( [ m_name '-vo.xyz' ], m_vop, ' ' );

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

