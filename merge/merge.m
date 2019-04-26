
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
        mkdir( [ m_path '/output/8_models_derive/' ] );

        % initialise index %
        m_index = 1;

        % initialise index %
        m_parse = 1;

        % merging process %
        while ( m_index < ( length( m_list ) - 2 ) )

            % compose segment path %
            m_segment = [ m_path '/output/8_models_derive/' num2str( m_parse ) ];

            % create directory %
            mkdir( m_segment );

            % create directory %
            mkdir( [ m_segment '/image' ] );

            % merge segment %
            [ m_index, m_vom, m_vop ] = merge_segment( m_path, m_list, m_index, m_segment );

            % export merged model %
            merge_export( m_segment, m_vom, m_vop );

            % update index %
            m_parse = m_parse + 1;

        end

    end

    function [ m_index, m_vom, m_vop ] = merge_segment( m_path, m_list, m_index, m_export )

        % push initial index %
        m_push = m_index;

        % cumulative matrix %
        m_rot = eye(3);

        % cumulative position %
        m_pos = zeros(1,3);

        % cumulative scale %
        m_scl = 1;

        % initialise odometry %
        m_vop = [];

        % initialise model %
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
                error( 'triplet not found' );

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
            m_model = dlmread( [ m_path '/output/6_sparse_3/' m_name ] );

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

            % check bootstrap state %
            if ( m_index == m_push )

                % create links %
                merge_link( m_export, m_list(m_index  ).name );
                merge_link( m_export, m_list(m_index+1).name );
                merge_link( m_export, m_list(m_index+2).name );

                % update index %
                m_index = m_index + 1;

                % initialise position %
                m_vop = [ m_p1; m_p2; m_p3 ];

                % initialise model %
                m_vom = m_model;

            else

                % compute position agreement %
                m_check = norm( m_p2 - m_vop(end,:) ) / norm( m_t23 );

                % check position agreement %
                if ( m_check < 0.1 )

                    % create link %
                    merge_link( m_export, m_list(m_index+2).name );

                    % update index %
                    m_index = m_index + 1;

                    % update position %
                    m_vop = [ m_vop; m_p3 ];

                    % update model %
                    m_vom = [ m_vom; m_model ];

                else

                    % abort incremental merge %
                    return;

                end

            end

            % update cumulative matrix %
            m_rot = m_rot * ( m_r12' );

            % update cumulative position %
            m_pos = m_p2;

            % update cumulative scale %
            m_scl = norm( m_t23 );

        end

    end

    function merge_link( m_path, m_image )

        % create link %
        fclose( fopen( [ m_path '/image/' m_image ], 'w' ) );

    end

    function merge_link_rt( m_path, m_image, m_r, m_t )

        % compose exportation matrix %
        m_transform = [ m_r, m_t ];

        % export link with transformation %
        dlmwrite( [ m_path '/image/' m_image ], m_transform, ' ' );

    end

    function merge_export( m_path, m_vom, m_vop )

        % export model data %
        dlmwrite( [ m_path '/model.xyz' ], m_vom, ' ' );

        % export visual odometry %
        dlmwrite( [ m_path '/path.xyz' ], m_vop, ' ' );

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

