
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

    function merge_sparse( m_path )

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
            [ m_index, m_vom, m_vop ] = merge_sparse_segment( m_path, m_list, m_index, m_segment );

            % export merged model %
            merge_sparse_export( m_segment, m_vom, m_vop );

            % update index %
            m_parse = m_parse + 1;

        end

    end

    function [ m_index, m_vom, m_vop ] = merge_sparse_segment( m_path, m_list, m_index, m_export )

        % push initial index %
        m_push = m_index;

        % cumulative matrix %
        m_rot = eye(3);

        % cumulative position %
        m_pos = zeros(3,1);

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

            % read estimated sparse  %
            m_model = dlmread( [ m_path '/output/6_sparse_3/' m_name ] );

            % read estimated pose %
            m_data = dlmread( [ m_path '/output/5_pose_3/' m_name ] );

            % extarct rotation 1-2 %
            m_r12 = m_data(1:3,1:3);

            % extract translation 1-2 %
            m_t12 = m_data(1:3,4);

            % extarct rotation 2-3 %
            m_r23 = m_data(1:3,5:7);

            % extract translation 2-3 %
            m_t23 = m_data(1:3,8);

            % compute scale factor %
            m_factor = m_scl / norm( m_t12 );

            % scale and transform translation %
            m_t12 = m_rot * ( m_t12 * m_factor );

            % scale translation %
            m_t23 = m_rot * ( m_t23 * m_factor );

            % scale and transform sparse %
            m_model = ( m_rot * ( m_model * m_factor )' )' + m_pos';

            % bootstrap state %
            if ( m_index == m_push )

                % update index %
                m_index = m_index + 1;

                % predict position %
                m_predict = m_pos - m_r12' * m_t12 - m_r12' * m_r23' * m_t23;

            else

                % compute position %
                m_check = m_pos - m_r12' * m_t12;

                % apply consistency check %
                if ( ( norm( m_check - m_predict ) / norm( m_t23 ) ) > 0.1 )

                    % abort incremental merge %
                    return;

                end

                % update index %
                m_index = m_index + 1;

                % predict position %
                m_predict = m_check - m_r12' * m_r23' * m_t23;

            end

            % create image link %
            merge_sparse_link( m_export, m_list(m_index).name, m_rot, m_pos, m_factor );

            % update odometry %
            m_vop = [ m_vop; m_pos' ];

            % update sparse model %
            m_vom = [ m_vom; m_model ];

            % update absolute rotation %
            m_rot = m_rot * m_r12';

            % update absolute translation %
            m_pos = m_pos - m_r12' * m_t12;

            % update absolute scale %
            m_scl = norm( m_t23 );

        end

    end

    function merge_sparse_link( m_path, m_image, m_r, m_t, m_f )

        % compose exportation matrix %
        m_transform = [ m_r, m_t, [ m_f; m_f; m_f ] ];

        % export link with transformation %
        dlmwrite( [ m_path '/image/' m_image ], m_transform, ' ' );

    end

    function merge_sparse_export( m_path, m_vom, m_vop )

        % export sparse model %
        dlmwrite( [ m_path '/model.xyz' ], m_vom, ' ' );

        % export sparse odometry %
        dlmwrite( [ m_path '/path.xyz' ], m_vop, ' ' );

    end


