
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

    function merge_dense( m_path )

        % create connected segment listing %
        m_list = dir( [ m_path '/output/9_geodesy_derive/*' ] );

        % create directory %
        mkdir( [ m_path '/output/11_dense_derive' ] );

        % parsing listing %
        for m_parse = 1 : size( m_list, 1 )

            % select directory only %
            if ( exist( [ m_path '/output/9_geodesy_derive/' m_list(m_parse).name ], 'dir' ) == 7 )

                % merge dense segment %
                merge_dense_segment( m_path, m_list(m_parse).name );

            end

        end

    end

    function merge_dense_segment( m_path, m_index )

        % create image listing %
        m_list = dir( [ m_path '/output/8_models_derive/' m_index '/image/*' ] );

        % parsing image listing %
        for m_parse = 1 : size( m_list, 1 ) - 2

            % compose triplet name %
            m_name = [ m_list(m_parse).name '_' m_list(m_parse+1).name '_' m_list(m_parse+2).name ];

            % display information %
            fprintf( 2, 'processing : %s ...\n', m_name );

            % read triplet absolute transformation %
            m_transform = dlmread( [ m_path '/output/8_models_derive/' m_index '/image/' m_list(m_parse).name ] );

            % extract absolute rotation %
            m_r = m_transform(1:3,1:3);

            % extract absolute position %
            m_t = m_transform(1:3,4);

            % extract absolute factor %
            m_f = m_transform(1,5);

            % read dense triplet %
            m_dense = dlmread( [ m_path '/output/10_dense_3_derive/' m_name '.xyz' ] );

            % apply scale factor %
            m_dense(:,1:3) = m_dense(:,1:3) * m_f;

            % apply aboslute rotation %
            m_dense(:,1:3) = merge_dense_rotation( m_dense(:,1:3), m_r );

            % apply aboslute tranlsation %
            m_dense(:,1:3) = merge_dense_position( m_dense(:,1:3), m_t );

            % export aligned dense triplet %
            merge_dense_append( [ m_path '/output/11_dense_derive/' m_index '.xyz' ], m_dense );

        end

    end

    function m_points = merge_dense_rotation( m_points, m_rotation )

        % parsing point list %
        for m_i = 1 : size( m_points, 1 )

            % apply rotation %
            m_points(m_i,:) = ( m_rotation * m_points(m_i,:)' )';

        end

    end

    function m_points = merge_dense_position( m_points, m_position )

        % parsing point list %
        for m_i = 1 : size( m_points, 1 )

            % apply translation %
            m_points(m_i,:) = m_points(m_i,:) + m_position';

        end

    end

    function merge_dense_append( m_path, m_point )

        % create stream %
        m_stream = fopen( m_path, 'a' );

        % parsing point list %
        for m_i = 1 : size( m_point, 1 )

            % export dense points %
            fprintf( m_stream, '%f %f %f %i %i %i\n', m_point(m_i,:) );

        end

        % delete stream %
        fclose( m_stream );

    end

