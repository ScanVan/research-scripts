
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

            % process dense triplet %
            merge_dense_process( [ m_path '/output/10_dense_3_derive/' m_name '.xyz' ], [ m_path '/output/11_dense_derive/' m_index '.xyz' ], m_r, m_t, m_f );

        end

    end

    function merge_dense_process( m_input, m_output, m_r, m_t, m_f )

        % check dense portion %
        if ( dir( m_input ).bytes == 0 )

            % interrupt merge of dense portion %
            return;

        end

        % create output stream %
        m_stream = fopen( m_output, 'a' );

        % read dense portion %
        m_dense = dlmread( m_input );

        % apply scale factor %
        m_dense(:,1:3) = m_dense(:,1:3) * m_f;

        % parsing dense portion elements %
        for m_i = 1 : size( m_dense, 1 )

            % apply absolute rotation %
            m_dense(m_i,1:3) = ( m_r  * m_dense(m_i,1:3)' )';

            % apply absolute translation %
            m_dense(m_i,1:3) = ( m_t' + m_dense(m_i,1:3)  );

            % export element %
            fprintf( m_stream, '%.14e %.14e %.14e %i %i %i\n', m_dense(m_i,:) );

        end

        % delete output stream %
        fclose( m_stream );

    end

