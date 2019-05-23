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

    function match( m_path )

        % create listing %
        m_list = dir( [ m_path '/output/2_matches_moving_index/*' ] );

        % extract listing size %
        m_size = size( m_list, 1 );

        % parsing listing %
        for m_i = 1 : m_size

            % display message %
            fprintf( 2, 'Importing %s ...\n', m_list(m_i).name );

            % import matches index %
            m_mindex{m_i} = dlmread( [ m_path '/output/2_matches_moving_index/' m_list(m_i).name ] );

            % import subsequent features %
            m_feature{m_i} = dlmread( [ m_path '/output/1_features/' strsplit( m_list(m_i).name, '_' ){1} ] );

        end

        % parsing listing %
        for m_i = 1 : m_size

            % display message %
            fprintf( 2, 'Process %s ...\n', strsplit( m_list(m_i).name, '_' ){1} );

            % parsing index %
            for m_j = 1 : size( m_mindex{m_i}, 1 )

                % check value of feature index %
                if ( m_mindex{m_i}(m_j,1) >= 0 )

                    % bootstrap link array %
                    m_link = [ m_mindex{m_i}(m_j,1), m_mindex{m_i}(m_j,2) ];

                    % invalidate feature index %
                    m_mindex{m_i}(m_j,1) = -1;

                    % compute link array - recursive %
                    m_link = match_detect( m_mindex, m_i, m_j, m_link );

                    % export match %
                    %match_export( m_path, strsplit( m_list(m_i).name, '_' ){1}, m_feature, m_link, m_i );

                end

            end

        end

    end

    function m_link = match_detect( m_mindex, m_base, m_local, m_link )

        % check listing %
        if ( ( m_base + 1 ) >= size( m_mindex, 2 ) )

            % interrupt search %
            return;

        end

        % extract link %
        m_idx = m_mindex{m_base}(m_local,2);

        % detect index %
        m_jdx = find( m_mindex{m_base+1}(:,1) == m_idx, 1, 'first' );

        % check search results %
        if ( length( m_jdx ) > 0 )

            % push feature on link array %
            m_link = [ m_link, m_mindex{m_base+1}(m_jdx,2) ];

            % invalidate feature index %
            m_mindex{m_base+1}(m_jdx,1) = -1;

            % compute link array - recursive %
            m_link = match_detect( m_mindex, m_base + 1, m_jdx, m_link );

        end

    end

    %function match_export( m_path, m_basename, m_feature, m_link, m_base )
    %    m_dpath = [ m_path '/output/3_hybrid' ];
    %    mkdir( m_dpath );
    %    m_dpath = sprintf( '%s/output/3_hybrid/%s', m_path, m_basename );
    %    mkdir( m_dpath );
    %    m_dpath = sprintf( '%s/output/3_hybrid/%s/%i', m_path, m_basename, length( m_link ) );
    %    m_stream = fopen( m_dpath, 'a' );
    %    for m_i =  1 : length( m_link )
    %        fprintf( m_stream, '%f %f ', m_feature{m_base + m_i - 1}(m_link(m_i)+1,1), m_feature{m_base + m_i - 1}(m_link(m_i)+1,2) );
    %    end
    %    fprintf( m_stream, '\n' );
    %    fclose( m_stream );
    %end

