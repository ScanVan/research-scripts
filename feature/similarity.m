
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

    function similarity( s_path )

        % read features %
        [ s_feature, s_list ] = similarity_read( s_path );

        % extract listing length %
        s_size = size( s_list, 1 );

        % initialise correlation matrix %
        s_matrix = zeros( s_size, s_size );

        % parsing listing %
        for s_i = 1 : s_size

            % import features %
            s_rfeat = s_feature{s_i};

            % display information %
            fprintf( 2, 'Reference : %s\n', s_list(s_i).name );

            % parsing listing %
            for s_j = s_i + 1 : s_size

                % import features %
                s_cfeat = s_feature{s_j};

                % similarity measure %
                s_matrix( s_j, s_i ) = similarity_distance( s_rfeat, s_cfeat, 64 );


            end

        end

        %%% research %%% display matrix
        %figure;
        %imagesc( [1:s_size], [1:s_size], s_matrix' );
        %%% research %%%

    end

    function [ s_feature, s_list ] = similarity_read( s_path )

        % create listing %
        s_list = dir( [ s_path '/output/1_features/*' ] );

        % parsing listing %
        for s_i = 1 : size( s_list, 1 )

            % display information %
            fprintf( 2, 'Reading %s ...\n', s_list(s_i).name );

            % import features %
            s_feature{s_i} = dlmread( [ s_path '/output/1_features/' s_list(s_i).name ] );

        end

    end

    function s_measure = similarity_distance( s_rfeat, s_cfeat, s_prec )

        % initialise value %
        s_measure = 0.0;

        % create random selection %
        s_rand = randperm( size( s_rfeat, 1 ) )(1:s_prec);

        % parsing reference features %
        for s_i = 1 : s_prec

            % compute distances component %
            s_distx = s_cfeat(:,1) - s_rfeat(s_rand(s_i),1);
            s_disty = s_cfeat(:,2) - s_rfeat(s_rand(s_i),2);

            % compute distance %
            s_dist = min( sqrt( s_distx .* s_distx + s_disty .* s_disty ) );

            % detect maximal distance %
            if ( s_dist > s_measure ); s_measure = s_dist; end

        end

    end

