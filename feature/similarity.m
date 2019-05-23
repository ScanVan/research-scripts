
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

    function similarity( s_path, s_sample )

        % read features %
        [ s_feature, s_list ] = similarity_read( s_path );

        % extract listing length %
        s_size = size( s_list, 1 );

        % initialise correlation matrix %
        s_matrix = zeros( s_size, s_size );

        % parsing listing %
        for s_i = 1 : s_size

            % display information %
            fprintf( 2, 'Reference : %s\n', s_list(s_i).name );

            % parsing listing %
            for s_j = s_i + 1 : s_size

                % similarity measure %
                s_matrix( s_j, s_i ) = similarity_distance( s_feature{s_i}, s_feature{s_j}, s_sample );

            end

        end

        % create figure %
        figure;

        % figure configuration %
        hold on;
        grid on;
        box  on;

        % display experimental similarity matrix %
        imagesc( [1:s_size], [1:s_size], s_matrix' );

        % axis labels %
        xlabel( 'image index' );
        ylabel( 'image index' );

        % reverse y-axis %
        set( gca, 'ydir', 'reverse' );

        % axis configuration %
        axis( [ 1.5, (s_size - 0.5), 1.5, (s_size - 0.5) ], 'square' );

    end

    function [ s_feature, s_list ] = similarity_read( s_path )

        % create listing %
        s_list = dir( [ s_path '/output/1_features/*' ] );

        s_list = s_list(1:700);

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

        % parsing reference features %
        for s_i = 1 : s_prec

            % compute distances component %
            s_dist = s_cfeat(:,:) - s_rfeat(s_i,:);

            % compute distance and search extremum %
            s_measure = max( s_measure, min( sqrt( s_dist(:,1) .* s_dist(:,1) + s_dist(:,2) .* s_dist(:,2) ) ) );

        end

    end

