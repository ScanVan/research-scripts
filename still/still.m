
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

    function still( s_path, s_size, s_threshold )

        % create listing %
        s_list = dir( [ s_path '/output/1_features/*' ] );

        % initialise array %
        s_crit = zeros( size( s_list, 1 ), 1 );

        % initialise array %
        s_select = zeros( size( s_list, 1 ), 1 );

        % parsing listing %
        for s_i = 2 : size( s_list, 1 )

            % display information %
            fprintf( 2, 'Checking %s ...\n', s_list(s_i).name );

            % compose matches file path %
            s_name = [ s_list(s_i-1).name '_' s_list(s_i).name ];

            % read matches %
            s_match = dlmread( [ s_path '/output/2_matches/' s_name ] ) / s_size;

            % compute distance components %
            s_comp_x = s_match(:,1) - s_match(:,3);
            s_comp_y = s_match(:,2) - s_match(:,4);

            % compute square components %
            s_comp_x = s_comp_x .* s_comp_x;
            s_comp_y = s_comp_y .* s_comp_y;

            % compute distances %
            s_dist = sqrt( s_comp_x + s_comp_y );

            % compute criterion %
            s_crit(s_i) = mean( s_dist(:) ) * std( s_dist(:) );

            % apply cirterion %
            if ( s_crit(s_i) > s_threshold )

                % set image as active in the dataset (not part of still area) %
                s_select(s_i) = 1;

            end

        end

        % display criterion %
        still_show( s_crit, s_select, s_threshold, 'export.png' );

    end

    function still_show( s_crit, s_select, s_threshold, s_export )

        % extract size %
        s_size = size( s_select, 1 );

        % create figure %
        figure;

        % create subplot %
        subplot( 6, 1, [ 1:5 ] );

        % figure configuration %
        hold on;
        grid on;
        box  on;

        % display criterion %
        plot( [1:s_size], s_crit, '-r' );

        % display threshold %
        plot( [1,s_size], [1,1] * s_threshold, 'b-' );

        % axis configuration %
        xlabel( 'Image index' );
        ylabel( 'Criterion' );

        % axis configuration %
        axis( [ 1, s_size, 0, 0.002 ] );

        % create subplot %
        subplot( 6, 1, 6 );

        % display selection function %
        s_area = area( [1:s_size], s_select );

        % area plot configuration %
        set( s_area(1), 'FaceColor', 'r' );

        % area plot configuration %
        set( s_area(1), 'EdgeColor', 'none' );

        % axis configuration %
        axis( [ 1, s_size, 0, 1 ] );


        % check exportation %
        if ( exist( 's_export', 'var' ) )

            % export figure %
            print( '-dpng', '-r300', s_export );

        end

    end

