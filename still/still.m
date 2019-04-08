
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

    function still( s_path, s_size )

        % create listing %
        s_list = dir( [ s_path '/output/1_features/*' ] );

        % parsing listing %
        for s_i = 2 : size( s_list, 1 )

            % display image name %
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

        end

        % display possible criterion %
        still_show( s_crit, 'Composite normalised criterion', 'export.png' );

    end

    function still_show( s_crit, s_title, s_export )

        % create figure %
        figure;

        % figure configuration %
        hold on;
        grid on;
        box  on;

        % display criterion %
        plot( [2:size(s_crit,2)+1], s_crit, '-r' );

        % axis configuration %
        xlabel( 'Image index' );
        ylabel( 'Criterion' );

        % ad-hoc axis scaling (avoid spikes) %
        ylim( [ 0, 0.002 ] );

        % figure title %
        title( s_title );

        % check exportation %
        if ( exist( 's_export', 'var' ) )

            % export figure %
            print( '-dpng', '-r300', s_export );

        end

    end

