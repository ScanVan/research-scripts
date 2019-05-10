
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

    function still( s_path, s_threshold, s_width, s_height )

        % create listing %
        s_list = dir( [ s_path '/output/1_features/*' ] );

        % initialise array %
        s_crit = zeros( size( s_list, 1 ), 1 );

        % initialise array %
        s_select = zeros( size( s_list, 1 ), 1 );

        % criterion accumulator %
        s_accum = 0.0;

        % creating output directory %
        mkdir( [ s_path '/output/2_1_selected_nh' ] );

        % parsing listing %
        for s_i = 2 : size( s_list, 1 )

            % display information %
            fprintf( 2, 'Checking %s ...\n', s_list(s_i).name );

            % compose matches file path %
            s_name = [ s_list(s_i-1).name '_' s_list(s_i).name ];

            % read matches %
            s_match = dlmread( [ s_path '/output/2_matches/' s_name ] );

            % compute distances %
            s_dist = still_distance( s_match, s_width, s_height );

            % compute criterion %
            s_crit(s_i) = mean( s_dist(:) ) * std( s_dist(:) );

            % accumulate criterion %
            s_accum = s_accum + s_crit(s_i);

            % apply cirterion %
            if ( s_accum > s_threshold )

                % set image as active in the dataset (not part of still area) %
                s_select(s_i) = 1;

                % create a link for the selected image %
                fclose( fopen( [ s_path '/output/2_1_selected_nh/' s_list(s_i-1).name '_' s_list(s_i).name ], 'w' ) );

                % reset accumulator %
                s_accum = 0.0;

            end

        end

        % display criterion %
        still_show( s_crit, s_select, s_threshold, 'export.png' );

    end

    function s_dist = still_distance( s_match, s_width, s_height )

        % compute features vectors %
        s_vect_a = still_cartesian( s_match(:,1:2), s_width, s_height );

        % compute features vectors %
        s_vect_b = still_cartesian( s_match(:,3:4), s_width, s_height );

        % compute dot products %
        s_dist = acos( dot( s_vect_a', s_vect_b' )' );

    end

    function s_point = still_cartesian( s_match, s_width, s_height )

        % coordinates re-normalisation %
        s_match(:,1) = ( ( s_match(:,1) - 1 ) / s_width ) * 2.0 * pi;

        % coordinates re-normalisation %
        s_match(:,2) = ( ( s_match(:,2) / s_height ) - 0.5 ) * pi;

        % initialise memory %
        s_point = zeros( size( s_match, 1 ), 3 );

        % coordinates conversion %
        s_point( :, 1 ) = cos( s_match( :, 2 ) ) .* cos( s_match( :, 1 ) );
        s_point( :, 2 ) = cos( s_match( :, 2 ) ) .* sin( s_match( :, 1 ) );
        s_point( :, 3 ) = sin( s_match( :, 2 ) );

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
        xlim( [ 1, s_size ] );

        % create subplot %
        subplot( 6, 1, 6 );

        % display threshold-based image selection mapping %
        imagesc( s_select' ); colormap( [ [ 0, 0, 0 ]; [ 1, 0, 0 ] ] );

        % axis configuration %
        axis( [ 0, s_size - 1, 0, 1 ] + 0.5 );


        % check exportation %
        if ( exist( 's_export', 'var' ) )

            % export figure %
            print( '-dpng', '-r300', s_export );

        end

    end

