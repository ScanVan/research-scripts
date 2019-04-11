
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

    function analyse( a_path )

        % decompose path %
        [ a_dir, a_name ] = fileparts( a_path );

        % display summary %
        fprintf( 2, 'Dataset analysis : %s\n', a_name );

        % compute image list %
        a_img_list = analyse_image_list( a_path );

        % display summary %
        fprintf( 2, 'Image(s) : %i\n', length( a_img_list ) );

        % compute segment list %
        a_seg_list = analyse_segment_list( a_path );

        % compute connected list %
        a_con_list = analyse_connected_list( a_path );

        % display summary %
        fprintf( 2, 'Connected segement(s) : %i over %i\n', length( a_con_list ), length( a_seg_list ) );

        % compute calibrated image %
        a_img_state = analyse_image_state( a_path, a_img_list, a_con_list );

        % compute calibrated image ratio %
        a_ratio = sum( a_img_state ) / length( a_img_list );

        % display summary %
        fprintf( 2, 'Calibrated image(s) %i over %i (%0.2f%%)\n', sum( a_img_state ), length( a_img_list ), a_ratio * 100.0 );

        % compute triplet list %
        a_tri_list = analyse_triplet_list( a_path );

        % compute triplet state %
        a_tri_state = analyse_triplet_state( a_tri_list, a_img_state );

        % compute calibrated triplet ratio %
        a_ratio = sum( a_tri_state ) / length( a_tri_list );

        % display summary %
        fprintf( 2, 'Calibrated triplet(s) %i over %i (%0.2f%%)\n', sum( a_tri_state ), length( a_tri_list ), a_ratio * 100.0 );

        % extract triplet estimation information %
        [ a_tri_iter, a_tri_ifeat, a_tri_ffeat ] = analyse_triplet_information( a_path, a_tri_list );

        % display plot %
        analyse_triplet_plot( a_tri_state, a_tri_iter, 'Iteration', 'triplet_iteration.png' );

        % display plot %
        analyse_triplet_plot( a_tri_state, a_tri_ifeat, 'Features (initial)', 'triplet_ifeatures.png' );

        % display plot %
        analyse_triplet_plot( a_tri_state, a_tri_ffeat, 'Featrues (final)', 'triplet_ffeatures.png' );

    end

    function a_img_list = analyse_image_list( a_path )

        % comptue image listing %
        a_img_list = dir( [ a_path '/output/1_features/*' ] );

    end

    function a_seg_list = analyse_segment_list( a_path )

        % compute segment listing %
        a_seg_list = dir( [ a_path '/output/8_models_derive/' ] );

        % remove dot and dot dot elements %
        a_seg_list= a_seg_list( ~ ismember( { a_seg_list.name }, { '.' , '..' } ) );

    end

    function a_con_list = analyse_connected_list( a_path )

        % compute connected parts listing %
        a_con_list = dir( [ a_path '/output/9_geodesy_derive/' ] );

        % remove dot and dot dot elements %
        a_con_list= a_con_list( ~ ismember( { a_con_list.name }, { '.' , '..' } ) );

    end

    function a_img_state = analyse_image_state( a_path, a_img_list, a_con_list )

        % initialise image state %
        a_img_state = zeros( length( a_img_list ), 1 );

        % parsing image listing %
        for a_i = 1 : length( a_img_list )

            % initialise index %
            a_j = 1;

            % parsing connected segment %
            while ( a_j <= length( a_con_list ) )

                % check image existance %
                if ( exist( [ a_path '/output/8_models_derive/' a_con_list(a_j).name '/image/' a_img_list(a_i).name ], 'file' ) == 2 )

                    % update image state %
                    a_img_state(a_i) = 1;

                    % cancel search %
                    a_j = length( a_con_list );

                end

                % update index %
                a_j = a_j + 1;

            end

        end

    end

    function a_tri_list = analyse_triplet_list( a_path )

        % compute triplet list %
        a_tri_list = dir( [ a_path '/output/3_triplets/*' ] );

    end

    function a_tri_state = analyse_triplet_state( a_tri_list, a_img_state )

        % initialise triplet list %
        a_tri_state = zeros( length( a_tri_list ), 1 );

        % initialise index %
        a_k = 1;

        % parsing triplet list %
        for a_i = 1 : length( a_tri_list )

            % check triplet state %
            if ( sum( a_img_state(a_k:a_k+2) ) == 3 )

                % update triplet state %
                a_tri_state(a_i) = 1;

            end

            % update index %
            a_k = a_k + 1;

        end

    end

    function [ a_tri_iter, a_tri_ifeat, a_tri_ffeat ] = analyse_triplet_information( a_path, a_tri_list )

        % initialise arrays %
        a_tri_iter = zeros( length( a_tri_list ), 1 );

        % initialise arrays %
        a_tri_ifeat = zeros( length( a_tri_list ), 1 );
        a_tri_ffeat = zeros( length( a_tri_list ), 1 );

        % parsing triplet list %
        for a_i = 1 : length( a_tri_list )

            % read pose estimation result file %
            a_content = textread( [ a_path '/output/5_pose_3/' a_tri_list(a_i).name ], '%s' );

            % extract information %
            a_tri_iter(a_i) = str2num( cell2mat( a_content(29) ) );

            % extract information %
            a_tri_ifeat(a_i) = str2num( cell2mat( a_content(35) ) );

            % extract information %
            a_tri_ffeat(a_i) = str2num( cell2mat( a_content(40) ) );

        end

    end

    function analyse_triplet_plot( a_tri_state, a_tri_information, a_label, a_export )

        % create figure %
        figure;

        % subplot configuration %
        subplot(6,1,[1:5]);

        % figure configuration %
        hold on;
        grid on;
        box  on;

        % display information %
        plot( [1:length(a_tri_state)], a_tri_information, '-r' );

        % axis label %
        xlabel( 'Triplets' );
        ylabel( a_label );

        % axis configuration %
        axis( [ 1 length(a_tri_state) 0 max(a_tri_information ) ] );

        % subplot configuration %
        subplot(6,1,6);

        % figure configuration %
        hold on;
        grid on;
        box  on;

        % display states %
        imagesc( a_tri_state' ); colormap( [ [ 0, 0, 0 ]; [ 1, 0, 0 ] ] );

        % axis configuration %
        axis( [ 1 length(a_tri_state) 0 1 ] + 0.5 );

        % check exportation %
        if ( exist( 'a_export', 'var' ) )

            % export figure %
            print( '-dpng', '-r300', a_export );

        end

    end

