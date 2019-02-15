
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

    function duplet( d_match, d_width, d_height )

        % import triplet matches %
        d_data = dlmread( d_match );

        % compute initial direction %
        d_1_d = duplet_cartesian( d_data( :, 1:2 ), d_width, d_height );
        d_2_d = duplet_cartesian( d_data( :, 3:4 ), d_width, d_height );

        % compute initial radius %
        d_1_r = ones( size( d_1_d, 1 ), 1 );
        d_2_r = ones( size( d_2_d, 1 ), 1 );

        % iteration counter %
        d_iter = 1;

        % iteration flag %
        d_flag = true;

        % initialise error %
        d_c_err = +0;
        d_p_err = -1;

        % iteration loop %
        while ( d_flag == true );

            % compute features %
            d_1_f = duplet_feature( d_1_d, d_1_r );

            d_2_f = duplet_feature( d_2_d, d_2_r );

            % compute pose estimation %
            [ d_r12, d_t12 ] = duplet_register( d_1_f, d_2_f );

            % compute common frame %
            [ d_1_p_, d_1_d_, d_2_p_, d_2_d_ ] = duplet_frame_on_sphere_1( d_1_d, d_2_d, d_r12, d_t12 );

            % compute radius correction and error %
            [ d_1_r, d_2_r, d_1_e, d_2_e ] = duplet_radius( d_1_p_, d_1_d_, d_2_p_, d_2_d_ );

            % compute error %
            d_c_err = duplet_error( d_1_e, d_2_e, d_t12 );

            % check error variation %
            if ( abs( d_c_err - d_p_err ) < 1e-8 )

                % abort iteration %
                d_flag = false;

            else

                % push error %
                d_p_err = d_c_err;

                % matches filtering %
                %[ d_1_d, d_1_r, d_2_d, d_2_r ] = duplet_consistency( d_1_d, d_1_r, d_2_d, d_2_r, 5.0, d_1_e, d_2_e );

                % stability filtering %
                [ d_1_d, d_1_r, d_2_d, d_2_r ] = duplet_filter( d_1_d, d_1_r, d_2_d, d_2_r, 3.0 ); % 2.5

                % compute triplet characteristic scale %
                d_norm = norm( d_t12 );

                % normalisation of radius %
                d_1_r = d_1_r / d_norm;
                d_2_r = d_2_r / d_norm;

            end

            % display information %
            fprintf( 2, 'Iteration %03i : d_norm : %g : with %i features : mean radius : (%g %g)\n', d_iter, norm( d_t12 ), size( d_1_d, 1 ), mean( d_1_r ), mean( d_2_r ) );

            % update iteration %
            d_iter = d_iter + 1;

        end

        % display results %
        display( d_r12 );

        % display results %
        display( d_t12 );

        % export results %
        duplet_export( d_1_d, d_1_r, d_2_d, d_2_r, d_r12, d_t12, d_match );

    end

    function d_point = duplet_cartesian( d_match, d_width, d_height )

        % coordinates re-normalisation %
        d_match(:,1) = ( ( d_match(:,1) - 1 ) / d_width ) * 2.0 * pi;

        % coordinates re-normalisation %
        d_match(:,2) = ( 0.5 - ( d_match(:,2) / d_height ) ) * pi;

        % initialise memory %
        d_point = zeros( size( d_match, 1 ), 3 );

        % parsing matches %
        for d_i = 1 : size( d_match, 1 )

            % coordinates conversion %
            d_point( d_i, 1 ) = cos( d_match( d_i, 2 ) ) * cos( d_match( d_i, 1 ) );
            d_point( d_i, 2 ) = cos( d_match( d_i, 2 ) ) * sin( d_match( d_i, 1 ) );
            d_point( d_i, 3 ) = sin( d_match( d_i, 2 ) );

        end

    end

    function d_f = duplet_feature( d_d, d_r )

        % compute feature position %
        d_f(:,1) = d_r(:) .* d_d(:,1);
        d_f(:,2) = d_r(:) .* d_d(:,2);
        d_f(:,3) = d_r(:) .* d_d(:,3);

    end

    function d_rotate = duplet_rotate( d_p, d_r )

        % initialise memory %
        d_rotate = zeros( size( d_p ) );

        % parsing points %
        for d_i = 1 : size( d_p, 1 )

            % apply rotation %
            d_rotate( d_i, : ) = ( d_r * d_p( d_i, : )' )';

        end

    end

    function d_corr = duplet_correlation( d_1_f, d_1_c, d_2_f, d_2_c )

        % initialise memory %
        d_corr = zeros( 3, 3 );

        % parsing features %
        for d_i = 1 : size( d_1_f, 1 )

            % accumulate correlation matrix %
            d_corr = d_corr + ( d_1_f - d_1_c )' * ( d_2_f - d_2_c );

        end

    end

    function [ d_r12, d_t12 ] = duplet_register( d_1_f, d_2_f )

        % compute centroids %
        d_1_c = sum( d_1_f ) / size( d_1_f, 1 );
        d_2_c = sum( d_2_f ) / size( d_2_f, 1 );

        % singular values decomposition %
        [ d_u d_s d_v ] = svd( duplet_correlation( d_1_f, d_1_c, d_2_f, d_2_c ) );

        % compute rotation matrix %
        d_r12 = d_v * d_u';

        % check rotation %
        if ( det( d_r12 ) < 0 )

            % display information %
            warning( 'inversion on 12' );

            % inversion correction %
            d_v(:,3) = -d_v(:,3); d_r12 = d_v * d_u';

        end

        % compute translation %
        d_t12 = d_2_c' - d_r12 * d_1_c';

    end

    function [ d_1_p_, d_1_d_, d_2_p_, d_2_d_ ] = duplet_frame_on_sphere_1( d_1_d, d_2_d, d_r12, d_t12 )

        % common frame - second camera - direction %
        d_1_d_ = d_1_d;

        % common frame - second camera - direction %
        d_2_d_ = duplet_rotate( d_2_d, d_r12' );

        % common frame - second camera - centers %
        d_1_p_ = + zeros( 3, 1 );

        % common frame - second camera - centers %
        d_2_p_ = - d_r12' * d_t12;

    end

    function [ d_1_p_, d_1_d_, d_2_p_, d_2_d_ ] = duplet_frame_on_sphere_2( d_1_d, d_2_d, d_r12, d_t12 )

        % common frame - second camera - direction %
        d_1_d_ = duplet_rotate( d_1_d, d_r12  );

        % common frame - second camera - direction %
        d_2_d_ = d_2_d;

        % common frame - second camera - centers %
        d_1_p_ = + d_t12;

        % common frame - second camera - centers %
        d_2_p_ = + zeros( 3, 1 );

    end

    function [ d_1_r, d_2_r, d_inter ] = duplet_intersect( d_1_p, d_1_d, d_2_p, d_2_d )

        % intermediate computation %
        d_0_d = d_2_p - d_1_p;

        % intermediate computation %
        d_aa = dot( d_1_d, d_1_d );
        d_bb = dot( d_2_d, d_2_d );
        d_ab = dot( d_1_d, d_2_d );
        d_ac = dot( d_1_d, d_0_d );
        d_bc = dot( d_2_d, d_0_d );

        % intermediate computation %
        d_dn = d_aa * d_bb - d_ab * d_ab;

        % compute radius %
        d_1_r = ( - d_ab * d_bc + d_ac * d_bb ) / d_dn;
        d_2_r = ( + d_ab * d_ac - d_bc * d_aa ) / d_dn;

        % compute intersection %
        d_inter = 0.5 * ( d_1_p + d_1_d * d_1_r + d_2_p + d_2_d * d_2_r );

    end

    function [ d_1_r, d_2_r, d_1_e, d_2_e ] = duplet_radius( d_1_p, d_1_d, d_2_p, d_2_d )

        % initialise memory %
        d_1_e = zeros( size( d_1_d, 1 ), 1 );
        d_2_e = zeros( size( d_2_d, 1 ), 1 );

        % parsing features %
        for d_i = 1 : size( d_1_d, 1 )

            % compute optimised intersection %
            [ d_1_r(d_i), d_2_r(d_i), d_inter ] = duplet_intersect( d_1_p', d_1_d(d_i,:), d_2_p', d_2_d(d_i,:) );

            % compute intersection disparity %
            d_1_e(d_i) = norm( d_1_p' + d_1_d(d_i,:) * d_1_r(d_i) - d_inter );
            d_2_e(d_i) = norm( d_2_p' + d_2_d(d_i,:) * d_2_r(d_i) - d_inter );

        end

    end

    function d_error = duplet_error( d_1_e, d_2_e, d_t12 )

        % compute error %
        d_error = max( d_1_e .+ d_2_e ) / norm( d_t12 );

    end

    function [ d_1_d_, d_1_r_, d_2_d_, d_2_r_ ] = duplet_consistency( d_1_d, d_1_r, d_2_d, d_2_r, d_tol, d_1_e, d_2_e )

        % compute statistics %
        d_1_s = std( d_1_e ) * d_tol;

        % indexation %
        d_index = 0;

        % parsing features %
        for d_i = 1 : size( d_1_d, 1 )

            if ( d_1_e(d_i) < d_1_s )

                % update index %
                d_index = d_index + 1;

                % select feature %
                d_1_r_( d_index ) = d_1_r( d_i );
                d_2_r_( d_index ) = d_2_r( d_i );

                % select feature %
                d_1_d_( d_index, : ) = d_1_d( d_i, : );
                d_2_d_( d_index, : ) = d_2_d( d_i, : );

            end

        end

    end

    function [ d_1_d_, d_1_r_, d_2_d_, d_2_r_ ] = duplet_filter( d_1_d, d_1_r, d_2_d, d_2_r, d_tol )

        % statistics %
        d_1_m = mean( d_1_r );
        d_2_m = mean( d_2_r );

        % statistics %
        d_1_s = std( d_1_r ) * d_tol;
        d_2_s = std( d_2_r ) * d_tol;

        % indexation parameter %
        d_index = 0;

        % filtering process %
        for d_i = 1 : size( d_1_d, 1 )

            if ( abs( d_1_r(d_i) - d_1_m ) <= d_1_s )
            if ( abs( d_2_r(d_i) - d_2_m ) <= d_2_s )

                % update index %
                d_index = d_index + 1;

                % select feature %
                d_1_r_( d_index ) = d_1_r( d_i );
                d_2_r_( d_index ) = d_2_r( d_i );

                % select feature %
                d_1_d_( d_index, : ) = d_1_d( d_i, : );
                d_2_d_( d_index, : ) = d_2_d( d_i, : );

            end
            end

        end

    end

    function duplet_export( d_1_d, d_1_r, d_2_d, d_2_r, d_r12, d_t12, d_output )

        % compute common frame %
        [ d_1_p_, d_1_d_, d_2_p_, d_2_d_ ] = duplet_frame_on_sphere_1( d_1_d, d_2_d, d_r12, d_t12 );

        % compute features %
        d_1_f = d_1_p_' + duplet_feature( d_1_d_, d_1_r );
        d_2_f = d_2_p_' + duplet_feature( d_2_d_, d_2_r );

        % create output stream %
        d_f = fopen( [ d_output '.disparity.xyz' ], 'w' );

        % parsing points %
        for d_i = 1 : size( d_1_d, 1 )

            % export feature position %
            fprintf( d_f, '%g %g %g 255 0 0\n', d_1_f(d_i,:) );

            % export feature position %
            fprintf( d_f, '%g %g %g 0 255 0\n', d_2_f(d_i,:) );

        end

        % delete output stream %
        fclose( d_f );

        % create output stream %
        d_f = fopen( [ d_output '.sparse.xyz' ], 'w' );

        % parsing points %
        for d_i = 1 : size( d_1_d, 1 )

            % compute best intesection %
            [ d_r1, d_r2, d_p ] = duplet_intersect( d_1_p_', d_1_d_(d_i,:), d_2_p_', d_2_d_(d_i,:) );

            % export point %
            fprintf( d_f, '%g %g %g 192 192 192\n', d_p(1), d_p(2), d_p(3) );

        end

        % delete output stream %
        fclose( d_f );

    end

