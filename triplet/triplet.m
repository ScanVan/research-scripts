
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

    function [ t_r12, t_r23, t_t12, t_t23, t_sparse ] = triplet( t_input, t_width, t_height, t_path )

        % display information %
        fprintf( 2, 'Pose estimation with %i common matches\n', size( t_input, 1 ) );

        % extract features directions %
        t_1_d = t_input(:,1:3);
        t_2_d = t_input(:,4:6);
        t_3_d = t_input(:,7:9);

        % compute initial radius %
        t_1_r = ones( size( t_1_d, 1 ), 1 );
        t_2_r = ones( size( t_2_d, 1 ), 1 );
        t_3_r = ones( size( t_3_d, 1 ), 1 );

        % iteration counter %
        t_iter = 1;

        % iteration flag %
        t_flag = true;

        % initialise error %
        t_c_err = +0;
        t_p_err = -1;

        % iteration loop %
        while ( t_flag == true );

            % compute features %
            t_1_f = triplet_feature( t_1_d, t_1_r );
            t_2_f = triplet_feature( t_2_d, t_2_r );
            t_3_f = triplet_feature( t_3_d, t_3_r );

            % compute pose estimation %
            [ t_r12, t_t12, t_r23, t_t23 ] = triplet_register( t_1_f, t_2_f, t_3_f );

            % compute common frame %
            [ t_1_p_, t_1_d_, t_2_p_, t_2_d_, t_3_p_, t_3_d_ ] = triplet_frame_on_sphere_1( t_1_d, t_2_d, t_3_d, t_r12, t_r23, t_t12, t_t23 );

            % compute radius correction and error %
            [ t_1_r, t_2_r, t_3_r, t_1_e, t_2_e, t_3_e ] = triplet_radius( t_1_p_, t_1_d_, t_2_p_, t_2_d_, t_3_p_, t_3_d_ );

            % compute error %
            t_c_err = triplet_error( t_1_e, t_2_e, t_3_e, t_t12, t_t23 );

            % check error variation %
            if ( abs( t_c_err - t_p_err ) < 1e-8 )

                % abort iteration %
                t_flag = false;

            else

                % push error %
                t_p_err = t_c_err;

                % compute triplet amplitude %
                t_amp = triplet_amplitude( t_1_p_, t_2_p_, t_3_p_ );

                % matches filtering %
                [ t_1_d, t_1_r, t_2_d, t_2_r, t_3_d, t_3_r ] = triplet_consistency( t_1_d, t_1_r, t_2_d, t_2_r, t_3_d, t_3_r, 5.0, t_1_e, t_2_e, t_3_e );

                % stability filtering %
                [ t_1_d, t_1_r, t_2_d, t_2_r, t_3_d, t_3_r ] = triplet_filter( t_1_d, t_1_r, t_2_d, t_2_r, t_3_d, t_3_r, t_amp, 1, 50 );

                % compute triplet characteristic scale %
                t_norm = norm( t_t12 ) + norm( t_t23 );

                % normalisation of radius %
                t_1_r = t_1_r / t_norm;
                t_2_r = t_2_r / t_norm;
                t_3_r = t_3_r / t_norm;

            end

            % display information %
            fprintf( 2, 'Iteration %03i : t_norm : %g, %g : with %i features : mean radius : (%g %g %g)\n', t_iter, norm( t_t12 ), norm( t_t23 ), size( t_1_d, 1 ), mean( t_1_r ), mean( t_2_r ), mean( t_3_r ) );

            % update iteration %
            t_iter = t_iter + 1;

        end

        % compute sparse model %
        t_sparse = triplet_model( t_1_d, t_1_r, t_2_d, t_2_r, t_3_d, t_3_r, t_r12, t_t12, t_r23, t_t23 );

    end

    function t_f = triplet_feature( t_d, t_r )

        % compute feature position %
        t_f(:,1) = t_d(:,1) .* t_r(:);
        t_f(:,2) = t_d(:,2) .* t_r(:);
        t_f(:,3) = t_d(:,3) .* t_r(:);

    end

    function t_rotate = triplet_rotate( t_p, t_r )

        % apply rotation %
        t_rotate = ( t_r * t_p' )';

    end

    function t_corr = triplet_correlation( t_1_f, t_1_c, t_2_f, t_2_c )

        % compute correlation matrix %
        t_corr = ( t_1_f - t_1_c )' * ( t_2_f - t_2_c );

    end

    function [ t_r12, t_t12, t_r23, t_t23 ] = triplet_register( t_1_f, t_2_f, t_3_f )

        % compute centroids %
        t_1_c = sum( t_1_f ) / size( t_1_f, 1 );
        t_2_c = sum( t_2_f ) / size( t_2_f, 1 );
        t_3_c = sum( t_3_f ) / size( t_3_f, 1 );

        % singular values decomposition %
        [ t_u t_s t_v ] = svd( triplet_correlation( t_1_f, t_1_c, t_2_f, t_2_c ) );

        % compute rotation matrix %
        t_r12 = t_v * t_u';

        % check rotation %
        if ( det( t_r12 ) < 0 )

            % display information %
            warning( 'inversion on 12' );

            % inversion correction %
            t_v(:,3) = -t_v(:,3); t_r12 = t_v * t_u';

        end

        % singular values decomposition %
        [ t_u t_s t_v ] = svd( triplet_correlation( t_2_f, t_2_c, t_3_f, t_3_c ) );

        % compute rotation matrix %
        t_r23 = t_v * t_u';

        % check rotation %
        if ( det( t_r23 ) < 0 )

            % display information %
            warning( 'inversion on 23' );

            % inversion correction %
            t_v(:,3) = -t_v(:,3); t_r23 = t_v * t_u';

        end

        % compute translation %
        t_t12 = t_2_c' - t_r12 * t_1_c';
        t_t23 = t_3_c' - t_r23 * t_2_c';

    end

    function [ t_1_p_, t_1_d_, t_2_p_, t_2_d_, t_3_p_, t_3_d_ ] = triplet_frame_on_sphere_1( t_1_d, t_2_d, t_3_d, t_r12, t_r23, t_t12, t_t23 )

        % common frame - second camera - direction %
        t_1_d_ = t_1_d;

        % common frame - second camera - direction %
        t_2_d_ = triplet_rotate( t_2_d, t_r12' );

        % common frame - second camera - direction %
        t_3_d_ = triplet_rotate( triplet_rotate( t_3_d, t_r23' ), t_r12' );

        % common frame - second camera - centers %
        t_1_p_ = + zeros( 3, 1 );

        % common frame - second camera - centers %
        t_2_p_ = - t_r12' * t_t12;

        % common frame - second camera - centers %
        t_3_p_ = + t_2_p_ - t_r12' * t_r23' * t_t23;

    end

    function [ t_1_p_, t_1_d_, t_2_p_, t_2_d_, t_3_p_, t_3_d_ ] = triplet_frame_on_sphere_2( t_1_d, t_2_d, t_3_d, t_r12, t_r23, t_t12, t_t23 )

        % common frame - second camera - direction %
        t_1_d_ = triplet_rotate( t_1_d, t_r12  );

        % common frame - second camera - direction %
        t_2_d_ = t_2_d;

        % common frame - second camera - direction %
        t_3_d_ = triplet_rotate( t_3_d, t_r23' );

        % common frame - second camera - centers %
        t_1_p_ = + t_t12;

        % common frame - second camera - centers %
        t_2_p_ = + zeros( 3, 1 );

        % common frame - second camera - centers %
        t_3_p_ = - t_r23' * t_t23;

    end

    function t_inter = triplet_intersect( t_1_p, t_1_d, t_2_p, t_2_d, t_3_p, t_3_d )

        % intermediate computation %
        t_w1 = eye(3,3) - t_1_d' * t_1_d;
        t_q1 = t_w1 * t_1_p;

        % intermediate computation %
        t_w2 = eye(3,3) - t_2_d' * t_2_d;
        t_q2 = t_w2 * t_2_p;

        % intermediate computation %
        t_w3 = eye(3,3) - t_3_d' * t_3_d;
        t_q3 = t_w3 * t_3_p;

        % compute best intersection point %
        t_inter = ( inv( t_w1 + t_w2 + t_w3 ) * ( t_q1 + t_q2 + t_q3 ) )';

    end

    function [ t_1_r, t_2_r, t_3_r, t_1_e, t_2_e, t_3_e ] = triplet_radius( t_1_p, t_1_d, t_2_p, t_2_d, t_3_p, t_3_d )

        % initialise memory %
        t_1_r = zeros( size( t_1_d, 1 ), 1 );
        t_2_r = zeros( size( t_2_d, 1 ), 1 );
        t_3_r = zeros( size( t_3_d, 1 ), 1 );

        % initialise memory %
        t_1_e = zeros( size( t_1_d, 1 ), 1 );
        t_2_e = zeros( size( t_2_d, 1 ), 1 );
        t_3_e = zeros( size( t_3_d, 1 ), 1 );

        % parsing features %
        for t_i = 1 : size( t_1_d, 1 )

            % compute optimised intersection %
            t_inter = triplet_intersect( t_1_p, t_1_d(t_i,:), t_2_p, t_2_d(t_i,:), t_3_p, t_3_d(t_i,:) );

            % compute radius correction %
            t_1_r(t_i) = dot( t_1_d(t_i,:), t_inter - t_1_p' );
            t_2_r(t_i) = dot( t_2_d(t_i,:), t_inter - t_2_p' );
            t_3_r(t_i) = dot( t_3_d(t_i,:), t_inter - t_3_p' );

            % compute intersection disparity %
            t_1_e(t_i) = norm( t_1_p' + t_1_d(t_i,:) * t_1_r(t_i) - t_inter );
            t_2_e(t_i) = norm( t_2_p' + t_2_d(t_i,:) * t_2_r(t_i) - t_inter );
            t_3_e(t_i) = norm( t_3_p' + t_3_d(t_i,:) * t_3_r(t_i) - t_inter );

        end

    end

    function t_error = triplet_error( t_1_e, t_2_e, t_3_e, t_t12, t_t23 )

        % compute error %
        t_error = max( [ max( t_1_e ), max( t_2_e ), max( t_3_e ) ] ) / min( norm( t_t12 ), norm( t_t23 ) );

    end

    function t_amplitude = triplet_amplitude( t_1_p, t_2_p, t_3_p )

        % compute and return amplitude %
        t_amplitude = max( [ norm( t_1_p - t_2_p ), norm( t_1_p - t_3_p ), norm( t_2_p - t_3_p ) ] );

    end

    function [ t_1_d_, t_1_r_, t_2_d_, t_2_r_, t_3_d_, t_3_r_ ] = triplet_consistency( t_1_d, t_1_r, t_2_d, t_2_r, t_3_d, t_3_r, t_tol, t_1_e, t_2_e, t_3_e )

        % compute statistics %
        t_1_s = std( t_1_e ) * t_tol;
        t_2_s = std( t_2_e ) * t_tol;
        t_3_s = std( t_3_e ) * t_tol;

        % indexation %
        t_index = 0;

        % parsing features %
        for t_i = 1 : size( t_1_d, 1 )

            if ( t_1_e(t_i) < t_1_s )
            if ( t_2_e(t_i) < t_2_s )
            if ( t_3_e(t_i) < t_3_s )

                % update index %
                t_index = t_index + 1;

                % select feature %
                t_1_r_( t_index ) = t_1_r( t_i );
                t_2_r_( t_index ) = t_2_r( t_i );
                t_3_r_( t_index ) = t_3_r( t_i );

                % select feature %
                t_1_d_( t_index, : ) = t_1_d( t_i, : );
                t_2_d_( t_index, : ) = t_2_d( t_i, : );
                t_3_d_( t_index, : ) = t_3_d( t_i, : );

            end
            end
            end

        end

    end

    function [ t_1_d_, t_1_r_, t_2_d_, t_2_r_, t_3_d_, t_3_r_ ] = triplet_filter( t_1_d, t_1_r, t_2_d, t_2_r, t_3_d, t_3_r, t_value, t_min, t_max )

        % range definition %
        t_min = t_min * t_value;
        t_max = t_max * t_value;

        % indexation parameter %
        t_index = 0;

        % parsing features %
        for t_i = 1 : size( t_1_d, 1 )

            % filtering condition %
            if ( t_1_r( t_i ) >= t_min )
            if ( t_1_r( t_i ) <= t_max )

                % update index %
                t_index = t_index + 1;

                % select feature %
                t_1_r_( t_index ) = t_1_r( t_i );
                t_2_r_( t_index ) = t_2_r( t_i );
                t_3_r_( t_index ) = t_3_r( t_i );

                % select feature %
                t_1_d_( t_index, : ) = t_1_d( t_i, : );
                t_2_d_( t_index, : ) = t_2_d( t_i, : );
                t_3_d_( t_index, : ) = t_3_d( t_i, : );

            end
            end

        end

    end

    function t_sparse = triplet_model( t_1_d, t_1_r, t_2_d, t_2_r, t_3_d, t_3_r, t_r12, t_t12, t_r23, t_t23 )

        % compute common frame %
        [ t_1_p_, t_1_d_, t_2_p_, t_2_d_, t_3_p_, t_3_d_ ] = triplet_frame_on_sphere_1( t_1_d, t_2_d, t_3_d, t_r12, t_r23, t_t12, t_t23 );

        % compute features %
        t_1_f = t_1_p_' + triplet_feature( t_1_d_, t_1_r );
        t_2_f = t_2_p_' + triplet_feature( t_2_d_, t_2_r );
        t_3_f = t_3_p_' + triplet_feature( t_3_d_, t_3_r );

        % initialise memory %
        t_sparse = zeros( size( t_1_d, 1 ), 3 );

        % parsing features %
        for t_i = 1 : size( t_1_d, 1 )

            % compute best intesection %
            t_sparse(t_i,:) = triplet_intersect( t_1_p_, t_1_d_(t_i,:), t_2_p_, t_2_d_(t_i,:), t_3_p_, t_3_d_(t_i,:) );

        end

    end

    %function triplet_export( t_1_d, t_1_r, t_2_d, t_2_r, t_3_d, t_3_r, t_r12, t_t12, t_r23, t_t23, t_output )
    %
    %   % compute common frame %
    %    [ t_1_p_, t_1_d_, t_2_p_, t_2_d_, t_3_p_, t_3_d_ ] = triplet_frame_on_sphere_1( t_1_d, t_2_d, t_3_d, t_r12, t_r23, t_t12, t_t23 );
    %
    %   % compute features %
    %    t_1_f = t_1_p_' + triplet_feature( t_1_d_, t_1_r );
    %    t_2_f = t_2_p_' + triplet_feature( t_2_d_, t_2_r );
    %    t_3_f = t_3_p_' + triplet_feature( t_3_d_, t_3_r );
    %
    %   % create output stream %
    %    t_f = fopen( [ t_output '.disparity.xyz' ], 'w' );
    %
    %   % parsing points %
    %    for t_i = 1 : size( t_1_d, 1 )
    %
    %       % export feature position %
    %        fprintf( t_f, '%g %g %g 255 0 0\n', t_1_f(t_i,:) );
    %
    %       % export feature position %
    %        fprintf( t_f, '%g %g %g 0 255 0\n', t_2_f(t_i,:) );
    %
    %       % export feature position %
    %        fprintf( t_f, '%g %g %g 0 0 255\n', t_3_f(t_i,:) );
    %
    %    end
    %
    %   % delete output stream %
    %    fclose( t_f );
    %
    %   % create output stream %
    %    t_f = fopen( [ t_output '.sparse.xyz' ], 'w' );
    %
    %   % parsing points %
    %    for t_i = 1 : size( t_1_d, 1 )
    %
    %       % compute best intesection %
    %        t_p = triplet_intersect( t_1_p_, t_1_d_(t_i,:), t_2_p_, t_2_d_(t_i,:), t_3_p_, t_3_d_(t_i,:) );
    %
    %       % export point %
    %        fprintf( t_f, '%g %g %g 192 192 192\n', t_p(1), t_p(2), t_p(3) );
    %
    %    end
    %
    %    % delete output stream %
    %    fclose( t_f );
    %
    %end

