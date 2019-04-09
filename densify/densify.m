
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

    function densify( d_path, d_estimation, d_f21u_path, d_f21v_path, d_f23u_path, d_f23v_path, d_mask_path )

        % display information %
        fprintf( 2, 'Import pose estimation\n' );

        % import pose estimation %
        [ d_r12, d_t12, d_r23, d_t23 ] = densify_read_estimation( [ d_path '/' d_estimation ] );

        % display information %
        fprintf( 2, 'Import flow\n' );

        % import flow %
        d_f21u = dlmread( [ d_path '/' d_f21u_path ] );
        d_f21v = dlmread( [ d_path '/' d_f21v_path ] );

        % import flow %
        d_f23u = dlmread( [ d_path '/' d_f23u_path ] );
        d_f23v = dlmread( [ d_path '/' d_f23v_path ] );

        % display information %
        fprintf( 2, 'Import mask\n' );

        % import mask %
        d_mask = imread( [ d_path '/' d_mask_path ] );

        % display information %
        fprintf( 2, 'Extract matches\n' );

        % create matches set %
        d_match = densify_match( d_f21u, d_f21v, d_f23u, d_f23v, d_mask );

        %figure; hold on; plot( d_match(:,1), d_match(:,2), '.r' ); set( gca, 'ydir', 'reverse' );

    end

    function [ d_rot12, d_tra12, d_rot23, d_tra23 ] = densify_read_estimation( d_estimation )

        % read pose estimation file %
        d_pose = dlmread( d_estimation );

        % extract rotation matrix %
        d_rot12 = d_pose(1:3,1:3);

        % extract translation vector %
        d_tra12 = d_pose(1:3,4);

        % extract rotation matrix %
        d_rot23 = d_pose(1:3,5:7);

        % extract translation vector %
        d_tra23 = d_pose(1:3,8);

    end

    function d_match = densify_match( d_f21u, d_f21v, d_f23u, d_f23v, d_mask )

        % initialise memory %
        d_match = zeros( prod( size( d_f21u ) ), 6 );

        % initialise index %
        d_k = 1;

        % parsing image %
        for d_y = 1 : size( d_f21u, 1 )

            % parsing image %
            for d_x = 1 : size( d_f21u, 2 )

                % check mask %
                if ( d_mask( d_y, d_x ) > 0 )

                    % compute match coordinates %
                    d_match(d_k,1) = d_x + d_f21u(d_y,d_x);
                    d_match(d_k,2) = d_y + d_f21v(d_y,d_x);

                    % compute match coordinates %
                    d_match(d_k,3) = d_x;
                    d_match(d_k,4) = d_y;

                    % compute match coordinates %
                    d_match(d_k,5) = d_x + d_f23u(d_y,d_x);
                    d_match(d_k,6) = d_y + d_f23v(d_y,d_x);

                    % update index %
                    d_k = d_k + 1;

                end

            end

        end

    end

