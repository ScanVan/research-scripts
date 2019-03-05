
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

    function odometry( o_path, o_width, o_height )

        % environment settings %
        addpath( '../duplet' );
        addpath( '../triplet' );
        addpath( '../merge' );
        addpath( '../core' );

        % stage : pose estimation %
        odometry_stage_pose_estimation( o_path, o_width, o_height );

        % stage : triplet incremental merge %
        odometry_stage_incremental_merge( o_path );

    end

    function odometry_stage_pose_estimation( o_path, o_width, o_height )

        % create listing %
        o_list = dir( [ o_path '/output/3_triplets/*' ] );

        % create directory %
        mkdir( [ o_path '/output/5_pose_3' ] );

        % create directory %
        mkdir( [ o_path '/output/6_sparse_3' ] );

        % parsing listing %
        for o_i = 1 : size( o_list, 1 )

            % read matched features %
            o_match = dlmread( [ o_path '/output/3_triplets/' o_list(o_i).name ] );

            % coordinates conversion : mapping to geocentric %
            o_match(:,1:2) = core_equirectangular_to_geocentric( o_match(:,1:2), o_width, o_height );
            o_match(:,3:4) = core_equirectangular_to_geocentric( o_match(:,3:4), o_width, o_height );
            o_match(:,5:6) = core_equirectangular_to_geocentric( o_match(:,5:6), o_width, o_height );

            % initialise memory %
            o_input = zeros( size( o_match, 1 ), 9 );

            % coordinates conversion : geocentric to unit sphere %
            o_input(:,1:3) = core_geocentric_to_cartesian( o_match(:,1:2) );
            o_input(:,4:6) = core_geocentric_to_cartesian( o_match(:,3:4) );
            o_input(:,7:9) = core_geocentric_to_cartesian( o_match(:,5:6) );

            % pose estimation %
            [ o_r12, o_r23, o_t12, o_t23, o_sparse ] = triplet( o_input, o_width, o_height );

            % compose compact result %
            o_pose = [ o_r12, o_t12, o_r23, o_t23 ];

            % export pose estimation %
            dlmwrite( [ o_path '/output/5_pose_3/' o_list(o_i).name ], o_pose, ' ' );

            % export sparse model %
            dlmwrite( [ o_path '/output/6_sparse_3/' o_list(o_i).name '.xyz' ], o_sparse, ' ' );

        end

    end

    function odometry_stage_incremental_merge( o_path )

        % incremental merge %
        merge( o_path );

    end
