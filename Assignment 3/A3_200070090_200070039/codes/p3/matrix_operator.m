classdef matrix_operator
    properties
        theta = []; % array of angles
        m = 0; % number of observations
        n = 0; % number of unknowns
        sq_size = 0; % for DCT and Radon implementation
        is_transpose = 0; % tracks if the matrix is transpose of the original one
        n_angles = 0; % number of angles
    end
    methods
        % constructor
        function obj = matrix_operator(m, n, theta, n_angles)
            obj.theta = theta;
            obj.m = m;
            obj.n = n;
            obj.sq_size = round(sqrt(n));
            obj.n_angles = n_angles;
        end
        % multiplier
        function y = mtimes(obj, x)
            if obj.is_transpose == 0
                y = AMatrix(x, obj.theta, obj.sq_size);
            else
                y = ATMatrix(x, obj.theta, obj.sq_size, obj.n_angles);
            end
        end
        % transpose
        function new_obj = ctranspose(obj)
            if obj.is_transpose == 0
                obj.is_transpose = 1;
            else
                obj.is_transpose = 0;
            end
            new_obj = obj;
        end
        % size
        function [ret_m, ret_n] = size(obj)
            if obj.is_transpose == 0
                ret_m = obj.m;
                ret_n = obj.n;
            else
                ret_m = obj.n;
                ret_n = obj.m;
            end
        end
    end
end