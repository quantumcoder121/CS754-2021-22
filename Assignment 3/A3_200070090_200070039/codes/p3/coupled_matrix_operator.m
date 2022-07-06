classdef coupled_matrix_operator
    properties
        theta1 = []; % first array of angles
        theta2 = []; % second array of angles
        m = 0; % number of observations per array
        n = 0; % number of unknowns per array
        sq_size = 0; % for DCT and Radon implementation
        is_transpose = 0; % tracks if the matrix is transpose of the original one
        n_angles = 0; % number of angles
    end
    methods
        % constructor
        function obj = coupled_matrix_operator(m, n, theta1, theta2, n_angles)
            obj.theta1 = theta1;
            obj.theta2 = theta2;
            obj.m = m;
            obj.n = n;
            obj.sq_size = round(sqrt(n));
            obj.n_angles = n_angles;
        end
        % multiplier
        function y = mtimes(obj, x)
            if obj.is_transpose == 0
                x1 = x(1:obj.n);
                x2 = x(obj.n + 1:end);
                y1 = AMatrix(x1, obj.theta1, obj.sq_size);
                y2 = AMatrix(x1, obj.theta2, obj.sq_size) + AMatrix(x2, obj.theta2, obj.sq_size);
                y1 = reshape(y1, 1, []);
                y2 = reshape(y2, 1, []);
                y = horzcat(y1, y2);
                y = reshape(y, [], 1);
            else
                x1 = x(1:obj.m);
                x2 = x(obj.m + 1:end);
                y1 = ATMatrix(x1, obj.theta1, obj.sq_size, obj.n_angles) + ATMatrix(x2, obj.theta2, obj.sq_size, obj.n_angles);
                y2 = ATMatrix(x2, obj.theta2, obj.sq_size, obj.n_angles);
                y1 = reshape(y1, 1, []);
                y2 = reshape(y2, 1, []);
                y = horzcat(y1, y2);
                y = reshape(y, [], 1);
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
                ret_m = 2 * obj.m;
                ret_n = 2 * obj.n;
            else
                ret_m = 2 * obj.n;
                ret_n = 2 * obj.m;
            end
        end
    end
end