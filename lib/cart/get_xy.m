% This module takes the basic geometry parameters like min/max
% r-coordinate, min/max z-coordinate, and delta_r, delta_z and outputs the
% following:
%
%                 r_i : array of integer radial coordinates
%                 r_h : array of half-integer radial coordinates
%                 z_i : array of integer z coordinates
%                 z_i : array of half-integer z coordinates


function [x_i, x_h, y_i, y_h] = get_xy(x_min,x_max,delta_x,...
                                       y_min,y_max,delta_y)



                                   
myeps = 1e-6 * eps;

x_i = x_min : delta_x + myeps : x_max + delta_x;

x_h = x_min + delta_x : delta_x + myeps : x_max + 2 * delta_x;



y_i = y_min : delta_y + myeps : y_max + delta_y;

y_h = y_min + delta_y : delta_y + myeps : y_max + 2 * delta_y;

