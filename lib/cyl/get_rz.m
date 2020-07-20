% This module takes the basic geometry parameters like min/max
% r-coordinate, min/max z-coordinate, and delta_r, delta_z and outputs the
% following:
%
%                 r_i : array of integer radial coordinates
%                 r_h : array of half-integer radial coordinates
%                 z_i : array of integer z coordinates
%                 z_i : array of half-integer z coordinates


function [r_i, r_h, z_i, z_h] = get_rz(r_min,r_max,delta_r,...
                                       z_min,z_max,delta_z)



                                   
myeps = 1e-6 * eps;

r_i = r_min : delta_r + myeps : r_max + delta_r;

r_h = r_min + delta_r : delta_r + myeps : r_max + 2 * delta_r;



z_i = z_min : delta_z + myeps : z_max + delta_z;

z_h = z_min + delta_z : delta_z + myeps : z_max + 2 * delta_z;

