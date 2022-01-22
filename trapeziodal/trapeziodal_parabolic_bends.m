clear all;
close all;
clc;


boundary_conditions(1).time = 0
boundary_conditions(1).position = 0;
boundary_conditions(1).velocity = 0;
boundary_conditions(1).acceleration = 0;
boundary_conditions(1).jerk = 0;

  
boundary_conditions(2).time = 6;
boundary_conditions(2).position = 30;
boundary_conditions(2).velocity = 0;
boundary_conditions(2).acceleration = 0;
boundary_conditions(2).jerk = 0;

constraints.accelTime = 2;
constraints.constantVel = 10;
constraints.maxAccel = 60;
constraints.maxVel = 10;


[q, v, a, j, t, coeff] = trapeziodal_composite(boundary_conditions, constraints, 'p');

#Finally plot the three profiles
plot_profile(q,v,a,j,t);