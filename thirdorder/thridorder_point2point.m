q0 = 0;  %Initial position mm
v0 = 0;  %Initial velocity mm/sec
v1 = 0;
t0 = 0;  %Start time in second
t1 = 2;  %End time in second
q1 = 50; %End position

clear all;
close all;

boundary_conditions(1).time = 0
boundary_conditions(1).position = 0;
boundary_conditions(1).velocity = 0;

  
boundary_conditions(2).time = 8;
boundary_conditions(2).position = 10;
boundary_conditions(2).velocity = 0;


#[q, a, v, t] = thirdorder_trajectory(q0,v0,t0,q1,v1,t1);

[q,v, a, j, t, coeff] = trajectory_generator(boundary_conditions, 3);

#Finally plot the three profiles
plot_profile(q,v,a,j,t);
