clear all;
close all;

boundary_conditions(1).time = 0
boundary_conditions(1).position = 0;
boundary_conditions(1).velocity = 0;
boundary_conditions(1).acceleration = 0;
  
boundary_conditions(2).time = 8;
boundary_conditions(2).position = 10;
boundary_conditions(2).velocity = 0;
boundary_conditions(2).acceleration = 0;

[q, v, a, j, t, coeff] = trajectory_generator(boundary_conditions, 5);

#Finally plot the three profiles
plot_profile(q,v,a,j, t);
