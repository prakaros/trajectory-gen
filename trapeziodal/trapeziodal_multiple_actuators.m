clear all;
close all;
clc;

#All actuators start at same time and end at same time, the end time will be 
#computed from actuator with longest stroke.
actuator(1).startTime = 0;
actuator(1).endTime = 6;
actuator(1).startPosition = 0;
actuator(1).endPosition = 50;
actuator(1).velocity = 0;
actuator(1).acceleration = 0;
color_id{1} = '+r';

actuator(2).startTime = 0;
actuator(2).endTime = 9;
actuator(2).startPosition = 0;
actuator(2).endPosition = 20;
actuator(2).velocity = 0;
actuator(2).acceleration = 0;
color_id{2} = '-b';

#Note that contraints for maximum acceleration and maximum velocity is impose on 
#actuator with largest stroke
constraints.maxAccel = 20;
constraints.maxVel = 20;

[trajectory_list] = trapeziodal_actuators(actuator, constraints);

#Finally plot the three profiles
figure;
for i=1:length(trajectory_list);
  plot_profile(trajectory_list(i).q,trajectory_list(i).v,trajectory_list(i).a,trajectory_list(i).j,trajectory_list(i).t, color_id(i));
endfor
hold off;
