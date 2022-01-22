clear all;
clc;
[wayPoints] = get_dlm_waypoints(true);
path = wayPoints;
display(wayPoints);
num_points = length(path);
display(num_points);

T = path(num_points).tT - path(1).tT;
h = path(num_points).q - path(1).q;
t_sample = 0.01;    %Sample time in seconds
N = T/t_sample;
display(N);
q = zeros(1,N);  %Position
f = [1:N];
prev_path_length = 0;
path_length = 0;
previous_path = path(1);
for i=1:num_points-1
#[q_path, v_path, a_path, j_path, t_path] = thirdorder_trajectory(previous_path.q, previous_path.v, previous_path.tT, path(i+1).q, path(i+1).v, path(i+1).tT);
bCond(1).time = previous_path.tT;
bCond(1).position = previous_path.q;
bCond(1).velocity = previous_path.v;

bCond(2).time = path(i+1).tT;
bCond(2).position = path(i+1).q;
bCond(2).velocity = path(i+1).v;

[q_path, v_path, a_path, j_path, t_path] = trajectory_generator(bCond, 3);

subpath_length = length(t_path);
path_length = subpath_length + path_length;
display(subpath_length);
display(path_length);
q(1 +  prev_path_length : path_length) = q_path;
a(1 +  prev_path_length : path_length) = a_path;
v(1 +  prev_path_length : path_length) = v_path;
j(1 +  prev_path_length : path_length) = j_path;
t(1 +  prev_path_length : path_length) = t_path;
prev_path_length = path_length;
previous_path = path(i+1);
endfor

#Finally plot the three profiles
plot_profile(q,v,a,j,t);
