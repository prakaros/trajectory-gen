q0 = 2;  %Initial position mm
v0 = 4;  %Initial velocity mm/sec
t0 = 0;  %Start time in second
t1 = 2;  %End time in second
q1 = 100; %End position
clear all;
close all;
clc;

h = q1 - q0;
T = t1-t0;
t_sample = 0.01;    %Sample time in seconds
N = T/t_sample;
q = zeros(1,N);  %Position
f = [1:N];

%First order polynomical q(t) = a0 + a1*(t - t0) + a2*(t-t0)^2;
%Degree 2, Three coefficient, Three boundary conditions
%n degree polynomial, (n+1) boundary conditions
% t=t0,  a0 = q0, a1 = v0, t = t1 q = q1
t = zeros(1,N);
t(1,1) = t0;
v = zeros(1,N); %Velocity
a = zeros(1,N);  %Acceleration
a2 = (h-v0*T)/(T*T);

for i=1:N-1;
  tT = t(1,i) - t0;
  q(1,i) = q0 + v0*tT + a2*tT*tT;
  t(1, i+1) = t(1, i) + t_sample;
endfor

tT = t(1,N) - t0;
q(1,N) = q0 + v0*tT + a2*tT*tT; %Last sample position
t_inv = 1/t_sample;   %Convert div operation into multiplication

%Calculate velocity from acceleration
for i=1:N-1;
  v(1,i) = q(1,i+1) - q(1,i);
  v(1,i) = v(1,i)*t_inv;
endfor
v(1,N) = v(1,N-1);  % Remove jump

%Calculate acceleration from velocity
for i=1:N-1;
  a(1,i) = v(1,i+1) - v(1,i);
  a(1,i) = a(1,i)*t_inv;
endfor
a(1,N-1) = a(1,N-2);  % Remove jump
a(1,N) = a(1,N-1);  % Remove jump

#Finally plot the three profiles
plot_profile(q,v,a,t);
