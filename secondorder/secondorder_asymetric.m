q0 = 0;  %Initial position mm
v0 = 10;  %Initial velocity mm/sec
v1 = 16;
t0 = 0;  %Start time in second
t1 = 8;  %End time in second
tf = 4;
Ta = tf-t0;
Td = t1-tf;

q1 = 10; %End position
h = q1 - q0;
T = Ta+Td;
t_sample = 0.01;    %Sample time in seconds
N = T/t_sample;
q = zeros(1,N);  %Position
f = [1:N];

%Second order polynomical q(t) = a0 + a1*(t - t0) + a2*(t-t0)^2; t = [t0,tf]
%qb(t) = a3 + a4*(t-tf) + a5*(t-tf)^2 From t = [tf,t1]
%Degree 2, Three coefficient, Three boundary conditions
% For each segment we need seperate three conditions hence total 6 conditions
% The two segments are symmetric in this case with tf as mid-point
%n degree polynomial, (n+1) boundary conditions.
% t=t0,  a0 = q0, a1 = v0, t = t1 q = q1, v=v1, vb(t) = va(t) at t=tf
% Note that above we have different boundary conditions than in sample
% secondorder_traingle_trajectory.m
t = zeros(1,N);
t(1,1) = t0;
v = zeros(1,N); %Velocity
a = zeros(1,N);  %Acceleration
a2 = (2*h - v0*(T+Ta) - v1*Td)/(2*T*Ta);
Na = Ta/t_sample;

for i=1:Na;
  q(1,i) = q0 + v0*t(1, i) + a2*t(1, i)*t(1, i);
  t(1, i+1) = t(1, i) + t_sample;
endfor
a3 = (2*q1*Ta + Td*(2*q0 + Ta*(v0-v1)))/(2*T);

a4 = (2*h - v0*Ta - v1*Td)/T;
a5 = -(2*h - v0*Ta - v1*(T + Td))/(2*T*Td);

display(a5);

for i=Na+1:N-1;
  tT = t(1, i) - tf;
  q(1,i) = a3 + a4*tT + a5*tT*tT;
  t(1, i+1) = t(1, i) + t_sample;
endfor
tT = t(1, N) - tf;
q(1,N) = a3 + a4*tT + a5*tT*tT; %Last sample position
t_inv = 1/t_sample;   %Convert div operation into multiplication

%Calculate velocity from acceleration
for i=1:N-1;
  v(1,i) = q(1,i+1) - q(1,i);
  v(1,i) = v(1,i)*t_inv;
endfor
%v(1,N-1) = v(1,N-2);  % Remove jump
v(1,N) = v(1,N-1);  % Remove jump

%Calculate acceleration from velocity
for i=1:N-1;
  a(1,i) = v(1,i+1) - v(1,i);
  a(1,i) = a(1,i)*t_inv;
endfor
%a(1,N-2) = a(1,N-3);
a(1,N-1) = a(1,N-2);  % Remove jump
a(1,N) = a(1,N-1);  % Remove jump

#Finally plot the three profiles
plot_profile(q,v,a,t);
