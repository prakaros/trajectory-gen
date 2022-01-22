q0 = 2;  %Initial position mm
v0 = 4;  %Initial velocity mm/sec
t0 = 0;  %Start time in second
t1 = 2;  %End time in second
t_sample = 0.05;    %Sample time in seconds
T = t1-t0;
N = T/t_sample;
q = zeros(1,N);  %Position
f = [1:N];

%First order polynomical q(t) = a0 + a1*(t - t0);
%Degree 1, two coefficient, two boundary conditions
% t=t0,  a0 = q0, a1 = v0
t = zeros(1,N);
t(1,1) = t0;
v = zeros(1,N); %Velocity
a = zeros(1,N);  %Acceleration
for i=1:N-1;
  q(1,i) = q0 + v0*t(1, i);
  t(1, i+1) = t(1, i) + t_sample;
endfor
q(1,N) = q0 + v0*t(1,N);  %Last sample position
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
a(1,N) = a(1,N-1);  % Remove jump

#Finally plot the three profiles
plot_profile(q,v,a,t);
