function [q, v, a, j, t, coeff] = thirdorder_trajectory(q0, v0, t0, q1, v1, t1)
T = t1 - t0;
h = q1 - q0;
t_sample = 0.01;    %Sample time in seconds
N = T/t_sample;
q = zeros(1,N);  %Position
f = [1:N];

%Third order polynomical 
%q(t) = a0 + a1*(t - t0) + a2*(t-t0)^2 + a3*(t-t0)^3; t = [t0,tf]
%Degree 3, Four coefficient, Four boundary conditions
%n degree polynomial, (n+1) boundary conditions.
% t=t0,  a0 = q0, a1 = v0, t = t1 q = q1, v=v1

t = zeros(1,N);
t(1,1) = t0;
v = zeros(1,N); %Velocity
a = zeros(1,N);  %Acceleration
coeff(1) = q0;
coeff(2) = v0;
coeff(3) = a2 = (3*h - T*(2*v0 + v1))/(T*T);
coeff(4) = a3 = (-2*h + T*(v0 + v1))/(T*T*T);

for i=1:N-1;
  tT = t(1,i) - t0;
  q(1,i) = q0 + v0*tT + a2*tT*tT + a3*tT*tT*tT;
  t(1, i+1) = t(1, i) + t_sample;
endfor
tT = t(1,N) - t0;
q(1,N) = q0 + v0*tT + a2*tT*tT + a3*tT*tT*tT;

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

for i=1:N-1;
  j(1,i) = a(1,i+1) - a(1,i);
  j(1,i) = j(1,i)*t_inv;
endfor

j(1,N-2) = j(1,N-3);  % Remove jump
j(1,N-1) = j(1,N-2);  % Remove jump
j(1,N) = j(1,N-1);  % Remove jump
endfunction
