q0 = 0;  %Initial position mm
v0 = 10;  %Initial velocity mm/sec
v1 = 14;
t0 = 0;  %Start time in second
t1 = 8;  %End time in second
tf = (t0+t1)/2; %Symmetric point
q1 = 10; %End position
h = q1 - q0;
T = t1-t0;
qf = (q0 + q1)/2;
t_sample = 0.01;    %Sample time in seconds
N = T/t_sample;
q = zeros(1,N);  %Position
f = [1:N];

%Second order polynomical q(t) = a0 + a1*(t - t0) + a2*(t-t0)^2; t = [t0,tf]
%qb(t) = a3 + a4*(t-tf) + a5*(t-tf)^2 From t = [tf,t1]
%Degree 2, Three coefficient, Three boundary conditions
% For each segment we need seperate three conditions hence total 6 conditions
% The two segments are symmetric in this case with tf as mid-point
%n degree polynomial, (n+1) boundary conditions
% t=t0,  a0 = q0, a1 = v0, t = t1 q = q1, tf = (t0 + t1)/2
%Note that if v0 != v1 then we have discontinuous velocity profile
t = zeros(1,N);
t(1,1) = t0;
v = zeros(1,N); %Velocity
a = zeros(1,N);  %Acceleration
j = zeros(1,N);  %Jerk
a2 = 2*(h-v0*T)/(T*T);

for i=1:N/2;
  q(1,i) = q0 + v0*t(1, i) + a2*t(1, i)*t(1, i);
  t(1, i+1) = t(1, i) + t_sample;
endfor
a4 = 2*h/T - v1;
display(N);
display(a4);
a5 = 2*(v1*T-h)/(T*T);
display(a5);
display(t(1,N/2+1));
va =  v0 + 4*(h-v0*T)*(tf-t0)/(T*T);
vb = (2*h)/T -v1;
display(va);
display(vb);
for i=N/2+1:N-1;
  tT = t(1, i) - tf;
  q(1,i) = qf + a4*tT + a5*tT*tT;
  t(1, i+1) = t(1, i) + t_sample;
endfor
tT = t(1, N) - tf;
q(1,N) = qf + a4*tT + a5*tT*tT; %Last sample position
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
plot_profile(q,v,a,j,t);
