function [q, v, a, j, t, coeff] = trajectory_generator(boundary_conditions, order)
  check_trajectory_order(order);
  check_boundary_conditions(boundary_conditions, order);
  
  if(order == 3)
    [q, v, a, j, t, coeff] = thirdorder_trajectory(boundary_conditions);
  endif
  
  if(order == 5)
    [q, v, a, j, t, coeff] = fifthorder_trajectory(boundary_conditions);
  endif
  
  if(order == 7)
    [q, v, a, j, t, coeff] = seventhorder_trajectory(boundary_conditions);
  endif
  
  if(order == 'c')
    [q, v, a, j, t, coeff] = circular_trajectory(boundary_conditions);
  endif
  
  if(order == 'cy')
    [q, v, a, j, t, coeff] = cycloid_trajectory(boundary_conditions);
  endif
  
  if(order == 'e')
    [q, v, a, j, t, coeff] = elliptical_trajectory(boundary_conditions);
  endif
  
  if(order == 'tr')
    [q, v, a, j, t, coeff] = trapeziod_with_circular_bends_trajectory(boundary_conditions);
  endif
  
endfunction

function [q, v, a, j, t, coeff] = thirdorder_trajectory(boundary_conditions)
  num_points = length(boundary_conditions);
  display(num_points);
  t0 = boundary_conditions(1).time;
  q0 = boundary_conditions(1).position;
  v0 = boundary_conditions(1).velocity;
 
  t1 = boundary_conditions(num_points).time;
  q1 = boundary_conditions(num_points).position;
  v1 = boundary_conditions(num_points).velocity;
 
  [T, h, N, q, v, a, j, t, t_sample] = initial_parameter(t0, t1, q0, q1);
  
%Third order polynomical 
%q(t) = a0 + a1*(t - t0) + a2*(t-t0)^2 + a3*(t-t0)^3; t = [t0,tf]
%Degree 3, Four coefficients, Four boundary conditions
%n degree polynomial, (n+1) boundary conditions.
% t=t0,  q(t0) = q0, v(t0) = v0, t = t1 q(t1) = q1, v(t1) = v1
% a0 = q0, a1 = v0, a2 = (3*h - T*(2*v0 + v1))/(T^2)
% a3 = (-2*h + T*(v0 + v1))/(T^3)

coeff(1) = q0;
coeff(2) = v0;
coeff(3) = a2 = (3*h - T*(2*v0 + v1))/(T^2);
coeff(4) = a3 = (-2*h + T*(v0 + v1))/(T^3);

for i=1:N-1;
  tT = t(1,i) - t0;
  q(1,i) = q0 + v0*tT + a2*tT^2 + a3*tT^3;
  t(1, i+1) = t(1, i) + t_sample;
endfor
tT = t(1,N) - t0;
q(1,N) = q0 + v0*tT + a2*tT^2 + a3*tT^3;

[v, a, j] = calculate_velacceljerk(q, t_sample);
endfunction

function [q, v, a, j, t, coeff] = fifthorder_trajectory(boundary_conditions)
  num_points = length(boundary_conditions);
  display(num_points);
  t0 = boundary_conditions(1).time;
  q0 = boundary_conditions(1).position;
  v0 = boundary_conditions(1).velocity;
  accel0 = boundary_conditions(1).acceleration;
  
  t1 = boundary_conditions(num_points).time;
  q1 = boundary_conditions(num_points).position;
  v1 = boundary_conditions(num_points).velocity;
  accel1 = boundary_conditions(num_points).acceleration;
 
  [T, h, N, q, v, a, j, t, t_sample] = initial_parameter(t0, t1, q0, q1);
%Fifth order polynomical 
%q(t) = a0 + a1*(t - t0) + a2*(t-t0)^2 + a3*(t-t0)^3 + a4*(t-t0)^4 + a5*(t-t0)^5; t = [t0,tf]
%Degree 5, Six coefficient, Six boundary conditions
%n degree polynomial, (n+1) boundary conditions.
% t=t0,  q(t0) = q0, v(t0) = v0, a(t0) = a0, t = t1 q(t1) = q1, v(t1) = v1, a(t1) = a1
% a0 = q0, a1 = v0, a2 = a0/2, a3 = (20*h - T*(8*v1 + 12*v0) - (3*accel0 - accel1)*T^2)/(2*T^3)
% a4 = (-30*h + T*(14*v1 + 16*v0) + (3*accel0 - 2*accel1 )*T^2)/(2*T^4)
% a5 = (12*h - 6*T*(v1 + v0) + (accel1 - accel0 )*T^2)/(2*T^5)

coeff(1) = q0;
coeff(2) = v0;
coeff(3) = a2 = accel0*0.5;
coeff(4) = a3 = (20*h - T*(8*v1 + 12*v0) - T*T*(3*accel0 - accel1 ))/(2*T^3);
coeff(5) = a4 = (-30*h + T*(14*v1 + 16*v0) + T*T*(3*accel0 - 2*accel1 ))/(2*T^4);
coeff(6) = a5 = (12*h - 6*T*(v1 + v0) + T*T*(accel1 - accel0 ))/(2*T^5);

for i=1:N-1;
  tT = t(1,i) - t0;
  q(1,i) = q0 + v0*tT + a2*tT^2 + a3*tT^3 + a4*tT^4 + a5*tT^5;
  t(1, i+1) = t(1, i) + t_sample;
endfor
tT = t(1,N) - t0;
q(1,N) = q0 + v0*tT + a2*tT^2 + a3*tT^3 + a4*tT^4 + a5*tT^5;

[v, a, j] = calculate_velacceljerk(q, t_sample);
endfunction

function [q, v, a, j, t, coeff] = seventhorder_trajectory(boundary_conditions)
  num_points = length(boundary_conditions);
  display(num_points);
  t0 = boundary_conditions(1).time;
  q0 = boundary_conditions(1).position;
  v0 = boundary_conditions(1).velocity;
  accel0 = boundary_conditions(1).acceleration;
  j0 = boundary_conditions(1).jerk;
  
  t1 = boundary_conditions(num_points).time;
  q1 = boundary_conditions(num_points).position;
  v1 = boundary_conditions(num_points).velocity;
  accel1 = boundary_conditions(num_points).acceleration;
  j1 = boundary_conditions(num_points).jerk;
 
  [T, h, N, q, v, a, j, t, t_sample] = initial_parameter(t0, t1, q0, q1);
%Seventh order polynomical 
%q(t) = a0 + a1*(t - t0) + a2*(t-t0)^2 + a3*(t-t0)^3 + a4*(t-t0)^4 + a5*(t-t0)^5 + a6*(t-t0)^6 + a7*(t-t0)^7; t = [t0,tf]
%Degree 7, Eight coefficients, Eight boundary conditions
%n degree polynomial, (n+1) boundary conditions.
% t=t0,  q(t0) = q0, v(t0) = v0, a(t0) = a0, j(t0) = j0, t = t1 q(t1) = q1, v(t1) = v1, a(t1) = a1, j(t1) = j1
% a0 = q0, a1 = v0, a2 = a0/2, a3 = j0/6, 
% a4 = (210*h - T*(T*(30*accel0 - 15*accel1) + (4*j0 + j1)*T^2 + 120*v0 + 90*v1))/(6*T^4)
% a5 = (-168*h + T*(T*(20*accel0 - 14*accel1) + (2*j0 + j1)*T^2 + 90*v0 + 78*v1))/(2*T^5)
% a6 = (420*h - T*(T*(45*accel0 - 39*accel1) + (4*j0 + 3*j1)*T^2 + 216*v0 + 204*v1))/(6*T^6)
% a7 = (-120*h + T*(T*(12*accel0 - 12*accel1) + (j0 + j1)*T^2 + 60*v0 + 60*v1))/(6*T^7)
coeff(1) = q0;
coeff(2) = v0;
coeff(3) = a2 = accel0*0.5;
coeff(4) = a3 = j0/6;
coeff(5) = a4 = (210*h - T*(T*(30*accel0 - 15*accel1) + (4*j0 + j1)*T^2 + 120*v0 + 90*v1))/(6*T^4);
coeff(6) = a5 = (-168*h + T*(T*(20*accel0 - 14*accel1) + (2*j0 + j1)*T^2 + 90*v0 + 78*v1))/(2*T^5);
coeff(7) = a6 = (420*h - T*(T*(45*accel0 - 39*accel1) + (4*j0 + 3*j1)*T^2 + 216*v0 + 204*v1))/(6*T^6);
coeff(8) = a7 = (-120*h + T*(T*(12*accel0 - 12*accel1) + (j0 + j1)*T^2 + 60*v0 + 60*v1))/(6*T^7);

for i=1:N-1;
  tT = t(1,i) - t0;
  q(1,i) = q0 + v0*tT + a2*tT^2 + a3*tT^3 + a4*tT^4 + a5*tT^5 +  a6*tT^6 +  a7*tT^7;
  t(1, i+1) = t(1, i) + t_sample;
endfor
tT = t(1,N) - t0;
q(1,N) = q0 + v0*tT + a2*tT^2 + a3*tT^3 + a4*tT^4 + a5*tT^5 +  a6*tT^6 +  a7*tT^7;

[v, a, j] = calculate_velacceljerk(q, t_sample);
endfunction

function [q, v, a, j, t, coeff] = circular_trajectory(boundary_conditions)
  num_points = length(boundary_conditions);
  display(num_points);
  t0 = boundary_conditions(1).time;
  q0 = boundary_conditions(1).position;

  t1 = boundary_conditions(num_points).time;
  q1 = boundary_conditions(num_points).position;

 
  [T, h, N, q, v, a, j, t, t_sample] = initial_parameter(t0, t1, q0, q1);
  
%circular trajectory 
%q(t) = q0 + h*0.5*(1-cos(pi*(t-t0)/T)) t = [t0, t1], h = q1-q0


coeff(1) = q0;
coeff(2) = h;
coeff(3) = T;

for i=1:N-1;
  tT = t(1,i) - t0;
  q(1,i) = q0 + h*0.5*(1-cos((pi*tT)/T));
  t(1, i+1) = t(1, i) + t_sample;
endfor
tT = t(1,N) - t0;
q(1,N) = q0 + h*0.5*(1-cos((pi*tT)/T));

[v, a, j] = calculate_velacceljerk(q, t_sample);
endfunction


function [q, v, a, j, t, coeff] = cycloid_trajectory(boundary_conditions)
  num_points = length(boundary_conditions);
  display(num_points);
  t0 = boundary_conditions(1).time;
  q0 = boundary_conditions(1).position;

  t1 = boundary_conditions(num_points).time;
  q1 = boundary_conditions(num_points).position;

 
  [T, h, N, q, v, a, j, t, t_sample] = initial_parameter(t0, t1, q0, q1);
  
%circular trajectory 
%q(t) = q0 + h*0.5*((t-t0)/T - (sin(2*pi*(t-t0)/T)/(2*pi))), t = [t0, t1], h = q1-q0


coeff(1) = q0;
coeff(2) = h;
coeff(3) = T;

for i=1:N-1;
  tT = t(1,i) - t0;
  q(1,i) = q0 + h*0.5*(tT/T-(sin(2*pi*tT/T)/(2*pi)));
  t(1, i+1) = t(1, i) + t_sample;
endfor
tT = t(1,N) - t0;
q(1,N) = q0 + h*0.5*(tT/T-(sin(2*pi*tT/T)/(2*pi)));

[v, a, j] = calculate_velacceljerk(q, t_sample);
endfunction

function [q, v, a, j, t, coeff] = elliptical_trajectory(boundary_conditions)
  num_points = length(boundary_conditions);
  display(num_points);
  t0 = boundary_conditions(1).time;
  q0 = boundary_conditions(1).position;

  t1 = boundary_conditions(num_points).time;
  q1 = boundary_conditions(num_points).position;

 
  [T, h, N, q, v, a, j, t, t_sample] = initial_parameter(t0, t1, q0, q1);
  
%circular trajectory 
%q(t) = q0 + h*0.5*(1 -  cos(pi*(t-t0)/T)/(sqrt(1-alpha*sin(pi*(t-t0)/T)^2)));, t = [t0, t1], h = q1-q0
%alpha = n^2-1/n^2; where n = a/b (Major Axis/Minor Axis)
majorAxis = 4;
minorAxis = 2;
n = majorAxis/minorAxis;
alpha = (n^2-1)/n^2;
coeff(1) = q0;
coeff(2) = h;
coeff(3) = T;

for i=1:N-1;
  tT = t(1,i) - t0;
  q(1,i) = q0 + h*0.5*(1 -  cos(pi*(tT)/T)/(sqrt(1-alpha*sin(pi*(tT)/T)^2)));
  t(1, i+1) = t(1, i) + t_sample;
endfor
tT = t(1,N) - t0;
q(1,N) = q0 + h*0.5*(1 -  cos(pi*(tT)/T)/(sqrt(1-alpha*sin(pi*(tT)/T)^2)));

[v, a, j] = calculate_velacceljerk(q, t_sample);
endfunction

function [q, v, a, j, t, coeff] = trapeziod_with_circular_bends_trajectory(boundary_conditions)
  num_points = length(boundary_conditions);
  display(num_points);
  t0 = boundary_conditions(1).time;
  q0 = boundary_conditions(1).position;

  t1 = boundary_conditions(num_points).time;
  q1 = boundary_conditions(num_points).position;

 
  [T, h, N, q, v, a, j, t, t_sample] = initial_parameter(t0, t1, q0, q1);
  
%Trapeziodal trajecttory with circular bends
% Three part of trajectory, 
% 1. Start part with circular trajectory 2. Constant velocity part 3. Stop part with circular trajectory 
% t = (t0, t0+Ta) where t0 is start time, Ta is circular trajectory time
% qa(t) = h*(1-sqrt(1-(t-t0)^2/(h^2)) + q0; 
% qc(t) = sqrt(h^2 - (t1-t)^2) + q0;  This is end part with circular trajectory
% qb(t) = a0 + a1(t-t0); Constant velocity mid part of trajectory
% Boundary conditions qa(t0+Ta) = qb(t0+Ta), va(t0+Ta) = vb(t0 + Ta)
% a1 = tan(alpha) , a1 = h*(cos(alpha) - 1)/cos(alpha) + q0
% alpha = acos((2*h^2 + T*sqrt(T^2 - 3*h^2))/(h^2 + T^2))
% T > sqrt(3)*h
%TO DO:: Fix discontinuity at t=t0 + Ta (Small dip in position)

invAngle = (2*(h*h) + T*sqrt(T*T - 3*(h*h)))/(h*h + T*T)
alpha = acos(invAngle)

Ta =1;  #Acceleration and Deceleration period
Nc = Ta/t_sample
Na = Nc
Nb = N - (Nc + Na)

#Acceleration phase
for i=1:Na;
  tT = t(1,i) - t0;
  tmp = (tT^2) / (h^2);
  q(1,i) = q0 + h*(1-sqrt(1-tmp)); 
  t(1, i+1) = t(1, i) + t_sample;
endfor
b1 = q(1,Na)

tf = ((cos(alpha)- 1)*h)/cos(alpha);
coeff(1) = a0 = tf + q0;
coeff(2) = a1 = tan(alpha);
b2 = a0 + a1*(t(1,Na)-t0)
#Here b1 != b2 are not equal which is wrong
for i=Na+1:Nb+Na;
  tT = t(1,i) - t0;
  q(1,i) = a0 + a1*tT;
  t(1, i+1) = t(1, i) + t_sample;
endfor
q(97:102)
for i=Nb+Na+1:N-1;
  tT = t1 - t(1,i);
  q(1,i) = q0 + sqrt(h*h - tT*tT);
  t(1, i+1) = t(1, i) + t_sample;
endfor

tT = t1 - t(1,N);
q(1,N) = q0 + sqrt(h*h - tT*tT);

[v, a, j] = velacceljerk(q, t_sample);
endfunction


function check_trajectory_order(order)
  if(order != 3 && order != 5 && order != 7 && order != 'c' && order != 'cy' && order != 'e'  && order != 'tr')
    error("Invalid supported order only order 3, 5 ,7, 'c', 'e', 'trpCb' ");
  endif
endfunction

function [num_points, res] = check_boundary_conditions(boundary_conditions, order)
  num_points = length(boundary_conditions);
  if(num_points < 2)
    error("Insufficitent boundary conditions");
  endif
  
  res = true;
  for i=1:num_points;
      result = false;
      [res_f] = check_boundary_conditions_field(boundary_conditions(i), order);
      res = res && res_f;
   endfor
   if(res == false)
      error("Boundary conditions field check failed ");
   endif
    [res_l] = check_boundary_conditions_limits(boundary_conditions, order);
    res = res && res_l;
endfunction

function [res] = check_boundary_conditions_field(boundary_condition, order)
    res = isfield(boundary_condition, 'time');
    res = res && isfield(boundary_condition, 'position');
    res = res && isfield(boundary_condition, 'velocity');
    if(order > 4)
       res = res && isfield(boundary_condition, 'acceleration');
    endif   
    if(order > 6)
      res = res && isfield(boundary_condition, 'jerk');
    endif
    
    if(res == 0)
      error("Boundary condition does not have complete fields");
    endif
endfunction

function [res] = check_boundary_conditions_limits(boundary_conditions, order)
  num_points = length(boundary_conditions);
  res = true;
  if(boundary_conditions(1).time >= boundary_conditions(num_points).time)
    error("End time must be greater that start time ");
    res = false;
  endif
  
  if(boundary_conditions(1).position == boundary_conditions(num_points).position)
    error("End position and Start position  must be different");
    res = false;
  endif
  
endfunction

function [v, a, j] = calculate_velacceljerk(q, tSample)
  if(tSample <= 0)
    error("Sample time must be positive");
  endif
  N = length(q);
  t_inv = 1/tSample;
  
  %Calculate velocity from position
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

function [T, h, N, q,v,a,j, t, t_sample] = initial_parameter(t0, t1, q0, q1)
  T = t1 - t0;
  h = q1 - q0;
  t_sample = 0.01;    %Sample time in seconds
  N = T/t_sample;
  q = zeros(1,N);  %Position
  
  t = zeros(1,N);
  t(1,1) = t0;
  v = zeros(1,N); %Velocity
  a = zeros(1,N);  %Acceleration
  j = zeros(1,N);  %Acceleration  
endfunction
