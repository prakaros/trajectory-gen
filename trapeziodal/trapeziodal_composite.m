function [q, v, a, j, t, coeff] = trapeziodal_composite(boundary_conditions, constraints, order)
    [q, v, a, j, t, coeff] = trapeziod_linear_parabolic(boundary_conditions, constraints);
endfunction

function [q, v, a, j, t, coeff] = trapeziod_linear_parabolic(boundary_conditions, constraints)
  num_points = length(boundary_conditions);
  display(num_points);
  [res, aRes, linearSeg, nVv, newTa, newT] = validate_constraints(boundary_conditions, constraints);
  
  t0 = boundary_conditions(1).time;
  q0 = boundary_conditions(1).position;

  t1 = boundary_conditions(num_points).time;
  q1 = boundary_conditions(num_points).position;
   
   Ta = newTa;
   Vv = constraints.constantVel;
   t1 = newT - t0;
   if(aRes == false)
    Vv = nVv
   endif
   
   useMaxAccel = false;
   if(constraints.maxAccel != 'Na')
      useMaxAccel = true;
   endif
   coeff = 0;
  [T, h, N, q, v, a, j, t, t_sample] = initialize_param(t0, t1, q0, q1);
  
%Linear trapeziodal trajecttory with parabolic bends
% Three part of trajectory, 
% 1. Acceleration phase 2. Constant velocity phase 3. Deceleration phase
% a) t = (t0, t0+Ta) where t0 is start time, Ta is acceleration phase time
% Initial conditions qa(t0) = q0, va(t0) = v0, va(t0+Ta) = Vv
% qa(t) = a0 + a1(t-t0) + a2*(t-t0)^2
% c) t = (t1-Ta, t1) where t1 is end time, Ta is deceleration phase time
% Final condition qc(t1) = q1, vc(t1) = v1, vc(t1-Ta) = Vv
% qc(t) = c0 + c1*(t1-t) + c2*(t1-t)^2
% b) t = (t0+Ta, t1-Ta) total time t1-Ta-t0-Ta = t1-t0-2Ta = T-2Ta
% Constant velocity during this phase Vv
% qb(t) = b0 + b1(t-t0-Ta)
% For simplicity t0 = 0, v0 = 0
% qa(t) = a0 + a1*t + a1*t^2, qb(t) = b0 + b1*t, qc(t) = c0 + c1*t + c2*t^2
% Solving with boundary conditions
% qa(t) = q0 + (Vv/2*Ta)*(t-t0)^2
% qb(t) = q0 + Vv(t-t0-Ta/2)
% qc(t) = q0 - (Vv/2*Ta)*(t1-t)^2

Tb = (t1-Ta) - (t0 + Ta)  %Constant velocity
Tc = Ta
Na = int64(Ta/t_sample)
Nb = 0;
if(linearSeg == true)
  Nb = int64(Tb/t_sample)
endif
Nc = int64(Tc/t_sample)

if(useMaxAccel == true)
  a2 = 0.5*constraints.maxAccel;
else
  a2 = Vv/(2*Ta);
endif
NStart = 1;
NEnd = Na;
for i=NStart:NEnd
  tT = t(1, i) - t0;
  q(1,i) = q0 + a2*tT*tT;
  t(1,i+1) = t(1,i) + t_sample;
endfor


NStart = Na+1;
NEnd = Na + Nb;

if(useMaxAccel == true)
  Vv = constraints.maxAccel*Ta;
endif

for i = NStart: NEnd
  tT = t(1,i) - t0 - (Ta/2);
  q(1,i) = q0 + Vv*tT;
  t(1,i+1) = t(1,i) + t_sample;
endfor
q(1,NEnd)

if(useMaxAccel == true)
  c2 = 0.5*constraints.maxAccel;
else
  c2 = Vv/(2*Ta);
endif

NStart = Na+Nb
NEnd = Na+ Nb + Nc -1;
for i = NStart : NEnd
  tT = t1 - t(1,i);
  q(1,i) = q1 - c2*(tT*tT);
  t(1,i+1) = t(1,i) + t_sample;
endfor
q(1,NStart)
tT = t1 - t(1,N);
q(1,N) = q1 - c2*(tT*tT);

[v, a, j] = velacceljerk(q, t_sample);
endfunction

function [res, accelRes, linearSeg, newVv, newTa, newT] = validate_constraints(boundary_conditions, constraints)
  num_points = length(boundary_conditions);
  res = true;
  accelRes = true;
  linearSeg = true;
  if(boundary_conditions(1).time >= boundary_conditions(num_points).time)
    error("End time must be greater that start time ");
    res = false;
  endif
  
  if(boundary_conditions(1).position == boundary_conditions(num_points).position)
    error("End position and Start position  must be different");
    res = false;
  endif
  
  %1. Constant velocity time interval = (t1-Ta)-(t0+Ta) = t1-t0-2Ta = T-2Ta >= 0
  %Ta <= T/2
  newT = T =  boundary_conditions(num_points).time - boundary_conditions(1).time;
  newTa = Ta = constraints.accelTime;
  newVv = Vv  = constraints.constantVel;
  if(constraints.accelTime > T/2)
    error("AccelTime should be <= TotalTrajectory Time / 2");
    res = false;
  endif
  
  %2. Aa*Ta = (Qm - Qa) /(Tm - Ta) Where Aa = Constant acceleration  
  %Tm = (t1-t0)/2 = T/2, Qm = q1+q0/2 = q0 + h/2, Qa = q(r0 + Ta)
  % Qa = Q0 + 0.5*Aa*Ta^2 , we get Aa*Ta*Ta - Aa*(t1-t0)*Ta + q1-q0 = 0
  % 3 Vv = (q1-q0)/(t1-t0-Ta) = h/(T-Ta);
  % Aa = Vv/Ta;
  Aa = Vv/Ta
  t1 = boundary_conditions(num_points).time;
  t0 = boundary_conditions(1).time;
  q1 = boundary_conditions(num_points).position;
  q0 =  boundary_conditions(1).position;
  h = q1-q0;
  m1 = Aa*Ta*Ta - Aa*T*Ta + h
  if(m1 != 0)
    display("Constant accel constaint not satisfied");
    display("Recomputing new Constant Accel and constant Vel");
    Aanew = h/(T*Ta - Ta*Ta);
    newVv = Aanew*Ta;
    m2 = Aanew*Ta*Ta - Aanew*T*Ta + h
    res = false;
    accelRes = false;
  endif
  
  %If maxAcceleration and max Veloocity are given then
  % Ta = Vmax / Amax, Vmax = h/T-Ta; T = (h*Amax + Vmax^2)/(Amax*Vmax)
  if (constraints.maxAccel != 'Na' && constraints.maxVel != 'Na')
    display('Calculating T and T for given Amac and Vmax');
    newTa = constraints.maxVel  / constraints.maxAccel 
    newT = (h*constraints.maxAccel + constraints.maxVel^2)/(constraints.maxAccel*constraints.maxVel)
    if( h <(constraints.maxVel^2/constraints.maxAccel))
      display('No linear segment');
      linearSeg = false;
      newTa = sqrt(h/constraints.maxAccel );
      newT = 2*newTa;
    endif
   endif
endfunction
