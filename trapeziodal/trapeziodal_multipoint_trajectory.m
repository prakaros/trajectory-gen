function [q, v, a, j, t] = trapeziodal_multipoint_trajectory(boundary_condition, constraints)
  num_actuators = length(boundary_condition);
  display(num_points);
  [res, aRes, linearSeg, nVv, newTa, newT] = validate_constraints(boundary_conditions, constraints);
  
  t0 = boundary_conditions(1).time;
  q0 = boundary_conditions(1).position;

  t1 = boundary_conditions(num_points).time;
  q1 = boundary_conditions(num_points).position;
  
  #Compute total duration and acceleration duration from Vmax and Amax applicable to
  #to actuator with longest stroke
  [Ta, T, linear_segment_present] = get_durations(boundary_conditions, constraints, 1)
  
  #Now from above T and Ta compute Vi and Ai for all actuators
  [result_list] = get_actuators_vel_accel(boundary_conditions, T, Ta);
  
  [q, v, a, j, t] = compute_trajectory(result_list(1), linear_segment_present, Ta);

endfunction


function [accel_time, total_duration, linear_segment] = get_durations(actuator_list, constraints, baseline_actuator)
  num_actuators = length(actuator_list)
  h = actuator_list(baseline_actuator).endPosition - actuator_list(baseline_actuator).startPosition;
  linear_segment = true;
  accel_time = 0;
  total_duration = 0;
  %If max acceleration and max veloocity are given then
  %Ta = Vmax / Amax, Vmax = h/T-Ta; T = (h*Amax + Vmax^2)/(Amax*Vmax)
  % If ha > Vmax^2 / Amx^2, then linear segment is present
  if (constraints.maxAccel != 'Na' && constraints.maxVel != 'Na')
    display('Calculating T and Ta for given Amax and Vmax');
    newTa = constraints.maxVel  / constraints.maxAccel 
    newT = (h*constraints.maxAccel + constraints.maxVel^2)/(constraints.maxAccel*constraints.maxVel)
    if( h <(constraints.maxVel^2/constraints.maxAccel))
      display('No linear segment');
      linear_segment = false;
      newTa = sqrt(h/constraints.maxAccel );
      newT = 2*newTa;
    endif
    accel_time = newTa;
    total_duration = newT;
   endif 
endfunction

function [result_list] = get_actuators_vel_accel(actuator_list,total_duration, accel_duration)
  num_actuators = length(actuator_list)
  #Relation used Vi = hi / (T - Ta) and Ai = hi / Ta*(T-Ta);
  T = total_duration;
  Ta = accel_duration;
  velDeno = T - Ta;
  accelDeno = Ta*velDeno;
  startTime = 0;
  for i=1:num_actuators;
     hi = actuator_list(i).endPosition - actuator_list(i).startPosition;
     result_list(i).startPosition = actuator_list(i).startPosition;
     result_list(i).endPosition = actuator_list(i).endPosition;
     result_list(i).velocity = hi/velDeno;
     result_list(i).acceleration = hi/accelDeno;
     result_list(i).start_time = startTime;
     result_list(i).end_time = T - startTime;
  endfor
endfunction

function [q, v, a, j, t] = compute_trajectory(actuator, linear_segment, Ta)
  t0 = actuator.start_time;
  q0 = actuator.startPosition

  t1 = actuator.end_time;
  q1 = actuator.endPosition
  [T, h, N, q, v, a, j, t, t_sample] = initialize_param(t0, t1, q0, q1);
  
  Tb = (t1-Ta) - (t0 + Ta)  %Constant velocity
  Tc = Ta
  Na = int64(Ta/t_sample)
  Nb = 0;
  if(linear_segment == true)
    Nb = int64(Tb/t_sample)
  endif
  Nc = int64(Tc/t_sample)
  a2 = 0.5*actuator.acceleration;
  NStart = 1;
  NEnd = Na;
  for i=NStart:NEnd
    tT = t(1, i) - t0;
    q(1,i) = q0 + a2*tT*tT;
    t(1,i+1) = t(1,i) + t_sample;
  endfor
  
  NStart = Na+1;
  NEnd = Na + Nb;
  
  Vv = actuator.acceleration*Ta;
  for i = NStart: NEnd
    tT = t(1,i) - t0 - (Ta/2);
    q(1,i) = q0 + Vv*tT;
    t(1,i+1) = t(1,i) + t_sample;
  endfor
  
  c2 = 0.5*actuator.acceleration;
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
