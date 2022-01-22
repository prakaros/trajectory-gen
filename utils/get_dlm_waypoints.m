function [wayPoints] = get_dlm_waypoints(assign_intermediate_vel)
  noVal = -999999;
  printf('Reading way points from wayPoints.txt');
  [data] = dlmread('dlmWayPoints.txt', "emptyvalue", noVal);
  num_points = length(data);
  if(num_points == 0)
    error("No data in dlmWayPoints.txt");
  endif
  display(num_points)
  display(data)
  display(data(2,3))
    
  for i=1:num_points;
    wayPoints(i).tT = data(i,1);
    wayPoints(i).q  = data(i,2);
    wayPoints(i).v  = data(i,3);
    
    check_time_position(wayPoints(i), noVal);
    
    if(i > 1)
      check_time_monotonicity(wayPoints(i), wayPoints(i-1));
    endif
    
    if(i == 1 || i == num_points)
      check_initial_final_waypoint(wayPoints(i), noVal);
    endif
     
    if((i > 2) && (wayPoints(i-1).v == noVal) && (assign_intermediate_vel == true))
      printf("WayPoint %d has no velocity assigning \n", (i-1));
      k = i-1;
      dk = intermediate_velocity(wayPoints(k), wayPoints(k-1));
      dk_1 = dk = intermediate_velocity(wayPoints(k+1), wayPoints(k));
      if(sign(dk) != sign(dk_1))
        wayPoints(k).v = 0;
      else
        wayPoints(k).v = (dk + dk_1)*0.5;
      endif
    endif
    
   endfor
     
endfunction

# if sign(dk) == sign(d(k-1)) then v = 0;
# if sign(dk) != sign(d(k-1)) then v = (d(k) - d(k-1))/2;
# d(k) = q(k) - q(k-1)/t(k) - t(k-1)
function [dk] = intermediate_velocity(wayPointCurrent, wayPointPrev)
     dk = (wayPointCurrent.q - wayPointPrev.q)/(wayPointCurrent.tT - wayPointPrev.tT);
endfunction

#Time and position of waypoint cannot be empty
function check_time_position(currentWayPoint, noVal)
  if(currentWayPoint.tT == noVal || currentWayPoint.q == noVal)
      error("Time and Position cannot be empty")
   endif
endfunction

#Time should be monotonic, Always increasing
function check_time_monotonicity(currentWayPoint, previouWayPoint)
    if(currentWayPoint.tT < previouWayPoint.tT)
      error("Time must be monotonic");
    endif
endfunction

#Initial and final waypoint should have complete boundary condition(t,q,v)
function check_initial_final_waypoint(currentWayPoint, noVal)
    if(currentWayPoint.tT == noVal || currentWayPoint.q == noVal  || currentWayPoint.v == noVal )
      error("Initial or Final waypoint cannot have no value");
    endif
endfunction
