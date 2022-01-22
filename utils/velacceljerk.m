function [v, a, j] = velacceljerk(q, tSample)
  if(tSample <= 0)
    error("Sample time must be positive");
  endif
  N = length(q);
  t_inv = 1/tSample;
  
  %Calculate velocity from position
  qPrev = q(1,1);
  for i=1:N-1;
    v(1,i) = q(1,i) - qPrev;
    v(1,i) = v(1,i)*t_inv;
    qPrev = q(1,i);
  endfor
  %v(1,N-1) = v(1,N-2);  % Remove jump
  v(1,N) = v(1,N-1);  % Remove jump

  %Calculate acceleration from velocity
  vPrev = v(1,1);
  for i=1:N-1;
   a(1,i) = v(1,i) - vPrev;
    a(1,i) = a(1,i)*t_inv;
    vPrev = v(1,i);
  endfor
  a(1,N-1) = a(1,N-2);  % Remove jump
  a(1,N) = a(1,N-1);  % Remove jump

   aPrev = a(1,1);
  for i=1:N-1;
    j(1,i) = a(1,i) - aPrev;
    j(1,i) = j(1,i)*t_inv;
    aPrev = a(1,i);
  endfor

  j(1,N-2) = j(1,N-3);  % Remove jump
  j(1,N-1) = j(1,N-2);  % Remove jump
  j(1,N) = j(1,N-1);  % Remove jump
endfunction