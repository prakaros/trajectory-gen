function [T, h, N, q,v,a,j, t, t_sample] = initialize_param(t0, t1, q0, q1)
  T = t1 - t0;
  h = q1 - q0;
  t_sample = 0.01;    %Sample time in seconds
  N = int64(T/t_sample)
  q = zeros(1,N);  %Position
  
  t = zeros(1,N);
  t(1,1) = t0;
  v = zeros(1,N); %Velocity
  a = zeros(1,N);  %Acceleration
  j = zeros(1,N);  %Acceleration  
endfunction