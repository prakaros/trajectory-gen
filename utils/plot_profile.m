function plot_profile(q,v,a,j,t, mkclr='+k')
  num_plots = 4;
  plot_single(q,t,num_plots,1,'Position', 'Time', 'Position', mkclr);
  plot_single(v,t,num_plots,2,'Velocity', 'Time', 'Velocity', mkclr);
  plot_single(a,t,num_plots,3,'Acceleration', 'Time', 'Acceleration', mkclr);
  plot_single(j,t,num_plots,4,'Jerk', 'Time', 'Jerk', mkclr);
endfunction
