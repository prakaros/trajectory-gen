function plot_single(yAxis,xAxis,numOfSubPlots,subPlotIndex,yLabel,xLabel,ttl, mkrclr = '+k')
  subplot(numOfSubPlots,1,subPlotIndex)
  plot(xAxis, yAxis, mkrclr);
  hold on;
  %stem(f,q,'b','filled');
  ylabel(yLabel);
  xlabel(xLabel);
  title(ttl);
#  hold off;
endfunction
