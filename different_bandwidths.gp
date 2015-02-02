set xrange [0:1]
set multiplot layout 2,2
set title '$\Delta f = \unit[64]{Hz}$'
plot 'filtered_64.data'  u 1:2 with lines notitle
set title '$\Delta f = \unit[128]{Hz}$'
plot 'filtered_128.data' u 1:2 with lines notitle 
set title '$\Delta f = \unit[256]{Hz}$'
plot 'filtered_256.data' u 1:2 with lines notitle 
set title '$\Delta f = \unit[512]{Hz}$'
plot 'filtered_512.data' u 1:2 with lines notitle

unset multiplot
reset
