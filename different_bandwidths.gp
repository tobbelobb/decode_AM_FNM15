set xrange [0:0.065]
set multiplot 

set key right 


set xtics (0,\
           0.0234375,\
           0.0468750,\
           0.0703125,\
           0.093750) format ""
set ytics (-1,1,1) format "% .0f"
set y2tics (-1,1,1) format ""

# TOP LEFT
set size 0.55,0.56
set origin 0,0.44
plot 'filtered_64.data' u 1:2 with lines lc rgb"white" title '$\Delta f = \unit[64]{Hz}$',\
     'filtered_64.data' u 1:2 with lines lc rgb"black" notitle

set ytics (-1,1,1) format ""
set y2tics (-1,1,1) format "% .0f"

# TOP RIGHT
set size 0.53,0.56
set origin 0.47,0.44
plot 'filtered_128.data' u 1:2 axes x1y2 with lines lc rgb"white" title '$\Delta f = \unit[128]{Hz}$',\
     'filtered_128.data' u 1:2 axes x1y2 with lines lc rgb"black" notitle

set xtics ('$0$' 0,\
           '$24$' 0.0234375,\
           '$48$' 0.0468750,\
           '$72$' 0.0703125,\
           '$96$' 0.093750)
set xlabel '$1024\cdot t/\unit{s}$'

set ytics (-1,1,1) format "% .0f"
set y2tics (-1,1,1) format ""

# BOTTOM LEFT
set size 0.55,0.54
set origin 0,0
plot 'filtered_256.data' u 1:2 with lines lc rgb"white" title '$\Delta f = \unit[256]{Hz}$',\
     'filtered_256.data' u 1:2 with lines lc rgb"black" notitle

set ytics (-1,1,1) format ""
set y2tics (-1,1,1) format "% .0f"

# BOTTOM RIGHT
set size 0.53,0.54
set origin 0.47,0
plot 'filtered_512.data' u 1:2 with lines lc rgb"white" title '$\Delta f = \unit[512]{Hz}$',\
     'filtered_512.data' u 1:2 with lines lc rgb"black" notitle

unset multiplot
reset
