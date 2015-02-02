set xrange [0:2048]
set yrange [-0.0:1.35]
set ytics nomirror
set xtics nomirror out
set key right 
set multiplot 

#set xtics (0,896,1024,1152) format ""
set xtics (0,1024) format ""
set ytics (-1,0,1)
set ytics format "% .0f"

# TOP LEFT
set size 0.55,0.56
set origin 0,0.44
plot\
     'fft_before_filtering.data' u 1:(abs($3)/2048) with lines lc rgb"white" title 'No filter',\
     'fft_before_filtering.data' u 1:(abs($3)/2048) with lines lc rgb"black" notitle 'Imag part',\
     'fft_before_filtering.data' u 1:(abs($2)/2048) with lines lc rgb"red" notitle 'Real part'

set ytics format ""

# TOP RIGHT
set size 0.53,0.56
set origin 0.47,0.44
plot\
     'filtered_fft_domain_64.data' u 1:($3/2048) with lines lc rgb"white" title '$\Delta f = \unit[64]{Hz}$',\
     'filtered_fft_domain_64.data' u 1:(abs($3)/2048) with lines lc rgb"black" notitle 'Imag part',\
     'filtered_fft_domain_64.data' u 1:(abs($2)/2048) with lines lc rgb"red" notitle 'Real part',\
               'Scale_fct_64.data' u 1:2 with lines notitle '$\Delta f = \unit[64]{Hz}$'

#set xtics (0, 896, 1024, 1152)
set xtics (0, 1024)
set xtics format "% .0f"
set ytics format "% .0f"

# BOTTOM LEFT
set size 0.55,0.54
set origin 0,0
plot\
     'filtered_fft_domain_128.data' u 1:(abs($3)/2048) with lines lc rgb"white" title '$\Delta f = \unit[128]{Hz}$',\
     'filtered_fft_domain_128.data' u 1:(abs($3)/2048) with lines lc rgb"black" notitle 'Imag part',\
     'filtered_fft_domain_128.data' u 1:(abs($2)/2048) with lines lc rgb"red" notitle 'Real part',\
               'Scale_fct_128.data' u 1:2 with lines notitle '$\Delta f = \unit[128]{Hz}$'

set ytics format ""

# BOTTOM RIGHT
set size 0.53,0.54
set origin 0.47,0
plot\
     'filtered_fft_domain_256.data' u 1:(abs($3)/2048) with lines lc rgb"white" title '$\Delta f = \unit[256]{Hz}$',\
     'filtered_fft_domain_256.data' u 1:(abs($3)/2048) with lines lc rgb"black" notitle 'Imag part',\
     'filtered_fft_domain_256.data' u 1:(abs($2)/2048) with lines lc rgb"red" notitle 'Real part',\
               'Scale_fct_256.data' u 1:2 with lines notitle '$\Delta f = \unit[256]{Hz}$'

unset multiplot
reset
