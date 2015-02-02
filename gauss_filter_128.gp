set xrange [-3000:3000]
set xlabel '$f$'

plot 'Scale_fct_128.data'   u 1:2 with lines title 'filter $\Delta f = 128$',\
     'Filtered_by_128.data' u 1:($2/2) with lines title '$\Re{H(f)}$',\
     'Filtered_by_128.data' u 1:($3/2) with lines title '$\Im{H(f)}$'

reset
