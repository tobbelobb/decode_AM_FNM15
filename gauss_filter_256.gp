set xrange [-3000:3000]
set xlabel '$f$'

plot 'Scale_fct_256.data'   u 1:2 with lines title 'filter $\Delta f = 256$',\
     'Filtered_by_256.data' u 1:($2/2) with lines title '$\Re{H(f)}$',\
     'Filtered_by_256.data' u 1:($3/2) with lines title '$\Im{H(f)}$'
