set yrange [-1:40]

plot 'warmup0_orig.data'   u 1:2 with lines title 'orig real',\
     'warmup0_orig.data'   u 1:3 with lines title 'orig imag',\
     'warmup0_transf.data' u ($1/1024):2 with lines title 'transf real',\
     'warmup0_transf.data' u ($1/1024):3 with lines title 'transf imag'

reset
