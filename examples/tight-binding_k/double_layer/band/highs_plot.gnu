set xtics( \
'g'   0.00000000E+00, \
'M'   0.36275987E+01, \
'K'   0.57219938E+01, \
'g'   0.99107840E+01\
)
set grid x linetype -1
set ylabel "E-E_F"
set nokey
plot 'highs_plot.dat'
pause -1 "Hit any key to continue"
