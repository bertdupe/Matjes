set xtics( \
'g'   0.00000000E+00, \
'X'   0.31415927E+01, \
'M'   0.62831853E+01, \
'g'   0.10726068E+02\
)
set grid x linetype -1
set ylabel "E-E_F"
set nokey
plot 'highs_plot.dat'
pause -1 "Hit any key to continue"
