
set style data lines

pl \
'maskq.7.5' u 1:3, \
'maskq.10'  u 1:3, \
'maskq.15'  u 1:3, \
'maskq.20'  u 1:3

pause -1

pl \
'maskq.7.5' u 1:4, \
'maskq.10'  u 1:4, \
'maskq.15'  u 1:4, \
'maskq.20'  u 1:4


