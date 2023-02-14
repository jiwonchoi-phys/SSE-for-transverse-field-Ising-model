set term pdf
set output "L4pure.pdf"
set log x
set xlabel "b"
set ylabel "U"

plot 'L4avg' u 2:(1-$6/$5**2/3) t "SSE",\
     'l4u' u 1:2 t "World-line"
