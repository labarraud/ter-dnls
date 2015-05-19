set terminal postscript eps enhanced color 
set output "graph/plot.eps"
set xlabel "x"
set ylabel "Phi"
set key on inside right top
plot 'result/Euler_e.dat' pointtype 1, 'result/Euler_i.dat' pointtype 8, 'result/RK4.dat' pointtype 2, 'result/Adam4.dat' pointtype 19, 'result/Adam4cor.dat' pointtype 4,  'result/exact.dat' with line linecolor -1 lw 2,  

set terminal postscript eps enhanced color 
set output "graph/int.eps"
set xlabel "x"
set ylabel "y"
set key on inside right top
plot 'result/intphi2.dat' with line linecolor 1 lw 2,  'result/exactint.dat' with line linecolor -1 lw 2

set terminal postscript eps enhanced color 
set output "graph/u0.eps"
set zlabel "x"
set xlabel "real(u)"
set ylabel "imaginaire(u)"			
set key on inside right top
splot 'result/u0.dat' with line linecolor 1 lw 2, 'result/u1.dat' with line linecolor 2 lw 2

set terminal postscript eps enhanced color 
set output "graph/modu0.eps"
set xlabel "x"
set ylabel "modu0"
set key on inside right top
plot 'result/modu0.dat' with line linecolor 1 lw 2, 'result/modu1.dat' with line linecolor 2 lw 2

set terminal postscript eps enhanced color 
set output "graph/error_n2.eps"
set logscale
set xlabel "step h"
set ylabel "error with euclidean norm"
set key on inside right bottom

Eeuler_n2(x) = ae2*x + be2
fit Eeuler_n2(x) 'result/error_euler_norme2.dat' u (log($1)):(log($2)) via ae2, be2
title_Eeuler_n2(ae2, be2) = sprintf('Eeulern2(x) = %.2fx + %.2f', ae2, be2)

Eeuleri_n2(x) = aei2*x + bei2
fit Eeuleri_n2(x) 'result/error_euleri_norme2.dat' u (log($1)):(log($2)) via aei2, bei2
title_Eeuleri_n2(aei2, bei2) = sprintf('Eeulerin2(x) = %.2fx + %.2f', aei2, bei2)

ERK4_n2(x) = aRK42*x + bRK42
fit ERK4_n2(x) 'result/error_RK4_norme2.dat' u (log($1)):(log($2)) via aRK42, bRK42
title_ERK4_n2(aRK42, bRK42) = sprintf('ERK4n2(x) = %.2fx + %.2f', aRK42, bRK42)

EAdam4_n2(x) = aAdam42*x + bAdam42
fit EAdam4_n2(x) 'result/error_Adam4_norme2.dat' u (log($1)):(log($2)) via aAdam42, bAdam42
title_EAdam4_n2(aAdam42, bAdam42) = sprintf('EAdam4n2(x) = %.2fx + %.2f', aAdam42, bAdam42)

EAdam4cor_n2(x) = aAdam4cor2*x + bAdam4cor2
fit EAdam4cor_n2(x) 'result/error_Adam4cor_norme2.dat' u (log($1)):(log($2)) via aAdam4cor2, bAdam4cor2
title_EAdam4cor_n2(aAdam4cor2, bAdam4cor2) = sprintf('EAdam4corn2(x) = %.2fx + %.2f', aAdam4cor2, bAdam4cor2)

plot 'result/error_euler_norme2.dat', exp(be2)*x**ae2 t title_Eeuler_n2(ae2, be2),'result/error_euleri_norme2.dat', exp(bei2)*x**aei2 t title_Eeuleri_n2(aei2, bei2), 'result/error_RK4_norme2.dat', exp(bRK42)*x**aRK42 t title_ERK4_n2(aRK42, bRK42),  'result/error_Adam4_norme2.dat', exp(bAdam42)*x**aAdam42 t title_EAdam4_n2(aAdam42, bAdam42),'result/error_Adam4cor_norme2.dat', exp(bAdam4cor2)*x**aAdam4cor2 t title_EAdam4cor_n2(aAdam4cor2, bAdam4cor2)


set terminal postscript eps enhanced color 
set output "graph/error_nf.eps"
set logscale
set xlabel "step h"
set ylabel "error with infinity norm"
set key on inside right bottom

Eeuler_nf(x) = aef*x + bef
fit Eeuler_nf(x) 'result/error_euler_nf.dat' u (log($1)):(log($2)) via aef, bef
title_Eeuler_nf(aef, bef) = sprintf('Eeulernf(x) = %.2fx + %.2f', aef, bef)

Eeuleri_nf(x) = aeif*x + beif
fit Eeuleri_nf(x) 'result/error_euleri_nf.dat' u (log($1)):(log($2)) via aeif, beif
title_Eeuleri_nf(aeif, beif) = sprintf('Eeulerinf(x) = %.2fx + %.2f', aeif, beif)

ERK4_nf(x) = aRK4f*x + bRK4f
fit ERK4_nf(x) 'result/error_RK4_nf.dat' u (log($1)):(log($2)) via aRK4f, bRK4f
title_ERK4_nf(aRK4f, bRK4f) = sprintf('ERK4nf(x) = %.2fx + %.2f', aRK4f, bRK4f)

EAdam4_nf(x) = aAdam4f*x + bAdam4f
fit EAdam4_nf(x) 'result/error_Adam4_nf.dat' u (log($1)):(log($2)) via aAdam4f, bAdam4f
title_EAdam4_nf(aAdam4f, bAdam4f) = sprintf('EAdam4nf(x) = %.2fx + %.2f', aAdam4f, bAdam4f)

EAdam4cor_nf(x) = aAdam4corf*x + bAdam4corf
fit EAdam4cor_nf(x) 'result/error_Adam4cor_nf.dat' u (log($1)):(log($2)) via aAdam4corf, bAdam4corf
title_EAdam4cor_nf(aAdam4corf, bAdam4corf) = sprintf('EAdam4cornf(x) = %.2fx + %.2f', aAdam4corf, bAdam4corf)

plot 'result/error_euler_nf.dat', exp(bef)*x**aef t title_Eeuler_nf(aef, bef), 'result/error_euleri_nf.dat', exp(beif)*x**aeif t title_Eeuleri_nf(aeif, beif),'result/error_RK4_nf.dat', exp(bRK4f)*x**aRK4f t title_ERK4_nf(aRK4f, bRK4f), 'result/error_Adam4_nf.dat', exp(bAdam4f)*x**aAdam4f t title_EAdam4_nf(aAdam4f, bAdam4f),'result/error_Adam4cor_nf.dat', exp(bAdam4corf)*x**aAdam4corf t title_EAdam4cor_nf(aAdam4corf, bAdam4corf)

set terminal postscript eps enhanced color 
set output "graph/error.eps"
set logscale
set xlabel "Pas h"
set ylabel "Error"
set key on inside right bottom

plot 'result/error_euler_norme2.dat', exp(be2)*x**ae2 t title_Eeuler_n2(ae2, be2),'result/error_euleri_norme2.dat', exp(bei2)*x**aei2 t title_Eeuleri_n2(aei2, bei2),  'result/error_RK4_norme2.dat', exp(bRK42)*x**aRK42 t title_ERK4_n2(aRK42, bRK42),  'result/error_Adam4_norme2.dat', exp(bAdam42)*x**aAdam42 t title_EAdam4_n2(aAdam42, bAdam42), 'result/error_euler_nf.dat', exp(bef)*x**aef t title_Eeuler_nf(aef, bef), 'result/error_euleri_nf.dat', exp(beif)*x**aeif t title_Eeuleri_nf(aeif, beif),'result/error_RK4_nf.dat', exp(bRK4f)*x**aRK4f t title_ERK4_nf(aRK4f, bRK4f), 'result/error_Adam4_nf.dat', exp(bAdam4f)*x**aAdam4f t title_EAdam4_nf(aAdam4f, bAdam4f),'result/error_Adam4cor_norme2.dat', exp(bAdam4cor2)*x**aAdam4cor2 t title_EAdam4cor_n2(aAdam4cor2, bAdam4cor2),'result/error_Adam4cor_nf.dat', exp(bAdam4corf)*x**aAdam4corf t title_EAdam4cor_nf(aAdam4corf, bAdam4corf)

set terminal postscript eps enhanced color 
set output "graph/test.eps"
test


