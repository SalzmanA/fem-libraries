set terminal tikz scale 0.65,0.54 fontscale 0.5
set output "curves_nl_time_all.tex"
#set title "Full simulation"
set xlabel "number of cores"
set logscale x
set logscale y
set grid
#set ylabel "Elapse time in s"
set key right top
set xrange [1.:128.1]
set yrange [10:2100.]
plot "curve_nl_time.txt" using 1:2:xtic(1) with lp pt 4 lc "green" title "MFEM full nl" ,"curve_nl_time.txt" using 1:18:xtic(1) with lp pt 5 lc "green" title "MFEM full nl (AD)" ,"curve_nl_time.txt" using 1:52 with lp pt 4 lc "forest-green" title "MFEM full nl fcst" ,"curve_nl_time.txt" using 1:68 with lp pt 5 lc "dark-green" title "MFEM full nl fcst (AD)" ,[1:128] 2000/x lc "black" notitle,"curve_nl_time.txt" using 1:34 with lp pt 1 lc "orange" title "FEniCSx full nl","curve_nl_time.txt" using 1:84 with lp pt 6 lc "orange" title "FEniCSx full nl ctrl apply"
set output "curves_nl_time_nloop.tex"
#set title "Non-linear loop"
#set ylabel "Elapse time in s"
set key right top vertical inside
set xrange [1.:128.1]
set yrange [1.:430]
plot "curve_nl_time.txt" using 1:($15/15):xtic(1) with lp pt 4 lc "green" title "MFEM full nl" ,"curve_nl_time.txt" using 1:($31/15):xtic(1) with lp pt 5 lc "green" title "MFEM full nl (AD)","curve_nl_time.txt" using 1:($65/15) with lp pt 4 lc "forest-green" title "MFEM full nl fcst" ,"curve_nl_time.txt" using 1:($81/15) with lp pt 5 lc "dark-green" title "MFEM full nl fcst (AD)" ,[1:128] 200/x lc "black" notitle,"curve_nl_time.txt" using 1:($49/11) with lp pt 1 lc "orange" title "FEniCSx full nl","curve_nl_time.txt" using 1:($99/11) with lp pt 6 lc "orange" title "FEniCSx full nl ctrl apply"
set output "curves_nl_time_elem_vect.tex"
#set title "Elementary creation vect"
#set ylabel "Elapse time in s"
set key right top inside vertical
set xrange [1.:128.1]
set yrange [0.1:100]
plot "curve_nl_time.txt" using 1:11 with lp pt 4 lc "green" title "MFEM full nl","curve_nl_time.txt" using 1:27:xtic(1) with lp pt 5 lc "green" title "MFEM full nl (AD)","curve_nl_time.txt" using 1:61 with lp pt 4 lc "dark-green" title "MFEM full nl fcst","curve_nl_time.txt" using 1:77 with lp pt 5 lc "dark-green" title "MFEM full nl fcst (AD)", [1:128] 30/x lc "black" notitle,"curve_nl_time.txt" using 1:45 with lp pt 1 lc "orange" title "FEniCSx full nl","curve_nl_time.txt" using 1:95 with lp pt 6 lc "orange" title "FEniCSx full nl ctrl apply"
set output "curves_nl_time_elem_mat.tex"
#set title "Elementary creation matrix"
#unset ylabel
set xrange [1.:128.1]
set yrange [0.35:150]
plot "curve_nl_time.txt" using 1:12 with lp pt 4 lc "green" notitle,"curve_nl_time.txt" using 1:28:xtic(1) with lp pt 5 lc "green" notitle, "curve_nl_time.txt" using 1:62 with lp pt 4 lc "dark-green" notitle,"curve_nl_time.txt" using 1:78 with lp pt 5 lc "dark-green" notitle, [1:128] 100/x lc "black" notitle,"curve_nl_time.txt" using 1:46 with lp pt 1 lc "orange" notitle,"curve_nl_time.txt" using 1:96 with lp pt 6 lc "orange" notitle
set output "curves_nl_time_elem_vect_r.tex"
#set title "Elementary creation vect"
#unset ylabel
#set ylabel "Elapse time in s"
#set key left bottom inside vertical
set xrange [1.:128.1]
set yrange [0.15:1.5]
set ytics ("%.2f (4)" 0.25,"%.3f (1.17)" 0.8547,1.,1.3,1.4)
plot "curve_nl_time.txt" using 1:($27/$11):xtic(1) with lp pt 5 lc "green" notitle, "curve_nl_time.txt" using 1:($61/$11) with lp pt 4 lc "dark-green" notitle,"curve_nl_time.txt" using 1:($77/$11) with lp pt 5 lc "dark-green" notitle, "curve_nl_time.txt" using 1:($45/$11) with lp pt 1 lc "orange" notitle, "curve_nl_time.txt" using 1:($95/$11) with lp pt 6 lc "orange" notitle
set output "curves_nl_time_elem_mat_r.tex"
#set title "Elementary creation matrix"
#unset ylabel
set xrange [1.:128.1]
set yrange [0.68:2.535]
set ytics ("%.3f (1.17)" 0.8547,1.,2.4)
#unset logscale y
plot "curve_nl_time.txt" using 1:($28/$12):xtic(1) with lp pt 5 lc "green" notitle, "curve_nl_time.txt" using 1:($62/$12) with lp pt 4 lc "dark-green" notitle,"curve_nl_time.txt" using 1:($78/$12) with lp pt 5 lc "dark-green" notitle, "curve_nl_time.txt" using 1:($46/$11) with lp pt 1 lc "orange" notitle, "curve_nl_time.txt" using 1:($96/$11) with lp pt 6 lc "orange" notitle
