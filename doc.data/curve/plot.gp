set terminal unknown
plot "curve_time.txt" using 1:2 every ::0::0
FENICSX1p=GPVAL_DATA_Y_MIN
plot "curve_time.txt" using 1:5 every ::0::0
REF1p=GPVAL_DATA_Y_MIN
set terminal tikz scale 0.65,0.54 fontscale 0.5
set output "curves_time_all.tex"
#set title "Full simulation"
set xlabel "number of cores"
set logscale x
set logscale y
set grid
#set ylabel "Elapse time in s"
set key right top
set xrange [1.:128.1]
set yrange [7.5:1800.]
plot "curve_time.txt" using 1:2 with lp pt 1 lc "orange" title "FEniCSx (handmade)" ,"curve_time.txt" using 1:20 with lp pt 2 lc "orange" title "FEniCSx (sympy)" ,"curve_time.txt" using 1:38 with lp pt 3 lc "orange" title "FEniCSx (sympy sym)" ,"curve_time.txt" using 1:56 with lp pt 0 lc "orange" title "FEniCSx (ufl)" ,"curve_time_py.txt" using 1:2 with lp pt 3 lc "blue" title "FEniCSx-py (sympy sym)" ,"curve_time.txt" using 1:74:xtic(1) with lp pt 4 lc "green" title "MFEM" ,"curve_time.txt" using 1:90:xtic(1) with lp pt 5 lc "green" title "MFEM (AD)" , "curve_time.txt" using 1:202 with lp pt 4 lc "forest-green" title "MFEM Fcst" ,"curve_time.txt" using 1:218 with lp pt 5 lc "forest-green" title "MFEM Fcst (AD)" ,[1:128] FENICSX1p/x lc "black" notitle
set output "curves_time_misc.tex"
#set title "Small elapse times"
#unset ylabel
set key center top vertical
set xrange [1.:128.1]
set yrange [0.00001:10.1]
plot "curve_time.txt" using 1:3 with lp pt 1 lc "orange" title "Initialize" ,"curve_time_py.txt" using 1:3 with lp pt 3 lc "blue" title " " ,"curve_time.txt" using 1:75:xtic(1) with lp pt 4 lc "green" title " " , "curve_time.txt" using 1:8 with lp pt 1 lc "dark-orange" title "Material constant" , "curve_time_py.txt" using 1:8 with lp pt 1 lc "dark-blue" title " " ,"curve_time.txt" using 1:80:xtic(1) with lp pt 4 lc "dark-green" title " " 
set output "curves_time_obj.tex"
#set title "space,form solver elapse times"
#unset ylabel
set key right bottom vertical outside
set xrange [1.:128.1]
set yrange [0.00001:3]
plot "curve_time.txt" using 1:6:xtic(1) with lp pt 1 lc "orange" title "space" ,"curve_time_py.txt" using 1:6 with lp pt 3 lc "blue" title " " ,"curve_time.txt" using 1:78 with lp pt 4 lc "green" title " ", "curve_time.txt" using 1:15 with lp pt 1 lc "dark-orange" title "form" ,"curve_time_py.txt" using 1:11 with lp pt 3 lc "dark-blue" title " " ,"curve_time.txt" using 1:85 with lp pt 4 lc "dark-green" title " ", "curve_time.txt" using 1:16 with lp pt 1 lc "red" title "solver" ,"curve_time_py.txt" using 1:26 with lp pt 3 lc "cyan" title " " ,"curve_time.txt" using 1:86 with lp pt 4 lc "dark-red" title " " 
set output "curves_time_BC.tex"
#set title "BC elapse times"
#unset ylabel
set key right top vertical inside
set xrange [1.:128.1]
set yrange [0.0005:3]
plot "curve_time.txt" using 1:9:xtic(1) with lp pt 1 lc "orange" title "Dirichlet" ,"curve_time_py.txt" using 1:9 with lp pt 3 lc "blue" title " " ,"curve_time.txt" using 1:81 with lp pt 4 lc "green" title " " , "curve_time.txt" using 1:10 with lp pt 1 lc "dark-orange" title "Neuman" , "curve_time_py.txt" using 1:10 with lp pt 3 lc "dark-blue" title " " ,"curve_time.txt" using 1:82 with lp pt 4 lc "dark-green" title " " 
set output "curves_time_mesh.tex"
#set title "mesh treatment"
#set ylabel "Elapse time in s"
set key center top right inside
set xrange [1.:128.1]
set yrange [0.015:25]
plot "curve_time.txt" using 1:4 with lp pt 1 lc "orange" title "Read" ,"curve_time_py.txt" using 1:4 with lp pt 3 lc "blue" title " " ,"curve_time.txt" using 1:76:xtic(1) with lp pt 4 lc "green" title " " , "curve_time.txt" using 1:5 with lp pt 1 lc "dark-orange" title "Refine" ,"curve_time_py.txt" using 1:5 with lp pt 3 lc "dark-blue" title " " ,"curve_time.txt" using 1:77 with lp pt 4 lc "dark-green" title " " ,[1:128] REF1p/x lc "black" notitle
set output "curves_time_dam.tex"
#unset ylabel
set key right top vertical inside
set xrange [1.:128.1]
set yrange [0.04:250]
plot "curve_time.txt" using 1:7 with lp pt 1 lc "orange" title "FEniCSx" ,"curve_time_py.txt" using 1:7 with lp pt 3 lc "blue" title "FEniCSx-py", "curve_time.txt" using 1:79:xtic(1) with lp pt 4 lc "green" title "MFEM" ,[1:128] REF1p*1.0/x lc "black" notitle
set output "curves_time_nl.tex"
#set title "Non-linear loop"
#set ylabel "Elapse time in s"
set key right top vertical inside
set xrange [1.:128.1]
set yrange [1.:430]
plot "curve_time.txt" using 1:($17/7) with lp pt 1 lc "orange" title "FEniCSx (handmade)" ,"curve_time.txt" using 1:($35/7) with lp pt 2 lc "orange" title "FEniCSx (sympy)" ,"curve_time.txt" using 1:($53/7) with lp pt 3 lc "orange" title "FEniCSx (sympy sym)" ,"curve_time.txt" using 1:($71/7) with lp pt 0 lc "orange" title "FEniCSx (ufl)" ,"curve_time_py.txt" using 1:($13/7) with lp pt 3 lc "blue" title "FEniCSx-py (sympy sym)" ,"curve_time.txt" using 1:($87/5):xtic(1) with lp pt 4 lc "green" title "MFEM" ,"curve_time.txt" using 1:($103/5):xtic(1) with lp pt 5 lc "green" title "MFEM (AD)" ,"curve_time.txt" using 1:($183/5) with lp pt 4 lc "dark-green" title "MFEM h1" ,"curve_time.txt" using 1:($199/5):xtic(1) with lp pt 5 lc "dark-green" title "MFEM h1 (AD)" ,"curve_time.txt" using 1:($215/5) with lp pt 4 lc "forest-green" title "MFEM Fcst" ,"curve_time.txt" using 1:($231/5):xtic(1) with lp pt 5 lc "forest-green" title "MFEM Fcst (AD)" ,[1:128] FENICSX1p/x/7 lc "black" notitle
set output "curves_time_ass.tex"
#set title "Elementary creation and assemby"
#unset ylabel
set key left bottom vertical inside
set xrange [1.:128.1]
set yrange [0.15:50]
set logscale y2
set ytics nomirror
set y2tics nomirror (22,30,37)
set grid y2
set ylabel "Elapse time in s"
set y2label "ratio in \\%"
plot "curve_time.txt" using 1:11:xtic(1) with lp pt 1 lc "orange" title "Vector crea+ass" ,"curve_time.txt" using 1:115 with lp pt 1 lc "red" title "Vector crea no prof+ass " , "curve_time.txt" using 1:12 with lp pt 1 lc "blue" title "Matrix crea+ass", "curve_time.txt" using 1:116 with lp pt 1 lc "cyan" title "Matrix crea no prof+ass",[1:128] REF1p/x lc "black" notitle,"curve_time.txt" using 1:(100.*$13/$11) with lp pt 1 lc "dark-orange" title "Vector $100\\times\\frac{crea}{crea+ass}$"  axes x1y2, "curve_time.txt" using 1:(100.*$14/$12) with lp pt 1 lc "skyblue" title "Matrix $100\\times\\frac{crea}{crea+ass}$" axes x1y2
unset ylabel
unset y2label
unset grid
set grid
set ytics mirror
set output "curves_time_elem_vect.tex"
#set title "Elementary creation vect"
#set ylabel "Elapse time in s"
set key right top inside vertical
set xrange [1.:128.1]
set yrange [0.05:45]
plot "curve_time.txt" using 1:13 with lp pt 1 lc "orange" title "FEniCSx (handmade)" ,"curve_time.txt" using 1:31 with lp pt 2 lc "orange" title "FEniCSx (sympy)" ,"curve_time.txt" using 1:49 with lp pt 3 lc "orange" title "FEniCSx (sympy sym)" ,"curve_time.txt" using 1:67 with lp pt 0 lc "orange" title "FEniCSx (ufl)" ,"curve_time.txt" using 1:195:xtic(1) with lp pt 5 lc "dark-green" title "MFEM h1 (AD)", "curve_time.txt" using 1:179 with lp pt 4 lc "dark-green" title "MFEM h1", "curve_time.txt" using 1:227 with lp pt 5 lc "forest-green" title "MFEM Fcst (AD)", "curve_time.txt" using 1:211 with lp pt 4 lc "forest-green" title "MFEM Fcst","curve_time.txt" using 1:99:xtic(1) with lp pt 5 lc "green" title "MFEM (AD)","curve_time.txt" using 1:83 with lp pt 4 lc "green" title "MFEM", [1:128] REF1p/x lc "black" notitle
set output "curves_time_elem_mat.tex"
#set title "Elementary creation matrix"
#unset ylabel
set xrange [1.:128.1]
set yrange [0.065:25]
plot "curve_time.txt" using 1:14 with lp pt 1 lc "orange" notitle ,"curve_time.txt" using 1:32 with lp pt 2 lc "orange" notitle ,"curve_time.txt" using 1:50 with lp pt 3 lc "orange" notitle ,"curve_time.txt" using 1:68 with lp pt 0 lc "orange" notitle ,"curve_time.txt" using 1:84 with lp pt 4 lc "green" notitle,"curve_time.txt" using 1:100:xtic(1) with lp pt 5 lc "green" notitle,"curve_time.txt" using 1:180 with lp pt 4 lc "dark-green" notitle,"curve_time.txt" using 1:196:xtic(1) with lp pt 5 lc "dark-green" notitle, [1:128] REF1p/x lc "black" notitle,"curve_time.txt" using 1:212 with lp pt 4 lc "forest-green" notitle,"curve_time.txt" using 1:228 with lp pt 5 lc "forest-green" notitle
set output "curves_time_output.tex"
#set title "Outputs"
set key right top horizontal inside
set xrange [1.:128.1]
set yrange [0.2:32]
plot "curve_time.txt" using 1:18 with lp pt 1 lc "orange" notitle ,"curve_time_py.txt" using 1:14 with lp pt 3 lc "blue" notitle ,"curve_time.txt" using 1:88:xtic(1) with lp pt 4 lc "green" notitle ,[1:128] 20.82/x lc "black" notitle
set output "curves_time_stress_strain.tex"
#set title "Outputs"
set key right top horizontal inside
set xrange [1.:128.1]
set yrange [0.015:10]
plot "curve_time.txt" using 1:19 with lp pt 1 lc "orange" notitle ,"curve_time.txt" using 1:37 with lp pt 2 lc "orange" notitle ,"curve_time.txt" using 1:55 with lp pt 3 lc "orange" notitle ,"curve_time.txt" using 1:73 with lp pt 0 lc "orange" notitle ,"curve_time_py.txt" using 1:15 with lp pt 3 lc "blue" notitle ,"curve_time.txt" using 1:89 with lp pt 4 lc "green" notitle,"curve_time.txt" using 1:105:xtic(1) with lp pt 5 lc "green" notitle, [1:128] 5/x lc "black" notitle
set output "curves_time_all_r.tex"
#set title "Full simulation"
set xlabel "number of cores"
#set ylabel "ratio vs MFEM"
set key center top vertical
set xrange [1.:128.1]
set yrange [0.9:2.3]
set ytics (0.9,1.,1.2,1.4,1.6,2.,2.25)
unset logscale y
plot "curve_time.txt" using 1:($2/$74) with lp pt 1 lc "orange" notitle,"curve_time.txt" using 1:($20/$74) with lp pt 2 lc "orange" notitle ,"curve_time.txt" using 1:($38/$74) with lp pt 3 lc "orange" notitle,"curve_time.txt" using 1:($56/$74) with lp pt 0 lc "orange" notitle ,"< paste curve_time.txt curve_time_py.txt" using 1:($235/$74) with lp pt 3 lc "blue" notitle,"curve_time.txt" using 1:($90/$74):xtic(1) with lp pt 5 lc "green" notitle , "curve_time.txt" using 1:($202/$74) with lp pt 4 lc "forest-green" notitle, "curve_time.txt" using 1:($218/$74) with lp pt 5 lc "forest-green" notitle
set output "curves_time_mesh_r.tex"
#set title "mesh treatment"
#set ylabel "ratio vs MFEM"
set key center top right inside
set xrange [1.:128.1]
set yrange [0.15:15]
set logscale y
set ytics ("%.1f (5)" 0.2,"%.1f (2)" 0.5,1.,1.5,8.,12.,14.)
plot "curve_time.txt" using 1:($4/$76):xtic(1) with lp pt 1 lc "orange" title "Read" ,"< paste curve_time.txt curve_time_py.txt" using 1:($237/$76) with lp pt 3 lc "blue" title " " ,"curve_time.txt" using 1:($5/$77) with lp pt 1 lc "dark-orange" title "Refine", "< paste curve_time.txt curve_time_py.txt" using 1:($238/$77) with lp pt 3 lc "dark-blue" title " " 
set output "curves_time_dam_r.tex"
set key right top horizontal inside
set xrange [1.:128.1]
set yrange [0.1:23.]
set logscale y
set ytics ("%.3f (6)" 0.16667,"%.1f (2)" 0.5,1.,7.,22.)
plot "curve_time.txt" using 1:($7/$79):xtic(1) with lp pt 1 lc "orange" notitle,"< paste curve_time.txt curve_time_py.txt" using 1:($240/$79) with lp pt 3 lc "blue" notitle
set output "curves_time_elem_vect_r.tex"
#set title "Elementary creation vect"
#unset ylabel
#set ylabel "Elapse time in s"
#set key left bottom inside vertical
set xrange [1.:128.1]
set yrange [0.15:3.]
set ytics ("%.3f (2.2)" 0.4545,"%.3f (1.6)" 0.625,"%.3f (1.2)" 0.833,1.,2.,2.7)
plot "curve_time.txt" using 1:($13/$83) with lp pt 1 lc "orange" notitle ,"curve_time.txt" using 1:($31/$83) with lp pt 2 lc "orange" notitle  ,"curve_time.txt" using 1:($49/$83) with lp pt 3 lc "orange" notitle,"curve_time.txt" using 1:($67/$83) with lp pt 0 lc "orange" notitle  ,"curve_time.txt" using 1:($99/$83):xtic(1) with lp pt 5 lc "green" notitle,"curve_time.txt" using 1:($179/$83):xtic(1) with lp pt 4 lc "dark-green" notitle,"curve_time.txt" using 1:($195/$83) with lp pt 5 lc "dark-green" notitle,"curve_time.txt" using 1:($211/$83) with lp pt 4 lc "forest-green" notitle,"curve_time.txt" using 1:($227/$83) with lp pt 5 lc "forest-green" notitle
set output "curves_time_elem_mat_r.tex"
#set title "Elementary creation matrix"
#unset ylabel
set xrange [1.:128.1]
set yrange [0.88:1.6]
set ytics ("%.3f (1.06)" 0.94339,1.,1.06,1.2,1.45)
plot "curve_time.txt" using 1:($14/$84) with lp pt 1 lc "orange" notitle ,"curve_time.txt" using 1:($32/$84) with lp pt 2 lc "orange" notitle ,"curve_time.txt" using 1:($50/$84) with lp pt 3 lc "orange" notitle ,"curve_time.txt" using 1:($68/$84) with lp pt 0 lc "orange" notitle ,"curve_time.txt" using 1:($100/$84):xtic(1) with lp pt 5 lc "green" notitle,"curve_time.txt" using 1:($180/$84):xtic(1) with lp pt 4 lc "dark-green" notitle,"curve_time.txt" using 1:($196/$84) with lp pt 5 lc "dark-green" notitle,"curve_time.txt" using 1:($212/$84) with lp pt 4 lc "forest-green" notitle,"curve_time.txt" using 1:($228/$84) with lp pt 5 lc "forest-green" notitle
set output "curves_time_stress_strain_r.tex"
#set title "Outputs"
set xrange [1.:128.1]
set yrange [0.29:1.7]
set logscale y
set ytics ("%.2f (3.)" 0.333,"%.2f (2.7)" 0.37,"%.1f (2.)" 0.5,1.,1.65)
plot "curve_time.txt" using 1:($19/$89) with lp pt 1 lc "orange" notitle ,"curve_time.txt" using 1:($37/$89) with lp pt 2 lc "orange" notitle ,"curve_time.txt" using 1:($55/$89) with lp pt 3 lc "orange" notitle ,"curve_time.txt" using 1:($73/$89) with lp pt 0 lc "orange" notitle ,"< paste curve_time.txt curve_time_py.txt" using 1:($248/$89) with lp pt 3 lc "blue" notitle ,"curve_time.txt" using 1:($105/$89):xtic(1) with lp pt 5 lc "green" notitle
set output "curves_time_nl_r.tex"
#set title "Non-linear loop"
#set ylabel "Elapse time in s"
set key center bottom vertical inside
set xrange [1.:128.1]
set yrange [0.65:1.24]
set logscale y
set ytics ("%.3f (1.35)" 0.7407,"%.3f (1.17)" 0.8547,"%.3f (1.1)" 0.90909,"%.3f (1.05)" 0.95238,1.,1.05,1.1,1.15,1.2)
plot "curve_time.txt" using 1:($17*5/$87/7) with lp pt 1 lc "orange" notitle  ,"curve_time.txt" using 1:($35*5/$87/7) with lp pt 2 lc "orange" notitle  ,"curve_time.txt" using 1:($53*5/$87/7) with lp pt 3 lc "orange" notitle  ,"curve_time.txt" using 1:($71*5/$87/7) with lp pt 0 lc "orange" notitle  ,"< paste curve_time.txt curve_time_py.txt" using 1:($246*5/$87/7) with lp pt 3 lc "blue" notitle  ,"curve_time.txt" using 1:($103/$87):xtic(1) with lp pt 5 lc "green" notitle,"curve_time.txt" using 1:($183/$87) with lp pt 4 lc "dark-green" notitle,"curve_time.txt" using 1:($199/$87) with lp pt 5 lc "dark-green" notitle,"curve_time.txt" using 1:($215/$87) with lp pt 4 lc "forest-green" notitle,"curve_time.txt" using 1:($231/$87) with lp pt 5 lc "forest-green" notitle
set output "curves_time_output_r.tex"
#set title "Outputs"
set key right top horizontal inside
set xrange [1.:128.1]
set yrange [0.075:12]
set logscale y
set ytics ("%.3f (11.3)" 0.088,"%.1f (5)" 0.2,"%.1f (2)" 0.5,1.,4.,10.)
plot "curve_time.txt" using 1:($18/$88):xtic(1) with lp pt 1 lc "orange" notitle,"< paste curve_time.txt curve_time_py.txt" using 1:($247/$88) with lp pt 3 lc "blue" notitle
unset ytics
set output "curves_time_misc_r.tex"
#set title "Small elapse times"
#unset ylabel
set key center top vertical
set xrange [1.:128.1]
set yrange [0.04:2090.]
set ytics
plot "curve_time.txt" using 1:($3/$75):xtic(1) with lp pt 1 lc "orange" title "Initialize" ,"< paste curve_time.txt curve_time_py.txt" using 1:($236/$75):xtic(1) with lp pt 3 lc "blue" title " " , "curve_time.txt" using 1:($8/$80) with lp pt 1 lc "dark-orange" title "Material constant" , "< paste curve_time.txt curve_time_py.txt" using 1:($177/$80) with lp pt 3 lc "dark-blue" title " " 
set output "curves_time_obj_r.tex"
#set title "space,form solver elapse times"
#unset ylabel
set key center bottom vertical
set xrange [1.:128.1]
set yrange [0.1:11000]
plot "curve_time.txt" using 1:($6/$78):xtic(1) with lp pt 1 lc "orange" notitle ,"< paste curve_time.txt curve_time_py.txt" using 1:($239/$78) with lp pt 3 lc "blue" notitle , "curve_time.txt" using 1:($15/$85) with lp pt 1 lc "dark-orange" notitle,"< paste curve_time.txt curve_time_py.txt" using 1:($244/$85) with lp pt 3 lc "dark-blue" notitle , "curve_time.txt" using 1:($16/$86) with lp pt 1 lc "red" notitle , "< paste curve_time.txt curve_time_py.txt" using 1:($245/$86) with lp pt 3 lc "cyan" notitle
set output "curves_time_BC_r.tex"
#set title "BC elapse times"
#unset ylabel
set key left bottom vertical
set xrange [1.:128.1]
set yrange [0.5:28]
plot "curve_time.txt" using 1:($9/$81):xtic(1) with lp pt 1 lc "orange" title "Dirichlet" ,"< paste curve_time.txt curve_time_py.txt" using 1:($242/$81) with lp pt 3 lc "blue" title " " , "curve_time.txt" using 1:($10/$82) with lp pt 1 lc "dark-orange" title "Neuman","< paste curve_time.txt curve_time_py.txt" using 1:($243/$82) with lp pt 3 lc "dark-blue" title " "  
exit
