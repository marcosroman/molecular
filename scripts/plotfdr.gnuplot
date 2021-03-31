#opcion 1:
#
#plot "d0.5fdr.data" u (column(1)/8.703555e+00):2 w lp lw 2, "d0.8fdr.data" u (column(1)/7.441435e+00):2 w lp lw 2, "d1fdr.data" u (column(1)/6.908016e+00):2 w lp lw 2, "d1.2fdr.data" u (column(1)/6.500692e+00):2 w lp lw 2, "d2fdr.data" u (column(1)/5.482896e+00):2 w lp lw 2

#opcion 2:
#->se grafica distancia normalizada vs fraccion de particulas

# normalized (both in x and y) radial distribution (frd) plots 
folder="examples/data/"

set terminal png enhanced size 1600,1200 font "arial,32"
set output "examples/fdr-radialdistfunct.png"
set title "Radial distribution function"
set ylabel "% of particles"
set xlabel "(Normalized) Length"
plot folder."d0.5fdr.data" u (column(1)/8.703555e+00):(column(2)/60.4113) w lp ps 3 lw 3 t "d=0.5", folder."d1fdr.data" u (column(1)/6.908016e+00):(column(2)/58.8892) w lp ps 3 lw 3 t "d=1", folder."d2fdr.data" u (column(1)/5.482896e+00):(column(2)/58.8892) w lp ps 3 lw 3 t "d=2"


