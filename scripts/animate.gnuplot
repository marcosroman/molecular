idxmax=pasos-pasosterm
idxstep=10

set xrange [0:lado_caja]
set yrange [0:lado_caja]
set zrange [0:lado_caja]
set border 4095
set ticslevel 0
set view 70,40
unset tics
set terminal png enhanced
# i need indices to increase by 
i=0
do for [idx=0:idxmax:idxstep] {
	set output outputfolder."/".i.".png"
	set title sprintf("density=%g\n(t=%g)",densidad,(0+idx)*dt)
	splot gpltposdatafile i idx pt 'o' notitle
	i=i+1
}
