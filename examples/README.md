Let's try to illustrate what 'molecular' does with a few examples.

We'll run some simulations with different parameters; this should result in different phases of the simulated systems (composed of bouncy but slighly pairwise-attractive particles in a cube with periodic boundary conditions).

From the main folder we run:

```
for d in 0.5 0.8 1 1.2 2; do
	./mol -m5 -t1.6 -d$d -h0.001 -p6000 -s1000 -e50 -wexamples/data/d$d --fdr -b100 --pos --info
done
```

We can visualize the effect of the parameters in the dynamics of the system by using the script (again from the main folder)
[This requires gnuplot (interactive plotting program) and ffmpeg (video converter)]

```
for d in 0.5 0.8 1 1.2 2; do
	./scripts/pos2mpg.sh examples/data d$d
	mv d$d.mpg examples
done
```

Also, the radial distribution function can be plotted (in this example, for 3 different densities while keeping all other parameters constant) can be plotted

```
gnuplot scripts/plotfdr.gnuplot
```

![radial distrubution function](examples/fdr-radialdistfunct.png)
