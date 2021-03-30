#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <getopt.h>

#define KI   0
#define VI   1
#define TI   2
#define PI   3
#define KP   4
#define VP   5
#define TP   6
#define PP   7
#define G    8
#define D    9
#define PRMS 10
#define TERM 11
#define INFO 12
#define SAVE 13
#define SPOS 13
#define NSWITCHES	15
#define NVARS 8

#define X 0
#define Y 1
#define Z 2

typedef double vector[3];

static int flag;

void imprimir_uso(char *av0) {
  fprintf(stderr,"Uso: %s\n\n      parametros requeridos:\n\t-m\t\tdefine cantidad de particulas N=4*m^3\n\t-t, --temp\ttemperatura\n\t-d, --dens\tdensidad\n\t-h, --dt\tpaso temporal\n\t-p, --pasos\tcantidad total de pasos a simular\n\t-s, --pasosterm\tcantidad de pasos de termalizacion [s<=p]\n\t-e, --cadapterm\tcada cuantos pasos de termalizacion ajustamos velocidades [cadapterm<=pterm]\n\n      evolucion temporal del sistema:\n\t--prk(i|p)\timprimir energia cinetica instantanea|promedio\n\t--prv(i|p)\timprimir energia potencial instantanea|promedio\n\t--prt(i|p)\timprimir temperatura instantanea|promedio\n\t--prp(i|p)\timprimir presion instantanea|promedio\n\t--prall(i|p)\timprimir todos los valores instantaneos|promedio\n\t--prall\t\timprimir todos los valores\n\t--prterm\timprimir etapa de termalizacion\n\n      funcion de distribucion radial:\n\t--fdr\t\tcalcular distribucion de funcion radial (flag)\n\t-b, --bins\tcant. de bins p/ la fdr\n\n      difusion:\n\t--dif\t\tcalcular difusion (flag)\n\t-1, --cadad\tcada cuantos pasos de evolucion se calculan distancias\n\t-2, --it0\tdefine cada cuantas veces que se calculan distancias\n\t\t\t se genera un nuevo origen de tiempos (cada it0*cadad pasos de evolucion)\n\t-3, --t0max\tcantidad de tiempos iniciales distintos a partir de los cuales se promedian distancias\n\t-4, --tmax\tcantidad de registros temporales para la curva distancia vs. tiempo\n\n      opcionales:\n\t--info\tincluir info de parametros y que se imprime\n\t-c, --rcorte\tradio de corte\n\t-w, --guardar\tguardar datos de {evolucion,fdr,difusion} (en [prefijo]{run,fdr,dif}.data)\n\t--pos\tguardar posiciones\n\n",av0);
}

void serror(char *s) {
	fprintf(stderr,"%s\n",s);
	exit(EXIT_FAILURE);
}

void argumentos(
    int ac, char *av[], int *s, int *m, int *p, int *pt, int *ct,
    double *d, double *T,  double *dt, double *rc,
		int *bins,
		int *it0, int *t0max, int *tmax, int *cadad,
		char *fnameprefix) {
	// switches|flags para verificar que parametros fueron dados
  int ssm=0, sst=0, ssd=0, ssp=0, ssdt=0, ssct=0, sspt=0, src=0; // (parametros)
	int sfdr=0, sbins=0; // (funcion de dist. radial)
	int sdif=0, sit0=0, st0max=0, stmax=0, scadad=0; // (difusion)
	int c,j;

	for (j=0;j<NSWITCHES;j++) s[j]=0; // switches=0

	// barremos sobre los argumentos hasta que se acaben
	while (1) {
		static struct option long_options[] = {
			{"m", required_argument, 0, 'm'},
			{"temp", required_argument, 0, 't'},
			{"dens", required_argument, 0, 'd'},
			{"dt", required_argument, 0, 'h'},
			{"pasos", required_argument, 0, 'p'},
			{"pasosterm", required_argument, 0, 's'},
			{"cadapterm", required_argument, 0, 'e'},
			{"guardar", optional_argument, 0, 'w'},
			{"pos", optional_argument, 0, 'x'},
			{"info", no_argument, &flag, INFO},
			{"rcorte", optional_argument, 0, 'c'},
			{"prterm", no_argument, &flag, TERM},
			{"pralli", no_argument, &flag, -1},
			{"prki", no_argument, &flag, KI},
			{"prvi", no_argument, &flag, VI},
			{"prti", no_argument, &flag, TI},
			{"prpi", no_argument, &flag, PI},
			{"prallp", no_argument, &flag, -2},
			{"prkp", no_argument, &flag, KP},
			{"prvp", no_argument, &flag, VP},
			{"prtp", no_argument, &flag, TP},
			{"prpp", no_argument, &flag, PP},
			{"prall",  no_argument, &flag, -3},
			{"fdr",  no_argument, &flag, -111},
			{"bins", optional_argument, 0, 'b'},
			{"dif", no_argument, &flag, -222},
			{"it0", optional_argument, 0,'1'},
			{"t0max", optional_argument, 0,'2'},
			{"tmax", optional_argument, 0,'3'},
			{"cadad", optional_argument, 0,'4'},
			{0, 0, 0, 0}
		};
		int option_index = 0;

		c = getopt_long (ac, av, "m:t:d:h:p:s:e:c::w::b::1::2::3::4::x::",
										 long_options, &option_index);

		if (c == -1) break; // fin de opciones

		switch (c) {
			case 0: // flags
				if(flag>=0) {
					s[flag]=1;
				} 
				if(flag==-1) {
					s[KI]=1; s[VI]=1;
					s[TI]=1; s[PI]=1;
				}
				if (flag==-2) {
					s[KP]=1; s[VP]=1;
					s[TP]=1; s[PP]=1;
				}
				if (flag==-3) {
					s[KI]=1; s[VI]=1;
					s[TI]=1; s[PI]=1;
					s[KP]=1; s[VP]=1;
					s[TP]=1; s[PP]=1;
				} 
				if (flag==-111) {
					// calculamos fdr
					s[G]=1;
					sfdr=1;
				}
				if (flag==-222) {
					// calculamos difusion
					s[D]=1;
					sdif=1;
				}
				break;
			case 'm': // (define nro de particulas)
				*m=atoi(optarg);
				ssm=1;
				break;
			case 't': // temperatura
				*T=atof(optarg);
				sst=1;
				break;
			case 'd': // densidad
				*d=atof(optarg);
				ssd=1;
				break;
			case 'h': // deltat, paso temporal
				*dt=atof(optarg);
				ssdt=1;
				break;
			case 'p': // pasos
				*p=atoi(optarg);
				ssp=1;
				break;
			case 's': // pasos de termalizacion
				*pt=atoi(optarg);
				sspt=1;
				break;
			case 'e': // intervalo de termalizacion
				*ct=atoi(optarg);
				ssct=1;
				break;
			case 'c': // radio de corte (rc)
				*rc=atof(optarg);
				src=1;
				break;
			case 'w': // guardar a archivo
				// si hay argumento, uso eso como prefijo para nombre de archivo
				if(optarg!=NULL)
					strcpy(fnameprefix,optarg);
				else
					strcpy(fnameprefix,"");
				s[SAVE]=1;
				break;
			case 'x': // guardar a archivo
				s[SPOS]=1;
				break;
			case 'b': // cant. de bins (necesario solo si se pide calcular fdr)
				if (optarg!=NULL) {
					*bins=atoi(optarg);
					sbins=1;
				} else serror("Indicar bins de la forma -b<nro_bins>, sin espacios, o --bins=<nro_bins>");
				break;
			case '1': // 
				if (optarg!=NULL) {
					*cadad=atoi(optarg);
					scadad=1;
				} else serror("Indicar cadad de la forma -4<valor>, sin espacios, o --cadad=<valor>");
				break;
			case '2': // 
				if (optarg!=NULL) {
					*it0=atoi(optarg);
					sit0=1;
				} else serror("Indicar it0 de la forma -1<valor>, sin espacios, o --it0=<valor>");
				break;
			case '3': // 
				if (optarg!=NULL) {
					*tmax=atoi(optarg);
					stmax=1;
				} else serror("Indicar tmax de la forma -3<valor>, sin espacios, o --tmax=<valor>");
				break;
			case '4': // 
				if (optarg!=NULL) {
					*t0max=atoi(optarg);
					st0max=1;
				} else serror("Indicar t0max de la forma -2<valor>, sin espacios, o --t0max=<valor>");
				break;
			case '?': // ?
				imprimir_uso(av[0]);
				serror("Opcion desconocida");
				break;
			default:
				abort();
		}
	}

	// verificamos parametros necesarios
  if(sst&&ssd&&ssm&&ssp&&ssdt&&ssct&&sspt) {
		// verificamos que pasos>=pasosterm
		if((*pt)>(*p)) serror("Se debe satisfacer pasos > pasosterm");

		// verificamos que se pide imprimir o calcular algo
		c=0;
		for (j=0;j<TERM;j++) c+=s[j];
		c+=s[SPOS]; // imprimir posiciones cuenta como 'algo'
		// contador hasta TERM para no contar el switch TERM (impr. term.)
		if(c==0) {
			imprimir_uso(av[0]);
			serror("Sin flags, nada para calcular|imprimir\n\n");
		}

		// si no se escogio radio de corte, ponemos rcorte=0.5*(lado de caja)/sqrt(3)
		// que es la max distancia entre particulas (en caja con cond. periodicas)
		// (lado de caja depende del nro de particulas y la densidad deseada)
		if (!src)
			(*rc) = (0.5/sqrt(3))*pow(4/(*d),1/3.)*(*m);

		// verificamos que --fdr & --bins juntos
		if(sfdr) {
			if (!sbins) serror("Seleccionar bins");
			else if (!((*bins)>0)) serror("Seleccionar nro. de bins (-b|--bins) > 0");
		}	else {
			if (sbins) serror("Flag --fdr no seleccionado");
		}

		// difusion: checkear aqui it0,tmax,etc
		if(sdif) {
			if(!(sit0&&stmax&&st0max&&scadad)) 
				serror("Indicar valores it0,t0max,tmax,cadad p/ calcular difusion\n");
			else {
				if((*cadad)>0 && (*it0)>0 && (*t0max)>0 && (*t0max)>0) {
					if((*cadad)*(*it0)*(*t0max)>(*p)-(*pt)) 
					serror("Se debe satisfacer la inecuacion cadad*it0*t0max <= pasos-pasosterm");
					else if( (*t0max)*(*it0)<(*tmax))
						serror("Se debe cumplir t0max*it0 > tmax");
				} else serror("Se debe cumplir que cadad,t0max,it0,tmax>0");
			}
		} else if(sit0||stmax||st0max||scadad) {
			serror("Usar flag --dif p/ calcular difusion (y especificar it0,t0max,tmax,cadad)\n");
		}
	} else {
		imprimir_uso(av[0]);
		serror("Parametros insuficientes");
  }
}

void iniciar_ctes(int *m, int *n, double *d, double *lcel, double *lcaj) {
  (*n) = 4*((*m)*(*m)*(*m));
  (*lcel) = pow(4/(*d),1/3.);
  (*lcaj) = (*lcel)*(*m);
}    

void dependencias(int *s, int *cs) {
  int k;

  for(k=0;k<NSWITCHES;k++) {
    cs[k]=s[k];
  }

  // si se calcula algun promedio, cs[PRMS]=1
  if(s[KP]||s[VP]||s[TP]||s[PP]) cs[PRMS]=1;
	// si se pidio calcular promedios, nos aseguramos de calcular valores inst.
  if(cs[PRMS]) {
    if(s[PP]) {
      cs[PP]=1; cs[PI]=1; cs[KP]=1; cs[PP]=1;
    }
    if(s[TP]) {
      cs[TP]=1; cs[TI]=1; cs[KP]=1; cs[KI]=1;
    }
    if(s[KP]) {
      cs[KP]=1; cs[KI]=1;
    }
  }
}

void iniciar_valores_inst(int *s, double *ki, double *vi, double *w) {
  if(s[KI]) *ki=0;
  if(s[VI]) *vi=0;
  if(s[PI]) *w=0;
}

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

// ran1, de Numerical Recipes in C
double ran1(long *idum) {
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

// gasdev, de Numerical Recipes in C
double gasdev(long *idum) {
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if(*idum<0) iset=0;
	if(iset==0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

void iniciar_posiciones(int N, vector *p, double l) {
	int m, i, j, k, q;
	vector pos;
  vector b[4] = {
		{0.0,0.0,0.0},
		{0.5,0.5,0.0},
		{0.5,0.0,0.5},
		{0.0,0.5,0.5}
	};

  for(i=0;i<N;i++) 
    for(j=0;j<N;j++) 
      for(k=0;k<N;k++) {
        pos[X]=i; pos[Y]=j; pos[Z]=k;
        for(q=0;q<4;q++,p++) 
          for(m=0;m<3;m++) 
						(*p)[m]=(pos[m]+b[q][m])*l;
			}
}

void iniciar_velocidades(int N, vector *v, double T) {
	int i, k;
	static long int s=0;
	double rt=sqrt(T);

	if(!s) {
		srand(time(NULL));
		s=-(long)rand();
	}

	for(i=0;i<N;i++,v++) {
		for(k=0;k<3;k++) (*v)[k] = gasdev(&s)*rt;
	}
}

void imprimir_info(int *s, int m, int p, int pt, int ct,
    double d, double T,  double dt, double rc, double l, FILE *fout) {
	int i;
	char prtlabels[40] = "";
	char *varlabel[] = {
	"ki", "vi", "ti", "pi",
	"kp", "vp", "tp", "pp"
	};

	// imprimo parametros
	fprintf(fout,"# m = %d\n", m);
	fprintf(fout,"# N = %d\n", 4*m*m*m);
	fprintf(fout,"# temperatura = %g\n", T);
	fprintf(fout,"# densidad = %g\n", d);
	fprintf(fout,"# lado_caja = %g\n", l);
	fprintf(fout,"# rcorte = %g\n", rc);
	fprintf(fout,"# dt = %g\n", dt);
	fprintf(fout,"# pasos = %d\n", p);
	fprintf(fout,"# pasosterm = %d\n", pt);
	fprintf(fout,"# cadapterm = %d\n", ct);
	// imprimo que se imprime de la simulacion
	for (i=0; i<NVARS; i++) {
		if(s[i]) strcat(prtlabels,", ");
		if(s[i]) strcat(prtlabels,varlabel[i]);
	}
	fprintf(fout,"# simulacion :: imprimiendo: p, t%s\n",prtlabels);
}

// distancia (minima) entre dos particulas
// (teniendo en cuenta condiciones periodicas de contorno en caja de lado l)
double min_dist2(double *p_i, double *p_j, double l) {
	int k;
	double s=0;
	double d;

	for(k=0;k<3;k++) {
		d=p_i[k]-p_j[k];

		if(fabs(d)>0.5*l)
			d-=(d>=0?l:-l);

		s+=d*d;
	}

	return s;
}

double dist(double *p_i, double *p_j) {
	int k;
	double s=0;

	for(k=0;k<3;k++,p_i++,p_j++)
		s+=(*p_i-*p_j)*(*p_i-*p_j);

	return sqrt(s);
}

void sumar_fuerza(
		double *p_i, double *p_j,
		double *f_i, double *f_j,
		double *w_inst, double *energia_p_inst,	double l, 
		int s_fuerza, int s_presion_inst, int s_ep_inst,
		double rcorte) {
	double rc2 = rcorte*rcorte;
	int k;
	double r2, rm2, rm6, f;
	vector d;

	// calculamos distancia minima entre el par
	// (considerando condiciones periodicas de contorno)
	for(k=0;k<3;k++) {
		d[k]=p_i[k]-p_j[k];
		if(fabs(d[k])>0.5*l)
			d[k]-=(d[k]>=0?l:-l);
	}
	r2=0;
	for(k=0;k<3;k++)
		r2+=d[k]*d[k];

	if(r2>=rc2) return;
	if(r2==0) r2=1e-6;

	rm2=1/r2;
	rm6=rm2*rm2*rm2;
	f=48*(rm6*(rm6-0.5));

	if(s_fuerza) {
		if(s_presion_inst) *w_inst+=f;
	  if(s_ep_inst) *energia_p_inst+=4*rm6*(rm6-1);
	}

	f*=rm2;

	for(k=0;k<3;k++,f_i++,f_j++) {
	  (*f_i) += f*d[k];
		(*f_j) -= f*d[k];
	} 
}

int verificar_condcontorno(int N, vector *p, double l) {
	int i, k;

	for(i=0;i<N;i++,p++) {
		for(k=0;k<3;k++) {
			if((*p)[k]>l) (*p)[k]-=l;
			else if ((*p)[k]<0) (*p)[k]+=l;

			if((*p)[k]>l||(*p)[k]<0) return 1;
		}
	}

	return 0;
}

double calcular_energia_k(int N, vector *v) {
	int i,k;
	double s=0;

	for(i=0;i<N;i++,v++)
		for(k=0;k<3;k++) s+=(*v)[k]*(*v)[k];

	return 0.5*s;
}

void termalizar(int N, vector *v, double T, double T_inst) {
	int i,k;
	double a=sqrt(T/T_inst);

	for(i=0;i<N;i++,v++)
		for(k=0;k<3;k++) (*v)[k]*=a;
}

// funcion de distribucion radial
// esta funcion tiene 2 usos, dependiendo del valor de b; si b>0,
// se utiliza para inicializar las variables y hacer memoria;
// de otra forma, se utiliza para imprimir los valores de la fdr
void fdistradial(int N, vector *p, double l, int b, int *switches, char *fnameprefix) {
	static double *distribucion;
	static long int n=0;
	static int gset=0, bins=0;
	static double hl, dr;
	double v, nn;
	int i,j;
	FILE *outputfile=stdout;
	char outputfname[120]="";

	if(b>0) {
		if(!gset) {
			hl=l*0.5;
			bins=b;
			//distribucion=(double*)malloc(sizeof(double)*(bins+1));
			distribucion=(double*)calloc(bins+1,sizeof(double));
			if(distribucion==NULL) serror("Memoria insuficiente. Reduzca numero de bins.");
			dr=hl/((double)bins);
			gset=1;
		}
		else {
			n++;
			// recorremos pares de particulas
			for(i=0;i<N-1;i++) {
				for(j=i+1;j<N;j++) {
					//printf("saving... ");
					//bintosaveto=(int)floor(sqrt(min_dist2(p[i],p[j],l))/dr);
					distribucion[(int)floor(sqrt(min_dist2(p[i],p[j],l))/dr)]+=2.;
					//distribucion[bintosaveto]+=2.;
					//printf("...saved\n");
				}
			}
		}
	} else {
		if(switches[SAVE]) { // si se pide imprimir en archivo...
			strcpy(outputfname,fnameprefix); strcat(outputfname,"fdr.data");
			outputfile=fopen(outputfname,"w");
			if(outputfile==NULL)
				serror("Error el abrir archivo para guardar fdr.");
		}
		
		if(switches[INFO]) // imprimir info si se pide
			fprintf(outputfile,"# funcion de distribucion radial (bins=%d)\n",bins);

		for(i=1;i<=bins;i++) {
			v=(pow(i+1,3)-pow(i,3))*pow(dr,3);
			nn=(4/3.)*PI*v*l;
			
			fprintf(outputfile,"%e %e\n",dr*(i+0.5),(distribucion[i-1])/(n*nn*N));
		}

    free(distribucion);
		// cerrar archivo (si se abrio)
		if(switches[SAVE]) fclose(outputfile);
	}
}

void difusion(int N, vector *x,
		int it0, int t0max, int tmax, int cadad, double dtime,
		int *switches, char *fnameprefix) {
	static int
		dset=0, // una vez que se inicializa la funcion, dset=1
		ntel=0, // contador de llamadas a esta funcion (una vez iniciada)
		t0=0,
		*time0;
	static double *ntime, *r2t;
	static vector **x0;
	int i, k, l, tt0, delt;
	FILE *outputfile=stdout;
	char outputfname[120]="";

	// inicializar variables para calcular difusion vs tiempo
	// ntel=0,r2t[0:tmax-1]=0,ntime[0:tmax-1],dtime=dt*cadad
	//	(se muestrea cada cadad pasos) 
	
	// if (x!=NULL) { if(!dset) iniciar(); else muestrear() } else imprimir();
	if(x!=NULL) {
		if(!dset) { // si es la primera vez que se llama a esta funcion...
			// hacemos memoria para vectores double 0:tmax-1; ntime, r2t
			ntime=(double *)malloc(sizeof(double)*(tmax+1));
			r2t=(double *)malloc(sizeof(double)*(tmax+1));	
			time0=(int *)malloc(sizeof(int)*(t0max+1));
			// ahora para las posiciones (y velocidades) hasta cierto tiempo limitado
			x0=malloc(sizeof(vector*)*(N+1));
			for(i=0;i<N;i++) {
				x0[i]=malloc(sizeof(vector)*(t0max+1));
				if(x0[i]==NULL) serror("Error al alocar memoria. Disminuir t0max");
			}

			if(ntime==NULL||r2t==NULL)
				serror("No se ha podido alocar memoria. Reducir tmax o algun otro parametro.");
			dset=1;
		}
		else /* muestreamos */ {
			// nuevo origen de tiempos, cada cadad*it0 pasos
			if(ntel%it0==0) {
				tt0=t0%t0max;
				time0[tt0]=ntel;
				// posiciones 
				for(i=0;i<N;i++) {
					for(k=0;k<3;k++) {
						x0[i][tt0][k]=x[i][k];
					}
				}
				t0++;
			}
			// se actualizan, hasta l<min(t0,t0max)
			for(l=0;l<((t0<t0max)?t0:t0max);l++) {
				delt=ntel-time0[l];
				if(delt<tmax) {
					ntime[delt]++;
					for(i=0;i<N;i++) {
						r2t[delt]+=dist(x[i],x0[i][l]);
					}
				}
			}
			ntel++;
		}
	}
	else {
		if(switches[SAVE]) { // si se pide imprimir en archivo...
			strcpy(outputfname,fnameprefix); strcat(outputfname,"dif.data");
			outputfile=fopen(outputfname,"w");
			if(outputfile==NULL)
				serror("Error el abrir archivo para guardar dif.");
		}
		
		if(switches[INFO]) // imprimir info si se pide
			fprintf(outputfile,
					"# difusion vs. t (it0=%d, t0max=%d, tmax=%d, cadad=%d)\n",
					it0,t0max,tmax,cadad);

		for(i=0;i<tmax;i++) {
			r2t[i]/=cadad*ntime[i];
			fprintf(outputfile,"%e\t%e\n",dtime*(i+0.5), r2t[i]);
		}
		// cerrar archivo
		free(ntime);
		free(r2t);
		free(time0);
		for(i=0;i<N;i++) {
			free(x0[i]);
		}
		free(x0);
	}
}

void guardar_posiciones(int n, double t, vector *p, FILE *f) {
	int i,k;
	//FILE *f = fopen("datapos","w"); 
	//FILE *pipe = popen("gnuplot -p","w");

	fprintf(f,"%e",t);
	for(i=0;i<n;i++) {
		for(k=0;k<3;k++) {
			fprintf(f," %e",p[i][k]);
		}
	}
	fprintf(f,"\n");
	//fclose(f);

	//fprintf(pipe, "splot \"datapos\" pt 7\n");
	//fclose(pipe);
}

