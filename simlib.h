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
#define MODO 10
#define PRMS 11

#define MD   1
#define MC   2

#define X 0
#define Y 1
#define Z 2

#define ERRMSG_DT "Reducir dt."

//temporalmente... (porque lo que se quiere luego es que sea un parametro!!!)
#define R_CORTE 2.5

typedef double vector[3];

void imprimir_uso(char *av0) {
	//while((c = getopt(ac,av,"n:t:d:p:h:i:q:c:m:z::"))!=-1) {
  fprintf(stderr, "Uso: %s\n\t-m (define nro. de particulas N=4*m^3)\n\t-temp|-t (temperatura)\n\t-dens|-d (densidad)\n\t-dt|-h (delta tiempo)\n\t-pasos|-p (pasos)\n\t-pasosterm|-pt (pasos de termalizacion) [pt<=p]\n\t-cadapterm|-e (cada cuantos pasos de termalizacion ajustamos velocidades) [cadapterm<=pterm]\n\t-c (radio de corte)\n",av0);//\t-m c|d (modo de simulacion: mc|md)\n\t-z* (dr, opcional)\n\t-b* (nro bins para dist radial g(r),opcional)\n",av0);
}

void parametros(
    int ac, char *av[],
    int *s, int *m, int *p, int *pt, int *ct,
    double *d, double *T,  double *dt, double *rc,
    double *dr) {
	// switches|flags para verificar que parametros fueron dados
  int ssn=0, sst=0, ssd=0, ssp=0, ssdt=0, ssct=0, sspt=0, ssrc=0, ssdr=0;
  char st[30];
	int c;

	while (1) {
		static struct option long_options[] =
			{
				/* These options don’t set a flag.
					 We distinguish them by their indices. */
				{"m",  required_argument, 0, 'm'},
				{"temp",  required_argument, 0, 't'},
				{"dens",    required_argument, 0, 'd'},
				{"dt",    required_argument, 0, 'h'},
				{"pasos",    required_argument, 0, 'p'},
				{"pasosterm",    required_argument, 0, 's'},
				{"cadapterm",    required_argument, 0, 'e'},
				{"rcorte",    optional_argument, 0, 'c'},
				//{"file",    optional_argument, 0, 'f'},
				{0, 0, 0, 0}
			};
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (ac, av, "m:t:d:h:p:s:e:c:",
										 long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1) break;

		switch (c) {
			/*
			case 0:
				if (long_options[option_index].flag != 0)
					break;
				printf ("option %s", long_options[option_index].name);
				if (optarg)
					printf (" with arg %s", optarg);
				printf ("\n");
				break;
			*/

			case 'm': // (define nro de particulas)
				*m=atoi(optarg);
				ssn=1;
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
				ssrc=1;
				break;
				/*
			case 'm': // modo (moleculardynamics o montecarlo)
				if(strcmp(optarg,"d")==0) {
					s[MODO]=MD;
				} else if (strcmp(optarg,"c")==0) {
          s[MODO]=MC;
				} else {
					fprintf(stderr,"Modo desconocido");
					abort();
				}
				break;
			case 'z': //dr
				*dr=atof(optarg);
				ssdr=1;
				break;
				*/
			case '?': // ?
				printf("what\n");
				imprimir_uso(av[0]);
				exit(1);
				break;
			default:
				imprimir_uso(av[0]);
				printf("\n");
				abort();
		}
	}

	// por el momento solo modo MD
	s[MODO]=MD;
	ssrc=1;// por ahora, testing

  if(sst&&ssd&&ssn&&ssp&&ssdt&&ssct&&sspt&&ssrc) {
    if(s[MODO]==MD) {
      sprintf(st, "Dinamica Molecular");
    } else if(s[MODO]==MC&&ssdr) {
      sprintf(st, "Monte Carlo");
    } else {printf("Indicar dr p/ MC\n");exit(1);} // parece q hace falta

    printf("########################################\n");
    printf("# Simulando en modo %s\n", st);
    printf("# ======================================\n");
    printf("# m = %d\n# (%d particulas)\n#\n", *m, 4*(*m)*(*m)*(*m));
    printf("# Temperatura = %g\n", *T);
    printf("# Densidad = %g\n", *d);
    printf("# Pasos totales = %d\n", *p);
    printf("# Paso de tiempo (dt) = %g\n", *dt);
    printf("# Pasos de termalizacion = %d\n", *pt);
    printf("# Intervalo de termalizacion = %d\n", *ct);
    printf("# Radio de corte=%g\n", *rc);
    printf("########################################\n");
  } else {
    fprintf(stderr, "Parametros insuficientes\n\n");
		imprimir_uso(av[0]);
    exit(1);
  }
}

void iniciar_ctes(int *m, int *n, double *d, double *lcel, double *lcaj) {
  (*n) = 4*((*m)*(*m)*(*m));
  (*lcel) = pow(4/(*d),1/3.);
  (*lcaj) = (*lcel)*(*m);
}    

void interruptores(int ac, char *av[], int *s) {
	int j;
	/*
  int i=1,j,k=0;
  char *cads[] = {
    "ki", "vi", "ti", "pi",
    "kp", "vp", "tp", "pp",
     "g", "d"
  };
	*/

  for(j=0;j<10;j++) s[j]=1; // todas on (para probar nomas)
	/*
  for(j=0;j<10;j++) s[j]=0;

  while(i<ac&&av[i][0]!='-') {
    for(j=0;j<10;j++) {
      if(!strcmp(av[i],cads[j])) {
        s[j]=1;
        k++;
      }
    }
    if(k!=i) {
      fprintf(stderr, "No se reconoce la opcion ingresada.\nLas opciones disponibles son: ");
      for(j=0;j<10;j++) fprintf(stderr, "%s ", cads[j]);
      fprintf(stderr,"\n");
      exit(EXIT_FAILURE);
    }
    i++;
  }
  if(k==0) {
    fprintf(stderr, "Debe escoger al menos una opcion: ");
    for(j=0;j<10;j++) fprintf(stderr, "%s ", cads[j]);
    fprintf(stderr,"\n");
    exit(EXIT_FAILURE);
  }
	*/
}

void dependencias(int *vs, int *s) {
  int k;

  for(k=0;k<12;k++) {
    if(k==MODO) continue;
    else s[k]=0;
  }
  
  if(vs[KP]||vs[VP]||vs[TP]||vs[PP]||vs[G]) s[PRMS]=1;
  if(s[PRMS]) {
    if(vs[PP]) {
      s[PP]=1;
      s[PI]=1;
      s[KP]=1;
      s[PP]=1;
    }
    if(vs[TP]) {
      s[TP]=1;
      s[TI]=1;
      if(s[MODO]==MD) {
        s[KP]=1;
        s[KI]=1;
      }
    }
    if(vs[KP]) {
      s[KP]=1;
      s[KI]=1;
    }
  }
  for(k=0;k<KP;k++) {
    if(vs[k]) s[k]=1;
  }
  for(k=8;k<MODO;k++) {
    if(vs[k]) s[k]=1;
  }
  
  if((s[MODO]==MC)&&(vs[KI]||vs[KP]||vs[G]||vs[TI]||vs[TP])) {
    fprintf(stderr, "Se escogio una opcion no compatible con el modo MC.\n");
    exit(EXIT_FAILURE);
  }
}

void iniciar_ctes_g(int ac, char *av[], int *b) {
  int i=1, c=0;

	//testing
	*b=10;

	/*
  while(i<ac) {
    if (!strcmp(av[i],"-b")) {
      if(i<ac-1&&av[i+1][0]!='-') {
        *b=atoi(av[++i]);
        c++;
      } else {
        printf("Especifique el valor para %s.\n",av[i]);
        exit(1);
      }
    break;
    }
    i++;
  }

  if(c!=1) {
    fprintf(stderr,"Especificar -b\n");
    exit(1);
  }
	*/
}

void iniciar_ctes_d(int ac, char *av[], int *it0, int *t0max, int *tmax, int *cadad) {
  int i=1, c=0;

	//testing
	*it0=4;
	*t0max=100;
	*tmax=300;
	*cadad=1;

	/*
  while(i<ac) {
    // it0
    if (!strcmp(av[i],"-t0")) {
      if(i<ac-1&&av[i+1][0]!='-') {
        *it0=atoi(av[++i]);
        c++;
      } else {
        printf("Especifique el valor para %s.\n",av[i]);
        exit(1);
      }
    // t0max
    } else if (!strcmp(av[i],"-t0m")) {
      if(i<ac-1&&av[i+1][0]!='-') {
        *t0max=atoi(av[++i]);
        c++;
      } else {
        printf("Especifique el valor para %s.\n",av[i]);
        exit(1);
      }
    // tmax
    } else if (!strcmp(av[i],"-tm")) {
      if(i<ac-1&&av[i+1][0]!='-') {
        *tmax=atoi(av[++i]);
        c++;
      } else {
        printf("Especifique el valor para %s.\n",av[i]);
        exit(1);
      }
    // cadad
    }  else if (!strcmp(av[i],"-cd")) {
      if(i<ac-1&&av[i+1][0]!='-') {
        *cadad=atoi(av[++i]);
        c++;
      } else {
        printf("Especifique el valor para %s.\n",av[i]);
        exit(1);
      }
    }
    i++;
  }

  if(c!=4) {
    fprintf(stderr, "Falta especificar alguno de los siguentes: -t0, -t0m, -tm, -cd\n");
    exit(1);
  }
	*/

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

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

void serror(char *s) {
	fprintf(stderr,"%s\n",s);
	exit(EXIT_FAILURE);
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

void cargar_posiciones(int N, vector *p, double l) {
	FILE *f;
	char fn[100], c;
	int i,n=0;

	printf("Archivo de posiciones: ('.' para posiciones.dat)\n");
	scanf("%s",fn);
	if(strcmp(fn,".")==0) sprintf(fn,"posiciones.dat");

	if((f=fopen(fn,"r"))==NULL) serror("No se puede abrir el archivo.");

	do {
		c=getc(f);
		if(c=='\n') n++;
	} while(c!=EOF);
	if(n!=N) serror("Numero de lineas y de particulas no concuerdan.");
	rewind(f);

	for(i=0;i<n;i++,p++) {
		fscanf(f,"%lf\t%lf\t%lf\n",&((*p)[X]),&((*p)[Y]),&((*p)[Z]));
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

void sumar_fuerza(
		double *p_i, double *p_j,
		double *f_i, double *f_j,
		double *w_inst, double *energia_p_inst,	double l, 
		int s_fuerza, int s_presion_inst, int s_ep_inst
		) {
	static double rc2 = R_CORTE*R_CORTE;
	int k;
	double r2=0, rm2, rm6, f;
	vector d;

	for(k=0;k<3;k++,p_i++,p_j++) {
		d[k]=*p_i-*p_j;

		if(fabs(d[k])>0.5*l)
			d[k]-=(d[k]>=0?l:-l);
	}

	for(k=0;k<3;k++) r2+=d[k]*d[k];
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

double min_dist2(double *p_i, double *p_j, double l) {
	int k;
	double s=0;
	vector d;

	for(k=0;k<3;k++,p_i++,p_j++) {
		d[k]=*p_i-*p_j;

		if(fabs(d[k])>0.5*l)
			d[k]-=(d[k]>=0?l:-l);

		s+=d[k]*d[k];
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

double pasitoMC(int N, int n, vector *p, double l, double T, double *dr, double *a) {
	int i, k;
	static int acep=0, reje=0, setdr=0;
	static double rcorte2=R_CORTE*R_CORTE, ds=0;
	double r2, rm6;
	vector *pn=&(p[n]),*p0=p;
	vector pp;
	double u=0, u0=0, du=0;
	static long int s=0;

	if(!s) srand(time(NULL));
	s=-(long)rand();

	if(!setdr) {
		ds=*dr;
		setdr=1;
	}
	else if(*dr!=ds) {
		acep=0;
		reje=0;
		ds=*dr;
	}

	for(i=0;i<N;i++,p++) {
		if(i!=n) {
		 	r2=min_dist2(*p,*pn,l);
			
			if(r2<rcorte2) {
				rm6=pow(r2,-3);
				u0+=4*rm6*(rm6-1);
			}
		}
	}

	for(k=0;k<3;k++) {
		pp[k]=(*pn)[k]+ds*(ran1(&s)-0.5);
		if(pp[k]>l) pp[k]-=l;
		else if (pp[k]<0) pp[k]+=l;
	}

	for(i=0,p=p0;i<N;i++,p++) {
		if(i!=n) {
		 	r2=min_dist2(*p,pp,l);
			if(r2<rcorte2) {
				rm6=pow(r2,-3);
				u+=4*rm6*(rm6-1);
			}
		}
	}

  du=u-u0;
	if(ran1(&s)<exp(-du/T)) {
		for(k=0;k<3;k++) (*pn)[k]=pp[k];
		acep++;
	} else {
		du=0;
		reje++;
	}

	*a=((double)acep)/((double)(acep+reje));

	return du;
}

double calcular_energia_p(int N, vector *p, double l) {
	double r2, rm6;
	double v=0;
	static double rc2=R_CORTE*R_CORTE;
	int i, j;

	for(i=0;i<N;i++) {
		for(j=i+1;j<N;j++) {
			r2=min_dist2(p[i],p[j],l);
			if(r2>=rc2) continue;
			rm6=pow(r2,-3);
			v+=rm6*(rm6-1);
		}
	}
	v*=4;

	return v;
}

double calcular_presion(int N, vector *p, double l, double T, double d) {
	int i, j;
	static double rc2=R_CORTE*R_CORTE;
	double r2, rm6, w=0;

	for(i=0;i<N;i++) {
		for(j=i+1;j<N;j++) {
			r2=min_dist2(p[i],p[j],l);
			if(r2<rc2) {
				rm6=pow(r2,-3);
				w+=rm6*(rm6-0.5);
			}
		}
	}
	w*=48;

	return d*(T+(w/(3*N)));
}

void fdistradial(int N, vector *p, double l, int b) {
	static double *distribucion;
	static long int n=0;
	static int gset=0, bins=0;
	static double hl, dr;
	double v, nn;
	int i,j;

	if(b>0) {
		if(!gset) {
			hl=l*0.5;
			bins=b;
			distribucion=(double*)malloc(sizeof(double)*(bins+1));
			if(distribucion==NULL) serror("Memoria insuficiente. Reduzca numero de bins.");

			gset=1;
			dr=hl/((double)bins);
		}
		else {
			n++;
			for(i=0;i<N-1;i++) {
				for(j=i+1;j<N;j++) {
					distribucion[(int)floor(sqrt(min_dist2(p[i],p[j],l))/dr)]+=2.;
				}
			}
		}
	} 
	else {
		printf("\n\n");
		for(i=1;i<=bins;i++) {
			v=(pow(i+1,3)-pow(i,3))*pow(dr,3);
			nn=(4/3.)*3.141593*v*l;
			
			printf("%e\t%e\n",dr*(i+0.5),(distribucion[i-1])/(n*nn*N));
		}
		// no entiendo por que me da segfault con estoooo!!!! si saco esta todo bien, supongo que no hace mucho dano al mundo
    //free(distribucion);
	}
}

/* parametros:
	 it0 (indica cada cuantos muestreos tengo un nuevo origen de tiempos),
	 tmax, indica la cantidad maxima de deltat (dentro de la rutina),
	 t0max, (cadad,timestep (al final))
*/
// quiero pasarle la posicion de las particulas (sin verificar condiciones de contorno, a partir de que finaliza la termalizacion)
// para el final necesito pasarle 'cadad'
void difusion(int N, vector *x, int it0, int t0max, int tmax, double dtime) {
	static int dset=0, ntel=0, t0=0, *time0;
	static double *ntime, *r2t;
	static vector **x0;
	int i, k, l, tt0, delt;

	// inicializar variables para coef de difusion...
	// ntel=0,r2t[0:tmax-1]=0,ntime[0:tmax-1],dtime=dt*cadad
	//	(se muestrea cada cadad pasos) 
	//
	if(x!=NULL) {
		if(!dset) {
			// alocar memoria para vectores double 0:tmax-1; ntime, vacf y r2t
			ntime=(double *)malloc(sizeof(double)*(tmax+1));
			r2t=(double *)malloc(sizeof(double)*(tmax+1));	
			time0=(int *)malloc(sizeof(int)*(t0max+1));
			// ahora para las posiciones (y velocidades) hasta cierto tiempo limitado
			x0=malloc(sizeof(vector*)*(N+1));
			for(i=0;i<N;i++) {
				x0[i]=malloc(sizeof(vector)*(t0max+1));
				if(x0[i]==NULL) serror("Error al alocar memoria. Disminuir t0max");
			}

			// VERIFICAR ALOCAMIENTO CORRECTO!!!!
			if(ntime==NULL||r2t==NULL)
				serror("No se ha podido alocar memoria. Reducir tmax o algun otro parametro.");
			dset=1;
		}
		// muestreamos cada 'cadad' pasos luego de la termalizacion
		// (en esta parte no existe referencia a 'cadad')
		else /**/ {
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
		printf("\n\n");
		for(i=0;i<tmax;i++) {
			// t0max<-cadad (p/ ahorrar unos cuantos bits)
			r2t[i]/=t0max*ntime[i];
			printf("%e\t%e\n",dtime*(i+0.5), r2t[i]);
		}
		free(ntime);
		free(r2t);
		free(time0);
		for(i=0;i<N;i++) {
			free(x0[i]);
		}
		free(x0);
	}
}

/*
void visualizar(vector *p) {
	int i,k;
	FILE *f = fopen("datapos","w"); 
	FILE *pipe = popen("gnuplot -p","w");

	for(i=0;i<NUM_PARTICULAS;i++,p++) {
		for(k=0;k<3;k++) {
			fprintf(f,"%e\t",(*p)[k]);
		}
		printf("\n");
	}
	fclose(f);

	fprintf(pipe, "splot \"datapos\" pt 7\n");
	fclose(pipe);
}*/


