#include "simlib.h"

int main(int argc, char *argv[]) {
	int i, j, k, p;  // contadores
	// switches (p/ determinar que variables se calculan e imprimen)
  int switches[NSWITCHES], // p/ variables se piden en los argumentos
			// p/ variables que hacen falta calcular (una vez resueltas dependencias)
	    calcsws[NSWITCHES],
			// contadores de variables
			cv=0, cvi=0;
	// mas switches
	int s_fuerza, s_termalizacion, s_c_difusion=0;
  double **pv, *pa[NVARS]; // p/ guardar variables calculadas
	// parametros
  int m=0, n; 
  double densidad, temperatura, lado_celda, lado_caja;
  double dt, rcorte;
	// pasos: nro de pasos totales; p_termalizacion<pasos
	// cadapterm: nro de pasos que se saltean, durante pasos de termalizacion,
	//            para ajustar velocidades (temperatura) nuevamente
  int p_termalizacion, pasos, cadapterm;
  int bins=0; // para funcion de dist radial
  int cadad, it0, t0max, tmax; // esto es para calcular difusion y/o g(r)
	// variables
	double pev;
	double energia_k_inst,   energia_k_prom=0,  energia_k_suma=0,
				 energia_p_inst,   energia_p_prom=0,  energia_p_suma=0,
				 w_inst,           w_prom=0,          w_suma=0,
				 presion_inst,     presion_prom,
				 temperatura_inst, temperatura_prom;
  vector *posicion, *velocidad, *fuerza;
  vector *x; double dx_k;
	// output file
	FILE *rundatafile=stdout, *infodatafile=stdout, *posdatafile;
	char fnameprefix[FNAMELEN]="",
			 rundatafname[FNAMELEN]="",
			 posdatafname[FNAMELEN]="",
			 infodatafname[FNAMELEN]="";
	char prtlabels[40] = "";
	char *varlabel[] = {
	"ki", "vi", "ti", "pi",
	"kp", "vp", "tp", "pp"
	};


	// pa es un array de punteros a las distintas variables
	// (usado para ordenarlas con un indice)
  pa[KI]=&energia_k_inst;   pa[KP]=&energia_k_prom;
  pa[VI]=&energia_p_inst;   pa[VP]=&energia_p_prom;
  pa[TI]=&temperatura_inst; pa[TP]=&temperatura_prom;
  pa[PI]=&presion_inst;     pa[PP]=&presion_prom;

	// lee los argumentos de entrada y guarda parametros especificados
  argumentos(argc, argv, switches, &m,
  	&pasos, &p_termalizacion, &cadapterm,
  	&densidad, &temperatura, &dt, &rcorte,
		&bins,&it0,&t0max,&tmax,&cadad,
		fnameprefix);

  iniciar_ctes(&m, &n, &densidad, &lado_celda, &lado_caja);
  
	// hacemos espacio en memoria
  posicion=malloc(sizeof(vector)*(n+1));
	velocidad=malloc(sizeof(vector)*(n+1));
	fuerza=malloc(sizeof(vector)*(n+1));

	// contamos valores a imprimir (cv) y valores instantaneos a imprimir (cvi)
  for(i=0;i<NVARS;i++) {
    cv+=switches[i];
    if(i<KP) cvi+=switches[i];
  }

	// array de punteros a las variables seleccionadas
  pv=malloc(sizeof(double*)*cv);
  for(i=0,j=0;i<NVARS;i++) {
    if(switches[i]) {
      pv[j]=pa[i];
      j++;
    }
  }

	// arregla dependencias (por ej., aun si solo pedimos imprimir promedios,
	// es de todas formas necesario calcular valores instantaneos,
	// aunque no se impriman)
	// (switches son los valores escogidos en argumentos, sw con dependencias resueltas)
  dependencias(switches,calcsws); 

	// iniciamos...
	iniciar_posiciones(m, posicion, lado_celda);
  iniciar_velocidades(n, velocidad, temperatura);
  // (las variables para promedios -si hacen falta- ya inician con valor nulo)
	// iniciamos tambien fdr y difusion, si fueron escogidas
	if(switches[G]) fdistradial(n,posicion,lado_caja,bins,switches,fnameprefix);
	if(switches[D]) s_c_difusion=0;
	
	// si se pidio imprimir la evolucion de alguna variable...
		// si se pidio guardar a un archivo...
	if(switches[SAVE]) {
		// abrimos archivo 'info'
		strcpy(infodatafname,fnameprefix); strcat(infodatafname,"info.data");
		infodatafile=fopen(infodatafname,"w");
		if(infodatafile==NULL)
			serror("Error al abrir archivo p/ guardar info \
					(parametros y condiciones) de simulacion.");
			
		// abrimos archivo 'run' (si se pidio imprimir alguna variable)
		if(cv>0) {
			strcpy(rundatafname,fnameprefix); strcat(rundatafname,"run.data");
			rundatafile=fopen(rundatafname,"w");
			if(rundatafile==NULL)
				serror("Error al abrir archivo p/ guardar evolucion vs. tiempo.");
		} // de otra forma, rundatafile=stdout
	}

	// si se pide imprimir info, imprimimos
	if(switches[INFO]) {
		imprimir_info(switches, m,	pasos, p_termalizacion, cadapterm,
			densidad, temperatura, dt, rcorte, lado_caja, infodatafile);
		// muestro que se imprime de la simulacion (si se imprime)
		if(cv>0) {
			for (i=0; i<NVARS; i++) {
				if(switches[i]) strcat(prtlabels,", ");
				if(switches[i]) strcat(prtlabels,varlabel[i]);
			}
			fprintf(rundatafile,"# simulacion (p, t%s)\n",prtlabels);
		}
	}
	// si se pidio guardar posiciones (no se imprime en pantalla)
	if(switches[SPOS]) {
		strcpy(posdatafname,fnameprefix); strcat(posdatafname,"pos.data");
			posdatafile=fopen(posdatafname,"w");
			if(posdatafile==NULL)
				serror("Error al abrir archivo p/ guardar posiciones vs tiempo.");
		// si se pide imprimir info, imprimimos
		if(switches[INFO])
			fprintf(posdatafile,
					"# posiciones vs t (t x0[0] x0[1] x0[2] ... xn[0] xn[1] xn[2])\n");
	}

	// bucle principal
	/////////////////////////////////////////////////////////
	// p: pasos (totales), pev: pasos de evolucion (luego de la termalizacion)
	for(p=0, pev=0, s_termalizacion=1; p<pasos; p++, pev+=!s_termalizacion) {
		if(switches[SPOS] && ((s_termalizacion&&switches[TERM])||!s_termalizacion))
			// guardamos posiciones, si se pide
			// (switches[TERM] define si se guardan posiciones durante termalizacion)
			guardar_posiciones(n,switches[TERM]?p*dt:pev*dt,posicion,posdatafile);

		// inicializacion (cero)
		energia_k_inst=0; energia_p_inst=0; w_inst=0;
		for(i=0;i<n;i++) for(k=0;k<3;k++) fuerza[i][k]=0;
	
		// calcular fuerzas
		// (ademas de energia_p_inst y w (factor p/ calcular presion), si se pide)
		s_fuerza=1;
		for(i=0;i<n-1;i++) for(j=i+1;j<n;j++)
			sumar_fuerza(posicion[i] , posicion[j],
					fuerza[i], fuerza[j],
					&w_inst, &energia_p_inst,
					lado_caja,s_fuerza,
					calcsws[PI],calcsws[VI],rcorte);

		// calcula posicion en t+h
		for(i=0;i<n;i++) for(k=0;k<3;k++) {
			dx_k=dt*(velocidad[i][k]+0.5*dt*fuerza[i][k]);
			posicion[i][k]+=dx_k;
			if(s_c_difusion) {
				x[i][k]+=dx_k;
			}
		}

		// suma fuerza en t+dt, para calcular aceleracion promedio
		s_fuerza=0;
		for(i=0;i<n-1;i++) for(j=i+1;j<n;j++)
				sumar_fuerza(posicion[i],  posicion[j],
						fuerza[i], fuerza[j],
						NULL, NULL,
						lado_caja,s_fuerza,
						calcsws[PI],calcsws[VI],rcorte);

		// verifica cond de contorno y termina si dt es muy grande
		if(verificar_condcontorno(n,posicion,lado_caja))
			serror("Reducir dt (paso temporal)");

		energia_k_inst = calcular_energia_k(n, velocidad);
		temperatura_inst = 2*energia_k_inst/(3*n); 
		presion_inst = densidad*(temperatura_inst+w_inst/(3*n));
	
		for(i=0;i<n;i++) for(k=0;k<3;k++) velocidad[i][k]+=0.5*dt*fuerza[i][k];

		if(s_termalizacion && p%cadapterm==0)
			termalizar(n, velocidad, temperatura, temperatura_inst);

		if(!s_termalizacion) {
			if(calcsws[PRMS]) {
				if(calcsws[KP]) {
					energia_k_suma+=energia_k_inst;
					energia_k_prom=energia_k_suma/pev;
				}
				if(calcsws[VP]) {
					energia_p_suma+=energia_p_inst;
					energia_p_prom=energia_p_suma/pev;
				}
				if(calcsws[PP]) {
					w_suma+=w_inst;
					w_prom=w_suma/pev;
					presion_prom = densidad*(2*energia_k_prom+w_prom)/(3*n);
				}
				if(calcsws[TP]) {
					temperatura_prom = 2*energia_k_prom/(3*n);
				}
			}
			if(switches[G]) {
				fdistradial(n,posicion,lado_caja,bins,switches,fnameprefix);
			}
    }

		// si ya se llevaron a cabo todos los pasos de la termalizacion...
		if(p==p_termalizacion) {
			s_termalizacion=0; // switch de termalizacion a cero

			// y si se pidio calcular difusion, es momento de empezar
			if(switches[D]) {
				s_c_difusion=1;

				x=malloc(sizeof(vector)*(n+1));
				if(x==NULL)
					serror("No se ha podido reservar memoria para nuevas posiciones.");
				for(i=0;i<n;i++) for(k=0;k<3;k++) x[i][k]=posicion[i][k];

				difusion(n,x,it0,t0max,tmax,cadad,0,switches,fnameprefix);
			}
		}

		// se calcula difusion cada 'cadad' pasos (una vez termalizado el sistema)
		// (s_c_difusion=1 recien luego de haber concluido la termalizacion)
		if(switches[D] && !s_termalizacion && (((int)pev)%cadad==0)) {
			difusion(n,x,it0,t0max,tmax,cadad,0,switches,fnameprefix);
		}

    // imprimir valores instantaneos y/o promedios seleccionados
    if(cv>0) { // cv>0 si se pidio calcular al menos algun valor
			// si switches[TERM]==0 y estamos en termalizacion, no se imprime nada
			if(s_termalizacion&&!switches[TERM]) continue;

      fprintf(rundatafile,"%d %g", switches[TERM]?p:(int)pev, switches[TERM]?p*dt:pev*dt);
      for(j=0;j<(s_termalizacion?cvi:cv);j++)
        fprintf(rundatafile," %e", *pv[j]);
      fprintf(rundatafile,"\n");
    }
  }
	// fin de loop principal
	/////////////////////////////////////////////////////////
		
	// cerrar archivos, si se abrio alguno
	if(switches[SAVE]) {
		fclose(infodatafile);
		fclose(rundatafile);
		fclose(posdatafile);
	}
	// imprimimos los resultados de fdr y dist (si se pidio computar)
	// y liberamos memoria para terminar
  if(switches[G]) {
		fdistradial(n,NULL,densidad,0,switches,fnameprefix); // indica finalizacion del muestreo
		// (pasamos densidad en lugar de lado de caja)
	}
	if(switches[D]) {
		difusion(n,NULL,it0,t0max,tmax,cadad,cadad*dt,switches,fnameprefix); // indica finalizacion del muestreo
		free(x);
	}
  free(posicion);
	free(velocidad);
	free(fuerza);
  free(pv);

	return 0;
}
