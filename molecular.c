/* molecular
 * =========
 *
 * simula la dinamica de un conjunto de particulas
 * interactuando de a pares mediante el potencial de Lennard-Jones
 *
 *     V_{LJ} (r) = 4*epsilon*((sigma/r)^12-(sigma/r)^6)
 *
 *     (donde r es la distancia entre un par de particulas)
 *
 *  el cual modela de forma simplificada una interaccion
 *    repulsiva a cortas distancias
 *    y levemente atractiva distancias mayores
 *
 *  +--------------------------------------------------------------------+
 *  |  +      *       +      +      +       +      +      +       +      |
 *  |         *                                                          |
 *  |         *                                                          |
 *  |         *   V_{LJ}(r)/epsilon                                      |
 *  |         *                                                          |
 *  |         *                                                          |
 *  |         *                                                          |
 *  |         *             *********************************************|
 *  |         *         *****                                            |
 *  |         *      ***                                                 |
 *  |         *     **                                                   |
 *  |          *   *                                                     |
 *  |          *  *                                                      |
 *  |          * *                                                       |
 *  |  +      + **    +      +      +       +      +      +       +      |   r
 *  +--------------------------------------------------------------------+ -----
 *    0.5     1      1.5     2     2.5      3     3.5     4      4.5     5 sigma
 *   
 *
 *  trabajo final de la materia 'fisica computacional'
 *    famaf-unc, 2011
 *
 */

#include "simlib.h"

int main(int argc, char *argv[]) {
	int i, j, k, p;
  int swv[10], sw[12], cv=0, cvi=0;
  double **pv, *pa[8];
  int m=0, n; 
  double densidad, temperatura, lado_celda, lado_caja;
  double dt, rcorte;
  double dr, acep;

	// pasos: nro de pasos totales; p_termalizacion<pasos
	// cadapterm: nro de pasos que se saltean, durante pasos de termalizacion,
	//            para ajustar velocidades (temperatura) nuevamente
  int p_termalizacion, pasos, cadapterm;

  int bins; // para funcion de dist radial
  int cadad=10, it0=12, t0max=20, tmax=200; // esto es para calcular difusion y/o g(r)

	// switches
	int s_fuerza, s_termalizacion, s_c_difusion=0;

	double t, pev;
	double energia_k_inst,   energia_k_prom,  energia_k_suma,
				 energia_p_inst,   energia_p_prom,  energia_p_suma,
				 w_inst,           w_prom,          w_suma,
				 presion_inst,     presion_prom,    presion_suma,
				 temperatura_inst, temperatura_prom;
  vector *posicion, *velocidad, *fuerza;
	
  vector *x;
	double dx_k;
 
	// pa es un array de punteros a las distintas variables
	// (usado para ordenarlas con un indice)
  pa[KI]=&energia_k_inst;
  pa[VI]=&energia_p_inst;
  pa[TI]=&temperatura_inst;
  pa[PI]=&presion_inst;
  pa[KP]=&energia_k_prom;
  pa[VP]=&energia_p_prom;
  pa[TP]=&temperatura_prom;
  pa[PP]=&presion_prom;

	// lee los argumentos de entrada y guarda parametros especificados
  parametros(argc, argv, sw, &m,
  	&pasos, &p_termalizacion, &cadapterm,
  	&densidad, &temperatura, &dt, &rcorte, &dr);

  iniciar_ctes(&m, &n, &densidad, &lado_celda, &lado_caja);
  
	// hacemos espacio en memoria para los vectores posicion de cada particula
	// (asi como velocidad y fuerza, si se selecciona el modo MD)
  posicion=malloc(sizeof(vector)*(n+1));
  if(sw[MODO]==MD) {
    velocidad=malloc(sizeof(vector)*(n+1));
    fuerza=malloc(sizeof(vector)*(n+1));
  }

	// lee la linea de entrada para saber que cosas calcular
  interruptores(argc, argv, swv);

	// que es esto??? ni idea (PENDIENTE ENTENDER, describir)
	// (G solo utilizo porque es el ultimo indice posible de swv,
	//  que guarda info de que cosas se calcula (swv=switch values, i guess))
	// (cv=calcular valores? cvi=cv instantaneos?)
	// [creo que lo que se hace es contar cuantas cosas tento que calcular
	// ...para hacer lugar en memoria]
  for(i=0;i<G;i++) {
    cv+=swv[i];
    if(i<KP) cvi+=swv[i];
  }

	// array de punteros a las variables seleccionadas
  pv=malloc(sizeof(double*)*cv);
  for(i=0,j=0;i<G;i++) {
    if(swv[i]) {
      pv[j]=pa[i];
      j++;
    }
  }

	// arregla dependencias (por ej., aun si solo pedimos imprimir promedios,
	// es de todas formas necesario calcular valores instantaneos)
  dependencias(swv,sw);
  
	// d por curva de difusion, g por funcion de distribucion radial
	// (saco mientras, porque me generaba problemas)
	/*
  if(sw[G]) iniciar_ctes_g(argc, argv, &bins);
  if(sw[D]) {
		iniciar_ctes_d(argc, argv, &it0, &t0max, &tmax, &cadad);
		if((t0max*it0)<tmax) {
			serror("t0max*it0<tmax");
		}
		if(cadad*it0*t0max>pasos-p_termalizacion) {
			serror("cadad*it0*t0max>pasos-p_termalizacion");
		}
	}
	*/

	// iniciamos...
	iniciar_posiciones(m, posicion, lado_celda);
  if(sw[MODO]==MD) iniciar_velocidades(n, velocidad, temperatura);
  if(sw[MODO]==MC) energia_p_inst=calcular_energia_p(n, posicion, lado_caja);
  // iniciamos variables para promedios
	if(sw[PRMS]) {
		if(sw[PP]) {presion_suma=0; w_suma=0; w_prom=0;}
		if(sw[KP]) {energia_k_suma=0; energia_k_prom=0;}
		if(sw[VP]) {energia_p_suma=0; energia_p_prom=0;}
		if(sw[G]) {fdistradial(n,posicion,lado_caja,bins);}
    if(sw[D]) {s_c_difusion=0;}
	}

	// bucle principal
	/////////////////////////////////////////////////////
	// PENDIENTE: que es pev? -> son los pasos luego de termalizar
	for(p=1,pev=0,t=0,s_termalizacion=1;
		p<=pasos;
		p++,t+=dt,pev+=!s_termalizacion) {

		// inicializacion (cero)
    if(sw[MODO]==MD)
			iniciar_valores_inst(sw,&energia_k_inst,&energia_p_inst,&w_inst);
    else if(sw[MODO]==MC) {
      // PENDIENTE. QUE SE HACE ACA?
    }

    if(sw[MODO]==MC) {
      for(i=0;i<n;i++) {
        energia_p_inst+=pasitoMC(n,i,posicion,lado_caja,temperatura,&dr,&acep);
      }
    }

    if(sw[MODO]==MD) {
      for(i=0;i<n;i++) {
        for(k=0;k<3;k++) fuerza[i][k]=0;
      }
    
      // calcula fuerzas
      /* ademas de Ep y W en caso de requerir el calculo de
         Ep y P instantaneos respectivamente */
      s_fuerza=1;
      for(i=0;i<n;i++) {
        for(j=i+1;j<n;j++) {
          sumar_fuerza(posicion[i] , posicion[j],
                       fuerza[i], fuerza[j],
                       &w_inst, &energia_p_inst,
                       lado_caja,s_fuerza,
                       sw[PI],sw[VI]);
        
        }
      }

      // calcula posicion en t+h
      for(i=0;i<n;i++) {
        for(k=0;k<3;k++) {
					dx_k=dt*(velocidad[i][k]+0.5*dt*fuerza[i][k]);
          posicion[i][k]+=dx_k;
          if(s_c_difusion) {
            x[i][k]+=dx_k;
          }
        }
      }

      // suma fuerza en t+dt, para calcular aceleracion promedio
      s_fuerza=0;
      for(i=0;i<n;i++) {
        for(j=i+1;j<n;j++) {
          sumar_fuerza(posicion[i],  posicion[j],
                       fuerza[i], fuerza[j],
                       NULL, NULL,
                       lado_caja,s_fuerza,
                       sw[PI],sw[VI]);
        }
      }

      // verifica cond de contorno y termina si dt es muy grande
      if(verificar_condcontorno(n,posicion,lado_caja)) {
        serror(ERRMSG_DT);
        return 1;
      }

      energia_k_inst = calcular_energia_k(n, velocidad);
      temperatura_inst = 2*energia_k_inst/(3*n); 
      presion_inst = densidad*(temperatura_inst+w_inst/(3*n));
    
      for(i=0;i<n;i++)
				for(k=0;k<3;k++)
          velocidad[i][k]+=0.5*dt*fuerza[i][k];

      if(s_termalizacion) if(p%cadapterm==0)
				termalizar(n, velocidad, temperatura, temperatura_inst);
    }

    if(sw[PRMS]) {
      if(!s_termalizacion) {
        if(sw[MODO]==MD) {
          if(sw[KP]) {
            energia_k_suma+=energia_k_inst;
            energia_k_prom=energia_k_suma/pev;
          }
          if(sw[VP]) {
            energia_p_suma+=energia_p_inst;
            energia_p_prom=energia_p_suma/pev;
          }
          if(sw[PP]) {
            w_suma+=w_inst;
            w_prom=w_suma/pev;
            presion_prom = densidad*(2*energia_k_prom+w_prom)/(3*n);
          }
          if(sw[TP]) {
            temperatura_prom = 2*energia_k_prom/(3*n);
          }
        }
        if(sw[MODO]==MC) {
          if(sw[VP]) {
            energia_p_suma+=energia_p_inst;
            energia_p_prom=energia_p_suma/pev;
          }
          if(sw[PP]) {
            presion_inst=calcular_presion(n,posicion,lado_caja,temperatura,densidad);
            presion_suma+=presion_inst;
            presion_prom=presion_suma/pev;
          }
        }

        if(sw[G]) {
          fdistradial(n,posicion,lado_caja,bins);
        }
      }
    }

		// si ya se llevaron a cabo todos los pasos de la termalizacion...
		if(p==p_termalizacion) {
			// switch de termalizacion a cero
			s_termalizacion=0;

			// y si se pidio calcular difusion, es momento de empezar
			/*
			if(sw[D]) {
				s_c_difusion=1;

				x=malloc(sizeof(vector)*(n+1));
				if(x==NULL)
					serror("No se ha podido alocar memoria para nuevas posiciones.");

				for(i=0;i<n;i++) {
					for(k=0;k<3;k++) x[i][k]=posicion[i][k];
				}
				difusion(n,x,it0,t0max,tmax,0);
			}
			*/
		}

		// se calcula difusion cada 'cadad' pasos
		// (s_c_difusion=1 recien luego de haber concluido la termalizacion)
		if(sw[D] && !s_termalizacion && (((int)pev)%cadad==0)) {
			difusion(n,x,it0,t0max,tmax,0);
		}

    // imprimir valores instantaneos y/o promedios seleccionados
    if(cv>0) { // cv>0 si se pidio calcular al menos algun valor
      printf("%5d %6.7g  ", p, t);
      if(sw[MODO]==MC) printf("%3.4g ", acep);
      for(j=0;j<(s_termalizacion?cvi:cv);j++) {
        printf("%e ", *pv[j]);
      }
      printf("\n");
    }
  }
	/////////////////////////////////////////////////////
	
	/*
  if(sw[G]) {
		// indica finalizacion del muestreo
		fdistradial(n,NULL,densidad,0);
	}
	if(sw[D]) {
		difusion(n,NULL,it0,cadad,tmax,cadad*dt);
		free(x);
	}
	*/

  free(posicion);
  if(sw[MODO]==MD) {
    free(velocidad);
    free(fuerza);
  }
  free(pv);

	return 0;
}
