molecular
=========

![simulaciones](./simtilevids.gif)

se simula la dinámica de un conjunto de partículas
(en una caja 3D, con condiciones periódicas de contorno)

mediante dinámica molecular (integración numerica de ecuaciones diferenciales)


las partículas interactuan de a pares mediante el

```
potencial de Lennard-Jones

$ V_{LJ} (r) = 4*\epsilon*((\sigma/r)^12-(\sigma/r)^6) $

(donde r es la distancia entre el par de particulas)
```

el cual modela de forma simplificada una interacción
	repulsiva a cortas distancias
	y levemente atractiva distancias mayores

```
+--------------------------------------------------------------------+
|  +      *       +      +      +       +      +      +       +      |
|         *                                                          |
|         *                                                          |
|         *   V_{LJ}(r)/epsilon                                      |
|         *                                                          |
|         *                                                          |
|         *                                                          |
|         *             *********************************************|
|         *         *****                                            |
|         *      ***                                                 |
|         *     **                                                   |
|          *   *                                                     |
|          *  *                                                      |
|          * *                                                       |
|  +      + **    +      +      +       +      +      +       +      |   r
+--------------------------------------------------------------------+ -----
0.5     1      1.5     2     2.5      3     3.5     4      4.5     5   sigma
```

Trabajo final de la materia Física Computacional @FaMaF-UNC, 2011



(2021: parsing de argumentos mejorado, entre otros)


