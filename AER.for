C     Last change:  ZDF  27 Apr 2015    9:48 am
C====                            PROGRAMA AMAE00.FOR                               ====
C====          PROGRAMA DEDICADO AL ANALISIS MATRICIAL DE ESTRUCTURAS              ====
C====                                                                              ====
C====          PROGRAMA ESCRITO POR EL PROF. ZEFERINO DAFONSECA                    ====
C==== 									           ====
C====          ELEMENTOS INCORPORADOS :                                            ====
C====          UNIDIMENSIONAL AXIAL (AXI) ; INDEL = 1                              ====
C====          ARMADURA PLANA       (A2D) ; INDEL = 2                              ====
C====          VIGA DE EJE RECTO    (VIG) ; INDEL = 3                              ====
C====          PORTICO PLANO        (P2D) ; INDEL = 4                              ====
C====          PARRILLA             (PAR) ; INDEL = 5                              ====
C====          ARMADURA ESPACIAL    (A3D) ; INDEL = 6                              ====
C====          PORTICO ESPACIAL     (P3D) ; INDEL = 7                              ====
C==== 									           ====
C====          NOTA1:                                                              ====
C====          EN EL CASO QUE EL DESPLAZAMIENTO DE UN NODO, EN UNA DIRECCION       ====
C====          DADA, SEA INCOGNITA, SE LE ASIGNARA EL VALOR DE 200.0; ES DECIR,    ====
C====          DIS(JJ) = 200.0                                                     ====
C====                                                                              ====
C====          NOTA2:                                                              ====
C====          (a) EL VECTOR DE CARGAS NODALES EQUIVALENTES CORRESPONDIENTE A      ====
C====              LA ACCION DE LAS CARGAS SOBRE LOS ELEMENTOS DEBE SER DADO,      ====
C====              AL PROGRAMA, POR EL USUARIO.                                    ====
C====          (b) NO SE INCLUYE LA POSIBLIDAD DE CARGA EXTERNA SOBRE UN APOYO.    ====
C====          (c) NO SE CONSIDERAN APOYOS INCLINADOS CON RELACION A LOS EJES      ====
C====              GLOBALES DE LA ESTRUCTURA.                                      ====
C======================================================================================
C======================================================================================
C====            EN LA UNIDAD 30 SE GUARDAN LAS PROPIEDADES GEOMETRICAS            ====
C====            Y FISICAS DEL MATERIAL DE LOS ELEMENTOS.                          ====
C====            EN LA UNIDAD 14 SE GUARDAN LOS COEFICIENTES DE LA MATRIZ          ====
C====            [K]*[R].                                                          ====
C====            EN LA UNIDAD 15 SE GUARDAN LAS COMPONENTES DEL VECTOR             ====
C====            DE CARGAS DE LOS ELEMENTOS CON RESPECTO AL SISTEMA LOCAL          ====
C====            DE REFERENCIA.                                                    ====
C======================================================================================
C======================================================================================
C====  VARIABLES PRINCIPALES:                                                      ====
C====                                                                              ====
C====  INCID   : VECTOR DE LAS INCIDENCIAS DE LOS ELEMENTOS.                       ====
C====  IMAT    : VECTOR DE LAS PROPIEDADES GEOMETRICAS Y FISICAS                   ====
C====            DE LOS ELEMENTOS.                                                 ====
C====  GLOBAL  : MATRIZ QUE CONTIENE LOS COEFICIENTES DE RIGIDEZ DE LA             ====
C====            ESTRUCTURA.                                                       ====
C====  COOR    : VECTOR DE COORDENADAS DE LOS NODOS DE LA ESTRUCTURA.              ====
C====  NQ      : VECTOR CUYOS ELEMENTOS IDENTIFICAN LOS ELEMENTOS CARGADOS.        ====
C====  NAPE    : VECTOR CUYOS ELEMENTOS IDENTIFICAN LOS NODOS CON APOYOS           ====
C====            ELASTICOS.                                                        ====
C====  DISA    : VECTOR CON LOS VALORES DE DE LA CONSTANTE DEL APOYO ELASTICO      ====
C====            EN LAS DIRECCIONES CORRESPONDIENTES.                              ====
C====  NAP     : VECTOR CUYOS ELEMENTOS IDENTIFICAN LOS NODOS CON DESPLAZAMIENTO   ====
C====            CONOCIDOS.                                                        ====
C====  DIS     : VECTOR CON LOS VALORES DE LOS DESPLAZAMIENTOS CONOCIDOS EN LAS    ====
C====            DIRESCCIONES CORRESPONDIENTES.                                    ====
C====  DES     : VECTOR CON LOS DESPLAZAMIENTOS DE LA ESTRUCTURA.                  ====
C====  CARGA   : VECTOR DE CARGAS DE LA ESTRUCTURA CON RELACION AL SISTEMA         ====
C===             GLOBAL DE REFERENCIA.                                             ====
C======================================================================================
      PROGRAM AER
      USE M1
      USE M2
      USE M3
      INCLUDE'AME'
c.....
          CALL DATG0
          CALL WDATG1
c.....   Lectura y Escritura de los datos generales de la estructura
	  CALL DATGE
c.....   Lectura y Escritura del Sistema de Cargas Actuante
	  IF(NNC.NE.0) CALL CARGN
          IF(NEC.NE.0) CALL CARGE
          CALL DATGE1
c.....   Montaje de la Matriz Global
	  CALL MONTA
c.....   Introduccion de las condiciones de contorno en la matriz global
c.....   y en el vector de cargas nodales equivalente
	  CALL CONTOR
c.....   Decomposicion de la matriz global segun el metodo de Cholesky
	  CALL FACTORA
c.....    Resolucion del sistema de ecuaciones resultante
	  CALL SOLUCION
c.....    Escritura del vector de desplazamientos
	  CALL ESCRIBE
c.....    Calculo y escritura de las acciones de extremo de miembro (siste-
c.....    ma local de referencia) y calculo de las reacciones
	  CALL ACCION
C....    Escritura de las reacciones
	  CALL ESCREA
      END
C=====================
      SUBROUTINE DATGE
      INCLUDE'AME'
c.....
      WRITE(OUP,19)
      READ(INP,10)TITULO
      WRITE(OUP,20)TITULO
      WRITE(OUP,19)
      WRITE (OUP,30)
      READ(INP,*)NN,NE,NGLN,NNPE,NCOPN,NNAE,NNDP,NNC,NEC,INDEL,NMAT
      WRITE(OUP,5766)NN,NE,NGLN,NNPE,NCOPN,NNAE,NNDP,NNC,NEC,INDEL,NMAT
c....
      NTEC = NN*NGLN                       !  Numero total de ecuaciones
      NGLL = NNPE*NGLN                     !  Numero de ecuaciones por elemento
       IF((INDEL.EQ.1).OR.(INDEL.EQ.3)) WRITE(OUP,601)
       IF((INDEL.EQ.2).OR.(INDEL.EQ.4)) WRITE(OUP,602)
       IF(INDEL.EQ.5) WRITE(OUP,602)
       IF((INDEL.EQ.6).OR.(INDEL.EQ.7)) WRITE(OUP,603)
       DO 171 JJ=1,NN
       LL=JJ*NCOPN
       L=LL-NCOPN+1
       READ(INP,*)N,(COOR(I),I=L,LL)       ! Lectura de las coordenadas de los nodos
       IF(JJ-N)170,171,170
 170  WRITE(OUP,210)JJ                                                   
      STOP
 171  CONTINUE
       DO 600 J=1,NN
       LL=J*NCOPN
       L=LL-NCOPN+1
       WRITE(OUP,110)J,(COOR(I),I=L,LL)
 600  CONTINUE
      WRITE(OUP,120)
       DO 172 K=1,NE
       LL=K*NNPE
       L=LL-NNPE+1
       READ(INP,*)N,(INCID(JJ),JJ=L,LL),IMAT(K)  ! Lectura de las incidencias de los elementos y
       IF(K-N)173,172,173                        ! tipo de material
 173  WRITE(OUP,200)K
      STOP
 172  CONTINUE                                                          
       DO 17 K=1,NE
       LL=K*NNPE
       L=LL-NNPE+1
       WRITE(OUP,140)K,(INCID(I),I=L,LL),IMAT(K)
 17   CONTINUE
c======                                       ! Consideracion de apoyos elasticos
       IF(NNAE.EQ.0) GO TO 402
       IF(INDEL.EQ.1) WRITE(OUP,249)
       IF((INDEL.EQ.2).OR.(INDEL.EQ.3)) WRITE(OUP,250)
       IF((INDEL.EQ.4).OR.(INDEL.EQ.5)) WRITE(OUP,251)
       IF(INDEL.EQ.6)WRITE(OUP,251)
       IF(INDEL.EQ.7)WRITE(OUP,252)
       DO 401 J = 1,NNAE
       LL = J*NGLN
       L = LL - NGLN + 1
       READ(INP,*)NAPE(J),(DISA(K),K=L,LL)    ! Lectura de los nodos con apoyos elasticos y sus valores
       WRITE(OUP,180)NAPE(J),(DISA(K),K=L,LL)
 401  CONTINUE
 402  CONTINUE                                ! Consideracion de los desplazamientos prescritos
c====
       IF(INDEL.EQ.1) WRITE(OUP,149)
       IF((INDEL.EQ.2).OR.(INDEL.EQ.3)) WRITE(OUP,150)
       IF((INDEL.EQ.4).OR.(INDEL.EQ.5)) WRITE(OUP,151)
       IF(INDEL.EQ.6)WRITE(OUP,151)
       IF(INDEL.EQ.7)WRITE(OUP,152)
       DO 41 J = 1,NNDP
       LL = J*NGLN
       L = LL-NGLN+1
       READ(INP,*)NAP(J),(DIS(K),K=L,LL)      ! Lectura de los nodos prescritos y los valores conocidos
       WRITE(OUP,180)NAP(J),(DIS(K),K=L,LL)
 41   CONTINUE
c.... Lectura de las propiedades geometricas y fisicas de los elementos
      DO 11 K=1,NMAT
       IF((INDEL.EQ.1).OR.(INDEL.EQ.2)) GO TO 700
       IF(INDEL.EQ.3) GO TO 710
       IF(INDEL.EQ.4) GO TO 720
       IF(INDEL.EQ.5) GO TO 730
       IF(INDEL.EQ.6) GO TO 700
       IF(INDEL.EQ.7) GO TO 740
c....  Propiedades de los elementos "AXIAL" , "A2D" y "A3D"
 700  READ(INP,*)E,AX
       WRITE(OUP,404)K,E,AX
       IH30 = K
       WRITE(30,REC=IH30)E,AX
       GO TO 11
c....  Propiedades del elemento "VIG"
 710  READ(INP,*)E,CIZ
       WRITE(OUP,405)K,E,CIZ
       IH30 = K
       WRITE(30,REC=IH30)E,CIZ
       GO TO 11
c....  Propiedades del elemento "P2D"
 720  READ(INP,*)E,AX,CIZ
       WRITE(OUP,407)K,E,AX,CIZ
       IH30 = K
      WRITE(30,REC=IH30)E,AX,CIZ
       GO TO 11
c....  Propiedades del elemento "PAR"
 730  READ(INP,*)E,G,CIX,CIY
       WRITE(OUP,406)K,E,G,CIX,CIY
       IH30 = K
       WRITE(30,REC=IH30)E,G,CIX,CIY
       GO TO 11
c....  Propiedades del elemento "P3D"
 740  READ(INP,*)E,G,AX,CIX,CIY,CIZ
       WRITE(OUP,408)K,E,G,AX,CIX,CIY,CIZ
       IH30 = K
       WRITE(30,REC=IH30)E,G,AX,CIX,CIY,CIZ
  11  CONTINUE
c==== Formatos de escritura
 10   FORMAT(20A4)
 19   FORMAT(1X,115('-'))
 20   FORMAT(1X,'-',2X,'PROBLEMA : ',20A4)
 30   FORMAT(/1X,'DATOS GENERALES DE LA ESTRUCTURA',
     1       /1X,'================================')
 5766 FORMAT(
     1/1X,'NUMERO DE NODOS............................(NN).............=
     1',I5,
     1/1X,'NUMERO DE ELEMENTOS........................(NE).............=
     2',I5,
     1/1X,'NUMERO DE GRADOS DE LIBERTAD/NODO..........(NGLN)...........=
     3',I5,
     1/1X,'NUMERO DE NODOS/ELEMENTO...................(NNPE)...........=
     4',I5,
     1/1X,'NUMERO DE COORDENADAS/NODO.................(NCOPN)..........=
     5',I5,
     1/1X,'NUMERO DE NODOS CON APOYOS ELASTICOS.......(NNAE)...........=
     6',I5,
     1/1X,'NUMERO DE NODOS PRESCRITOS.................(NNDP)...........=
     6',I5,
     1/1X,'NUMERO DE NODOS CARGADOS...................(NNC)............=
     7',I5,
     1/1X,'NUMERO DE ELEMENTOS CARGADOS...............(NEC)............=
     8',I5,
     1/1X,'ELEMENTO.....(AXI/A2D/VIG/P2D/PAR/A3D/P3D:1/2/3/4/5/6/7)....=
     9',I5,
     1/1X,'NUMERO DE MATERIALES.......................(NMAT)...........=
     9',I5)
c==============================================
 601  FORMAT(/1X,'COORDENADAS DE LOS NODOS',/1X,
     1            '========================',/,1X,'NODO',10X,'COOR. X')
 602  FORMAT(/1X,'COORDENADAS DE LOS NODOS',/1X,
     1           '========================',/,1X,'NODO',10X,'COOR. X',1
     21X,'COOR. Y')
 603  FORMAT(/1X,'COORDENADAS DE LOS NODOS',/1X,
     1            '========================',/,1X,'NODO',10X,'COOR. X',1
     21X,'COOR. Y'11X,'COOR. Z')
 110  FORMAT(I3,3X,3E18.5)
 210  FORMAT(//5X,'====  ERROR DE NUMERACION EN EL NODO ',I5,' ====')
c==================================================
 120  FORMAT(/1X,'INCIDENCIAS DE LOS ELEMENTOS',/1X,
     1           '============================',/1X,'ELEMENTO',5X,'NODOS
     2 CONCURRENTES', 4X,'TIPO DE MATERIAL')
 140  FORMAT(I5,5X,2I9,12X,I3)
 200  FORMAT(//5X,'===== ERROR DE NUMERACION EN EL ELEMENTO ',I5,'===='
     1)
c======================================
 249  FORMAT(/1X,'APOYOS ELASTICOS',/1X,
     1            '===============',/1X,'NODO',4X,'DIR.1')
 250  FORMAT(/1X,'APOYOS ELASTICOS',/1X,
     1          '================',/1X,'NODO',4X,'DIR.1',9X,'DIR.2')
 251  FORMAT(/1X,'APOYOS ELASTICOS',/1X,
     1          '================',/1X,'NODO',4X,'DIR.1',9X,'DIR.2',9X,
     1'DIR.3')
 252  FORMAT(/1X,'APOYOS ELASTICOS',/1X,
     1          '================',/1X,'NODO',4X,'DIR.1',9X,'DIR2.',9X,
     1'DIR.3',9X,'DIR.4',9X,'DIR.5',9X,'DIR.6')
c=============================================
 149  FORMAT(/1X,'CONDICIONES DE CONTORNO',/1X,
     1            '======================',/1X,'NODO',4X,'DIR.1')
 150  FORMAT(/1X,'CONDICIONES DE CONTORNO',/1X,
     1          '=======================',/1X,'NODO',4X,'DIR.1',9X,'DIR
     2.2')
 151  FORMAT(/1X,'CONDICIONES DE CONTORNO',/1X,
     1          '=======================',/1X,'NODO',4X,'DIR.1',9X,'DIR
     2.2',9X,'DIR.3')
 152  FORMAT(/1X,'CONDICIONES DE CONTORNO',/1X,
     1          '=======================',/1X,'NODO',4X,'DIR.1',9X,'DIR
     2.2',9X,'DIR.3',9X,'DIR.4',9X,'DIR.5',9X,'DIR.6')
 180  FORMAT(I3,6E14.5)
c================================================================================
c==== Propiedades geometricas y fisicas del elemento: Axial                 (AXI)
c==== Propiedades geometricas y fisicas del elemento: Armadura Plana        (A2D)
c==== Propiedades geometricas y fisicas del elemento: Armadura Plana        (A3D)
 404  FORMAT(/1X,'PROPIEDADES DEL MATERIAL NUMER0  =',I3,/1X,
     1           '===============================',
     2/1X,'MODULO DE ELASTICIDAD................(E)....=',E15.5
     3/1X,'AREA DE LA SECCION TRANSVERSAL.......(AX)...=',E15.5)           ! OK
c==== Propiedades geometricas y fisicas del elemento: Viga de Eje Recto     (VIG)
 405  FORMAT(/1X,'PROPIEDADES DEL MATERIAL NUMER0  =',I3,/1X,
     1           '===============================',
     2/1X,'MODULO DE ELASTICIDAD................(E).....=',E15.5,
     3/1X,'MOMENTO DE INERCIA DE LA SECCION....(CIZ)....=',E15.5)
c==== Propiedades geometricas y fisicas del elemento: Parrilla              (PAR)
 406  FORMAT(/1X,'PROPIEDADES DEL MATERIAL NUMER0  =',I3,/1X,
     1           '===============================',
     2/1X,'MODULO DE ELASTICIDAD...............(E)......=',E15.5,
     3/1X,'MODULO DE RIGIDEZ...................(G)......=',E15.5,
     4/1X,'MOMENTO DE INERCIA DE LA SECCION....(CIX)....=',E15.5,
     5/1X,'MOMENTO DE INERCIA DE LA SECCION....(CIY)....=',E15.5)
c==== Propiedades geometricas y fisicas del elemento: Portico Plano         (P2D)
 407  FORMAT(/1X,'PROPIEDADES DEL MATERIAL NUMER0  =',I3,/1X,
     1           '===============================',
     2/1X,'MODULO DE ELASTICIDAD...............(E)......=',E15.5,
     3/1X,'AREA DE LA SECCION TRANSVERSAL......(AX).....=',E15.5,
     4/1X,'MOMENTO DE INERCIA DE LA SECCION....(CIZ)....=',E15.5)
c==== Propiedades geometricas y fisicas del elemento: Portico Espacial      (P3D)
 408  FORMAT(/1X,'PROPIEDADES DEL MATERIAL NUMER0  =',I3,/1X,
     1           '===============================',
     2/1X,'MODULO DE ELASTICIDAD...............(E)......=',E15.5,
     3/1X,'MODULO DE RIGIDEZ...................(G)......=',E15.5,
     4/1X,'AREA DE LA SECCION TRANSVERSAL......(AX).....=',E15.5,
     5/1X,'MOMENTO DE INERCIA DE LA SECCION....(CIX)....=',E15.5,
     6/1X,'MOMENTO DE INERCIA DE LA SECCION....(CIY)....=',E15.5,
     7/1X,'MOMENTO DE INERCIA DE LA SECCION....(CIZ)....=',E15.5)
      PAUSE
      RETURN
      END
C=====================
      SUBROUTINE CARGN
      INCLUDE'AME'
      DIMENSION CAUX(12)
c.....
       IF(INDEL.EQ.1) WRITE(OUP,101)
       IF((INDEL.EQ.2).OR.(INDEL.EQ.3)) WRITE(OUP,102)
       IF((INDEL.EQ.4).OR.(INDEL.EQ.5)) WRITE(OUP,103)
       IF(INDEL.EQ.6) WRITE(OUP,103)
       IF(INDEL.EQ.7) WRITE(OUP,104)
       DO 205 J=1,NTEC
 205  CARGA(J) = 0.0
       DO 206 J = 1,NNC
       READ(INP,*) NP,(CAUX(K),K=1,NGLN)
       WRITE(OUP,215)NP,(CAUX(K),K=1,NGLN)
       DO 220 K=1,NGLN
       I=(NP-1)*NGLN + K
 220  CARGA(I) = CARGA(I) + CAUX(K)
 206  CONTINUE
      IF(INDEL.EQ.7) GO TO 315
       WRITE(OUP,20)
       WRITE(OUP,22)(CARGA(J),J=1,NTEC)
 315  CONTINUE
c==== Formatos de escritura
 101  FORMAT(/1X,'CARGAS EN LOS NODOS',/1X,
     1          '===================',/1X,'NODO',4X,'DIR.1')
 102  FORMAT(/1X,'CARGAS EN LOS NODOS',/1X,
     1          '===================',/1X,'NODO',4X,'DIR.1',9X,'DIR.2')
 103  FORMAT(/1X,'CARGAS EN LOS NODOS',/1X,
     1          '===================',/1X,'NODO',4X,'DIR.1',9X,'DIR.2',
     29X,'DIR.3')
 104  FORMAT(/1X,'CARGAS EN LOS NODOS',/1X,
     1          '===================',/1X,'NODO',4X,'DIR.1',9X,'DIR.2',
     39X,'DIR.3',9X,'DIR.4',9X,'DIR.5',9X,'DIR.6')
 215  FORMAT(I3,6E14.5)
  20  FORMAT(/1X,'VECTOR DE CARGAS NODALES',/1X,
     1           '========================')
  22  FORMAT(E15.5)
      RETURN
      END
C=====================
      SUBROUTINE CARGE
      INCLUDE 'AME'
      DIMENSION CAUX(12),JP(2),ROT(12,12),AUX(12)
c......
       IF(INDEL.EQ.3) WRITE(OUP,100)
       IF((INDEL.EQ.4).OR.(INDEL.EQ.5)) WRITE(OUP,101)
       IF(INDEL.EQ.6) WRITE(OUP,101)
       IF(INDEL.EQ.7) WRITE(OUP,102)
       IF(NNC) 2,3,2
  3    DO 4 J=1,NTEC
  4   CARGA(J) = 0.0
  2   CONTINUE
       DO 5 IK = 1,NEC
       READ(INP,*) NQ(IK),(CAUX(K),K=1,NGLL)
      IF(INDEL.EQ.3)  WRITE(OUP,7)NQ(IK),(CAUX(K),K=1,NGLL)
      IF((INDEL.EQ.4).OR.(INDEL.EQ.5)) WRITE(OUP,8)NQ(IK),(CAUX(K),
     1K=1,NGLL)
      IF(INDEL.EQ.6)  WRITE(OUP,8)NQ(IK),(CAUX(K),K=1,NGLL)
      IF(INDEL.EQ.7)  WRITE(OUP,9)NQ(IK),(CAUX(K),K=1,NGLL)
       IH = NQ(IK)
       IH15 = IK
       WRITE(15,REC=IH15) (CAUX(I),I=1,NGLL)
       CALL ROTAC(IH,XL,ROT)
       DO 80 I=1,NGLL
       AUX(I) = 0.0
       DO 80 J=1,NGLL
 80   AUX(I) = AUX(I) + ROT(J,I)*CAUX(J)
       DO 90 I=1,NNPE
       I1 =(IH-1)*NNPE + I
 90   JP(I) = INCID(I1)
      IF(INDEL.EQ.3) GO TO 50
      IF((INDEL.EQ.4).OR.(INDEL.EQ.5)) GO TO 60
      IF(INDEL.EQ.6) GO TO 60
      IF(INDEL.EQ.7) GO TO 70
  50  CONTINUE
       CARGA(2*JP(1)-1) = CARGA(2*JP(1)-1)  - CAUX(1)
       CARGA(2*JP(1))   = CARGA(2*JP(1))    - CAUX(2)
       CARGA(2*JP(2)-1) = CARGA(2*JP(2)-1)  - CAUX(3)
       CARGA(2*JP(2))   = CARGA(2*JP(2))    - CAUX(4)
        GO TO 5
  60  CONTINUE
       CARGA(3*JP(1)-2) = CARGA(3*JP(1)-2)  - AUX(1)
       CARGA(3*JP(1)-1) = CARGA(3*JP(1)-1)  - AUX(2)
       CARGA(3*JP(1))   = CARGA(3*JP(1))    - AUX(3)
       CARGA(3*JP(2)-2) = CARGA(3*JP(2)-2)  - AUX(4)
       CARGA(3*JP(2)-1) = CARGA(3*JP(2)-1)  - AUX(5)
       CARGA(3*JP(2))   = CARGA(3*JP(2))    - AUX(6)
        GO TO 5
  70  CONTINUE
       CARGA(6*JP(1)-5) = CARGA(6*JP(1)-5) - AUX(1)
       CARGA(6*JP(1)-4) = CARGA(6*JP(1)-4) - AUX(2)
       CARGA(6*JP(1)-3) = CARGA(6*JP(1)-3) - AUX(3)
       CARGA(6*JP(1)-2) = CARGA(6*JP(1)-2) - AUX(4)
       CARGA(6*JP(1)-1) = CARGA(6*JP(1)-1) - AUX(5)
       CARGA(6*JP(1))   = CARGA(6*JP(1))   - AUX(6)
       CARGA(6*JP(2)-5) = CARGA(6*JP(2)-5) - AUX(7)
       CARGA(6*JP(2)-4) = CARGA(6*JP(2)-4) - AUX(8)
       CARGA(6*JP(2)-3) = CARGA(6*JP(2)-3) - AUX(9)
       CARGA(6*JP(2)-2) = CARGA(6*JP(2)-2) - AUX(10)
       CARGA(6*JP(2)-1) = CARGA(6*JP(2)-1) - AUX(11)
       CARGA(6*JP(2))   = CARGA(6*JP(2))   - AUX(12)
  5   CONTINUE

c.... Formatos de escritura
 100  FORMAT(/1X,'CARGAS EN LOS ELEMENTOS',/1X,
     1           '=======================',/1X,
     1'ELEM.',5X,'CAR.1',9X,'CAR.2'/11X'CAR.3',9X,'CAR.4')
 101  FORMAT(/1X,'CARGAS EN LOS ELEMENTOS',/1X,
     1           '=======================',/1X,
     1'ELEM.',5X,'CAR.1',9X,'CAR.2',9X,'CAR.3'/11X,'CAR.4',9X,'CAR.5',9X
     2,'CAR.6')
 102  FORMAT(/1X,'CARGAS EN LOS ELEMENTOS',/1X,
     1           '=======================',/1X,
     1'ELEM.',5X,'CAR.1',9X,'CAR.2',9X,'CAR.3',9X,'CAR.4',9X,'CAR.5',9X,
     2'CAR.6'/11X,'CAR.7',9X,'CAR.8',9X,'CAR.9',9X,'CAR.10',8X,'CAR.11',
     38X,'CAR.12')
 22   FORMAT(E12.3)
  7   FORMAT(I3,2X,2E14.5,/5X,2E14.5,/)
  8   FORMAT(I3,2X,3E14.5,/5X,3E14.5,/)
  9   FORMAT(I3,2X,6E14.5,/5X,6E14.5,/)
      RETURN
      END
C==============================
      SUBROUTINE AXI(IH,RLOCAL)
      INCLUDE'AME'
c.....
      DIMENSION RLOCAL(12,12),JP(2)
       IH30 = IMAT(IH)
       READ(30,REC=IH30)E,AX
       DO 1 I=1,2
       DO 1 J=1,2
1     RLOCAL(I,J) = 0.0
       DO 2 J=1,2
       I = (IH-1)*NNPE + J
2     JP(J) = INCID(I)
       I = (JP(1) - 1)*NCOPN
       J = (JP(2) - 1)*NCOPN
       XCL = COOR(J+1) - COOR(I+1)
       XL = DSQRT(XCL**2)
       S1 = E*AX/XL
c.... Matriz de rigidez del elemento referida al sistema local/global
       RLOCAL(1,1) = S1
       RLOCAL(1,2) = -S1
       RLOCAL(2,1) = -S1
       RLOCAL(2,2) =  S1
       IH14 = IH
       WRITE(14,REC=IH14) ((RLOCAL(I,J),I=1,2),J=1,2)
      RETURN
      END
C==============================
      SUBROUTINE A2D(IH,RLOCAL)
      INCLUDE'AME'
      DIMENSION RLOCAL(12,12),REAC(12,12),ROT(12,12)
c.....
       IH30 = IMAT(IH)
       READ(30,REC=IH30)E,AX
       DO 1 I=1,4
       DO 1 J=1,4
       REAC(I,J) = 0.0
  1   RLOCAL(I,J) = 0.0
       CALL ROTAC(IH,XL,ROT)
c.... Matriz de rigidez del elemento referida al sistema local
       S1 = E*AX/XL
       RLOCAL(1,1) =  S1
       RLOCAL(1,3) = -S1
       RLOCAL(3,1) = RLOCAL(1,3)
       RLOCAL(3,3) =  S1
c.... Matriz de rigidez del elemento referida al sistema global
       DO 4 I=1,4
       DO 4 J=1,4
       REAC(I,J) = 0.0
       DO 4 K=1,4
  4   REAC(I,J) = REAC(I,J) + RLOCAL(I,K)*ROT(K,J)
       DO 5 I=1,4
       DO 5 J=1,4
  5   RLOCAL(I,J) = REAC(I,J)
       IH14 = IH
       WRITE(14,REC=IH14) ((RLOCAL(I,J),I=1,4),J=1,4)
       DO 6 I=1,4
       DO 6 J=1,4
  6   REAC(I,J) = 0.0
       DO 7 I=1,4
       DO 7 J=1,4
       REAC(I,J) = 0.0
       DO 7 K=1,4
  7   REAC(I,J) = REAC(I,J) + ROT(K,I)*RLOCAL(K,J)
       DO 8 I=1,4
       DO 8 J=1,4
  8   RLOCAL(I,J) = REAC(I,J)
      RETURN
      END
C==============================
      SUBROUTINE VIG(IH,RLOCAL)
      INCLUDE'AME'
      DIMENSION RLOCAL(12,12),JP(2)
c.....
       IH30 = IMAT(IH)
       READ(30,REC=IH30)E,CIZ
       DO 1 I=1,4
       DO 1 J=1,4
  1   RLOCAL(I,J) = 0.0
       DO 39 J=1,2
       I = (IH-1)*NNPE + J
 39   JP(J) = INCID(I)
       I = (JP(1) - 1)*NCOPN
       J = (JP(2) - 1)*NCOPN
       XCL = COOR(J+1) - COOR(I+1)
       XL = DSQRT(XCL**2.0)
       S1 = 12*E*CIZ/(XL**3.0)
       S2 = 6.0*E*CIZ/(XL**2.0)
       S3 = 4.0*E*CIZ/XL
c.... Matriz de rigidez del elemento referida al sistema local/global
       RLOCAL(1,1) =  S1
       RLOCAL(1,2) =  S2
       RLOCAL(1,3) = -S1
       RLOCAL(1,4) =  S2
       RLOCAL(2,2) =  S3
       RLOCAL(2,3) = -S2
       RLOCAL(2,4) =  S3/2.0
       RLOCAL(3,3) =  S1
       RLOCAL(3,4) = -S2
       RLOCAL(4,4) =  S3
       DO 3 I=1,4
       DO 3 J=1,4
  3   RLOCAL(J,I) = RLOCAL(I,J)
       IH14 = IH
       WRITE(14,REC=IH14) ((RLOCAL(I,J),I=1,4),J=1,4)
      RETURN
      END
C==============================
      SUBROUTINE P2D(IH,RLOCAL)
      INCLUDE'AME'
      DIMENSION RLOCAL(12,12),REAC(12,12),ROT(12,12)
c....
       IH30 = IMAT(IH)
       READ(30,REC=IH30)E,AX,CIZ
       DO 1 I = 1,6
       DO 1 J = 1,6
       ROT(I,J) = 0.0
  1   RLOCAL(I,J) = 0.0
       CALL ROTAC(IH,XL,ROT)
c....  Matriz de rigidez referida al sistema local
       S1 = E*AX/XL
       S2 = E*CIZ/XL**3.0
       RLOCAL(1,1) =  S1
       RLOCAL(2,2) =  12.*S2
       RLOCAL(3,2) =  6.0*S2*XL
       RLOCAL(3,3) =  4.0*S2*(XL**2.0)
       RLOCAL(4,1) = -S1
       RLOCAL(4,4) =  RLOCAL(1,1)
       RLOCAL(5,2) = -RLOCAL(2,2)
       RLOCAL(5,3) = -RLOCAL(3,2)
       RLOCAL(5,5) =  RLOCAL(2,2)
       RLOCAL(6,2) =  RLOCAL(3,2)
       RLOCAL(6,3) =  RLOCAL(3,3)/2.0
       RLOCAL(6,5) =  RLOCAL(5,3)
       RLOCAL(6,6) =  RLOCAL(3,3)
       DO 3 J1 = 1,6
       DO 3 J2 = J1,6
  3   RLOCAL(J1,J2) = RLOCAL(J2,J1)
       CALL ROTAC(IH,XL,ROT)
c.... Matriz de rigidez referida al sistema global
       DO 4 I=1,6
       DO 4 J=1,6
       REAC(I,J) = 0.0
       DO 4 K=1,6
   4  REAC(I,J) = REAC(I,J) + RLOCAL(I,K)*ROT(K,J)
       DO 5 I=1,6
       DO 5 J=1,6
   5  RLOCAL(I,J) = REAC(I,J)
       IH14 = IH
       WRITE(14,REC=IH14) ((RLOCAL(I,J),I=1,6),J=1,6)
       DO 6 I=1,6
       DO 6 J=1,6
   6  REAC(I,J) = 0.0
       DO 7 I=1,6
       DO 7 J=1,6
       REAC(I,J) = 0.0
       DO 7 K=1,6
   7  REAC(I,J) = REAC(I,J) + ROT(K,I)*RLOCAL(K,J)
       DO 8 I=1,6
       DO 8 J=1,6
  8   RLOCAL(I,J) = REAC(I,J)
      RETURN
      END
C==============================
      SUBROUTINE PAR(IH,RLOCAL)
      INCLUDE'AME'
      DIMENSION RLOCAL(12,12),REAC(12,12),ROT(12,12)
c....
      IH30 = IMAT(IH)
      READ(30,REC=IH30)E,G,CIX,CIY
       DO 1 I = 1,6
       DO 1 J = 1,6
       ROT(I,J) = 0.0
 1     RLOCAL(I,J) = 0.0
       CALL ROTAC(IH,XL,ROT)
       C2 = XL**2.0
       C3 = XL**3.0
c....  Matriz de rigidez referida al sistema local
       RLOCAL(1,1)  =   G*CIX/XL
       RLOCAL(2,2)  =   4.0*E*CIY/XL
       RLOCAL(3,2) =   -6.0*E*CIY/C2
       RLOCAL(3,3) =   12.0*E*CIY/C3
       RLOCAL(4,1) =  -RLOCAL(1,1)
       RLOCAL(4,4) =   RLOCAL(1,1)
       RLOCAL(5,2) =   2.0*E*CIY/XL
       RLOCAL(5,3) =   RLOCAL(3,2)
       RLOCAL(5,5) =   RLOCAL(2,2)
       RLOCAL(6,2) =  -RLOCAL(3,2)
       RLOCAL(6,3) =  -RLOCAL(3,3)
       RLOCAL(6,5) =  -RLOCAL(3,2)
       RLOCAL(6,6) =   RLOCAL(3,3)
       DO 3 J1 = 1,6
       DO 3 J2 = J1,6
  3   RLOCAL(J1,J2) = RLOCAL(J2,J1)
       CALL ROTAC(IH,XL,ROT)
c.... Matriz de rigidez referida al sistema global
       DO 4 I=1,6
       DO 4 J=1,6
       REAC(I,J) = 0.0
       DO 4 K=1,6
  4   REAC(I,J) = REAC(I,J) + RLOCAL(I,K)*ROT(K,J)
       DO 5 I=1,6
       DO 5 J=1,6
  5   RLOCAL(I,J) = REAC(I,J)
       IH14 = IH
       WRITE(14,REC=IH14) ((RLOCAL(I,J),I=1,6),J=1,6)
       DO 6 I=1,6
       DO 6 J=1,6
  6   REAC(I,J) = 0.0
       DO 7 I=1,6
       DO 7 J=1,6
       REAC(I,J) = 0.0
       DO 7 K=1,6
  7   REAC(I,J) = REAC(I,J) + ROT(K,I)*RLOCAL(K,J)
       DO 8 I=1,6
       DO 8 J=1,6
  8   RLOCAL(I,J) = REAC(I,J)
      RETURN
      END
C==============================
      SUBROUTINE A3D(IH,RLOCAL)
      INCLUDE'AME'
      DIMENSION RLOCAL(12,12),REAC(12,12),ROT(12,12)
c....
      IH30 = IMAT(IH)
      READ(30,REC=IH30)E,AX
       DO 1 I = 1,6
       DO 1 J = 1,6
       ROT(I,J) = 0.0
 1     RLOCAL(I,J) = 0.0
       CALL ROTAC(IH,XL,ROT)
c....  Matriz de rigidez referida al sistema local
       S1 = E*AX/XL
       RLOCAL(1,1) =  S1
       RLOCAL(1,4) = -S1
       RLOCAL(4,1) = -S1
       RLOCAL(4,4) =  S1
       CALL ROTAC(IH,XL,ROT)
c.... Matriz de rigidez referida al sistema global
       DO 4 I=1,6
       DO 4 J=1,6
       REAC(I,J) = 0.0
       DO 4 K=1,6
  4   REAC(I,J) = REAC(I,J) + RLOCAL(I,K)*ROT(K,J)
       DO 5 I=1,6
       DO 5 J=1,6
  5   RLOCAL(I,J) = REAC(I,J)
       IH14 = IH
       WRITE(14,REC=IH14) ((RLOCAL(I,J),I=1,6),J=1,6)
       DO 6 I=1,6
       DO 6 J=1,6
  6   REAC(I,J) = 0.0
       DO 7 I=1,6
       DO 7 J=1,6
       REAC(I,J) = 0.0
       DO 7 K=1,6
  7   REAC(I,J) = REAC(I,J) + ROT(K,I)*RLOCAL(K,J)
       DO 8 I=1,6
       DO 8 J=1,6
  8   RLOCAL(I,J) = REAC(I,J)
      RETURN
      END
C==============================
      SUBROUTINE P3D(IH,RLOCAL)
      INCLUDE'AME'
      DIMENSION RLOCAL(12,12),REAC(12,12),ROT(12,12)
c....
       IH30 = IMAT(IH)
       READ(30,REC=IH30)E,G,AX,CIX,CIY,CIZ
       DO 1 I = 1,12
       DO 1 J = 1,12
       ROT(I,J) = 0.0
 1     RLOCAL(I,J) = 0.0
       CALL ROTAC(IH,XL,ROT)
c....  Matriz de rigidez referida al sistema local
       S1 = E*CIZ/XL
       S2 = E*CIY/XL
       RLOCAL(1,1)   =   E*AX/XL
       RLOCAL(2,2)   =   12.0*S1 /XL**2.0
       RLOCAL(6,2)   =   6.0*S1/XL
       RLOCAL(3,3)   =   12.0*S2/XL**2.0
       RLOCAL(5,3)   =  -6.0*S2/XL
       RLOCAL(4,4)   =   G*CIX/XL
       RLOCAL(5,5)   =   4.0*S2
       RLOCAL(6,6)   =   4.0*S1
       RLOCAL(8,8)   =   RLOCAL(2,2)
       RLOCAL(7,1)   =  -RLOCAL(1,1)
       RLOCAL(7,7)   =   RLOCAL(1,1)
       RLOCAL(8,2)   =  -RLOCAL(2,2)
       RLOCAL(8,6)   =  -RLOCAL(6,2)
       RLOCAL(9,3)   =  -RLOCAL(3,3)
       RLOCAL(9,5)   =  -RLOCAL(5,3)
       RLOCAL(9,9)   =   RLOCAL(3,3)
       RLOCAL(10,4)  =  -RLOCAL(4,4)
       RLOCAL(10,10) =   RLOCAL(4,4)
       RLOCAL(11,3)  =   RLOCAL(5,3)
       RLOCAL(11,5)  =   RLOCAL(5,5)/2.0
       RLOCAL(11,9)  =  -RLOCAL(5,3)
       RLOCAL(11,11) =   RLOCAL(5,5)
       RLOCAL(12,2)  =   RLOCAL (6,2)
       RLOCAL(12,6)  =   RLOCAL(6,6)/2.0
       RLOCAL(12,8)  =  -RLOCAL(6,2)
       RLOCAL(12,12) =   RLOCAL(6,6)
       DO 3 J1 = 1,12
       DO 3 J2 = J1,12
  3   RLOCAL(J1,J2) = RLOCAL(J2,J1)
       CALL ROTAC(IH,XL,ROT)
c.... Matriz de rigidez referida al sistema global
       DO 4 I=1,12
       DO 4 J=1,12
       REAC(I,J) = 0.0
       DO 4 K=1,12
  4   REAC(I,J) = REAC(I,J) + RLOCAL(I,K)*ROT(K,J)
       DO 5 I=1,12
       DO 5 J=1,12
  5   RLOCAL(I,J) = REAC(I,J)
       IH14 = IH
       WRITE(14,REC=IH14) ((RLOCAL(I,J),I=1,12),J=1,12)
       DO 6 I=1,12
       DO 6 J=1,12
  6   REAC(I,J) = 0.0
       DO 7 I=1,12
       DO 7 J=1,12
       REAC(I,J) = 0.0
       DO 7 K=1,12
  7   REAC(I,J) = REAC(I,J) + ROT(K,I)*RLOCAL(K,J)
       DO 8 I=1,12
       DO 8 J=1,12
  8   RLOCAL(I,J) = REAC(I,J)
      RETURN
      END
C================================
      SUBROUTINE ROTAC(IH,XL,ROT)
      INCLUDE'AME'
      DIMENSION ROT(12,12),JP(2)
       DO 1 I=1,NGLL
       DO 1 J=1,NGLL
1     ROT(I,J) = 0.0
       DO 2 J=1,NNPE
       I = (IH-1)*NNPE + J
2     JP(J) = INCID(I)
       I = (JP(1) - 1)*NCOPN
       J = (JP(2) - 1)*NCOPN
       XCL = COOR(J + 1) - COOR(I + 1)
       YCL = COOR(J + 2) - COOR(I + 2)
      IF((INDEL.EQ.6).OR.(INDEL.EQ.7)) ZCL = COOR(J + 3) - COOR(I + 3)
       XL = DSQRT(XCL**2.0 + YCL**2.0)
      IF((INDEL.EQ.6).OR.(INDEL.EQ.7)) XL = DSQRT(XCL**2.0 + YCL**2.0
     1                                      + ZCL**2.0)
       CX = XCL/XL
       CY = YCL/XL
      IF((INDEL.EQ.6).OR.(INDEL.EQ.7)) CZ = ZCL/XL
      IF(INDEL.EQ.2) GO TO 10
      IF((INDEL.EQ.4).OR.(INDEL.EQ.5)) GO TO 20
      IF(INDEL.EQ.6) GO TO 30
      IF(INDEL.EQ.7) GO TO 40
c.... Matriz de rotacion del elemento: Armadura Plana
 10     ROT(1,1) =  CX
        ROT(1,2) =  CY
        ROT(2,1) = -CY
        ROT(2,2) =  CX
                     ROT(3,3) =  CX
                     ROT(3,4) =  CY
                     ROT(4,3) = -CY
                     ROT(4,4) =  CX
       GO TO 90
c.... Matriz de rotacion de los elementos: Parrilla y Portico Plano
 20    ROT(1,1) =  CX
       ROT(1,2) =  CY
       ROT(2,1) = -CY
       ROT(2,2) =  CX
       ROT(3,3) = 1.0
                     ROT(4,4) = CX
	             ROT(4,5) = CY
	             ROT(5,4) = -CY
	             ROT(5,5) = CX
	             ROT(6,6) = 1.0
        GO TO 90
c....   Matriz de rotacion del elemento: armadura espacial
 30    C1 = CX*CY
       C2 = CY*CZ
       CD = DSQRT(CX**2.0 + CZ**2.0)
      IF(CD - 0.001)120,130,130
c       (a) Caso elemento vertical
 120   ROT(1,2)  =  CY
       ROT(2,1)  = -CY
       ROT(3,3) =  1.0
                      ROT(4,5) =  CY
                      ROT(5,4) = -CY
                      ROT(6,6) =  1.0
       GO TO 115
c       (b) caso elemento con orientacion arbitraria
 130   ROT(1,1)   =   CX
       ROT(1,2)   =   CY
       ROT(1,3)   =   CZ
       ROT(2,1)   = -C1/CD
       ROT(2,2)   =  CD
       ROT(2,3)   = -C2/CD
       ROT(3,1)   = -CZ/CD
       ROT(3,3)   =  CX/CD
                          ROT(4,4)   =   ROT(1,1)
                          ROT(4,5)   =   ROT(1,2)
                          ROT(4,6)   =   ROT(1,3)
                          ROT(5,4)   =   ROT(2,1)
                          ROT(5,5)   =   ROT(2,2)
                          ROT(5,6)   =   ROT(2,3)
                          ROT(6,4)   =   ROT(3,1)
                          ROT(6,6)   =   ROT(3,3)
 115  CONTINUE
       GO TO 90
c.... Matriz de rotacion de los elementos: Portico Espacial
c       (a) Caso elemento vertical
  40  Q = DSQRT(CX**2.0 + CZ**2.0)
      IF(Q - 0.001)5, 6, 6
   5   ROT(1,2) =  CY
       ROT(2,1) = -CY
       ROT(3,3) = 1.0
                     ROT(4,5) =  CY
                     ROT(5,4) = -CY
                     ROT(6,6) = 1.0
                                   ROT(7,8) =  CY
                                   ROT(8,7) = -CY
                                   ROT(9,9) = 1.0
                                                 ROT(10,11) =  CY
                                                 ROT(11,10) = -CY
                                                 ROT(12,12) = 1.0
c       (b) caso elemento con orientacion arbitraria
	GO TO 8
  6    ROT(1,1) = CX
       ROT(1,2) = CY
       ROT(1,3) = CZ
       ROT(2,1) = -CX*CY/Q
       ROT(2,2) = Q
       ROT(2,3) = -CY*CZ/Q
       ROT(3,1) = -CZ/Q
       ROT(3,3) = CX/Q
                     ROT(4,4) = CX
                     ROT(4,5) = CY
                     ROT(4,6) = CZ
                     ROT(5,4) = -CX*CY/Q
                     ROT(5,5) = Q
                     ROT(5,6) = -CY*CZ/Q
                     ROT(6,4) = -CZ/Q
                     ROT(6,6) = CX/Q
                                   ROT(7,7) = CX
                                   ROT(7,8) = CY
                                   ROT(7,9) = CZ
                                   ROT(8,7) = -CX*CY/Q
                                   ROT(8,8) = Q
                                   ROT(8,9) = -CY*CZ/Q
                                   ROT(9,7) = -CZ/Q
                                   ROT(9,9) = CX/Q
                                                 ROT(10,10) = CX
                                                 ROT(10,11) = CY
                                                 ROT(10,12) = CZ
                                                 ROT(11,10) = -CX*CY/Q
                                                 ROT(11,11) = Q
                                                 ROT(11,12) = -CY*CZ/Q
                                                 ROT(12,10) = -CZ/Q
                                                 ROT(12,12) = CX/Q
  8   CONTINUE
 90   CONTINUE
      RETURN
      END
C=====================
      SUBROUTINE MONTA
      INCLUDE'AME'
      DIMENSION RLOCAL(12,12)
c====
c====   ESTA SUBRUTINA ENSAMBLA LAS MATRICES LOCALES DE LOS ELEMENTOS
c====   EN UNA MATRIZ CUADRADA.
c====
       DO 1 IH = 1,NE
      IF(INDEL.EQ.1) CALL AXI(IH,RLOCAL)
      IF(INDEL.EQ.2) CALL A2D(IH,RLOCAL)
      IF(INDEL.EQ.3) CALL VIG(IH,RLOCAL)
      IF(INDEL.EQ.4) CALL P2D(IH,RLOCAL)
      IF(INDEL.EQ.5) CALL PAR(IH,RLOCAL)
      IF(INDEL.EQ.6) CALL A3D(IH,RLOCAL)
      IF(INDEL.EQ.7) CALL P3D(IH,RLOCAL)
c====
       DO 2 I = 1,NNPE
       N1 = INCID((IH-1)*NNPE +I)
       I1 = (I-1)*NGLN
       J1 = (N1-1)*NGLN
       DO 2  J = 1,NNPE
       N2 = INCID((IH - 1)*NNPE + J)
       I2 = (J - 1)*NGLN
       J2 = (N2 - 1)*NGLN
       DO 2 K = 1,NGLN
       IL = I1 + K
       IG = J1 + K
       DO 2 L = 1,NGLN
       JL = I2 + L
       JG = J2 + L
       GLOBAL(IG,JG) =  GLOBAL(IG,JG) + RLOCAL(IL,JL)
  2   CONTINUE
  1   CONTINUE
      RETURN
      END
C======================
      SUBROUTINE CONTOR
      INCLUDE'AME'
c====
c====   ESTA  SUBRUTINA INTRODUCE  LAS CONDICIONES DE CONTORNO EN LA
c====   MATRIZ GLOBAL DEL SISTEMA Y CORRIGE EL TERMINO INDEPENDIENTE.
c====
       IF(NNAE.EQ.0) GO TO 101     ! Consideracion de apoyos elasticos
       DO 100 I = 1,NNAE
       J = NAPE(I)
       K1 = (J - 1)*NGLN
       DO 100 JK = 1,NGLN
       JJ = (I - 1)*NGLN + JK
       IG = K1 + JK
 100  GLOBAL(IG,IG) = GLOBAL(IG,IG) + DISA(JJ)
 101  CONTINUE
       DO 1 I = 1,NNDP
       J = NAP(I)
       K1 = (J - 1)*NGLN
       DO 1 JK = 1,NGLN
       JJ = (I - 1)*NGLN + JK
       IG = K1 + JK
       IF((DIS(JJ)).EQ.200.0) GO TO 1                  ! Direccion libre
       DO 2 JI = 1,NTEC
       CARGA(JI) = CARGA(JI) - GLOBAL(IG,JI)*DIS(JJ)   ! Correccion del vector de cargas
       GLOBAL(IG,JI) =  0.0
       GLOBAL(JI,IG) =  0.0                            ! Ceros en las filas de la matriz global
   2  CONTINUE                                         ! Ceros en las columnas de la matriz global
       GLOBAL(IG,IG) = 1.0                             ! Uno (1.0) en la diagonal
       CARGA(IG) = DIS(JJ)                             ! Desplazamiento conocido en el vector de cargas
   1  CONTINUE
      RETURN
      END
C=======================
      SUBROUTINE FACTORA
      INCLUDE'AME'
c....
c.... ESTA SUBRUTINA  FACTORIZA  LA MATRIZ GLOBAL  DEL SISTEMA SE-
c.... GUN EL METODO DE CHOLESKY.
c....
      IF(GLOBAL(1,1).LE.0.0) GO TO 10
       DO 1  J = 2,NTEC
       J1 = J -1
      IF(J1.EQ.1) GO TO 3
       DO 2 I = 2,J1
       SUM = GLOBAL(I,J)
       I1 = I - 1
       DO 4 K = 1,I1
  4   SUM = SUM - GLOBAL(K,I)*GLOBAL(K,J)
       GLOBAL(I,J) = SUM
  2   CONTINUE
  3   SUM = GLOBAL(J,J)
       DO 5 K = 1,J1
       TEMP = GLOBAL(K,J)/GLOBAL(K,K)
       SUM = SUM - TEMP*GLOBAL(K,J)
  5   GLOBAL(K,J) = TEMP
      IF (SUM.LE.0.0) GO TO 10
  1   GLOBAL(J,J) = SUM
        GO TO 20
 10   WRITE(OUP,12) SUM
12    FORMAT(//1X,'SUBRUTINA INAPROPIADA PARA RESOLVER EL SISTEMA',/1X,'
     1SUM ='D15.3,5X,'SE SUSPENDE EL PROCESAMIENTO')
      STOP
 20   CONTINUE
      RETURN
      END
C========================
      SUBROUTINE SOLUCION
      INCLUDE'AME'
c....
c.... ESTA SUBRUTINA RESUELVE EL SISTEMA DE ECUACIONES RESULTANTE
c.... MEDIANTE UNA SUSTITUCION REGRESIVA.
c....
       DO 1 I=1,NTEC
       SUM = CARGA(I)
       K1 = I- 1
      IF(I.EQ.1) GO TO 1
       DO 2 K = 1,K1
 2    SUM = SUM - GLOBAL(K,I)*DES(K)
 1    DES(I) = SUM
       DO 3 I=1,NTEC
       DES(I) = DES(I)/GLOBAL(I,I)
 3    CONTINUE
       DO 4 I1=1,NTEC
       I = NTEC - I1 + 1
       SUM = DES(I)
       K2 = I + 1
       IF(I.EQ.NTEC) GO TO 4
       DO 5 K = K2,NTEC
 5    SUM = SUM - GLOBAL(I,K)*DES(K)
 4    DES(I) = SUM
      RETURN
      END
C=======================
      SUBROUTINE ESCRIBE
      INCLUDE'AME'
c.....
       WRITE(OUP,10)
       IF(INDEL.EQ.1) THEN
       WRITE(OUP,20)
       ELSE IF((INDEL.EQ.2).OR.(INDEL.EQ.3)) THEN
       WRITE(OUP,30)
       ELSE IF((INDEL.EQ.4).OR.(INDEL.EQ.5))THEN
       WRITE (OUP,40)
       ELSE IF(INDEL.EQ.6)THEN
       WRITE (OUP,40)
       ELSE IF(INDEL.EQ.7)THEN
       WRITE (OUP,50)
       END IF
       DO 1 J = 1,NN
       LL = J*NGLN
       L = LL -NGLN + 1
   1  WRITE(OUP,60) J,(DES(K),K=L,LL)
c.... Formatos de escritura
  10  FORMAT(/1X,'DESPLAZAMIENTOS DE LOS NODOS',/1X,
     1           '============================')
  20  FORMAT(1X,'NODO',4X,'DIR.1')
  30  FORMAT(1X,'NODO',4X,'DIR.1',9X,'DIR.2')
  40  FORMAT(1X,'NODO',4X,'DIR.1',9X,'DIR.2',9X,'DIR.3')
  50  FORMAT(1X,'NODO',4X,'DIR.1',9X,'DIR.2',9X,'DIR.3',9X,'DIR.4',9X,'D
     1IR.5',9X,'DIR.6')
  60  FORMAT(I3,6E14.5)
      RETURN
      END
C======================
      SUBROUTINE ACCION
      INCLUDE'AME'
      DIMENSION RLOCAL(12,12),S(12),DJ(12),CAUX(12),JP(2)
c.....
       GO TO (300,310,320,330,340,350,360),INDEL
 300  WRITE(OUP,101)
       WRITE(OUP,201)
        GO TO  900
 310  WRITE(OUP,102)
       WRITE(OUP,202)
        GO TO 900
 320  WRITE(OUP,103)
       WRITE(OUP,203)
        GO TO 900
 330  WRITE(OUP,105)
       WRITE(OUP,205)
        GO TO 900
 340  WRITE(OUP,104)
       WRITE(OUP,204)
        GO TO 900
 350  WRITE(OUP,106)
       WRITE(OUP,206)
        GO TO 900
 360  WRITE(OUP,107)
       WRITE(OUP,207)
 900  CONTINUE
       DO 5 I =1,NTEC
  5   CARGA(I) = 0.0
       DO 1 IH=1,NE
       L1 = (IH-1)*NNPE + 1
       IH14 = IH
       READ(14,REC=IH14)((RLOCAL(I,J),I=1,NGLL),J=1,NGLL)
       J1 = 0
       DO 2 L=1,NNPE
       I = (IH-1)*NNPE + L
       N1 = INCID(I)
       JP(L) = N1
       DO 2 I=1,NGLN
       J1 = J1 + 1
       K1 = (N1-1)*NGLN + I
 2    DJ(J1) = DES(K1)
       DO 3 J=1,NGLL
       S(J) = 0.0
       DO 3 J1=1,NGLL
 3    S(J) = S(J) + RLOCAL(J,J1)*DJ(J1)
      IF(NEC.EQ.0) GO TO 46
       IE = 0
 45   IE = IE + 1
       IC = NQ(IE)
      IF(IH - IC) 50, 51, 50
  51  IL = IL + 1
       IH15 = IL
       READ(15,REC=IH15) (CAUX(I1),I1 = 1,NGLL)
       DO 43 I3=1,NGLL
  43  S(I3) = S(I3) + CAUX(I3)
  50  IF(IE - NEC)45, 46, 46
  46  CONTINUE
       CALL REACCI(IH,S,JP)
       GO TO (15,16,18,19,19,19,20),INDEL
  15  WRITE(OUP,30)IH,INCID(L1),(S(J1),J1=1,1),INCID(L1+1),(S(I),I=2,2)
        GO TO 1
  16  WRITE(OUP,31)IH,INCID(L1),(S(J1),J1=1,2),INCID(L1+1),(S(I),I=3,4)
        GO TO 1
  18  WRITE(OUP,31)IH,INCID(L1),(S(J1),J1=1,2),INCID(L1+1),(S(I),I=3,4)
        GO TO 1
  19  WRITE(OUP,32)IH,INCID(L1),(S(J1),J1=1,3),INCID(L1+1),(S(I),I=4,6)
        GO TO 1
  20  WRITE(OUP,33)IH,INCID(L1),(S(J1),J1=1,6),INCID(L1+1),(S(I),I=7,12)
1     CONTINUE
c==== FORMATOS DE ESCRITURA
 101  FORMAT(/1X,'ELEMENTO: AXIAL',/1X,
     1           '===============')
 201  FORMAT(1X,'ELEM',3X,'NODO',4X,' FX ')
 102  FORMAT(/1X,'ELEMENTO: ARMADURA PLANA',/1X,
     1           '========================')
 202  FORMAT(1X,'ELEM',3X,'NODO',4X,' FX ',9X,' FY ')

 103  FORMAT(/1X,'ELEMENTO: VIGA DE EJE RECTO',/1X,
     1           '===========================')
 203  FORMAT(1X,'ELEM',3X,'NODO',4X,' FY ',9X,' MZ ')

 104  FORMAT(/1X,'ELEMENTO: PARRILLA',/1X,
     1           '==================')
 204  FORMAT(1X,'ELEM',3X,'NODO',4X,' MX ',9X,' MY ',9X,' FZ ')

 105  FORMAT(/1X,'ELEMENTO: PORTICO PLANO',/1X,
     1           '=======================')
 205  FORMAT(1X,'ELEM',3X,'NODO',4X,' FX ',9X,' FY ',9X,' MZ ')

 106  FORMAT(/1X,'ELEMENTO: ARMADURA ESPACIAL',/1X,
     1           '===========================')
 206  FORMAT(1X,'ELEM',3X,'NODO',4X,' FX ',9X,' FY ',9X,' FZ ')
 107  FORMAT(/1X,'ELEMENTO: PORTICO ESPACIAL',/1X,
     1           '==========================')
 207  FORMAT(1X,'ELEM',3X,'NODO',4X,' FX ',9X,' FY ',9X,' FZ ',9X,' MX '
     1,9X,' MY',9X,'  MZ ')
 30   FORMAT(I3,I7,E13.4,/(I10,E13.4))
 31   FORMAT(I3,I7,2E13.4,/(I10,2E13.4))
 32   FORMAT(I3,I7,3E13.4,/(I10,3E13.4))
 33   FORMAT(I3,I7,6E13.4,/(I10,6E13.4))
      RETURN
      END
C===============================
      SUBROUTINE REACCI(IH,S,JP)
      INCLUDE'AME'
      DIMENSION JP(2),S(12),SR(6),ROT(12,12)
c.....
       II = JP(1)
       JJ = JP(2)
       IA = 0
 11   IA = IA + 1
       JN = NAP(IA)
      IF((JN.NE.II).AND.(JN.NE.JJ)) GO TO 36
      IF(JN - II)14,15,14
 15     GO TO(17,18,19,20,20,20,21),INDEL
 17   CARGA(JP(1)) = CARGA(JP(1)) + S(1)
        GO TO 36
 18   CALL ROTAC(IH,XL,ROT)
       DO 16 I3 = 1,NGLN
       SR(I3) = 0.0
       DO 16 J3 = 1,NGLN
 16   SR(I3) = SR(I3) + S(J3)*ROT(J3,I3)
       CARGA(2*JP(1)-1) = CARGA(2*JP(1)-1) + SR(1)
       CARGA(2*JP(1))   = CARGA(2*JP(1))   + SR(2)
        GO TO 36
 19   CARGA(2*JP(1) - 1)   =  CARGA(2*JP(1) - 1) + S(1)
      CARGA(2*JP(1))   =  CARGA(2*JP(1))     + S(2)
        GO TO 36
 20   CALL ROTAC(IH,XL,ROT)
       DO 160 I3 = 1,NGLN
       SR(I3) = 0.0
       DO 160 J3 = 1,NGLN
 160  SR(I3) = SR(I3) + S(J3)*ROT(J3,I3)
       CARGA(3*JP(1)-2) = CARGA(3*JP(1)-2) + SR(1)
       CARGA(3*JP(1)-1) = CARGA(3*JP(1)-1) + SR(2)
       CARGA(3*JP(1))   = CARGA(3*JP(1))   + SR(3)
         GO TO 36
  21  CALL ROTAC(IH,XL,ROT)
       DO 161 I3 = 1,NGLN
       SR(I3) = 0.0
       DO 161 J3 = 1,NGLN
 161  SR(I3) = SR(I3) + S(J3)*ROT(J3,I3)
       CARGA(6*JP(1)-5) = CARGA(6*JP(1)-5) + SR(1)
       CARGA(6*JP(1)-4) = CARGA(6*JP(1)-4) + SR(2)
       CARGA(6*JP(1)-3) = CARGA(6*JP(1)-3) + SR(3)
       CARGA(6*JP(1)-2) = CARGA(6*JP(1)-2) + SR(4)
       CARGA(6*JP(1)-1) = CARGA(6*JP(1)-1) + SR(5)
       CARGA(6*JP(1))   = CARGA(6*JP(1))   + SR(6)
	 GO TO 36
 14   CONTINUE
      IF(JN - JJ)36,25,36
25    GO TO(40,41,42,43,43,43,44),INDEL
40    CARGA(JP(2))   = CARGA(JP(2)) + S(2)
        GO TO 36
41    CALL ROTAC(IH,XL,ROT)
       KI = NGLN + 1
       DO 26 I3 = KI,NGLL
       L3 = I3 - 2
       SR(L3) = 0.0
       DO 26 J3=3,4
 26   SR(L3) = SR(L3) + S(J3)*ROT(J3,I3)
       CARGA(2*JP(2)-1) = CARGA(2*JP(2)-1)  + SR(1)
       CARGA(2*JP(2))   = CARGA(2*JP(2))    + SR(2)
          GO TO 36
 42   CARGA(2*JP(2) - 1 )   =  CARGA(2*JP(2) - 1) + S(3)
      CARGA(2*JP(2))        =  CARGA(2*JP(2))     + S(4)
          GO TO 36
 43   CALL ROTAC(IH,XL,ROT)
       KI = NGLN + 1
       DO 260 I3 = KI,NGLL
       L3 = I3 - NGLN
       SR(L3) = 0.0
       DO 260 J3 = KI,NGLL
 260  SR(L3) = SR(L3) + S(J3)*ROT(J3,I3)
       CARGA(3*JP(2)-2) = CARGA(3*JP(2)-2)  + SR(1)
       CARGA(3*JP(2)-1) = CARGA(3*JP(2)-1)  + SR(2)
       CARGA(3*JP(2))   = CARGA(3*JP(2))    + SR(3)
          GO TO 36
  44  CALL ROTAC(IH,XL,ROT)
       KI = NGLN + 1
       DO 261 I3 = KI,NGLL
       L3 = I3 - NGLN
       SR(L3) = 0.0
       DO 261 J3 = KI,NGLL
 261  SR(L3) = SR(L3) + S(J3)*ROT(J3,I3)
       CARGA(6*JP(2)-5) = CARGA(6*JP(2)-5)  + SR(1)
       CARGA(6*JP(2)-4) = CARGA(6*JP(2)-4)  + SR(2)
       CARGA(6*JP(2)-3) = CARGA(6*JP(2)-3)  + SR(3)
       CARGA(6*JP(2)-2) = CARGA(6*JP(2)-2)  + SR(4)
       CARGA(6*JP(2)-1) = CARGA(6*JP(2)-1)  + SR(5)
       CARGA(6*JP(2))   = CARGA(6*JP(2))    + SR(6)
 36   CONTINUE
	  IF(IA-NNDP)11,27,27
  27  CONTINUE
       RETURN
      END
C======================
      SUBROUTINE ESCREA
      INCLUDE'AME'
c.....
        WRITE(OUP,10)
        GO TO(501,502,502,504,504,504,505),INDEL
 501	WRITE(OUP,20)
        GO TO 500
 502	WRITE(OUP,30)
        GO TO 500
 504    WRITE(OUP,40)
        GO TO 500
 505    WRITE(OUP,50)
        GO TO 500
 500  CONTINUE
	DO 1 I=1,NNDP
	K = NAP(I)
	LL = K*NGLN
	L = LL - NGLN + 1
  1   WRITE(OUP,60)K,(CARGA(J),J=L,LL)
c==== Formatos de escritura
  10  FORMAT(/1X,'REACCIONES EN LOS APOYOS',/1X,
     1           '=======================')
  20  FORMAT(1X,'NODO',3X,'REAC 1')
  30  FORMAT(1X,'NODO',3X,'REAC 1',8X,'REAC 2')
  40  FORMAT(1X,'NODO',3X,'REAC 1',8X,'REAC 2',8X,'REAC 3')
  50  FORMAT(1X,'NODO',3X,'REAC 1',8X,'REAC 2',8X,'REAC 3',8X,'REAC 4',8
     1X,'REAC 5',8X,'REAC 6')
  60  FORMAT(I3,6E14.5)
       RETURN
      END
