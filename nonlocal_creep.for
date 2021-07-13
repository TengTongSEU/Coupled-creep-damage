C ======================================================================
C User Subroutine for Conventional Gradient Damage
C Type: 8 noded SERENDIPITY element for displacement and e_eqv
C ======================================================================
C Please cite the related paper as:
C Authors- S Sarkar, IV Singh, BK Mishra, AS Shedbale, LH Poh
C Title- A comparative study and ABAQUS implementation of conventional and 
C        localizing gradient enhanced damage models
C Journal- Finite Elements in Analysis and Design 160 (2019) 1-31.
C DOI-10.1016/j.finel.2019.04.001
C ======================================================================
C Material properties to be given through the input file (*.inp), are
C
C For U1 element:
C PROPS(1) = Young's modulus (E)
C PROPS(2) = Poisson's ratio (nu)
C PROPS(3) = Averaging lenght parameter squared (c)
C PROPS(4) = Threshold fracture strain (ki)
C PROPS(5) = Second parameter in exp damage law (beta)
C
C ---- Used variables ---------------------
C N_ELEM - number of elements used in the model divided
C          by 2 - (N+N_UMAT)/2 (to be changed for each model)
C NSDV - solution dependent variables for the user element
C          strain(e)(x3)---> SDV = 1,2,3
C          loc_e_eqv(x1)---> SDV = 4
C          gp_e_eqv(x1) ---> SDV = 5
C          k_gp(x1)     ---> SDV = 6
C          k0(x1)       ---> SDV = 7
C          damage(D)(x1)---> SDV = 8
C          ....total=8
C NSDV - overall solution dependent variables (NSDV + 2), where
C        the additional 2 variables are the: time and iteration number
C ======================================================================
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     &     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     &     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     &     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     &     PERIOD)
C     ==================================================================
      INCLUDE 'ABA_PARAM.INC'
C     ==================================================================
      PARAMETER(N_ELEM=10000,NSDV=120,NGP=9)
C     ==================================================================
C     Initialization for all the element types
C     ==================================================================
      DIMENSION RHS(MLVARX,1),AMATRX(NDOFEL,NDOFEL),
     &     SVARS(NSVARS),ENERGY(8),PROPS(NPROPS),COORDS(MCRD,NNODE),
     &     U(NDOFEL),DU(MLVARX,1),V(NDOFEL),A(NDOFEL),TIME(2),
     &     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     &     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     &     JPROPS(*)
C
      DIMENSION GP(2,NGP),GW(2,NGP),XI(2),AJACOB(2,2),dNdxi(NNODE,2),
     & dHdxi(4,2),dNdx(NNODE,2),dHdx(4,2),AJABOBINV(2,2),AN(1,NNODE),
     & AH(1,4),BU(3,2*NNODE),BE(2,4),DMAT(3,3)
C
      DIMENSION E(3,NGP), ALOC_E_EQV(1,NGP), AK_GP(1,NGP), AK0(1,NGP)
      DIMENSION GP_E_EQV(1,NGP), DEEQV_DE(3), D(1,NGP)
	DIMENSION DISP(2*NNODE,1), E_EQV(4,1)
      DIMENSION DISP_INC(2*NNODE,1)
      DIMENSION E_TOT_INC(3,NGP)
C
      DIMENSION AK_UU(2*NNODE,2*NNODE),AK_EE(4,4),AK_UE(2*NNODE,4)
      DIMENSION AK_EU(4,2*NNODE),F_EE(4,1)
      DIMENSION F_UU(2*NNODE,1)
      DIMENSION F_UU_T1(2*NNODE,3),F_UU_T2(2*NNODE,2*NNODE)
      DIMENSION F_UU_T3(2*NNODE,1)
      DIMENSION E_RHS(3,1)
      
      DIMENSION E_TOT(3,NGP), UEL_STRESS_EFF(3,NGP),
     +          UEL_STRESS_EFF_PRE(3,NGP),E_MPS(3,NGP),
     +          UEL_STRESS(3,NGP),E_SHRINK(3,NGP)
      
      DIMENSION YU_T_MU(30), YU_DELTA_Y_MU(30)
      DIMENSION YU_L_MU(30), YU_A_MU(30), YU_LAM_MU(30)
      DIMENSION YU_IN(30,3,NGP), YU_EN(30,3,NGP), YU_S(3,NGP)
      DIMENSION YU_STRAIN(3,NGP),YU_DELTA_STRESS(3,NGP)
      DIMENSION E_DRY(3,NGP)
      
      DIMENSION AK_UU_T1(3,2*NNODE),AK_UU_T2(2*NNODE,2*NNODE)
      
      DIMENSION YU_D(3,3)
C
      COMMON/KUSER/USRVAR(N_ELEM,NSDV,NGP)
C
C     ******************************************************************
C     Constructing element U1
C     ******************************************************************
C     ==================================================================
C     Nodal variables saved in local array
C     ==================================================================
      DISP = 0.0
      DISP_INC = 0.0D0
      E_EQV = 0.0
C      
      DO I=1,NNODE
           DISP(2*I-1,1)     = U(3*I-2)
           DISP(2*I,1)       = U(3*I-1)	
      END DO
	DO I=1,4
		 E_EQV(I,1)        = U(3*I)
      END DO
      DO I=1,NNODE
           DISP_INC(2*I-1,1) = DU(3*I-2,1)
           DISP_INC(2*I,1)   = DU(3*I-1,1)	
      END DO
C
C     ==================================================================
C     Saving time and increment no.
C     ==================================================================
      TIMEZ=USRVAR(JELEM,9,1)
      IF (TIMEZ.LT.TIME(2)) THEN
          USRVAR(JELEM,9,1)=TIME(2)
          USRVAR(JELEM,10,1)=0.0
      ELSE
          USRVAR(JELEM,10,1)=USRVAR(JELEM,10,1)+1.0
      ENDIF
      STEPITER=USRVAR(JELEM,10,1)
C     ==================================================================
C     Material parameters 
C     ==================================================================
      EMOD = PROPS(1)
      ENU = PROPS(2)
      C = PROPS(3)
      AKI = PROPS(4)
      beta = PROPS(5)
C       
      alpha = 0.99 
      AK = 10.
      c1 = (AK-1.)/(2.*AK*(1.-2.*ENU))
      c2 = ((AK-1.)/(1.-2.*ENU))**2.
      c3 = 2.*AK/(1.+ENU)**2.
C  c     
      YU_Q2       = 200.0D-6
      YU_Q3       =  20.0D-6
      YU_Q4       =  46.0D-6
C
C      YU_Q2       = 1.0D-16
C      YU_Q3       = 1.0D-16
C      YU_Q4       = 0.0D-6
C      
      TIME_PRIME  = 13.9D0
      TIME_ACT    = 14.6D0
C      
      NS_YU = 15 
      YU_N = 0.1D0
C
C     ==================================================================
C     Calculating materials constitutive matrix (plane strain)
C     ==================================================================
      DMAT = 0.
      DMAT(1,1) = EMOD*(1.-ENU)/((1.+ENU)*(1.-2.*ENU))
      DMAT(2,2) = EMOD*(1.-ENU)/((1.+ENU)*(1.-2.*ENU))
      DMAT(3,3) = EMOD*(0.5-ENU)/((1.+ENU)*(1.-2.*ENU))
      DMAT(1,2) = EMOD*ENU/((1.+ENU)*(1.-2.*ENU))
      DMAT(2,1) = EMOD*ENU/((1.+ENU)*(1.-2.*ENU))
C     ==================================================================
C     Calculating Paramters for creep
C     ==================================================================
C     Retardation Time
      YU_T_MU(1) = 1.0D-6
	DO I = 2, NS_YU
	   YU_T_MU(I) = YU_T_MU(I-1)*10.0D0
      ENDDO     
C     Dtime/Retardation time
	DO I = 1, NS_YU
         YU_DELTA_Y_MU(I) = DTIME/YU_T_MU(I)
      ENDDO      
C     A_MU      
	DO I = 1, NS_YU
C	    TEMX = 3*YU_T_MU(I)
C	    X1   = 0.336/TEMX**2.40/(10+TEMX**0.60)
C	    X2   = 0.528/TEMX**1.80/(10+TEMX**0.60)**2
C          X3   = 0.432/TEMX**1.20/(10+TEMX**0.60)**3
C	    X4   = 1.296/TEMX**0.60/(10+TEMX**0.60)**4                                                 
C	    X5   = TEMX**3/2
C	    X6   = (X1+X2+X3-X4)*X5
C          YU_GA = 2.35D0
C	    YU_A_MU(I) = X6*log(10.0d0)*YU_GA 

          TEMX = 3*YU_T_MU(I)
          X1 = YU_N*(YU_N-1.0D0)*(YU_N-2.0D0)*TEMX**(YU_N-3.0D0)
          X2 = (1.0D0 + TEMX**YU_N)
          X3 = YU_N*YU_N*(YU_N-1.0D0)*TEMX**(2.0D0*YU_N-3.0D0)
          X4 = (1.0D0 + TEMX**YU_N)**2.0D0
          X5 = YU_N*YU_N*(2.0D0*YU_N-2.0D0)*TEMX**(2.0D0*YU_N-3.0D0)
          X6 = (1.0D0 + TEMX**YU_N)**2.0D0
          X7 = 2.0D0*YU_N*YU_N*YU_N*TEMX**(3.0D0*YU_N-3.0D0)
          X8 = (1.0D0 + TEMX**YU_N)**3.0D0
          X9   = TEMX**3.0D0/2.0D0
          X10   = (X1/X2 - X3/X4 - X5/X6 + X7/X8)*X9
	    YU_A_MU(I) = X10*log(10.0d0)
C         TO REMDEDY THE DEFICIENCY IN THE EARLY-AGE QUICK CREEP, COEFFICIENT = 2  
          YU_A_MU(I) = 2.0D0 * YU_A_MU(I)
      ENDDO
C     LAMBDA
	DO I = 1, NS_YU
	   TEMP_YU      = EXP( 0.0D0-YU_DELTA_Y_MU(I) )
	   YU_LAM_MU(I) = (1.0D0 - TEMP_YU)/YU_DELTA_Y_MU(I)
      ENDDO      
C     CALCULATE THE INEALSTIC INCREMENT OF STRAIN DUE TO CREEP
C     YU_E: EFFECTIVE MODULUS CONSIDERED THE CREEP
C     YU_V: COEFFICIENT OF VISCOSITY
C     YU_F: EFFECT OF APPLIED STRESS. HERE WE CONSIDER IT EQUAL 1.0
C           LATER, WE WILL WRITE IT AS EQUATION OF PRINCIPLE STRESS
C     YU_MEDTIME: TIME AT MID-STEP OF LOG-TIME SCALE      
      IF(TIME(2).GT.TIME_PRIME)THEN
          YU_T0 = TIME_PRIME
	    YU_MEDTIME = YU_T0 + SQRT((TIME(2)-YU_T0)*(TIME(2)+DTIME-YU_T0))
C
          YU_E = 0.0D0
	    YU_V = SQRT(1.0D0/YU_MEDTIME) + YU_Q3/YU_Q2
          YU_F = 1.0D0
C         CALCULATE EFFECTIVE MODULUS     
          DO I = 1, NS_YU
	        YU_E = YU_E + ( 1.0D0-YU_LAM_MU(I) )*YU_Q2*YU_A_MU(I) 
          ENDDO 
          YU_Q1 = 1.0D0 / EMOD
	    YU_E  = 1.0D0 / YU_Q1
      ENDIF

     
C     ==================================================================
C     Initialisation
C     ==================================================================
      RHS=0.;AMATRX=0.;AK_UU=0.;AK_EE=0.;AK_UE=0.;AK_EU=0.;F_EE=0.
      F_UU=0.;F_UU_T1=0.;F_UU_T2=0.;F_UU_T3=0.
C     ==================================================================
C     Local coordinates and weights
C     ==================================================================
C     Three Point Gauss Quadrature
      GW(1,1) = 5./9.; GW(2,1) = 5./9.
      GW(1,2) = 8./9.; GW(2,2) = 5./9.
      GW(1,3) = 5./9.; GW(2,3) = 5./9.
      GW(1,4) = 5./9.; GW(2,4) = 8./9.
      GW(1,5) = 8./9.; GW(2,5) = 8./9.
      GW(1,6) = 5./9.; GW(2,6) = 8./9.
      GW(1,7) = 5./9.; GW(2,7) = 5./9.
      GW(1,8) = 8./9.; GW(2,8) = 5./9.
      GW(1,9) = 5./9.; GW(2,9) = 5./9.
      GP(1,1) = -SQRT(3./5.); GP(2,1) = -SQRT(3./5.)
      GP(1,2) =  0.000000 ;   GP(2,2) = -SQRT(3./5.)
      GP(1,3) =  SQRT(3./5.); GP(2,3) = -SQRT(3./5.)
      GP(1,4) = -SQRT(3./5.); GP(2,4) =  0.0000000
      GP(1,5) =  0.000000 ;   GP(2,5) =  0.0000000
      GP(1,6) =  SQRT(3./5.); GP(2,6) =  0.0000000
      GP(1,7) = -SQRT(3./5.); GP(2,7) =  SQRT(3./5.)
      GP(1,8) =  0.000000 ;   GP(2,8) =  SQRT(3./5.)
      GP(1,9) =  SQRT(3./5.); GP(2,9) =  SQRT(3./5.)
C
      DO INPT=1,NGP!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>GP LOOP BEGIN
C     ==================================================================
C     Shape Functions and derivatives
C     ==================================================================
C     Local coordinates of the integration point
      XI(1) = GP(1,INPT)
      XI(2) = GP(2,INPT)
C     Shape functions and local derivatives
      CALL SHAPEFUN(AN,dNdxi,AH,dHdxi,XI)
C     Jacobian
      AJACOB = 0.0
      DO I = 1,2
          DO J = 1,2
              DO K = 1,NNODE
                  AJACOB(I,J) = AJACOB(I,J) + COORDS(I,K)*dNdxi(K,J)
              END DO
          END DO
      END DO
C        
      DTM = 0.0
      DTM = AJACOB(1,1)*AJACOB(2,2)-AJACOB(1,2)*AJACOB(2,1)
      IF (DTM.LT.0.0) THEN
          WRITE(7,*) 'Negative Jacobian',DTM
          CALL XIT	
      END IF
C     Inverse of Jacobian
      AJABOBINV(1,1)=AJACOB(2,2)/DTM
      AJABOBINV(1,2)=-AJACOB(1,2)/DTM
      AJABOBINV(2,1)=-AJACOB(2,1)/DTM
      AJABOBINV(2,2)=AJACOB(1,1)/DTM       
C     Derivatives of shape functions Q8
      dNdx = 0.0
      DO K = 1,NNODE
          DO I = 1,2
              DO J = 1,2
                  dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*AJABOBINV(J,I)
              END DO
          END DO
      END DO
C     Derivatives of shape functions Q4
      dHdx = 0.0
      DO K = 1,4
          DO I = 1,2
              DO J = 1,2
                  dHdx(K,I) = dHdx(K,I) + dHdxi(K,J)*AJABOBINV(J,I)
              END DO
          END DO
      END DO
C
C     Calculating B matrix for disp and e_eqv
      BU = 0.0
      BU(1,1:2*NNODE:2) = dNdx(1:NNODE,1)
      BU(2,2:2*NNODE:2) = dNdx(1:NNODE,2)
      BU(3,1:2*NNODE:2) = dNdx(1:NNODE,2)
      BU(3,2:2*NNODE:2) = dNdx(1:NNODE,1)
C
      BE = 0.0
      BE(1,1:4)=dHdx(1:4,1)
      BE(2,1:4)=dHdx(1:4,2)
C
C     ==================================================================
C     Updating SVARS and USRVAR
C     ==================================================================
C     Transfering maximum of last LOAD increment K_GP->K0
	IF (STEPITER .LT. 1D-8) THEN
         AK0(1,INPT) = USRVAR(JELEM,6,INPT)
      ELSE
         AK0(1,INPT) = USRVAR(JELEM,7,INPT)
      ENDIF
C     Strains{1ST-3RD SDV}.............................................
      E_TOT(1,INPT)     = DOT_PRODUCT(BU(1,:),DISP(:,1))
      E_TOT(2,INPT)     = DOT_PRODUCT(BU(2,:),DISP(:,1))
      E_TOT(3,INPT)     = DOT_PRODUCT(BU(3,:),DISP(:,1))
C     Strain increments{1ST-3RD SDV} 
      E_TOT_INC(1,INPT) = DOT_PRODUCT(BU(1,:),DISP_INC(:,1))
      E_TOT_INC(2,INPT) = DOT_PRODUCT(BU(2,:),DISP_INC(:,1))
      E_TOT_INC(3,INPT) = DOT_PRODUCT(BU(3,:),DISP_INC(:,1))
C      
      IF(TIME(2).LE.TIME_PRIME) THEN
         E(1,INPT) = E_TOT(1,INPT)
         E(2,INPT) = E_TOT(2,INPT)
         E(3,INPT) = E_TOT(3,INPT)
      ELSE
C        EFFECTIVE STRESS AT PREVIOUS STEP          
         YU_S(1,INPT) = SVARS(NSDV*(INPT-1)+11)
         YU_S(2,INPT) = SVARS(NSDV*(INPT-1)+12)
         YU_S(3,INPT) = SVARS(NSDV*(INPT-1)+13)
C
C        ==================================================================
C        Calculating for shrinkage 
C        ==================================================================      
C        TIME HISTORY
C        AT THE END OF CURRENT STEP
         TIME_DRY = TIME(2) - TIME_ACT
C     
C        DRYING     
         IF(TIME_DRY .GT. 0.0D0)THEN
		    T_CURRENT  = TIME_DRY
		    T_PREVIOUS = T_CURRENT - DTIME
	        IF(T_PREVIOUS .LT. 0.0D0)THEN
	          T_PREVIOUS =  0.0D0
	       ENDIF
         ENDIF      
C
         IF(T_CURRENT.GT.0.0d0) THEN
C             ADJUST ULTIMATE STRAIN
              SHU_ULTIMATE = -920.0D-6
              T_F  = 40
C             SHRINKAGE INCREMENT     
          DELTA_SHRINKAGE = SHU_ULTIMATE*(T_CURRENT/(T_CURRENT + T_F)
     +	                   - T_PREVIOUS/(T_PREVIOUS + T_F))	      
          ENDIF
          

         E_SHRINK(1,INPT) = 1.5d0 * DELTA_SHRINKAGE / 1.2D0
         E_SHRINK(2,INPT) = 1.5d0 * DELTA_SHRINKAGE / 1.2D0
         E_SHRINK(3,INPT) = 0.0d0
         
C         E_SHRINK(1,INPT) = 0.0d0
C         E_SHRINK(2,INPT) = 0.0d0
C
C     ==================================================================
C     Calculating for BASIC CREEP (Q2, Q3) 
C     ==================================================================      
         
C        GAMMA CAUSED BY CREEP
         K1 = 1
         DO I = 1, NS_YU
	       DO J = 1, 3
                 YU_IN(I,J,INPT) = SVARS(NSDV*(INPT-1)+20+K1)
                 K1 = K1 + 1
             ENDDO
         ENDDO
C         
C        CALCULATE THE INELASTIC STRAIN INCREMENT
         DO I = 1,3
	       YU_STRAIN(I,INPT) = 0.0D0
	   ENDDO
	   DO I = 1, NS_YU
	        DO J = 1, 3
                  YU_STRAIN(J,INPT) = YU_STRAIN(J,INPT) 
     +		    + YU_IN(I,J,INPT) * (1.0D0 - EXP( 0.0D0-YU_DELTA_Y_MU(I) ))
	        ENDDO
         ENDDO 
C        CONSIDER COMBINED EFFECT OF Q2 AND Q3
         DO I = 1,3
	       YU_STRAIN(I,INPT) = YU_STRAIN(I,INPT) * YU_V
         ENDDO


C
C     ==================================================================
C     Calculating for DRYING CREEP (Q4) 
C     ==================================================================       
         DO I = 1, 3
	        E_DRY(I,INPT) = YU_Q4*DTIME/YU_MEDTIME *YU_S(I,INPT)
C             E_DRY(I,INPT) = 0.0d0
         ENDDO
C
C     ==================================================================
C     Calculating for total CREEP and shrinkage strains
C     ==================================================================                 
C        TOTAL MPS INELASTIC STRAIN BY  CREEP
         E_MPS(1,INPT) = SVARS(NSDV*(INPT-1)+17)
         E_MPS(2,INPT) = SVARS(NSDV*(INPT-1)+18)
         E_MPS(3,INPT) = SVARS(NSDV*(INPT-1)+19)
C
         E_MPS(1,INPT) = E_MPS(1,INPT) + YU_STRAIN(1,INPT) +  
     1                   E_SHRINK(1,INPT) +  
     1                   E_DRY(1,INPT)
C         
         E_MPS(2,INPT) = E_MPS(2,INPT) + YU_STRAIN(2,INPT) + 
     1                   E_SHRINK(2,INPT) +  
     1                   E_DRY(2,INPT)
C
         E_MPS(3,INPT) = E_MPS(3,INPT) + YU_STRAIN(3,INPT) + 
     1                   E_SHRINK(3,INPT) +  
     1                   E_DRY(3,INPT) 

C
         SVARS(NSDV*(INPT-1)+17)  = E_MPS(1,INPT)
         SVARS(NSDV*(INPT-1)+18)  = E_MPS(2,INPT)
         SVARS(NSDV*(INPT-1)+19)  = E_MPS(3,INPT)
         USRVAR(JELEM,17,INPT)    = E_MPS(1,INPT)
         USRVAR(JELEM,18,INPT)    = E_MPS(2,INPT)
         USRVAR(JELEM,19,INPT)    = E_MPS(3,INPT)
C
C     ==================================================================
C     update internal variables
C     ==================================================================               
C        ELASTIC STRAIN INCREMENT
C        DELTA EPSILON = DELTA TOTAL EPSILON - DELTA EPSILON''
	   DO I = 1, 3
	       YU_STRAIN(I,INPT) = E_TOT_INC(I,INPT) - YU_STRAIN(I,INPT) - 
     1                           E_SHRINK(I,INPT) -  E_DRY(I,INPT)
	   ENDDO
         YU_DELTA_STRESS = 0.0D0
C
C        CALCULATE THE STRERSS INCREMENT
         DO K1=1,3
             DO K2=1,3
                 YU_DELTA_STRESS(K1,INPT) = YU_DELTA_STRESS(K1,INPT) + 
     +		                    DMAT(K1,K2)*YU_STRAIN(K2,INPT)
             ENDDO
         ENDDO      
      
C        UPDATE GAMA_I
         YU_EN = 0.0D0
C         YU_Q2 = 1.0D0/EMOD
	   DO I = 1, NS_YU
	        DO J = 1, 3
              YU_EN(I,J,INPT) = YU_IN(I,J,INPT) * 
     +                          EXP( 0.0D0-YU_DELTA_Y_MU(I) )
     +         + YU_LAM_MU(I)*YU_DELTA_STRESS(J,INPT)*YU_Q2*YU_A_MU(I)
              ENDDO
	   ENDDO		            
          
          
C         UPDATE SOLUTION DEPENDENT VARIABLE 
C         (INITAL INELASTIC STRAIN IN EACH KELVIN CHAIN)
          K1 = 1
          DO I = 1,NS_YU
	        DO J = 1,3
	            SVARS(NSDV*(INPT-1)+20+K1) = YU_EN(I,J,INPT)
                  USRVAR(JELEM,20+K1,INPT)   = YU_EN(I,J,INPT) 
                  K1 = K1 + 1
              ENDDO
          ENDDO
C  
          E(1,INPT) = E_TOT(1,INPT) - E_MPS(1,INPT)
          E(2,INPT) = E_TOT(2,INPT) - E_MPS(2,INPT)
          E(3,INPT) = E_TOT(3,INPT) - E_MPS(3,INPT)
C          
       ENDIF
C      
       SVARS(NSDV*(INPT-1)+1) = E(1,INPT)
       SVARS(NSDV*(INPT-1)+2) = E(2,INPT)
       SVARS(NSDV*(INPT-1)+3) = E(3,INPT)
       USRVAR(JELEM,1,INPT)   = E(1,INPT)
       USRVAR(JELEM,2,INPT)   = E(2,INPT)
       USRVAR(JELEM,3,INPT)   = E(3,INPT)
C       
C     ALOC_E_EQV{4TH SDV}..............................................
      ALOC_E_EQV(1,INPT) = c1*(E(1,INPT) + E(2,INPT)) + (1./(2.*AK))*
     &      SQRT((c2*(E(1,INPT)+E(2,INPT))**2.)+(c3*(E(1,INPT)**2. + 
     &      E(2,INPT)**2 - E(1,INPT)*E(2,INPT) + 3.*(E(3,INPT)**2.))))
	SVARS(NSDV*(INPT-1)+4) = ALOC_E_EQV(1,INPT)
	USRVAR(JELEM,4,INPT)   = ALOC_E_EQV(1,INPT)
C 
C     GP_E_EQV(5TH SDV)................................................
      GP_E_EQV(1,INPT) = DOT_PRODUCT(AH(1,:),E_EQV(:,1))
      SVARS(NSDV*(INPT-1)+5) = GP_E_EQV(1,INPT)
      USRVAR(JELEM,5,INPT) = GP_E_EQV(1,INPT)
C 
C     AK_GP(6TH SDV)...................................................
	IF (GP_E_EQV(1,INPT) .LT. AK0(1,INPT)) THEN
		AK_GP(1,INPT) = AK0(1,INPT)
	ELSE
          AK_GP(1,INPT) = GP_E_EQV(1,INPT)
      END IF
	SVARS(NSDV*(INPT-1)+6) = AK_GP(1,INPT)
      USRVAR(JELEM,6,INPT) = AK_GP(1,INPT)
C 
C     AK0(7th SDV).....................................................
	SVARS(NSDV*(INPT-1)+7) = AK0(1,INPT)!AK0_OUT
      USRVAR(JELEM,7,INPT)   = AK0(1,INPT)!AK0_OUT
C     DAMAGE(8TH SDV)..................................................	
      IF (AK_GP(1,INPT) .LE. AKI) THEN
		D(1,INPT) = 0.
	ELSE
          D(1,INPT) = 1.-(AKI/AK_GP(1,INPT))*(1.-alpha+alpha*
     +                                  DEXP(-beta*(AK_GP(1,INPT)-AKI)))
      END IF
	SVARS(NSDV*(INPT-1)+8) = D(1,INPT)
      USRVAR(JELEM,8,INPT) = D(1,INPT)
C      
C     SOTRE THE STRESS
      UEL_STRESS_EFF(1,INPT) =    (DMAT(1,1) *E(1,INPT) +
     +                   DMAT(1,2) *E(2,INPT)  +  DMAT(1,3) *E(3,INPT) )
      UEL_STRESS_EFF(2,INPT) =    (DMAT(2,1) *E(1,INPT) +
     +                   DMAT(2,2) *E(2,INPT)  +  DMAT(2,3) *E(3,INPT) )
      UEL_STRESS_EFF(3,INPT) =    (DMAT(3,1) *E(1,INPT) +
     +                   DMAT(3,2) *E(2,INPT)  +  DMAT(3,3) *E(3,INPT) )
C 
      UEL_STRESS(1,INPT) =  (1.0d0 - D(1,INPT)) * UEL_STRESS_EFF(1,INPT)
      UEL_STRESS(2,INPT) =  (1.0d0 - D(1,INPT)) * UEL_STRESS_EFF(2,INPT)
      UEL_STRESS(3,INPT) =  (1.0d0 - D(1,INPT)) * UEL_STRESS_EFF(3,INPT)
C
      SVARS(NSDV*(INPT-1)+11)   = UEL_STRESS_EFF(1,INPT)
      SVARS(NSDV*(INPT-1)+12)   = UEL_STRESS_EFF(2,INPT)
      SVARS(NSDV*(INPT-1)+13)   = UEL_STRESS_EFF(3,INPT)
      USRVAR(JELEM,11,INPT)     = UEL_STRESS_EFF(1,INPT) 
      USRVAR(JELEM,12,INPT)     = UEL_STRESS_EFF(2,INPT) 
      USRVAR(JELEM,13,INPT)     = UEL_STRESS_EFF(3,INPT)
C
      SVARS(NSDV*(INPT-1)+14)   = UEL_STRESS(1,INPT)
      SVARS(NSDV*(INPT-1)+15)   = UEL_STRESS(2,INPT)
      SVARS(NSDV*(INPT-1)+16)   = UEL_STRESS(3,INPT)
      USRVAR(JELEM,14,INPT)     = UEL_STRESS(1,INPT) 
      USRVAR(JELEM,15,INPT)     = UEL_STRESS(2,INPT) 
      USRVAR(JELEM,16,INPT)     = UEL_STRESS(3,INPT) 
      
      USRVAR(JELEM,117,INPT)     = YU_DELTA_STRESS(1,INPT) 
      USRVAR(JELEM,118,INPT)     = YU_DELTA_STRESS(2,INPT) 
      USRVAR(JELEM,119,INPT)     = YU_DELTA_STRESS(3,INPT) 
C   
C
C     ==================================================================
C     Calculating element stiffness matrix
C     ==================================================================
      T_INT = GW(1,INPT)*GW(2,INPT)*DTM
C     K_UU..............................................................
C      DO K=1,2*NNODE
C          DO L=1,2*NNODE
C              DO I=1,3
C                  DO J=1,3
C                      AK_UU(K,L)=AK_UU(K,L)+
C     &                    BU(I,K)*DMAT(I,J)*BU(J,L)*(1.-D(1,INPT))*T_INT
C                  END DO
C              END DO
C          END DO
C      END DO
      
      AK_UU_T1 = MATMUL(DMAT, BU)
      AK_UU_T2 = MATMUL(TRANSPOSE(BU), AK_UU_T1)
C
      DO K=1,2*NNODE
          DO L=1,2*NNODE
              AK_UU(K,L) = AK_UU(K,L) + AK_UU_T2(K,L) *
     +            (1.-D(1,INPT))*T_INT          
          ENDDO
      ENDDO
C     F_UU..............................................................	   
      F_UU_T1 = MATMUL(TRANSPOSE(BU), DMAT)
      E_RHS(1,1) = E_TOT(1,INPT) - E_MPS(1,INPT)
      E_RHS(2,1) = E_TOT(2,INPT) - E_MPS(2,INPT)
      E_RHS(3,1) = E_TOT(3,INPT) - E_MPS(3,INPT)
      F_UU_T2 = MATMUL(F_UU_T1, E_RHS)
      
C      F_UU_T3 = MATMUL(F_UU_T2 , DISP)
      DO K=1,2*NNODE
          F_UU(K,1) = F_UU(K,1) +  F_UU_T2(K,1)* 
     +      (1.-D(1,INPT))*T_INT
      END DO      
       
      
      
C     K_EE..............................................................
      DO K=1,4
          DO L=1,4
              AK_EE(K,L) = AK_EE(K,L) + AH(1,K)*AH(1,L)*T_INT
              DO I=1,2
                  AK_EE(K,L) = AK_EE(K,L) + BE(I,K)*BE(I,L)*C*T_INT
              END DO
          END DO
      END DO
C     K_UE..............................................................
      !Derivative of damage wrt history parameter 
	IF (AK_GP(1,INPT) .LT. AKI) THEN
	    dD_dK = 0.
      ELSEIF (AK_GP(1,INPT) .GT. AK0(1,INPT)) THEN
          dD_dk = ((AKI/(AK_GP(1,INPT)**2.))*
     +        (1.-alpha+alpha*DEXP(-beta*(AK_GP(1,INPT)-AKI)))) +
     +        ((alpha*beta*AKI*DEXP(-beta*(AK_GP(1,INPT)-AKI)))/
     +        AK_GP(1,INPT))
       ELSE
          dD_dk = 0.
       END IF
C
       DO K = 1, 2*NNODE
          DO L = 1, 4
              DO I = 1, 3
                  DO J = 1, 3
                      AK_UE(K,L) = AK_UE(K,L) - 
     &                  BU(I,K)*DMAT(I,J)*E(J,INPT)*AH(1,L)*dD_dK*T_INT
                  END DO
              END DO
          END DO
       END DO
C      K_EU..............................................................
       Denominator = c2*(E(1,INPT) + E(2,INPT))**2. + 
     &     (c3*(E(1,INPT)**2. + E(2,INPT)**2. - 
     &     E(1,INPT)*E(2,INPT)+3.*(E(3,INPT)**2.)))
      IF (Denominator .LT. 1D-10) THEN
		DEEQV_DE = 0.0
	ELSE
C      Computing derivatives of principal strains
          DEEQV_DE(1) = c1 + (c3*(2.*E(1,INPT) - E(2,INPT))
     &     + 2.*c2*(E(1,INPT) + E(2,INPT))) /(4.*AK*SQRT(Denominator))
          DEEQV_DE(2) = c1 + (c3*(2.*E(2,INPT) - E(1,INPT))
     &     + 2.*c2*(E(1,INPT) + E(2,INPT))) /(4.*AK*SQRT(Denominator))
          DEEQV_DE(3) = (3.*c3*E(3,INPT)) /(2.*AK*SQRT(Denominator))
      END IF
C
      DO K=1,4
          DO L=1,2*NNODE
              DO J=1,3
                  AK_EU(K,L) = AK_EU(K,L) - 
     +                     AH(1,K)*DEEQV_DE(J)*BU(J,L)*T_INT
              END DO
          END DO
      END DO
C     F_EE..............................................................	   
      DO K=1,4
          F_EE(K,1) = F_EE(K,1) + AH(1,K)*ALOC_E_EQV(1,INPT)*T_INT
      END DO

      
      
       
      
      
C
      END DO!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>GP LOOP END
C
C     =================================================================
C     Assembly of K into AMATRX........................................
	     !Assemble K_UU
      DO I=1,NNODE
          DO J=1,NNODE
              AMATRX(3*I-2,3*J-2) = AK_UU(2*I-1,2*J-1)
              AMATRX(3*I-1,3*J-2) = AK_UU(2*I,2*J-1)
              AMATRX(3*I-2,3*J-1) = AK_UU(2*I-1,2*J)
              AMATRX(3*I-1,3*J-1) = AK_UU(2*I,2*J)
		END DO
	END DO
C	Assemble K_UE
	DO I=1,NNODE
		DO J=1,4
              AMATRX(3*I-2,3*J) = AK_UE(2*I-1,J)
              AMATRX(3*I-1,3*J) = AK_UE(2*I,J)
		END DO
	END DO
C	Assemble K_EU
	DO I=1,4
		DO J=1,NNODE
              AMATRX(3*I,3*J-2) = AK_EU(I,2*J-1)
              AMATRX(3*I,3*J-1) = AK_EU(I,2*J)
		END DO
	END DO
C     Assemble K_EE
	DO I=1,4
		DO J=1,4
              AMATRX(3*I,3*J) = AK_EE(I,J)
          END DO
      END DO
C     Conditioning
	AMATRX(15,15) = 1.0; AMATRX(18,18) = 1.0
	AMATRX(21,21) = 1.0; AMATRX(24,24) = 1.0
C
C     ==================================================================
C     Residual vector
C     ==================================================================
C     F_U
C      DO I=1,NNODE
C         DO J=1,2*NNODE
C              RHS(3*I-2,1) = RHS(3*I-2,1)-AK_UU(2*I-1,J)*(DISP(J,1))
C              RHS(3*I-1,1) = RHS(3*I-1,1)-AK_UU(2*I,J)  *(DISP(J,1))              
C         END DO
C	END DO
      DO I=1,NNODE
          RHS(3*I-2,1) = RHS(3*I-2,1) - F_UU(2*I-1,1)
          RHS(3*I-1,1) = RHS(3*I-1,1) - F_UU(2*I,1)
      ENDDO
C	F_E
      DO I=1,4
          DO J=1,4
              RHS(3*I,1) = RHS(3*I,1)-AK_EE(I,J)*E_EQV(J,1)
          END DO
          RHS(3*I,1) = RHS(3*I,1) + F_EE(I,1)
      END DO
C
      RETURN
      END
C***********************************************************************
C     ==================================================================
C     SHAPE FUNCTION SUBROUTINE
C     ==================================================================
      SUBROUTINE SHAPEFUN(AN,dNdxi,AH,dHdxi,XI)
      INCLUDE 'ABA_PARAM.INC'
      Real*8 AN(1,9),dNdxi(8,2)
						Real*8 AH(1,4),dHdxi(4,2)
      Real*8 XI(2)
C
C     Values of shape functions for Q8
						AN(1,1) = -0.25*(1.-XI(1))*(1.-XI(2))*(1.+XI(1)+XI(2))
      AN(1,2) = -0.25*(1.+XI(1))*(1.-XI(2))*(1.-XI(1)+XI(2))
      AN(1,3) = -0.25*(1.+XI(1))*(1.+XI(2))*(1.-XI(1)-XI(2))
      AN(1,4) = -0.25*(1.-XI(1))*(1.+XI(2))*(1.+XI(1)-XI(2))
      AN(1,5) =  0.5*(1.-XI(1))*(1.+XI(1))*(1.-XI(2))
      AN(1,6) =  0.5*(1.+XI(1))*(1.+XI(2))*(1.-XI(2))
      AN(1,7) =  0.5*(1.-XI(1))*(1.+XI(1))*(1.+XI(2))
      AN(1,8) =  0.5*(1.-XI(1))*(1.+XI(2))*(1.-XI(2))
C
C     Values of shape functions for Q4
      AH(1,1) = 0.25*(1.-XI(1))*(1.-XI(2))
      AH(1,2) = 0.25*(1.+XI(1))*(1.-XI(2))
      AH(1,3) = 0.25*(1.+XI(1))*(1.+XI(2))
      AH(1,4) = 0.25*(1.-XI(1))*(1.+XI(2))
C
C     Derivatives of shape functions Q8
      dNdxi = 0.0
      dNdxi(1,1) = -0.25*(-1.+XI(2))*(2.*XI(1)+XI(2))
      dNdxi(1,2) = -0.25*(-1.+XI(1))*(XI(1)+2.*XI(2))
      dNdxi(2,1) =  0.25*(-1.+XI(2))*(XI(2)-2.*XI(1))
      dNdxi(2,2) =  0.25*(1.+XI(1))*(2.*XI(2)-XI(1))
      dNdxi(3,1) =  0.25*(1.+XI(2))*(2.*XI(1)+XI(2))
      dNdxi(3,2) =  0.25*(1.+XI(1))*(XI(1)+2.*XI(2))
      dNdxi(4,1) = -0.25*(1.+XI(2))*(XI(2)-2.*XI(1))
      dNdxi(4,2) = -0.25*(-1.+XI(1))*(2.*XI(2)-XI(1))
						dNdxi(5,1) =  XI(1)*(-1.+XI(2))  
      dNdxi(5,2) =  0.5*(1.+XI(1))*(-1.+XI(1))
      dNdxi(6,1) = -0.5*(1.+XI(2))*(-1.+XI(2)) 
      dNdxi(6,2) = -XI(2)*(1.+XI(1))
      dNdxi(7,1) = -XI(1)*(1.+XI(2)) 
      dNdxi(7,2) = -0.5*(1.+XI(1))*(-1.+XI(1))
      dNdxi(8,1) =  0.5*(1.+XI(2))*(-1.+XI(2)) 
      dNdxi(8,2) =  XI(2)*(-1.+XI(1))
C
C     Derivatives of shape functions Q4
      dHdxi = 0.0
      dHdxi(1,1) = -0.25*(1.-XI(2))
      dHdxi(1,2) = -0.25*(1.-XI(1))
      dHdxi(2,1) =  0.25*(1.-XI(2))
      dHdxi(2,2) = -0.25*(1.+XI(1))
      dHdxi(3,1) =  0.25*(1.+XI(2))
      dHdxi(3,2) =  0.25*(1.+XI(1))
      dHdxi(4,1) = -0.25*(1.+XI(2))
      dHdxi(4,2) =  0.25*(1.-XI(1))
      RETURN
      END
C***********************************************************************      
C Subroutine UMAT  : 
C Dummy material
C
C ==============================================================
C !!! NOTE: N_ELEM has to be changed according to the UEL !!!!!
C ==============================================================
C
       SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
       CHARACTER*80 CMNAME
       DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS),
     3 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C 
       PARAMETER (ONE=1.0,TWO=2.0,THREE=3.0,SIX=6.0, HALF =0.5,
     1 N_ELEM=10000,NSDV=120)
       DATA NEWTON,TOLER/40,1.D-2/ 
C       
       COMMON/KUSER/USRVAR(N_ELEM,NSDV,9)
C 
C ----------------------------------------------------------- 
C          Material properties
C ----------------------------------------------------------- 
C          PROPS(1) - Young's modulus 
C          PROPS(2) - Poisson ratio 
C ----------------------------------------------------------- 
C
C	Elastic properties
C
       EMOD=PROPS(1)
       ENU=PROPS(2)
       EG=EMOD/(TWO*(ONE+ENU))
       EG2=EG*TWO
       ELAM=EG2*ENU/(ONE-TWO*ENU)
C
C	Stiffness tensor
C
       DO K1=1, NTENS
        DO K2=1, NTENS
         DDSDDE(K2, K1)=0.0
        END DO
       END DO
C
       DO K1=1, NDI
        DO K2=1, NDI
         DDSDDE(K2, K1)=ELAM
        END DO
        DDSDDE(K1, K1)=EG2+ELAM
       END DO 
C
       DO K1=NDI+1, NTENS
        DDSDDE(K1, K1)=EG
       END DO
C
C	Calculate Stresses
C
       DO K1=1, NTENS
        DO K2=1, NTENS
         STRESS(K2)=STRESS(K2)+DDSDDE(K2, K1)*DSTRAN(K1)
        END DO
       END DO 
C
        NELEMAN=NOEL-N_ELEM
C       
       DO I=1,NSTATV
        STATEV(I)=USRVAR(NELEMAN,I,NPT)
       END DO
C       
       RETURN
       END      