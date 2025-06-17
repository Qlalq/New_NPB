static void verify(double xcr[5], double xce[5], double xci,
		   char *class, boolean *verified) {

/*--------------------------------------------------------------------
c  verification routine                         
--------------------------------------------------------------------*/

  double xcrref[5],xceref[5],xciref, 
    xcrdif[5],xcedif[5],xcidif,
    epsilon, dtref;
  int m;

/*--------------------------------------------------------------------
c   tolerance level
--------------------------------------------------------------------*/
  epsilon = 1.0e-08;

  *class = 'U';
  *verified = TRUE;

  for (m = 0; m < 5; m++) {
    xcrref[m] = 1.0;
    xceref[m] = 1.0;
  }
  xciref = 1.0;

  if ( nx0 == 12 && ny0 == 12 && nz0 == 12 && itmax == 50)  {
    *class = 'S';
    dtref = 5.0e-1;

/*--------------------------------------------------------------------
c   Reference values of RMS-norms of residual, for the (12X12X12) grid,
c   after 50 time steps, with  DT = 5.0d-01
--------------------------------------------------------------------*/
    xcrref[0] = 1.6196343210976702e-02;
    xcrref[1] = 2.1976745164821318e-03;
    xcrref[2] = 1.5179927653399185e-03;
    xcrref[3] = 1.5029584435994323e-03;
    xcrref[4] = 3.4264073155896461e-02;

/*--------------------------------------------------------------------
c   Reference values of RMS-norms of solution error, for the (12X12X12) grid,
c   after 50 time steps, with  DT = 5.0d-01
--------------------------------------------------------------------*/
    xceref[0] = 6.4223319957960924e-04;
    xceref[1] = 8.4144342047347926e-05;
    xceref[2] = 5.8588269616485186e-05;
    xceref[3] = 5.8474222595157350e-05;
    xceref[4] = 1.3103347914111294e-03;

/*--------------------------------------------------------------------
c   Reference value of surface integral, for the (12X12X12) grid,
c   after 50 time steps, with DT = 5.0d-01
--------------------------------------------------------------------*/
    xciref = 7.8418928865937083;

  } else if ( nx0 == 33 && ny0 == 33 && nz0 == 33 && itmax == 300) {

    *class = 'W';   /* SPEC95fp size */
    dtref = 1.5e-3;
/*--------------------------------------------------------------------
c   Reference values of RMS-norms of residual, for the (33x33x33) grid,
c   after 300 time steps, with  DT = 1.5d-3
--------------------------------------------------------------------*/
    xcrref[0] =   0.1236511638192e+02;
    xcrref[1] =   0.1317228477799e+01;
    xcrref[2] =   0.2550120713095e+01;
    xcrref[3] =   0.2326187750252e+01;
    xcrref[4] =   0.2826799444189e+02;

/*--------------------------------------------------------------------
c   Reference values of RMS-norms of solution error, for the (33X33X33) grid,
--------------------------------------------------------------------*/
    xceref[0] =   0.4867877144216;
    xceref[1] =   0.5064652880982e-01;
    xceref[2] =   0.9281818101960e-01;
    xceref[3] =   0.8570126542733e-01;
    xceref[4] =   0.1084277417792e+01;

/*--------------------------------------------------------------------
c   Reference value of surface integral, for the (33X33X33) grid,
c   after 300 time steps, with  DT = 1.5d-3
--------------------------------------------------------------------*/
    xciref    =   0.1161399311023e+02;

  } else if ( nx0 == 64 && ny0 == 64 && nz0 == 64 && itmax == 250) {

    *class = 'A';
    dtref = 2.0e+0;
/*--------------------------------------------------------------------
c   Reference values of RMS-norms of residual, for the (64X64X64) grid,
c   after 250 time steps, with  DT = 2.0d+0.0
--------------------------------------------------------------------*/
    xcrref[0] = 7.7902107606689367e+02;
    xcrref[1] = 6.3402765259692870e+01;
    xcrref[2] = 1.9499249727292479e+02;
    xcrref[3] = 1.7845301160418537e+02;
    xcrref[4] = 1.8384760349464247e+03;

/*--------------------------------------------------------------------
c   Reference values of RMS-norms of solution error, for the (64X64X64) grid,
c   after 250 time steps, with  DT = 2.0d+0.0
--------------------------------------------------------------------*/
    xceref[0] = 2.9964085685471943e+01;
    xceref[1] = 2.8194576365003349;
    xceref[2] = 7.3473412698774742;
    xceref[3] = 6.7139225687777051;
    xceref[4] = 7.0715315688392578e+01;

/*--------------------------------------------------------------------
c   Reference value of surface integral, for the (64X64X64) grid,
c   after 250 time steps, with DT = 2.0d+0.0
--------------------------------------------------------------------*/
    xciref = 2.6030925604886277e+01;

    } else if ( nx0 == 102 && ny0 == 102 && nz0 == 102 && itmax == 250) {

      *class = 'B';
      dtref = 2.0e+0;

/*--------------------------------------------------------------------
c   Reference values of RMS-norms of residual, for the (102X102X102) grid,
c   after 250 time steps, with  DT = 2.0d+0.0
--------------------------------------------------------------------*/
      xcrref[0] = 3.5532672969982736e+03;
      xcrref[1] = 2.6214750795310692e+02;
      xcrref[2] = 8.8333721850952190e+02;
      xcrref[3] = 7.7812774739425265e+02;
      xcrref[4] = 7.3087969592545314e+03;

/*--------------------------------------------------------------------
c   Reference values of RMS-norms of solution error, for the (102X102X102) 
c   grid, after 250 time steps, with  DT = 2.0d+0.0
--------------------------------------------------------------------*/
      xceref[0] = 1.1401176380212709e+02;
      xceref[1] = 8.1098963655421574;
      xceref[2] = 2.8480597317698308e+01;
      xceref[3] = 2.5905394567832939e+01;
      xceref[4] = 2.6054907504857413e+02;

/*--------------------------------------------------------------------
c   Reference value of surface integral, for the (102X102X102) grid,
c   after 250 time steps, with DT = 2.0d+0.0
--------------------------------------------------------------------*/
      xciref = 4.7887162703308227e+01;

      } else if ( nx0 == 162 && ny0 == 162 && nz0 == 162 && itmax == 250) {

	*class = 'C';
	dtref = 2.0e+0;

/*--------------------------------------------------------------------
c   Reference values of RMS-norms of residual, for the (162X162X162) grid,
c   after 250 time steps, with  DT = 2.0d+0.0
--------------------------------------------------------------------*/
	xcrref[0] = 1.03766980323537846e+04;
	xcrref[1] = 8.92212458801008552e+02;
	xcrref[2] = 2.56238814582660871e+03;
	xcrref[3] = 2.19194343857831427e+03;
	xcrref[4] = 1.78078057261061185e+04;

/*--------------------------------------------------------------------
c   Reference values of RMS-norms of solution error, for the (162X162X162) 
c   grid, after 250 time steps, with  DT = 2.0d+0.0
--------------------------------------------------------------------*/
	xceref[0] = 2.15986399716949279e+02;
	xceref[1] = 1.55789559239863600e+01;
	xceref[2] = 5.41318863077207766e+01;
	xceref[3] = 4.82262643154045421e+01;
	xceref[4] = 4.55902910043250358e+02;

/*--------------------------------------------------------------------
c   Reference value of surface integral, for the (162X162X162) grid,
c   after 250 time steps, with DT = 2.0d+0.0
--------------------------------------------------------------------*/
	xciref = 6.66404553572181300e+01;
      } else {
	*verified = FALSE;
      }

/*--------------------------------------------------------------------
c    verification test for residuals if gridsize is either 12X12X12 or 
c    64X64X64 or 102X102X102 or 162X162X162
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c    Compute the difference of solution values and the known reference values.
--------------------------------------------------------------------*/
  for (m = 0; m < 5; m++) {
           
    xcrdif[m] = fabs((xcr[m]-xcrref[m])/xcrref[m]);
    xcedif[m] = fabs((xce[m]-xceref[m])/xceref[m]);
           
  }
  xcidif = fabs((xci - xciref)/xciref);

/*--------------------------------------------------------------------
c    Output the comparison of computed results to known cases.
--------------------------------------------------------------------*/

  if (*class != 'U') {
    printf("\n Verification being performed for class %1c\n", *class);
    printf(" Accuracy setting for epsilon = %20.13e\n", epsilon);
    if (fabs(dt-dtref) > epsilon) {  
      *verified = FALSE;
      *class = 'U';
      printf(" DT does not match the reference value of %15.8e\n", dtref);
    }
  } else {
    printf(" Unknown class\n");
  }

  if (*class != 'U') {
    printf(" Comparison of RMS-norms of residual\n");
  } else {
    printf(" RMS-norms of residual\n");
  }

  for (m = 0; m < 5; m++) {
    if (*class  ==  'U') {
      printf("          %2d  %20.13e\n", m, xcr[m]);
    } else if (xcrdif[m] > epsilon) {
      *verified = FALSE;
      printf(" FAILURE: %2d  %20.13e%20.13e%20.13e\n",
	     m,xcr[m],xcrref[m],xcrdif[m]);
    } else {
      printf("          %2d  %20.13e%20.13e%20.13e\n",
	     m,xcr[m],xcrref[m],xcrdif[m]);
    }
  }

  if (*class != 'U') {
    printf(" Comparison of RMS-norms of solution error\n");
  } else {
    printf(" RMS-norms of solution error\n");
  }
        
  for (m = 0; m < 5; m++) {
    if (*class  ==  'U') {
      printf("          %2d  %20.13e\n", m, xce[m]);
    } else if (xcedif[m] > epsilon) {
      *verified = FALSE;
      printf(" FAILURE: %2d  %20.13e%20.13e%20.13e\n",
	     m,xce[m],xceref[m],xcedif[m]);
    } else {
      printf("          %2d  %20.13e%20.13e%20.13e\n",
	     m,xce[m],xceref[m],xcedif[m]);
    }
  }
        
  if (*class != 'U') {
    printf(" Comparison of surface integral\n");
  } else {
    printf(" Surface integral\n");
  }

  if (*class  ==  'U') {
    printf("              %20.13e\n", xci);
  } else if (xcidif > epsilon) {
    *verified = FALSE;
    printf(" FAILURE:     %20.13e%20.13e%20.13e\n", 
	   xci, xciref, xcidif);
  } else {
    printf("              %20.13e%20.13e%20.13e\n",
	   xci, xciref, xcidif);
  }

  if (*class  ==  'U') {
    printf(" No reference values provided\n");
    printf(" No verification performed\n");
  } else if (*verified) {
    printf(" Verification Successful\n");
  } else {
    printf(" Verification failed\n");
  }
}