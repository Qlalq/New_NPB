static void verify(int no_time_steps, char *class, boolean *verified) {
  double xcrref[5],xceref[5],xcrdif[5],xcedif[5], 
    epsilon, xce[5], xcr[5], dtref;
  int m;
  epsilon = 1.0e-08;
  error_norm(xce);
  compute_rhs();
  rhs_norm(xcr);
  for (m = 0; m < 5; m++) {
    xcr[m] = xcr[m] / dt;
  }
  *class = 'U';
  *verified = TRUE;
  for (m = 0; m < 5; m++) {
    xcrref[m] = 1.0;
    xceref[m] = 1.0;
  }
  if (grid_points[0] == 12 &&
      grid_points[1] == 12 &&
      grid_points[2] == 12 &&
      no_time_steps == 60) {
    *class = 'S';
    dtref = 1.0e-2;
    xcrref[0] = 1.7034283709541311e-01;
    xcrref[1] = 1.2975252070034097e-02;
    xcrref[2] = 3.2527926989486055e-02;
    xcrref[3] = 2.6436421275166801e-02;
    xcrref[4] = 1.9211784131744430e-01;
    xceref[0] = 4.9976913345811579e-04;
    xceref[1] = 4.5195666782961927e-05;
    xceref[2] = 7.3973765172921357e-05;
    xceref[3] = 7.3821238632439731e-05;
    xceref[4] = 8.9269630987491446e-04;
    } else if (grid_points[0] == 24 &&
	       grid_points[1] == 24 &&
	       grid_points[2] == 24 &&
	       no_time_steps == 200) {
      *class = 'W';
      dtref = 0.8e-3;
      xcrref[0] = 0.1125590409344e+03;
      xcrref[1] = 0.1180007595731e+02;
      xcrref[2] = 0.2710329767846e+02;
      xcrref[3] = 0.2469174937669e+02;
      xcrref[4] = 0.2638427874317e+03;
      xceref[0] = 0.4419655736008e+01;
      xceref[1] = 0.4638531260002e+00;
      xceref[2] = 0.1011551749967e+01;
      xceref[3] = 0.9235878729944e+00;
      xceref[4] = 0.1018045837718e+02;
    } else if (grid_points[0] == 64 &&
	       grid_points[1] == 64 &&
	       grid_points[2] == 64 &&
	       no_time_steps == 200) {
      *class = 'A';
      dtref = 0.8e-3;
      xcrref[0] = 1.0806346714637264e+02;
      xcrref[1] = 1.1319730901220813e+01;
      xcrref[2] = 2.5974354511582465e+01;
      xcrref[3] = 2.3665622544678910e+01;
      xcrref[4] = 2.5278963211748344e+02;
      xceref[0] = 4.2348416040525025e+00;
      xceref[1] = 4.4390282496995698e-01;
      xceref[2] = 9.6692480136345650e-01;
      xceref[3] = 8.8302063039765474e-01;
      xceref[4] = 9.7379901770829278e+00;
    } else if (grid_points[0] == 102 &&
	       grid_points[1] == 102 &&
	       grid_points[2] == 102 &&
	       no_time_steps == 200) {
      *class = 'B';
      dtref = 3.0e-4;
      xcrref[0] = 1.4233597229287254e+03;
      xcrref[1] = 9.9330522590150238e+01;
      xcrref[2] = 3.5646025644535285e+02;
      xcrref[3] = 3.2485447959084092e+02;
      xcrref[4] = 3.2707541254659363e+03;
      xceref[0] = 5.2969847140936856e+01;
      xceref[1] = 4.4632896115670668e+00;
      xceref[2] = 1.3122573342210174e+01;
      xceref[3] = 1.2006925323559144e+01;
      xceref[4] = 1.2459576151035986e+02;
    } else if (grid_points[0] == 162 &&
	       grid_points[1] == 162 &&
	       grid_points[2] == 162 &&
	       no_time_steps == 200) {
      *class = 'C';
      dtref = 1.0e-4;
      xcrref[0] = 0.62398116551764615e+04;
      xcrref[1] = 0.50793239190423964e+03;
      xcrref[2] = 0.15423530093013596e+04;
      xcrref[3] = 0.13302387929291190e+04;
      xcrref[4] = 0.11604087428436455e+05;
      xceref[0] = 0.16462008369091265e+03;
      xceref[1] = 0.11497107903824313e+02;
      xceref[2] = 0.41207446207461508e+02;
      xceref[3] = 0.37087651059694167e+02;
      xceref[4] = 0.36211053051841265e+03;
    } else {
      *verified = FALSE;
    }
 #pragma omp parallel for default(shared) private(m)
  for (m = 0; m < 5; m++) {
    xcrdif[m] = fabs((xcr[m]-xcrref[m])/xcrref[m]);
    xcedif[m] = fabs((xce[m]-xceref[m])/xceref[m]);
  }
  if (*class != 'U') {
    printf(" Verification being performed for class %1c\n", *class);
    printf(" accuracy setting for epsilon = %20.13e\n", epsilon);
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
    if (*class == 'U') {
      printf("          %2d%20.13e\n", m, xcr[m]);
    } else if (xcrdif[m] > epsilon) {
      *verified = FALSE;
      printf(" FAILURE: %2d%20.13e%20.13e%20.13e\n",
	     m, xcr[m], xcrref[m], xcrdif[m]);
    } else {
      printf("          %2d%20.13e%20.13e%20.13e\n",
	     m, xcr[m], xcrref[m], xcrdif[m]);
    }
  }
  if (*class != 'U') {
    printf(" Comparison of RMS-norms of solution error\n");
  } else {
    printf(" RMS-norms of solution error\n");
  }
  for (m = 0; m < 5; m++) {
    if (*class == 'U') {
      printf("          %2d%20.13e\n", m, xce[m]);
    } else if (xcedif[m] > epsilon) {
      *verified = FALSE;
      printf(" FAILURE: %2d%20.13e%20.13e%20.13e\n",
	     m, xce[m], xceref[m], xcedif[m]);
    } else {
      printf("          %2d%20.13e%20.13e%20.13e\n",
	     m, xce[m], xceref[m], xcedif[m]);
    }
  }
  if (*class == 'U') {
    printf(" No reference values provided\n");
    printf(" No verification performed\n");
  } else if (*verified == TRUE) {
    printf(" Verification Successful\n");
  } else {
    printf(" Verification failed\n");
  }
}