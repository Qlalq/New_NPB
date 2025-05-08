static void exact_solution(double xi, double eta, double zeta,
			   double dtemp[5]) {
  int m;
  for (m = 0; m < 5; m++) {
    dtemp[m] =  ce[0][m] +
      xi*(ce[1][m] + xi*(ce[4][m] + 
			     xi*(ce[7][m] + xi*ce[10][m]))) +
      eta*(ce[2][m] + eta*(ce[5][m] + 
			       eta*(ce[8][m] + eta*ce[11][m])))+
      zeta*(ce[3][m] + zeta*(ce[6][m] +
				 zeta*(ce[9][m] + 
				       zeta*ce[12][m])));
  }
}