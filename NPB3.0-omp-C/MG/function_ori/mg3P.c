static void mg3P(double ****u, double ***v, double ****r, double a[4],
		 double c[4], int n1, int n2, int n3, int k) {

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c     multigrid V-cycle routine
c-------------------------------------------------------------------*/

    int j;

/*--------------------------------------------------------------------
c     down cycle.
c     restrict the residual from the find grid to the coarse
c-------------------------------------------------------------------*/

    for (k = lt; k >= lb+1; k--) {
	j = k-1;
	rprj3(r[k], m1[k], m2[k], m3[k],
	      r[j], m1[j], m2[j], m3[j], k);
    }

    k = lb;
/*--------------------------------------------------------------------
c     compute an approximate solution on the coarsest grid
c-------------------------------------------------------------------*/
    zero3(u[k], m1[k], m2[k], m3[k]);
    psinv(r[k], u[k], m1[k], m2[k], m3[k], c, k);

    for (k = lb+1; k <= lt-1; k++) {
	j = k-1;
/*--------------------------------------------------------------------
c        prolongate from level k-1  to k
c-------------------------------------------------------------------*/
	zero3(u[k], m1[k], m2[k], m3[k]);
	interp(u[j], m1[j], m2[j], m3[j],
	       u[k], m1[k], m2[k], m3[k], k);
/*--------------------------------------------------------------------
c        compute residual for level k
c-------------------------------------------------------------------*/
	resid(u[k], r[k], r[k], m1[k], m2[k], m3[k], a, k);
/*--------------------------------------------------------------------
c        apply smoother
c-------------------------------------------------------------------*/
	psinv(r[k], u[k], m1[k], m2[k], m3[k], c, k);
    }

    j = lt - 1;
    k = lt;
    interp(u[j], m1[j], m2[j], m3[j], u[lt], n1, n2, n3, k);
    resid(u[lt], v, r[lt], n1, n2, n3, a, k);
    psinv(r[lt], u[lt], n1, n2, n3, c, k);
}