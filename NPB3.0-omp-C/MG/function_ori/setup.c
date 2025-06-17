static void setup(int *n1, int *n2, int *n3, int lt) {

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

    int k;

    for ( k = lt-1; k >= 1; k--) {
	nx[k] = nx[k+1]/2;
	ny[k] = ny[k+1]/2;
	nz[k] = nz[k+1]/2;
    }

    for (k = 1; k <= lt; k++) {
	m1[k] = nx[k]+2;
	m2[k] = nz[k]+2;
	m3[k] = ny[k]+2;
    }

    is1 = 1;
    ie1 = nx[lt];
    *n1 = nx[lt]+2;
    is2 = 1;
    ie2 = ny[lt];
    *n2 = ny[lt]+2;
    is3 = 1;
    ie3 = nz[lt];
    *n3 = nz[lt]+2;

    if (debug_vec[1] >=  1 ) {
	printf(" in setup, \n");
	printf("  lt  nx  ny  nz  n1  n2  n3 is1 is2 is3 ie1 ie2 ie3\n");
	printf("%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d\n",
	       lt,nx[lt],ny[lt],nz[lt],*n1,*n2,*n3,is1,is2,is3,ie1,ie2,ie3);
    }
}