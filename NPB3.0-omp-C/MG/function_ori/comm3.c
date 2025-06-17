static void comm3(double ***u, int n1, int n2, int n3, int kk) {

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c     comm3 organizes the communication on all borders 
c-------------------------------------------------------------------*/

    int i1, i2, i3;

    /* axis = 1 */
#pragma omp parallel default(shared) private(i1,i2,i3)
{
#pragma omp for
    for ( i3 = 1; i3 < n3-1; i3++) {
	for ( i2 = 1; i2 < n2-1; i2++) {
	    u[i3][i2][n1-1] = u[i3][i2][1];
	    u[i3][i2][0] = u[i3][i2][n1-2];
	}
//    }

    /* axis = 2 */
//#pragma omp for
//    for ( i3 = 1; i3 < n3-1; i3++) {
	for ( i1 = 0; i1 < n1; i1++) {
	    u[i3][n2-1][i1] = u[i3][1][i1];
	    u[i3][0][i1] = u[i3][n2-2][i1];
	}
    }

    /* axis = 3 */
#pragma omp for nowait
    for ( i2 = 0; i2 < n2; i2++) {
	for ( i1 = 0; i1 < n1; i1++) {
	    u[n3-1][i2][i1] = u[1][i2][i1];
	    u[0][i2][i1] = u[n3-2][i2][i1];
	}
    }
}//end #pragma omp parallel
}