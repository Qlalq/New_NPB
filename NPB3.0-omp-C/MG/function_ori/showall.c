static void showall(double ***z, int n1, int n2, int n3) {

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

    int i1,i2,i3;
    int m1, m2, m3;

    m1 = min(n1,18);
    m2 = min(n2,14);
    m3 = min(n3,18);

    printf("\n");
    for (i3 = 0; i3 < m3; i3++) {
	for (i1 = 0; i1 < m1; i1++) {
	    for (i2 = 0; i2 < m2; i2++) {
		printf("%6.3f", z[i3][i2][i1]);
	    }
	    printf("\n");
	}
	printf(" - - - - - - - \n");
    }
    printf("\n");
}