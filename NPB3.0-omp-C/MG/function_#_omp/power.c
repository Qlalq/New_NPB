static double power( double a, int n ) {
    double aj;
    int nj;
    double rdummy;
    double power;
    power = 1.0;
    nj = n;
    aj = a;
    while (nj != 0) {
	if( (nj%2) == 1 ) rdummy =  randlc( &power, aj );
	rdummy = randlc( &aj, aj );
	nj = nj/2;
    }
    return (power);
}