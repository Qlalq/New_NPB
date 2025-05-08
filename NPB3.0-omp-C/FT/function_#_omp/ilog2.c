static int ilog2(int n) {
    int nn, lg;
    if (n == 1) {
	return 0;
    }
    lg = 1;
    nn = 2;
    while (nn < n) {
	nn = nn << 1;
	lg++;
    }
    return lg;
}