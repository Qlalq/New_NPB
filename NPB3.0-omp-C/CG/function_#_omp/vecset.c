static void vecset(
    int n,
    double v[],	
    int iv[],	
    int *nzv,
    int i,
    double val)
{
    int k;
    boolean set;
    set = FALSE;
    for (k = 1; k <= *nzv; k++) {
	if (iv[k] == i) {
            v[k] = val;
            set  = TRUE;
	}
    }
    if (set == FALSE) {
	*nzv = *nzv + 1;
	v[*nzv] = val;
	iv[*nzv] = i;
    }
}