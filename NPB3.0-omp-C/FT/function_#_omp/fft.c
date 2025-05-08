static void fft(int dir, dcomplex x1[NZ][NY][NX], dcomplex x2[NZ][NY][NX]) {
    dcomplex y0[NX][FFTBLOCKPAD];
    dcomplex y1[NX][FFTBLOCKPAD];
    if (dir == 1) {
        cffts1(1, dims[0], x1, x1, y0, y1);	
        cffts2(1, dims[1], x1, x1, y0, y1);	
        cffts3(1, dims[2], x1, x2, y0, y1);	
    } else {
	cffts3(-1, dims[2], x1, x1, y0, y1);	
    cffts2(-1, dims[1], x1, x1, y0, y1);	
    cffts1(-1, dims[0], x1, x2, y0, y1);	
    }
}