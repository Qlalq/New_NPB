static void fft(int dir, dcomplex x1[NZ][NY][NX], dcomplex x2[NZ][NY][NX]) {

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

    dcomplex y0[NX][FFTBLOCKPAD];
    dcomplex y1[NX][FFTBLOCKPAD];

/*--------------------------------------------------------------------
c note: args x1, x2 must be different arrays
c note: args for cfftsx are (direction, layout, xin, xout, scratch)
c       xin/xout may be the same and it can be somewhat faster
c       if they are
c-------------------------------------------------------------------*/

    if (dir == 1) {
        cffts1(1, dims[0], x1, x1, y0, y1);	/* x1 -> x1 */
        cffts2(1, dims[1], x1, x1, y0, y1);	/* x1 -> x1 */
        cffts3(1, dims[2], x1, x2, y0, y1);	/* x1 -> x2 */
    } else {
	cffts3(-1, dims[2], x1, x1, y0, y1);	/* x1 -> x1 */
    cffts2(-1, dims[1], x1, x1, y0, y1);	/* x1 -> x1 */
    cffts1(-1, dims[0], x1, x2, y0, y1);	/* x1 -> x2 */
    }
}