static void setcoeff(void) {
#pragma omp parallel default(none) \
    shared(nx0, ny0, nz0, ce, \
           dxi, deta, dzeta, \
           tx1, tx2, tx3, \
           ty1, ty2, ty3, \
           tz1, tz2, tz3, \
           ii1, ii2, ji1, ji2, ki1, ki2, \
           dx1, dx2, dx3, dx4, dx5, \
           dy1, dy2, dy3, dy4, dy5, \
           dz1, dz2, dz3, dz4, dz5, \
           dssp)
  {
#pragma omp sections
    {
#pragma omp section
      {
        dxi   = 1.0 / (nx0 - 1);
        deta  = 1.0 / (ny0 - 1);
        dzeta = 1.0 / (nz0 - 1);
      }

#pragma omp section
      {
        tx1 = 1.0 / (dxi * dxi);
        tx2 = 1.0 / (2.0 * dxi);
        tx3 = 1.0 / dxi;
        ty1 = 1.0 / (deta * deta);
        ty2 = 1.0 / (2.0 * deta);
        ty3 = 1.0 / deta;
        tz1 = 1.0 / (dzeta * dzeta);
        tz2 = 1.0 / (2.0 * dzeta);
        tz3 = 1.0 / dzeta;
      }

#pragma omp section
      {
        ii1 = 1;           ii2 = nx0 - 2;
        ji1 = 1;           ji2 = ny0 - 3;
        ki1 = 2;           ki2 = nz0 - 2;
      }

#pragma omp section
      {
        dx1 = dx2 = dx3 = dx4 = dx5 = 0.75;
        dy1 = dy2 = dy3 = dy4 = dy5 = 0.75;
        dz1 = dz2 = dz3 = dz4 = dz5 = 1.00;
      }

#pragma omp section
      {
        dssp = (max(dx1, max(dy1, dz1))) / 4.0;
      }

#pragma omp section
      {
        /* CE row 0 */
        ce[0][0]=2.0;   ce[0][1]=0.0;   ce[0][2]=0.0;   ce[0][3]=4.0;
        ce[0][4]=5.0;   ce[0][5]=3.0;   ce[0][6]=0.5;   ce[0][7]=0.02;
        ce[0][8]=0.01;  ce[0][9]=0.03;  ce[0][10]=0.5;  ce[0][11]=0.4;
        ce[0][12]=0.3;
      }

#pragma omp section
      {
        /* CE row 1 */
        ce[1][0]=1.0;   ce[1][1]=0.0;   ce[1][2]=0.0;   ce[1][3]=0.0;
        ce[1][4]=1.0;   ce[1][5]=2.0;   ce[1][6]=3.0;   ce[1][7]=0.01;
        ce[1][8]=0.03;  ce[1][9]=0.02;  ce[1][10]=0.4;  ce[1][11]=0.3;
        ce[1][12]=0.5;
      }

#pragma omp section
      {
        /* CE row 2 */
        ce[2][0]=2.0;   ce[2][1]=2.0;   ce[2][2]=0.0;   ce[2][3]=0.0;
        ce[2][4]=0.0;   ce[2][5]=2.0;   ce[2][6]=3.0;   ce[2][7]=0.04;
        ce[2][8]=0.03;  ce[2][9]=0.05;  ce[2][10]=0.3;  ce[2][11]=0.5;
        ce[2][12]=0.4;
      }

#pragma omp section
      {
        /* CE row 3 */
        ce[3][0]=2.0;   ce[3][1]=2.0;   ce[3][2]=0.0;   ce[3][3]=0.0;
        ce[3][4]=0.0;   ce[3][5]=2.0;   ce[3][6]=3.0;   ce[3][7]=0.03;
        ce[3][8]=0.05;  ce[3][9]=0.04;  ce[3][10]=0.2;  ce[3][11]=0.1;
        ce[3][12]=0.3;
      }

#pragma omp section
      {
        /* CE row 4 */
        ce[4][0]=5.0;   ce[4][1]=4.0;   ce[4][2]=3.0;   ce[4][3]=2.0;
        ce[4][4]=0.1;   ce[4][5]=0.4;   ce[4][6]=0.3;   ce[4][7]=0.05;
        ce[4][8]=0.04;  ce[4][9]=0.03;  ce[4][10]=0.1;  ce[4][11]=0.3;
        ce[4][12]=0.2;
      }
    }  // end sections
  }    // end parallel
}