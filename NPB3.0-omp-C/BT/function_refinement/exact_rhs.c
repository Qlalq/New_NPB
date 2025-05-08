static void exact_rhs(void) {
    // Original local variables
    double dtemp[5];
    double xi, eta, zeta, dtpp;
    int m, i, j, k, ip1, im1, jp1, jm1, km1, kp1;

    // Assuming grid_points, dnxm1, dnym1, dnzm1,
    // tx2, dx1tx1, xxcon1, dx2tx1, c2, cuf, q,
    // xxcon2, dx3tx1, dx4tx1, c1, xxcon3, xxcon4, xxcon5, dx5tx1, dssp,
    // ty2, dy1ty1, yycon2, dy2ty1, dy3ty1, yycon1, dy4ty1, dy5ty1,
    // tz2, dz1tz1, zzcon2, dz2tz1, dz3tz1, dz4tz1, zzcon1, dz5tz1
    // are accessible globally or from a higher scope.
    // Assuming forcing is a global or file-scope array.
    // The original ue, buf, cuf, q are treated as temporary storage
    // that needs to be privatized for parallelization.

    // 1. Initialization of forcing array
    // This loop nest is fully parallelizable across i, j, and k.
    // No dependencies between iterations.
    for (i = 0; i < grid_points[0]; i++) {
        for (j = 0; j < grid_points[1]; j++) {
            for (k = 0; k < grid_points[2]; k++) {
                for (m = 0; m < 5; m++) {
                    forcing[i][j][k][m] = 0.0;
                }
            }
        }
    }

    // The following three major blocks calculate contributions to forcing
    // from derivatives in x, y, and z directions respectively.
    // As written, they accumulate contributions to the same forcing array.
    // This structure implies sequential execution of the three blocks.
    // Parallelism is applied *within* each block by parallelizing the outer loops
    // and making the temporary arrays (ue, buf, cuf, q from original code) private.

    // Block 2: Calculate contributions from xi-direction (i)
    // Loop structure suggests parallelizing over j and k.
    // Temporaries (ue, buf, cuf, q) are used along the i-dimension for a fixed (j, k).
    // Declare them locally within the parallelizable loop scope.
    for (j = 1; j < grid_points[1]-1; j++) { // Candidate for outer parallelization
        eta = (double)j * dnym1;
        for (k = 1; k < grid_points[2]-1; k++) { // Candidate for outer parallelization
            zeta = (double)k * dnzm1;

            // Declare temporary arrays as VLAs (C99) or similar (e.g., thread-private).
            // These replace the global/static ue, buf, cuf, q arrays used in the i-loop.
            double ue_private[grid_points[0]][5];
            double buf_private[grid_points[0]][5];
            double cuf_private[grid_points[0]];
            double q_private[grid_points[0]];

            // Calculate intermediate terms along the i-direction for fixed (j, k)
            // This loop runs sequentially for each (j,k) pair handled by a thread.
            for (i = 0; i < grid_points[0]; i++) {
                xi = (double)i * dnxm1;
                exact_solution(xi, eta, zeta, dtemp); // dtemp is local
                for (m = 0; m < 5; m++) {
                    ue_private[i][m] = dtemp[m];
                }
                dtpp = 1.0 / dtemp[0];
                for (m = 1; m <= 4; m++) {
                    buf_private[i][m] = dtpp * dtemp[m];
                }
                cuf_private[i]   = buf_private[i][1] * buf_private[i][1]; // Note: Uses index 1 for i-dir
                buf_private[i][0] = cuf_private[i] + buf_private[i][2] * buf_private[i][2] +
                                   buf_private[i][3] * buf_private[i][3];
                q_private[i] = 0.5*(buf_private[i][1]*ue_private[i][1] + buf_private[i][2]*ue_private[i][2] +
                                   buf_private[i][3]*ue_private[i][3]);
            }

            // Update forcing based on i-derivatives for fixed (j, k)
            // These updates are independent for different i *for this fixed (j,k)*.
            // This loop runs sequentially for each (j,k) pair.
            for (i = 1; i < grid_points[0]-1; i++) {
                im1 = i-1;
                ip1 = i+1;
                forcing[i][j][k][0] = forcing[i][j][k][0] -
                                      tx2*(ue_private[ip1][1]-ue_private[im1][1])+
                                      dx1tx1*(ue_private[ip1][0]-2.0*ue_private[i][0]+ue_private[im1][0]);
                forcing[i][j][k][1] = forcing[i][j][k][1] -
                                      tx2 * ((ue_private[ip1][1]*buf_private[ip1][1]+c2*(ue_private[ip1][4]-q_private[ip1]))-
                                             (ue_private[im1][1]*buf_private[im1][1]+c2*(ue_private[im1][4]-q_private[im1])))+
                                      xxcon1*(buf_private[ip1][1]-2.0*buf_private[i][1]+buf_private[im1][1])+
                                      dx2tx1*( ue_private[ip1][1]-2.0* ue_private[i][1]+ ue_private[im1][1]);
                forcing[i][j][k][2] = forcing[i][j][k][2] -
                                      tx2 * (ue_private[ip1][2]*buf_private[ip1][1]-ue_private[im1][2]*buf_private[im1][1])+
                                      xxcon2*(buf_private[ip1][2]-2.0*buf_private[i][2]+buf_private[im1][2])+
                                      dx3tx1*( ue_private[ip1][2]-2.0* ue_private[i][2]+ ue_private[im1][2]);
                forcing[i][j][k][3] = forcing[i][j][k][3] -
                                      tx2*(ue_private[ip1][3]*buf_private[ip1][1]-ue_private[im1][3]*buf_private[im1][1])+
                                      xxcon2*(buf_private[ip1][3]-2.0*buf_private[i][3]+buf_private[im1][3])+
                                      dx4tx1*( ue_private[ip1][3]-2.0* ue_private[i][3]+ ue_private[im1][3]);
                forcing[i][j][k][4] = forcing[i][j][k][4] -
                                      tx2*(buf_private[ip1][1]*(c1*ue_private[ip1][4]-c2*q_private[ip1])-
                                           buf_private[im1][1]*(c1*ue_private[im1][4]-c2*q_private[im1]))+
                                      0.5*xxcon3*(buf_private[ip1][0]-2.0*buf_private[i][0]+buf_private[im1][0])+
                                      xxcon4*(cuf_private[ip1]-2.0*cuf_private[i]+cuf_private[im1])+
                                      xxcon5*(buf_private[ip1][4]-2.0*buf_private[i][4]+buf_private[im1][4])+
                                      dx5tx1*( ue_private[ip1][4]-2.0* ue_private[i][4]+ ue_private[im1][4]);
            }

            // Apply i-derivatives for boundary points
            // These loops also run sequentially for each (j,k) pair.
            for (m = 0; m < 5; m++) {
                i = 1;
                forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
                                     (5.0*ue_private[i][m] - 4.0*ue_private[i+1][m] +ue_private[i+2][m]);
                i = 2;
                forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
                                     (-4.0*ue_private[i-1][m] + 6.0*ue_private[i][m] -
                                      4.0*ue_private[i+1][m] +     ue_private[i+2][m]);
            }
            for (m = 0; m < 5; m++) {
                for (i = 1*3; i <= grid_points[0]-3*1-1; i++) {
                    forcing[i][j][k][m] = forcing[i][j][k][m] - dssp*
                                         (ue_private[i-2][m] - 4.0*ue_private[i-1][m] +
                                          6.0*ue_private[i][m] - 4.0*ue_private[i+1][m] + ue_private[i+2][m]);
                }
            }
            for (m = 0; m < 5; m++) {
                i = grid_points[0]-3;
                forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
                                     (ue_private[i-2][m] - 4.0*ue_private[i-1][m] +
                                      6.0*ue_private[i][m] - 4.0*ue_private[i+1][m]);
                i = grid_points[0]-2;
                forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
                                     (ue_private[i-2][m] - 4.0*ue_private[i-1][m] + 5.0*ue_private[i][m]);
            }
        }
    }

    // Block 3: Calculate contributions from eta-direction (j)
    // Loop structure suggests parallelizing over i and k.
    // Temporaries (ue, buf, cuf, q) are used along the j-dimension for a fixed (i, k).
    // Declare them locally within the parallelizable loop scope.
    for (i = 1; i < grid_points[0]-1; i++) { // Candidate for outer parallelization
        xi = (double)i * dnxm1;
        for (k = 1; k < grid_points[2]-1; k++) { // Candidate for outer parallelization
            zeta = (double)k * dnzm1;

            // Declare temporary arrays as VLAs (C99) or similar (e.g., thread-private).
            // These replace the global/static ue, buf, cuf, q arrays used in the j-loop.
            double ue_private[grid_points[1]][5];
            double buf_private[grid_points[1]][5];
            double cuf_private[grid_points[1]];
            double q_private[grid_points[1]];

            // Calculate intermediate terms along the j-direction for fixed (i, k)
            // This loop runs sequentially for each (i,k) pair handled by a thread.
            for (j = 0; j < grid_points[1]; j++) {
                eta = (double)j * dnym1;
                exact_solution(xi, eta, zeta, dtemp);
                for (m = 0; m < 5; m++) {
                    ue_private[j][m] = dtemp[m];
                }
                dtpp = 1.0/dtemp[0];
                for (m = 1; m <= 4; m++) {
                    buf_private[j][m] = dtpp * dtemp[m];
                }
                cuf_private[j]   = buf_private[j][2] * buf_private[j][2]; // Note: Uses index 2 for j-dir
                buf_private[j][0] = cuf_private[j] + buf_private[j][1] * buf_private[j][1] +
                                   buf_private[j][3] * buf_private[j][3];
                q_private[j] = 0.5*(buf_private[j][1]*ue_private[j][1] + buf_private[j][2]*ue_private[j][2] +
                                   buf_private[j][3]*ue_private[j][3]);
            }

            // Update forcing based on j-derivatives for fixed (i, k)
            // These updates are independent for different j *for this fixed (i,k)*.
            // This loop runs sequentially for each (i,k) pair.
            for (j = 1; j < grid_points[1]-1; j++) {
                jm1 = j-1;
                jp1 = j+1;
                forcing[i][j][k][0] = forcing[i][j][k][0] -
                                      ty2*( ue_private[jp1][2]-ue_private[jm1][2] )+
                                      dy1ty1*(ue_private[jp1][0]-2.0*ue_private[j][0]+ue_private[jm1][0]);
                forcing[i][j][k][1] = forcing[i][j][k][1] -
                                      ty2*(ue_private[jp1][1]*buf_private[jp1][2]-ue_private[jm1][1]*buf_private[jm1][2])+
                                      yycon2*(buf_private[jp1][1]-2.0*buf_private[j][1]+buf_private[jm1][1])+
                                      dy2ty1*( ue_private[jp1][1]-2.0* ue_private[j][1]+ ue_private[jm1][1]);
                forcing[i][j][k][2] = forcing[i][j][k][2] -
                                      ty2*((ue_private[jp1][2]*buf_private[jp1][2]+c2*(ue_private[jp1][4]-q_private[jp1]))-
                                           (ue_private[jm1][2]*buf_private[jm1][2]+c2*(ue_private[jm1][4]-q_private[jm1])))+
                                      yycon1*(buf_private[jp1][2]-2.0*buf_private[j][2]+buf_private[jm1][2])+
                                      dy3ty1*( ue_private[jp1][2]-2.0*ue_private[j][2] +ue_private[jm1][2]);
                forcing[i][j][k][3] = forcing[i][j][k][3] -
                                      ty2*(ue_private[jp1][3]*buf_private[jp1][2]-ue_private[jm1][3]*buf_private[jm1][2])+
                                      yycon2*(buf_private[jp1][3]-2.0*buf_private[j][3]+buf_private[jm1][3])+
                                      dy4ty1*( ue_private[jp1][3]-2.0*ue_private[j][3]+ ue_private[jm1][3]);
                forcing[i][j][k][4] = forcing[i][j][k][4] -
                                      ty2*(buf_private[jp1][2]*(c1*ue_private[jp1][4]-c2*q_private[jp1])-
                                           buf_private[jm1][2]*(c1*ue_private[jm1][4]-c2*q_private[jm1]))+
                                      0.5*yycon3*(buf_private[jp1][0]-2.0*buf_private[j][0]+
                                                  buf_private[jm1][0])+
                                      yycon4*(cuf_private[jp1]-2.0*cuf_private[j]+cuf_private[jm1])+
                                      yycon5*(buf_private[jp1][4]-2.0*buf_private[j][4]+buf_private[jm1][4])+
                                      dy5ty1*(ue_private[jp1][4]-2.0*ue_private[j][4]+ue_private[jm1][4]);
            }

            // Apply j-derivatives for boundary points
            // These loops also run sequentially for each (i,k) pair.
            for (m = 0; m < 5; m++) {
                j = 1;
                forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
                                     (5.0*ue_private[j][m] - 4.0*ue_private[j+1][m] +ue_private[j+2][m]);
                j = 2;
                forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
                                     (-4.0*ue_private[j-1][m] + 6.0*ue_private[j][m] -
                                      4.0*ue_private[j+1][m] +       ue_private[j+2][m]);
            }
            for (m = 0; m < 5; m++) {
                for (j = 1*3; j <= grid_points[1]-3*1-1; j++) {
                    forcing[i][j][k][m] = forcing[i][j][k][m] - dssp*
                                         (ue_private[j-2][m] - 4.0*ue_private[j-1][m] +
                                          6.0*ue_private[j][m] - 4.0*ue_private[j+1][m] + ue_private[j+2][m]);
                }
            }
            for (m = 0; m < 5; m++) {
                j = grid_points[1]-3;
                forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
                                     (ue_private[j-2][m] - 4.0*ue_private[j-1][m] +
                                      6.0*ue_private[j][m] - 4.0*ue_private[j+1][m]);
                j = grid_points[1]-2;
                forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
                                     (ue_private[j-2][m] - 4.0*ue_private[j-1][m] + 5.0*ue_private[j][m]);
            }
        }
    }

    // Block 4: Calculate contributions from zeta-direction (k)
    // Loop structure suggests parallelizing over i and j.
    // Temporaries (ue, buf, cuf, q) are used along the k-dimension for a fixed (i, j).
    // Declare them locally within the parallelizable loop scope.
    for (i = 1; i < grid_points[0]-1; i++) { // Candidate for outer parallelization
        xi = (double)i * dnxm1;
        for (j = 1; j < grid_points[1]-1; j++) { // Candidate for outer parallelization
            eta = (double)j * dnym1;

            // Declare temporary arrays as VLAs (C99) or similar (e.g., thread-private).
            // These replace the global/static ue, buf, cuf, q arrays used in the k-loop.
            double ue_private[grid_points[2]][5];
            double buf_private[grid_points[2]][5];
            double cuf_private[grid_points[2]];
            double q_private[grid_points[2]];

            // Calculate intermediate terms along the k-direction for fixed (i, j)
            // This loop runs sequentially for each (i,j) pair handled by a thread.
            for (k = 0; k < grid_points[2]; k++) {
                zeta = (double)k * dnzm1;
                exact_solution(xi, eta, zeta, dtemp);
                for (m = 0; m < 5; m++) {
                    ue_private[k][m] = dtemp[m];
                }
                dtpp = 1.0/dtemp[0];
                for (m = 1; m <= 4; m++) {
                    buf_private[k][m] = dtpp * dtemp[m];
                }
                cuf_private[k]   = buf_private[k][3] * buf_private[k][3]; // Note: Uses index 3 for k-dir
                buf_private[k][0] = cuf_private[k] + buf_private[k][1] * buf_private[k][1] +
                                   buf_private[k][2] * buf_private[k][2];
                q_private[k] = 0.5*(buf_private[k][1]*ue_private[k][1] + buf_private[k][2]*ue_private[k][2] +
                                   buf_private[k][3]*ue_private[k][3]);
            }

            // Update forcing based on k-derivatives for fixed (i, j)
            // These updates are independent for different k *for this fixed (i,j)*.
            // This loop runs sequentially for each (i,j) pair.
            for (k = 1; k < grid_points[2]-1; k++) {
                km1 = k-1;
                kp1 = k+1;
                forcing[i][j][k][0] = forcing[i][j][k][0] -
                                      tz2*( ue_private[kp1][3]-ue_private[km1][3] )+
                                      dz1tz1*(ue_private[kp1][0]-2.0*ue_private[k][0]+ue_private[km1][0]);
                forcing[i][j][k][1] = forcing[i][j][k][1] -
                                      tz2 * (ue_private[kp1][1]*buf_private[kp1][3]-ue_private[km1][1]*buf_private[km1][3])+
                                      zzcon2*(buf_private[kp1][1]-2.0*buf_private[k][1]+buf_private[km1][1])+
                                      dz2tz1*( ue_private[kp1][1]-2.0* ue_private[k][1]+ ue_private[km1][1]);
                forcing[i][j][k][2] = forcing[i][j][k][2] -
                                      tz2 * (ue_private[kp1][2]*buf_private[kp1][3]-ue_private[km1][2]*buf_private[km1][3])+
                                      zzcon2*(buf_private[kp1][2]-2.0*buf_private[k][2]+buf_private[km1][2])+
                                      dz3tz1*(ue_private[kp1][2]-2.0*ue_private[k][2]+ue_private[km1][2]);
                forcing[i][j][k][3] = forcing[i][j][k][3] -
                                      tz2 * ((ue_private[kp1][3]*buf_private[kp1][3]+c2*(ue_private[kp1][4]-q_private[kp1]))-
                                             (ue_private[km1][3]*buf_private[km1][3]+c2*(ue_private[km1][4]-q_private[km1])))+
                                      zzcon1*(buf_private[kp1][3]-2.0*buf_private[k][3]+buf_private[km1][3])+
                                      dz4tz1*( ue_private[kp1][3]-2.0*ue_private[k][3] +ue_private[km1][3]);
                forcing[i][j][k][4] = forcing[i][j][k][4] -
                                      tz2 * (buf_private[kp1][3]*(c1*ue_private[kp1][4]-c2*q_private[kp1])-
                                             buf_private[km1][3]*(c1*ue_private[km1][4]-c2*q_private[km1]))+
                                      0.5*zzcon3*(buf_private[kp1][0]-2.0*buf_private[k][0]
                                                  +buf_private[km1][0])+
                                      zzcon4*(cuf_private[kp1]-2.0*cuf_private[k]+cuf_private[km1])+
                                      zzcon5*(buf_private[kp1][4]-2.0*buf_private[k][4]+buf_private[km1][4])+
                                      dz5tz1*( ue_private[kp1][4]-2.0*ue_private[k][4]+ ue_private[km1][4]);
            }

            // Apply k-derivatives for boundary points
            // These loops also run sequentially for each (i,j) pair.
            for (m = 0; m < 5; m++) {
                k = 1;
                forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
                                     (5.0*ue_private[k][m] - 4.0*ue_private[k+1][m] +ue_private[k+2][m]);
                k = 2;
                forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
                                     (-4.0*ue_private[k-1][m] + 6.0*ue_private[k][m] -
                                      4.0*ue_private[k+1][m] +       ue_private[k+2][m]);
            }
            for (m = 0; m < 5; m++) {
                for (k = 1*3; k <= grid_points[2]-3*1-1; k++) {
                    forcing[i][j][k][m] = forcing[i][j][k][m] - dssp*
                                         (ue_private[k-2][m] - 4.0*ue_private[k-1][m] +
                                          6.0*ue_private[k][m] - 4.0*ue_private[k+1][m] + ue_private[k+2][m]);
                }
            }
            for (m = 0; m < 5; m++) {
                k = grid_points[2]-3;
                forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
                                     (ue_private[k-2][m] - 4.0*ue_private[k-1][m] +
                                      6.0*ue_private[k][m] - 4.0*ue_private[k+1][m]);
                k = grid_points[2]-2;
                forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
                                     (ue_private[k-2][m] - 4.0*ue_private[k-1][m] + 5.0*ue_private[k][m]);
            }
        }
    }

    // Final negation of the forcing array
    // This loop nest is fully parallelizable across i, j, k, and m.
    // No dependencies between iterations.
    for (i = 1; i < grid_points[0]-1; i++) { // Candidate for parallelization
        for (j = 1; j < grid_points[1]-1; j++) { // Candidate for parallelization
            for (k = 1; k < grid_points[2]-1; k++) { // Candidate for parallelization
                for (m = 0; m < 5; m++) { // Candidate for parallelization
                    forcing[i][j][k][m] = -1.0 * forcing[i][j][k][m];
                }
            }
        }
    }
}