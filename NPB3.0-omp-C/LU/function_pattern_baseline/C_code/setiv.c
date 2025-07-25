#include <omp.h>

static void setiv(void) {
    int i, j, k, m;
    int iglob, jglob;
    double xi, eta, zeta;
    double pxi, peta, pzeta;
    double ue_1jk[5], ue_nx0jk[5], ue_i1k[5],
           ue_iny0k[5], ue_ij1[5], ue_ijnz[5];

    /* Parallelize over j and k dimensions; each thread
       works on a distinct (j,k) slab. */
    #pragma omp parallel for collapse(2)                         \
        private(i, iglob, jglob, xi, eta, zeta,                  \
                pxi, peta, pzeta,                                \
                ue_1jk, ue_nx0jk, ue_i1k, ue_iny0k, ue_ij1, ue_ijnz) \
        schedule(static)
    for (j = 0; j < ny; j++) {
        for (k = 1; k < nz - 1; k++) {
            jglob = j;
            zeta  = (double)k / (nz - 1);

            /* skip physical boundaries in j */
            if (jglob != 0 && jglob != ny0 - 1) {
                eta = (double)jglob / (ny0 - 1);

                for (i = 0; i < nx; i++) {
                    iglob = i;
                    /* skip physical boundaries in i */
                    if (iglob != 0 && iglob != nx0 - 1) {
                        xi = (double)iglob / (nx0 - 1);

                        /* sample exact solution on the six faces */
                        exact(0,           jglob, k, ue_1jk);
                        exact(nx0 - 1,     jglob, k, ue_nx0jk);
                        exact(iglob,       0,     k, ue_i1k);
                        exact(iglob, ny0 - 1,     k, ue_iny0k);
                        exact(iglob,   jglob,   0, ue_ij1);
                        exact(iglob,   jglob, nz - 1, ue_ijnz);

                        /* vectorize over the 5-component state */
                        #pragma omp simd
                        for (m = 0; m < 5; m++) {
                            pxi   = (1.0 - xi) * ue_1jk[m]   + xi   * ue_nx0jk[m];
                            peta  = (1.0 - eta) * ue_i1k[m]   + eta  * ue_iny0k[m];
                            pzeta = (1.0 - zeta) * ue_ij1[m]   + zeta * ue_ijnz[m];

                            u[i][j][k][m] = pxi + peta + pzeta
                                          - pxi * peta
                                          - peta * pzeta
                                          - pzeta * pxi
                                          + pxi * peta * pzeta;
                        }
                    }
                }
            }
        }
    }
}