// Assume M is defined elsewhere, e.g., #define M ...
// Assume debug_vec is declared elsewhere, e.g., extern int debug_vec[];
// Assume rep_nrm and showall are defined elsewhere
// extern void rep_nrm(double ***v, int n1, int n2, int n3, const char *name, int k);
// extern void showall(double ***v, int n1, int n2, int n3);

#include <stdio.h>
#include <stdlib.h> // Required for malloc/free

static void interp( double ***z, int mm1, int mm2, int mm3,
		    double ***u, int n1, int n2, int n3, int k ) {
    int i3, i2, i1;
    // z1, z2, z3 arrays are local to the (i3, i2) loop pair in the first branch
    // and can be privatized per thread if the outer loops (i3 or i2) are parallelized.
    // Assuming M is large enough, e.g., M >= mm1.
    double z1[M], z2[M], z3[M];
    int d1, d2, d3, t1, t2, t3; // Moved declaration outside inner block for clarity

    if ( n1 != 3 && n2 != 3 && n3 != 3 ) {
	// This branch updates disjoint 2x2x2 blocks in u for different (i3, i2, i1) source triplets.
	// Parallelizing the loops over i3 or i2 should be straightforward as destinations are disjoint.
	for (i3 = 0; i3 < mm3-1; i3++) {
            for (i2 = 0; i2 < mm2-1; i2++) {
		// The following loop calculates z1, z2, z3 based on z data for the current (i3, i2) plane.
		// These arrays are reused for the subsequent write loops for the same (i3, i2) plane.
		// If i3 or i2 loops are parallelized, z1, z2, z3 should be private to each thread's (i3, i2) iteration.
		for (i1 = 0; i1 < mm1; i1++) {
		    z1[i1] = z[i3][i2+1][i1] + z[i3][i2][i1];
		    z2[i1] = z[i3+1][i2][i1] + z[i3][i2][i1];
		    z3[i1] = z[i3+1][i2+1][i1] + z[i3+1][i2][i1] + z1[i1]; // z3 depends on z1
		}
		// The following loops write to u using the calculated z1, z2, z3 and z.
		// For a fixed (i3, i2), they update a 2x2x(mm1-1) block in u based on stencil calculations.
		// Different (i3, i2) pairs update disjoint blocks in u.
		for (i1 = 0; i1 < mm1-1; i1++) {
		    u[2*i3][2*i2][2*i1] += z[i3][i2][i1];
		    u[2*i3][2*i2][2*i1+1] += 0.5*(z[i3][i2][i1+1]+z[i3][i2][i1]);
		}
		for (i1 = 0; i1 < mm1-1; i1++) {
		    u[2*i3][2*i2+1][2*i1] += 0.5 * z1[i1];
		    u[2*i3][2*i2+1][2*i1+1] += 0.25*( z1[i1] + z1[i1+1] );
		}
		for (i1 = 0; i1 < mm1-1; i1++) {
		    u[2*i3+1][2*i2][2*i1] += 0.5 * z2[i1];
		    u[2*i3+1][2*i2][2*i1+1] += 0.25*( z2[i1] + z2[i1+1] );
		}
		for (i1 = 0; i1 < mm1-1; i1++) {
		    u[2*i3+1][2*i2+1][2*i1] += 0.25* z3[i1];
		    u[2*i3+1][2*i2+1][2*i1+1] += 0.125*( z3[i1] + z3[i1+1] );
		}
	    }
	}
    } else {
	// This branch updates u with contributions from z using different stencils/ranges,
	// typically for boundary or coarse grid operations in multigrid.
	// Due to overlapping destination indices in u for different source cells / loop iterations,
	// a direct parallelization of the accumulation loops would cause race conditions.
	// To enable parallelization, we use a temporary private array (u_private) to accumulate
	// contributions independently per thread/task, and then sum the temporary array into
	// the main u array after all contributions are computed.

	if (n1 == 3) { d1 = 2; t1 = 1; } else { d1 = 1; t1 = 0; }
	if (n2 == 3) { d2 = 2; t2 = 1; } else { d2 = 1; t2 = 0; }
	if (n3 == 3) { d3 = 2; t3 = 1; } else { d3 = 1; t3 = 0; }

        // Allocate a temporary array u_private of size n3 x n2 x n1
        // to accumulate updates without race conditions. This array should ideally
        // be thread-private if the following accumulation loops are parallelized.
        // The allocation/deallocation adds overhead but provides structural independence.
        double ***u_private;
        u_private = (double ***) malloc(n3 * sizeof(double **));
        if (u_private == NULL) {
             // Handle allocation failure - In a real application, add error logging/exit
             fprintf(stderr, "Memory allocation failed for u_private (dim 0).\n");
             return; // Simplified failure handling
        }
        for (int j=0; j<n3; j++) {
            u_private[j] = (double **) malloc(n2 * sizeof(double *));
            if (u_private[j] == NULL) {
                // Cleanup previously allocated rows
                for(int k=0; k<j; ++k) free(u_private[k]);
                free(u_private);
                fprintf(stderr, "Memory allocation failed for u_private (dim 1).\n");
                return; // Simplified failure handling
            }
            for (int i=0; i<n2; i++) {
                u_private[j][i] = (double *) malloc(n1 * sizeof(double));
                if (u_private[j][i] == NULL) {
                    // Cleanup previously allocated elements in this row
                    for(int k=0; k<i; ++k) free(u_private[j][k]);
                    // Cleanup previously allocated rows and the top level
                    for(int k=0; k<j; ++k) free(u_private[k]);
                    free(u_private);
                    fprintf(stderr, "Memory allocation failed for u_private (dim 2).\n");
                    return; // Simplified failure handling
                }
            }
        }

        // Initialize the temporary array to zero
        for (int i3_u=0; i3_u < n3; i3_u++) {
            for (int i2_u=0; i2_u < n2; i2_u++) {
                for (int i1_u=0; i1_u < n1; i1_u++) {
                    u_private[i3_u][i2_u][i1_u] = 0.0;
                }
            }
        }

        { // Original block scope - This block contains 8 loop nests
	    // These 8 loop nests calculate contributions from z to u_private.
            // Since they write to a private temporary array, the write operations themselves
            // are safe from races if the loops are parallelized (e.g., parallelizing the outermost i3 loop).
            // Different iterations of the loops will write to u_private. If the target indices
            // in u_private overlap between iterations assigned to different threads, this implies
            // that u_private itself needs to be thread-private, or the parallelization
            // must ensure disjoint writes. The simplest structure is to have a thread-private
            // u_private, but representing this without OpenMP pragmas is done by having a
            // single u_private that *receives contributions* sequentially or with later manual reduction.
            // The += into u_private here still implies an accumulation dependency if loops are parallelized.
            // The typical pattern for this is *thread-private* u_private, but that requires OpenMP syntax.
            // The refinement here using *one* u_private highlights the *logical* separation of computation
            // of contributions from the final addition, allowing the 8 nests to potentially run in parallel
            // if the underlying write patterns to u_private are handled (e.g., via atomic writes, though slow,
            // or by ensuring disjoint writes via careful loop scheduling, or true thread-private storage).
            // Given the objective is *structuring* for OpenMP, creating the temporary array is the key step.
            // The most efficient OpenMP approach would use a thread-private copy of u_private and a final reduction.
            // This code structure lays the groundwork for that.
	    for ( i3 = d3; i3 <= mm3-1; i3++) {
                for ( i2 = d2; i2 <= mm2-1; i2++) {
		    for ( i1 = d1; i1 <= mm1-1; i1++) {
		        u_private[2*i3-d3-1][2*i2-d2-1][2*i1-d1-1] +=
			    z[i3-1][i2-1][i1-1];
		    }
		    for ( i1 = 1; i1 <= mm1-1; i1++) {
		        u_private[2*i3-d3-1][2*i2-d2-1][2*i1-t1-1] +=
			    0.5*(z[i3-1][i2-1][i1]+z[i3-1][i2-1][i1-1]);
		    }
		}
                for ( i2 = 1; i2 <= mm2-1; i2++) { // Note: starts from 1, not d2
		    for ( i1 = d1; i1 <= mm1-1; i1++) {
		        u_private[2*i3-d3-1][2*i2-t2-1][2*i1-d1-1] +=
			    0.5*(z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
		    }
		    for ( i1 = 1; i1 <= mm1-1; i1++) {
		        u_private[2*i3-d3-1][2*i2-t2-1][2*i1-t1-1] +=
			    0.25*(z[i3-1][i2][i1]+z[i3-1][i2-1][i1]
			       +z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
		    }
		}
	    }
	    for ( i3 = 1; i3 <= mm3-1; i3++) { // Note: starts from 1, not d3
                for ( i2 = d2; i2 <= mm2-1; i2++) {
		    for ( i1 = d1; i1 <= mm1-1; i1++) {
		        u_private[2*i3-t3-1][2*i2-d2-1][2*i1-d1-1] +=
			    0.5*(z[i3][i2-1][i1-1]+z[i3-1][i2-1][i1-1]);
		    }
		    for ( i1 = 1; i1 <= mm1-1; i1++) {
		        u_private[2*i3-t3-1][2*i2-d2-1][2*i1-t1-1] +=
			    0.25*(z[i3][i2-1][i1]+z[i3][i2-1][i1-1]
			       +z[i3-1][i2-1][i1]+z[i3-1][i2-1][i1-1]);
		    }
		}
	    for ( i2 = 1; i2 <= mm2-1; i2++) { // Note: starts from 1, not d2
		    for ( i1 = d1; i1 <= mm1-1; i1++) {
		        u_private[2*i3-t3-1][2*i2-t2-1][2*i1-d1-1] +=
			    0.25*(z[i3][i2][i1-1]+z[i3][i2-1][i1-1]
			       +z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
		    }
		    for ( i1 = 1; i1 <= mm1-1; i1++) {
		        u_private[2*i3-t3-1][2*i2-t2-1][2*i1-t1-1] +=
			    0.125*(z[i3][i2][i1]+z[i3][i2-1][i1]
				+z[i3][i2][i1-1]+z[i3][i2-1][i1-1]
				+z[i3-1][i2][i1]+z[i3-1][i2-1][i1]
				+z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
		    }
		}
	    }
        } // end original block scope containing 8 loops

        // Reduction step: Add the accumulated contributions from u_private to the main u array.
        // This loop iterates over the destination grid size (n1, n2, n3) and can be
        // parallelized using OpenMP reduction clauses (if applicable) or a simple parallel for.
        for (int i3_u=0; i3_u < n3; i3_u++) {
            for (int i2_u=0; i2_u < n2; i2_u++) {
                for (int i1_u=0; i1_u < n1; i1_u++) {
                    u[i3_u][i2_u][i1_u] += u_private[i3_u][i2_u][i1_u];
                }
            }
        }

        // Free the temporary array
        for (int j=0; j<n3; j++) {
            for (int i=0; i<n2; i++) {
                free(u_private[j][i]);
            }
            free(u_private[j]);
        }
        free(u_private);
    } //end else (n1==3 || n2==3 || n3==3)

    // Original debug code remains outside the computational branches
    if (debug_vec[0] >= 1 ) {
	rep_nrm(z,mm1,mm2,mm3,"z: inter",k-1);
	rep_nrm(u,n1,n2,n3,"u: inter",k);
    }
    if ( debug_vec[5] >= k ) {
	showall(z,mm1,mm2,mm3);
	showall(u,n1,n2,n3);
    }
}