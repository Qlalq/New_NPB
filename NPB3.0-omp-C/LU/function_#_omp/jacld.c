static void jacld(int k) {
  int i, j;
  double  r43;
  double  c1345;
  double  c34;
  double  tmp1, tmp2, tmp3;
  r43 = ( 4.0 / 3.0 );
  c1345 = C1 * C3 * C4 * C5;
  c34 = C3 * C4;
  for (i = ist; i <= iend; i++) {
    for (j = jst; j <= jend; j++) {
      tmp1 = 1.0 / u[i][j][k][0];
      tmp2 = tmp1 * tmp1;
      tmp3 = tmp1 * tmp2;
      d[i][j][0][0] =  1.0
	+ dt * 2.0 * (   tx1 * dx1
			 + ty1 * dy1
			 + tz1 * dz1 );
      d[i][j][0][1] =  0.0;
      d[i][j][0][2] =  0.0;
      d[i][j][0][3] =  0.0;
      d[i][j][0][4] =  0.0;
      d[i][j][1][0] =  dt * 2.0
	* (  tx1 * ( - r43 * c34 * tmp2 * u[i][j][k][1] )
	     + ty1 * ( -       c34 * tmp2 * u[i][j][k][1] )
	     + tz1 * ( -       c34 * tmp2 * u[i][j][k][1] ) );
      d[i][j][1][1] =  1.0
	+ dt * 2.0 
	* (  tx1 * r43 * c34 * tmp1
	     + ty1 *       c34 * tmp1
	     + tz1 *       c34 * tmp1 )
	+ dt * 2.0 * (   tx1 * dx2
			 + ty1 * dy2
			 + tz1 * dz2  );
      d[i][j][1][2] = 0.0;
      d[i][j][1][3] = 0.0;
      d[i][j][1][4] = 0.0;
      d[i][j][2][0] = dt * 2.0
	* (  tx1 * ( -       c34 * tmp2 * u[i][j][k][2] )
	     + ty1 * ( - r43 * c34 * tmp2 * u[i][j][k][2] )
	     + tz1 * ( -       c34 * tmp2 * u[i][j][k][2] ) );
      d[i][j][2][1] = 0.0;
      d[i][j][2][2] = 1.0
	+ dt * 2.0
	* (  tx1 *       c34 * tmp1
	     + ty1 * r43 * c34 * tmp1
	     + tz1 *       c34 * tmp1 )
	+ dt * 2.0 * (  tx1 * dx3
			+ ty1 * dy3
			+ tz1 * dz3 );
      d[i][j][2][3] = 0.0;
      d[i][j][2][4] = 0.0;
      d[i][j][3][0] = dt * 2.0
	* (  tx1 * ( -       c34 * tmp2 * u[i][j][k][3] )
	     + ty1 * ( -       c34 * tmp2 * u[i][j][k][3] )
	     + tz1 * ( - r43 * c34 * tmp2 * u[i][j][k][3] ) );
      d[i][j][3][1] = 0.0;
      d[i][j][3][2] = 0.0;
      d[i][j][3][3] = 1.0
	+ dt * 2.0
	* (  tx1 *       c34 * tmp1
	     + ty1 *       c34 * tmp1
	     + tz1 * r43 * c34 * tmp1 )
	+ dt * 2.0 * (  tx1 * dx4
			+ ty1 * dy4
			+ tz1 * dz4 );
      d[i][j][3][4] = 0.0;
      d[i][j][4][0] = dt * 2.0
	* ( tx1 * ( - ( r43*c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][1]) )
		    - ( c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][2]) )
		    - ( c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][3]) )
		    - ( c1345 ) * tmp2 * u[i][j][k][4] )
	    + ty1 * ( - ( c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][1]) )
		      - ( r43*c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][2]) )
		      - ( c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][3]) )
		      - ( c1345 ) * tmp2 * u[i][j][k][4] )
	    + tz1 * ( - ( c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][1]) )
		      - ( c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][2]) )
		      - ( r43*c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][3]) )
		      - ( c1345 ) * tmp2 * u[i][j][k][4] ) );
      d[i][j][4][1] = dt * 2.0
	* ( tx1 * ( r43*c34 - c1345 ) * tmp2 * u[i][j][k][1]
	    + ty1 * (     c34 - c1345 ) * tmp2 * u[i][j][k][1]
	    + tz1 * (     c34 - c1345 ) * tmp2 * u[i][j][k][1] );
      d[i][j][4][2] = dt * 2.0
	* ( tx1 * ( c34 - c1345 ) * tmp2 * u[i][j][k][2]
	    + ty1 * ( r43*c34 -c1345 ) * tmp2 * u[i][j][k][2]
	    + tz1 * ( c34 - c1345 ) * tmp2 * u[i][j][k][2] );
      d[i][j][4][3] = dt * 2.0
	* ( tx1 * ( c34 - c1345 ) * tmp2 * u[i][j][k][3]
	    + ty1 * ( c34 - c1345 ) * tmp2 * u[i][j][k][3]
	    + tz1 * ( r43*c34 - c1345 ) * tmp2 * u[i][j][k][3] );
      d[i][j][4][4] = 1.0
	+ dt * 2.0 * ( tx1 * c1345 * tmp1
		       + ty1 * c1345 * tmp1
		       + tz1 * c1345 * tmp1 )
        + dt * 2.0 * (  tx1 * dx5
			+  ty1 * dy5
			+  tz1 * dz5 );
      tmp1 = 1.0 / u[i][j][k-1][0];
      tmp2 = tmp1 * tmp1;
      tmp3 = tmp1 * tmp2;
      a[i][j][0][0] = - dt * tz1 * dz1;
      a[i][j][0][1] =   0.0;
      a[i][j][0][2] =   0.0;
      a[i][j][0][3] = - dt * tz2;
      a[i][j][0][4] =   0.0;
      a[i][j][1][0] = - dt * tz2
	* ( - ( u[i][j][k-1][1]*u[i][j][k-1][3] ) * tmp2 )
	- dt * tz1 * ( - c34 * tmp2 * u[i][j][k-1][1] );
      a[i][j][1][1] = - dt * tz2 * ( u[i][j][k-1][3] * tmp1 )
	- dt * tz1 * c34 * tmp1
	- dt * tz1 * dz2 ;
      a[i][j][1][2] = 0.0;
      a[i][j][1][3] = - dt * tz2 * ( u[i][j][k-1][1] * tmp1 );
      a[i][j][1][4] = 0.0;
      a[i][j][2][0] = - dt * tz2
	* ( - ( u[i][j][k-1][2]*u[i][j][k-1][3] ) * tmp2 )
	- dt * tz1 * ( - c34 * tmp2 * u[i][j][k-1][2] );
      a[i][j][2][1] = 0.0;
      a[i][j][2][2] = - dt * tz2 * ( u[i][j][k-1][3] * tmp1 )
	- dt * tz1 * ( c34 * tmp1 )
	- dt * tz1 * dz3;
      a[i][j][2][3] = - dt * tz2 * ( u[i][j][k-1][2] * tmp1 );
      a[i][j][2][4] = 0.0;
      a[i][j][3][0] = - dt * tz2
	* ( - ( u[i][j][k-1][3] * tmp1 ) *( u[i][j][k-1][3] * tmp1 )
	    + 0.50 * C2
	    * ( ( u[i][j][k-1][1] * u[i][j][k-1][1]
		  + u[i][j][k-1][2] * u[i][j][k-1][2]
		  + u[i][j][k-1][3] * u[i][j][k-1][3] ) * tmp2 ) )
	- dt * tz1 * ( - r43 * c34 * tmp2 * u[i][j][k-1][3] );
      a[i][j][3][1] = - dt * tz2
	* ( - C2 * ( u[i][j][k-1][1] * tmp1 ) );
      a[i][j][3][2] = - dt * tz2
	* ( - C2 * ( u[i][j][k-1][2] * tmp1 ) );
      a[i][j][3][3] = - dt * tz2 * ( 2.0 - C2 )
	* ( u[i][j][k-1][3] * tmp1 )
	- dt * tz1 * ( r43 * c34 * tmp1 )
	- dt * tz1 * dz4;
      a[i][j][3][4] = - dt * tz2 * C2;
      a[i][j][4][0] = - dt * tz2
	* ( ( C2 * (  u[i][j][k-1][1] * u[i][j][k-1][1]
                      + u[i][j][k-1][2] * u[i][j][k-1][2]
                      + u[i][j][k-1][3] * u[i][j][k-1][3] ) * tmp2
	      - C1 * ( u[i][j][k-1][4] * tmp1 ) )
	    * ( u[i][j][k-1][3] * tmp1 ) )
	- dt * tz1
	* ( - ( c34 - c1345 ) * tmp3 * (u[i][j][k-1][1]*u[i][j][k-1][1])
	    - ( c34 - c1345 ) * tmp3 * (u[i][j][k-1][2]*u[i][j][k-1][2])
	    - ( r43*c34 - c1345 )* tmp3 * (u[i][j][k-1][3]*u[i][j][k-1][3])
	    - c1345 * tmp2 * u[i][j][k-1][4] );
      a[i][j][4][1] = - dt * tz2
	* ( - C2 * ( u[i][j][k-1][1]*u[i][j][k-1][3] ) * tmp2 )
	- dt * tz1 * ( c34 - c1345 ) * tmp2 * u[i][j][k-1][1];
      a[i][j][4][2] = - dt * tz2
	* ( - C2 * ( u[i][j][k-1][2]*u[i][j][k-1][3] ) * tmp2 )
	- dt * tz1 * ( c34 - c1345 ) * tmp2 * u[i][j][k-1][2];
      a[i][j][4][3] = - dt * tz2
	* ( C1 * ( u[i][j][k-1][4] * tmp1 )
            - 0.50 * C2
            * ( (  u[i][j][k-1][1]*u[i][j][k-1][1]
		   + u[i][j][k-1][2]*u[i][j][k-1][2]
		   + 3.0*u[i][j][k-1][3]*u[i][j][k-1][3] ) * tmp2 ) )
	- dt * tz1 * ( r43*c34 - c1345 ) * tmp2 * u[i][j][k-1][3];
      a[i][j][4][4] = - dt * tz2
	* ( C1 * ( u[i][j][k-1][3] * tmp1 ) )
	- dt * tz1 * c1345 * tmp1
	- dt * tz1 * dz5;
      tmp1 = 1.0 / u[i][j-1][k][0];
      tmp2 = tmp1 * tmp1;
      tmp3 = tmp1 * tmp2;
      b[i][j][0][0] = - dt * ty1 * dy1;
      b[i][j][0][1] =   0.0;
      b[i][j][0][2] = - dt * ty2;
      b[i][j][0][3] =   0.0;
      b[i][j][0][4] =   0.0;
      b[i][j][1][0] = - dt * ty2
	* ( - ( u[i][j-1][k][1]*u[i][j-1][k][2] ) * tmp2 )
	- dt * ty1 * ( - c34 * tmp2 * u[i][j-1][k][1] );
      b[i][j][1][1] = - dt * ty2 * ( u[i][j-1][k][2] * tmp1 )
	- dt * ty1 * ( c34 * tmp1 )
	- dt * ty1 * dy2;
      b[i][j][1][2] = - dt * ty2 * ( u[i][j-1][k][1] * tmp1 );
      b[i][j][1][3] = 0.0;
      b[i][j][1][4] = 0.0;
      b[i][j][2][0] = - dt * ty2
	* ( - ( u[i][j-1][k][2] * tmp1 ) *( u[i][j-1][k][2] * tmp1 )
	    + 0.50 * C2 * ( (  u[i][j-1][k][1] * u[i][j-1][k][1]
			       + u[i][j-1][k][2] * u[i][j-1][k][2]
			       + u[i][j-1][k][3] * u[i][j-1][k][3] )
			    * tmp2 ) )
	- dt * ty1 * ( - r43 * c34 * tmp2 * u[i][j-1][k][2] );
      b[i][j][2][1] = - dt * ty2
	* ( - C2 * ( u[i][j-1][k][1] * tmp1 ) );
      b[i][j][2][2] = - dt * ty2 * ( ( 2.0 - C2 )
				  * ( u[i][j-1][k][2] * tmp1 ) )
	- dt * ty1 * ( r43 * c34 * tmp1 )
	- dt * ty1 * dy3;
      b[i][j][2][3] = - dt * ty2
	* ( - C2 * ( u[i][j-1][k][3] * tmp1 ) );
      b[i][j][2][4] = - dt * ty2 * C2;
      b[i][j][3][0] = - dt * ty2
	* ( - ( u[i][j-1][k][2]*u[i][j-1][k][3] ) * tmp2 )
	- dt * ty1 * ( - c34 * tmp2 * u[i][j-1][k][3] );
      b[i][j][3][1] = 0.0;
      b[i][j][3][2] = - dt * ty2 * ( u[i][j-1][k][3] * tmp1 );
      b[i][j][3][3] = - dt * ty2 * ( u[i][j-1][k][2] * tmp1 )
	- dt * ty1 * ( c34 * tmp1 )
	- dt * ty1 * dy4;
      b[i][j][3][4] = 0.0;
      b[i][j][4][0] = - dt * ty2
	* ( ( C2 * (  u[i][j-1][k][1] * u[i][j-1][k][1]
		      + u[i][j-1][k][2] * u[i][j-1][k][2]
		      + u[i][j-1][k][3] * u[i][j-1][k][3] ) * tmp2
	      - C1 * ( u[i][j-1][k][4] * tmp1 ) )
	    * ( u[i][j-1][k][2] * tmp1 ) )
	- dt * ty1
	* ( - (     c34 - c1345 )*tmp3*(pow2(u[i][j-1][k][1]))
	    - ( r43*c34 - c1345 )*tmp3*(pow2(u[i][j-1][k][2]))
	    - (     c34 - c1345 )*tmp3*(pow2(u[i][j-1][k][3]))
	    - c1345*tmp2*u[i][j-1][k][4] );
      b[i][j][4][1] = - dt * ty2
	* ( - C2 * ( u[i][j-1][k][1]*u[i][j-1][k][2] ) * tmp2 )
	- dt * ty1
	* ( c34 - c1345 ) * tmp2 * u[i][j-1][k][1];
      b[i][j][4][2] = - dt * ty2
	* ( C1 * ( u[i][j-1][k][4] * tmp1 )
	    - 0.50 * C2 
	    * ( (  u[i][j-1][k][1]*u[i][j-1][k][1]
                   + 3.0 * u[i][j-1][k][2]*u[i][j-1][k][2]
		   + u[i][j-1][k][3]*u[i][j-1][k][3] ) * tmp2 ) )
	- dt * ty1
	* ( r43*c34 - c1345 ) * tmp2 * u[i][j-1][k][2];
      b[i][j][4][3] = - dt * ty2
	* ( - C2 * ( u[i][j-1][k][2]*u[i][j-1][k][3] ) * tmp2 )
	- dt * ty1 * ( c34 - c1345 ) * tmp2 * u[i][j-1][k][3];
      b[i][j][4][4] = - dt * ty2
	* ( C1 * ( u[i][j-1][k][2] * tmp1 ) )
	- dt * ty1 * c1345 * tmp1
	- dt * ty1 * dy5;
      tmp1 = 1.0 / u[i-1][j][k][0];
      tmp2 = tmp1 * tmp1;
      tmp3 = tmp1 * tmp2;
      c[i][j][0][0] = - dt * tx1 * dx1;
      c[i][j][0][1] = - dt * tx2;
      c[i][j][0][2] =   0.0;
      c[i][j][0][3] =   0.0;
      c[i][j][0][4] =   0.0;
      c[i][j][1][0] = - dt * tx2
	* ( - ( u[i-1][j][k][1] * tmp1 ) *( u[i-1][j][k][1] * tmp1 )
	    + C2 * 0.50 * (  u[i-1][j][k][1] * u[i-1][j][k][1]
                             + u[i-1][j][k][2] * u[i-1][j][k][2]
                             + u[i-1][j][k][3] * u[i-1][j][k][3] ) * tmp2 )
	- dt * tx1 * ( - r43 * c34 * tmp2 * u[i-1][j][k][1] );
      c[i][j][1][1] = - dt * tx2
	* ( ( 2.0 - C2 ) * ( u[i-1][j][k][1] * tmp1 ) )
	- dt * tx1 * ( r43 * c34 * tmp1 )
	- dt * tx1 * dx2;
      c[i][j][1][2] = - dt * tx2
	* ( - C2 * ( u[i-1][j][k][2] * tmp1 ) );
      c[i][j][1][3] = - dt * tx2
	* ( - C2 * ( u[i-1][j][k][3] * tmp1 ) );
      c[i][j][1][4] = - dt * tx2 * C2;
      c[i][j][2][0] = - dt * tx2
	* ( - ( u[i-1][j][k][1] * u[i-1][j][k][2] ) * tmp2 )
	- dt * tx1 * ( - c34 * tmp2 * u[i-1][j][k][2] );
      c[i][j][2][1] = - dt * tx2 * ( u[i-1][j][k][2] * tmp1 );
      c[i][j][2][2] = - dt * tx2 * ( u[i-1][j][k][1] * tmp1 )
	- dt * tx1 * ( c34 * tmp1 )
	- dt * tx1 * dx3;
      c[i][j][2][3] = 0.0;
      c[i][j][2][4] = 0.0;
      c[i][j][3][0] = - dt * tx2
	* ( - ( u[i-1][j][k][1]*u[i-1][j][k][3] ) * tmp2 )
	- dt * tx1 * ( - c34 * tmp2 * u[i-1][j][k][3] );
      c[i][j][3][1] = - dt * tx2 * ( u[i-1][j][k][3] * tmp1 );
      c[i][j][3][2] = 0.0;
      c[i][j][3][3] = - dt * tx2 * ( u[i-1][j][k][1] * tmp1 )
	- dt * tx1 * ( c34 * tmp1 )
	- dt * tx1 * dx4;
      c[i][j][3][4] = 0.0;
      c[i][j][4][0] = - dt * tx2
	* ( ( C2 * (  u[i-1][j][k][1] * u[i-1][j][k][1]
		      + u[i-1][j][k][2] * u[i-1][j][k][2]
		      + u[i-1][j][k][3] * u[i-1][j][k][3] ) * tmp2
	      - C1 * ( u[i-1][j][k][4] * tmp1 ) )
	    * ( u[i-1][j][k][1] * tmp1 ) )
	- dt * tx1
	* ( - ( r43*c34 - c1345 ) * tmp3 * ( pow2(u[i-1][j][k][1]) )
	    - (     c34 - c1345 ) * tmp3 * ( pow2(u[i-1][j][k][2]) )
	    - (     c34 - c1345 ) * tmp3 * ( pow2(u[i-1][j][k][3]) )
	    - c1345 * tmp2 * u[i-1][j][k][4] );
      c[i][j][4][1] = - dt * tx2
	* ( C1 * ( u[i-1][j][k][4] * tmp1 )
	    - 0.50 * C2
	    * ( (  3.0*u[i-1][j][k][1]*u[i-1][j][k][1]
		   + u[i-1][j][k][2]*u[i-1][j][k][2]
		   + u[i-1][j][k][3]*u[i-1][j][k][3] ) * tmp2 ) )
	- dt * tx1
	* ( r43*c34 - c1345 ) * tmp2 * u[i-1][j][k][1];
      c[i][j][4][2] = - dt * tx2
	* ( - C2 * ( u[i-1][j][k][2]*u[i-1][j][k][1] ) * tmp2 )
	- dt * tx1
	* (  c34 - c1345 ) * tmp2 * u[i-1][j][k][2];
      c[i][j][4][3] = - dt * tx2
	* ( - C2 * ( u[i-1][j][k][3]*u[i-1][j][k][1] ) * tmp2 )
	- dt * tx1
	* (  c34 - c1345 ) * tmp2 * u[i-1][j][k][3];
      c[i][j][4][4] = - dt * tx2
	* ( C1 * ( u[i-1][j][k][1] * tmp1 ) )
	- dt * tx1 * c1345 * tmp1
	- dt * tx1 * dx5;
    }
  }
}