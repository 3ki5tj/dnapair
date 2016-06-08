#ifndef COM_H__
#define COM_H__

/* common routines */



/* load mass from .psf file */
static int loadmass(const char *fnpsf,
    double *mass, int n)
{
  FILE *fp;
  char s[512];
  int i;

  if ( (fp = fopen(fnpsf, "r")) == NULL ) {
    fprintf(stderr, "cannot read %s\n", fnpsf);
    return -1;
  }

  /* search the atom information
   * by looking for the "!NATOM" */
  while ( fgets(s, sizeof s, fp) ) {
    if ( strstr(s, "!NATOM") != NULL )
      break;
  }

  for ( i = 0; i < n; i++ ) {
    char tok[7][128];

    if ( fgets(s, sizeof s, fp) == NULL ) {
      fprintf(stderr, "%s: corrupted in scanning atom %d\n", fnpsf, i);
      fclose(fp);
      return -1;
    }

    sscanf(s, "%s%s%s%s%s%s%s%lf",
        tok[0], tok[1], tok[2], tok[3],
        tok[4], tok[5], tok[6], &mass[i]);
  }

  fclose(fp);
  return 0;
}


/* check the mass of two helices are the same */
static int checkmass(const double *mass, int np)
{
  int i, ns = np/2;

  for ( i = 0; i < ns; i++ ) {
    if ( fabs(mass[i] - mass[i+ns]) > 0.001 ) {
      fprintf(stderr, "mass %d != mass %d, %g vs. %g\n",
          i, i + ns, mass[i], mass[i+ns]);
      return -1;
    }
  }
  fprintf(stderr, "mass is ok!\n");
  return 0;
}



/* compute the total mass */
__inline static double getmtot(const double *m, int n)
{
  int i;
  double mtot;

  if ( m == NULL ) return n;
  mtot = 0;
  for ( i = 0; i < n; i++ )
    mtot += m[i];
  return mtot;
}



/* compute the center of mass */
static void calccom(double (*x)[3],
    const double *m, int n, double xc[3])
{
  int i, j;
  double mi, mtot = 0;

  for ( j = 0; j < 3; j++ ) {
    xc[j] = 0;
  }

  for ( i = 0; i < n; i++ ) {
    mi = (m != NULL) ? m[i] : 1.0;
    for ( j = 0; j < 3; j++ ) {
      xc[j] += x[i][j] * mi;
    }
    mtot += mi;
  }

  for ( j = 0; j < 3; j++ ) {
    xc[j] /= mtot;
  }

  printf("%d, mtot %g: %g %g %g\n", n, mtot, xc[0], xc[1], xc[2]);
}





#endif /* COM_H__ */
