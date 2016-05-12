#include "corr.h"


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
    //printf("%d: %g %g %g\n", i, x[i][0], x[i][1], x[i][2]); getchar();
  }

  for ( j = 0; j < 3; j++ ) {
    xc[j] /= mtot;
  }

  printf("%d, mtot %g: %g %g %g\n", n, mtot, xc[0], xc[1], xc[2]);
}



/* compute the radial force of a single frame */
static double calcrf(float (*f)[3], int np)
{
  int i, ns = np / 2;
  double fr;

  /* compute the mean force in frame ifr */
  fr = 0;
  for ( i = 0; i < np; i++ ) {
    if ( i < ns ) {
      fr -= f[i][0];
    } else {
      fr += f[i][0];
    }
  }
  fr /= 2;
  return fr;
}


/* compute the torque of a single frame */
static double calctorq(double (*x)[3], float (*f)[3],
    int np, double xc[2][3])
{
  int i, sig, sgn, ns = np / 2;
  double dx[3], torq;

  /* compute the mean force in frame ifr */
  torq = 0;
  for ( i = 0; i < np; i++ ) {
    if ( i < ns ) {
      /* flip the sign of the weight for the first half */
      sig = 0;
      sgn = -1;
    } else {
      sig = 1;
      sgn = 1;
    }
    dx[0] = x[i][0] - xc[sig][0];
    dx[1] = x[i][1] - xc[sig][1];
    torq += sgn * (-dx[1] * f[i][0] + dx[0] * f[i][1]);
  }
  torq /= 2;
  return torq;
}



/* compute the mean force */
__inline static int calcmf(xf_t *xf, const double *m,
    double *mf, double *var, int docorr)
{
  int k, ifr, nfr = xf->nfr, np = xf->np, ns = np / 2;
  double xc[2][3], f[2], std[2];
  corr_t *corr = NULL;

  if ( docorr ) {
    corr = corr_open(2, nfr);
  }

  /* compute the two centers of mass */
  calccom(xf->x, m, ns, xc[0]);
  calccom(xf->x + ns, m, ns, xc[1]);

  mf[0] = 0;
  mf[1] = 0;
  for ( ifr = 0; ifr < nfr; ifr++ ) {
    /* compute the mean force in frame ifr */
    f[0] = calcrf(xf->f + ifr * np, np);
    f[1] = calctorq(xf->x, xf->f + ifr * np, np, xc);
    if ( docorr ) {
      corr_add(corr, f);
    }
    for ( k = 0; k < 2; k++ ) {
      mf[k] += f[k];
      var[k] += f[k] * f[k];
    }
  }
  for ( k = 0; k < 2; k++ ) {
    mf[k] /= nfr;
    var[k] = var[k] / nfr - mf[k] * mf[k];
    std[k] = sqrt( var[k] );
  }
  printf("fr %g, %g | fa %g, %g | nfr %d\n",
      mf[0], std[0], mf[1], std[1], nfr);

  if ( docorr ) {
    corr_save(corr, 1, 10, 1e-2, 1, "corr.dat");
    corr_printfluc(corr, 1, NULL);
    corr_close(corr);
  }
  return 0;
}



