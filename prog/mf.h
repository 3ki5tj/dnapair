#include "com.h"



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


/* compute the torque of a single frame
 * the default (return value) contains the total torque
 * on the second DNA molecule
 * for reference, we also compute the averaged torque, *symmtorq
 * on the two DNA molecules */
static double calctorq(double (*x)[3], float (*f)[3],
    int np, double xc[2][3], double *symmtorq)
{
  int i, who, sgn, ns = np / 2;
  double dx[3], t, torq;

  /* compute the mean force in frame ifr */
  torq = 0;
  *symmtorq = 0;
  for ( i = 0; i < np; i++ ) {
    if ( i < ns ) {
      /* flip the sign of the weight for the first half */
      who = 0;
      sgn = -1;
    } else {
      who = 1;
      sgn = 1;
    }
    dx[0] = x[i][0] - xc[who][0];
    dx[1] = x[i][1] - xc[who][1];
    t = sgn * (-dx[1] * f[i][0] + dx[0] * f[i][1]);
    *symmtorq += t;
    if ( who == 1 ) {
      torq += t;
    }
  }
  *symmtorq /= 2;
  return torq;
}



#define MFCNT 3
/* load coordinates and compute mean force in the same time
 * sums[0]: 0 counts, 1 sum, 2 square sum of the radial force
 * sums[1]: those for the angular torque
 * */
__inline static int calcmf_inplace(xf_t *xf, const char *fn,
    const double *mass, double sums[MFCNT][3])
{
  FILE *fp;
  char s[256];
  int i, np = xf->np, ns = np / 2;
  clock_t starttime = clock();
  double nfr0 = sums[0][0], f[MFCNT], xc[2][3];
  int k, once = 0;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "calcmf_inplace: cannot open %s\n", fn);
    return -1;
  }

  for ( ; ; ) {
    if ( fgets(s, sizeof s, fp) == NULL
      || strncmp(s, "timestep", 8) != 0 ) {
      break;
    }

    /* read in the coordinates */
    for ( i = 0; i < np; i++ ) {
      char tok[4][128];

      if ( fgets(s, sizeof s, fp) == NULL ) {
        fprintf(stderr, "cannot read frame %d from %s\n",
            xf->nfr, fn);
        break;
      }
      /* for the coordinates, we only scan them as strings
       * without converting them to real numbers
       * This helps saving time */
      sscanf(s, "%s%s%s%s%f%f%f", tok[0], tok[1], tok[2], tok[3],
          &xf->f[i][0], &xf->f[i][1], &xf->f[i][2]);

      /* save the coordinates only for the first frame
       * because those of a later frame are the same */
      if ( xf->nfr == 0 ) {
        xf->x[i][0] = atof(tok[1]);
        xf->x[i][1] = atof(tok[2]);
        xf->x[i][2] = atof(tok[3]);
      }
    }

    /* compute the centers of mass */
    if ( !once ) {
      calccom(xf->x, mass, ns, xc[0]);
      calccom(xf->x + ns, mass, ns, xc[1]);
      once = 1;
    }

    /* compute the force and torque */
    f[0] = calcrf(xf->f, np);
    f[1] = calctorq(xf->x, xf->f, np, xc, &f[2]);

    for ( k = 0; k < MFCNT; k++ ) {
      sums[k][0] += 1;
      sums[k][1] += f[k];
      sums[k][2] += f[k] * f[k];
    }

    /* something wrong has happened */
    if ( i < np ) break;

    /* after the first frame,
     * no position will be read */
    xf->nfr = 1;
  }
  fprintf(stderr, "loaded %s in %.3f seconds, %g -> %g frames\n",
      fn, (double)(clock() - starttime) / CLOCKS_PER_SEC,
      nfr0, sums[0][0]);

  fclose(fp);
  return 0;
}





