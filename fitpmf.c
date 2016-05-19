#include "util.h"
#include "lu.h"
#include "mtrand.h"



char *fninp = "mfrt.log";
char *fnout = "mfrt.out";
double angsp = M_PI / 3;
int nerr = 100; /* number of runs to compute error */


typedef struct {
  double dis, ang;
  double fr, ft;
  double varfr, varft;
  double pmf;
  double pmferr;
  int nbcnt;
  int nb[4];
} rtpoint_t;



/* load data from input file,
 * in which each line has,
 * distance, angle, mean force, mean torque */
static rtpoint_t *loaddata(const char *fn, int *cnt)
{
  rtpoint_t *rt = NULL;
  FILE *fp;
  int i, n, err;
  char s[512];

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return NULL;
  }

  /* count the number of lines */
  n = 0;
  while ( fgets(s, sizeof s, fp) ) {
    if ( strncmp(s, "dis", 3) == 0 )
      n += 1;
  }
  rewind(fp);

  /* allocate the array */
  xnew(rt, n);

  /* read each line */
  for ( i = 0; i < n; i++ ) {
    double ang2, rmsd, frstd, frcnt, ftstd, ftcnt;

    fgets(s, sizeof s, fp);
    err = sscanf(s, "dis %lf, ang %lf/%lf, rmsd %lf | "
                    "f %lf %lf %lf | torq %lf %lf %lf",
        &rt[i].dis, &rt[i].ang, &ang2, &rmsd, &rt[i].fr, &frstd, &frcnt,
        &rt[i].ft, &ftstd, &ftcnt);
    if ( rt[i].ang < -0.05 ) {
      rt[i].ang += 2 * M_PI;
    }
    if ( err != 10 ) {
      fprintf(stderr, "error reading line %d/%d, err %d\n", i, n, err);
      exit(1);
    }
    /* convert to error measure */
    rt[i].varfr = frstd * frstd / frcnt;
    rt[i].varft = ftstd * ftstd / ftcnt;
  }

  *cnt = n;

  fclose(fp);
  return rt;
}



/* locate an entry with the appoximate value of dis and ang */
static int locate(rtpoint_t *rt, int cnt,
    double dis, double ang)
{
  int i;

  for ( i = 0; i < cnt; i++ ) {
    if ( fabs(rt[i].dis - dis) < 0.05
      && fabs( fmod(rt[i].ang - ang + 5 * M_PI, 2 * M_PI) - M_PI ) < 0.01 )
      return i;
  }
  //fprintf(stderr, "cannot match dis %g, ang %.0f\n", dis, ang*180/M_PI);
  return -1;
}



/* find neighbors */
static void findnb(rtpoint_t *rt, int cnt)
{
  int i, j1, j2, j3, j4;
  double dis, ang;

  for ( i = 0; i < cnt; i++ ) {
    dis = rt[i].dis;
    ang = rt[i].ang;
    j1 = locate(rt, cnt, dis - 0.5, ang);
    if ( j1 >= 0 ) {
      rt[i].nb[ rt[i].nbcnt ] = j1;
      rt[i].nbcnt += 1;
    }
    j2 = locate(rt, cnt, dis + 0.5, ang);
    if ( j2 >= 0 ) {
      rt[i].nb[ rt[i].nbcnt ] = j2;
      rt[i].nbcnt += 1;
    }
    j3 = locate(rt, cnt, dis, ang - angsp);
    if ( j3 >= 0 ) {
      rt[i].nb[ rt[i].nbcnt ] = j3;
      rt[i].nbcnt += 1;
    }
    j4 = locate(rt, cnt, dis, ang + angsp);
    if ( j4 >= 0 ) {
      rt[i].nb[ rt[i].nbcnt ] = j4;
      rt[i].nbcnt += 1;
    }
    //printf("point %d %g, %g has %d neighbors, %d, %d, %d, %d\n",
    //    i, dis, ang*180/M_PI, rt[i].nbcnt, j1, j2, j3, j4);
    //getchar();
  }
}



/* solve the PMF by minimizing the error */
static void solvepmf(rtpoint_t *rt, int n)
{
  int i, j, k;
  double *mat, *x, wt, dx, dy, mf;


  xnew(mat, n * n);
  xnew(x, n);

  for ( i = 0; i < n * n; i++ ) mat[i] = 0;
  for ( i = 0; i < n; i++ ) x[i] = 0;

  for ( i = 0; i < n; i++ ) {
    for ( k = 0; k < rt[i].nbcnt; k++ ) {
      j = rt[i].nb[k];
      dx = rt[i].dis - rt[j].dis;
      dy = fmod(rt[i].ang - rt[j].ang + 5 * M_PI, 2 * M_PI) - M_PI;
      wt = 1 / ( (rt[i].varfr + rt[j].varfr) * dx * dx
               + (rt[i].varft + rt[j].varft) * dy * dy );
      mf = 0.5 * (rt[i].fr + rt[j].fr) * dx
         + 0.5 * (rt[i].ft + rt[j].ft) * dy;
      mat[i*n + i] += wt;
      mat[i*n + j] -= wt;
      x[i] -= mf * wt;
    }
  }

  /* replace the last equation by Sum V_i = 0 */
  i = n - 1;
  for ( j = 0; j < n; j++ ) {
    mat[i * n + j] = 1;
  }
  x[i] = 0;

  lusolve(mat, x, n, 1e-14);
  for ( i = 0; i < n; i++ ) {
    rt[i].pmf = x[i];
  }

  free(mat);
  free(x);
}


static int savepmf(rtpoint_t *rt, int cnt, const char *fn)
{
  FILE *fp;
  int i, i0 = 0;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }

  fprintf(fp, "# %d %d\n", cnt, nerr);
  for ( i = 0; i < cnt; i++ ) {
    fprintf(fp, "%g %g %g %g %g %g\n",
        rt[i].dis, rt[i].ang, rt[i].pmf, rt[i].pmferr,
        rt[i].fr, rt[i].ft);
    if ( i == cnt - 1 || fabs(rt[i].dis - rt[i+1].dis) > 0.05 ) {
      fprintf(fp, "%g %g %g %g %g %g\n",
          rt[i0].dis, rt[i0].ang + 2*M_PI,
          rt[i0].pmf, rt[i0].pmferr,
          rt[i0].fr, rt[i0].ft);
      if ( i < cnt - 1 ) {
        i0 = i + 1;
        fprintf(fp, "\n");
      }
    }
  }
  fclose(fp);

  fprintf(stderr, "successfully saved PMF to %s\n", fn);

  return 0;
}



/* compute the mean force */
static rtpoint_t *fitpmf(int *cnt)
{
  rtpoint_t *rtarr;

  rtarr = loaddata(fninp, cnt);
  findnb(rtarr, *cnt);
  solvepmf(rtarr, *cnt);
  return rtarr;
}


/* vary the mean forces and recompute the pmf */
static rtpoint_t *varpmf(rtpoint_t *rtarr0, int cnt)
{
  rtpoint_t *rtarr;
  int i;

  xnew(rtarr, cnt);
  for ( i = 0; i < cnt; i++ ) {
    memcpy(rtarr + i, rtarr0 + i, sizeof(rtarr[0]));
    rtarr[i].fr += randgaus() * sqrt( rtarr[i].varfr );
    rtarr[i].ft += randgaus() * sqrt( rtarr[i].varft );
  }
  solvepmf(rtarr, cnt);
  return rtarr;
}


int main(int argc, char **argv)
{
  rtpoint_t *rtarr, *rtarr0;
  double *x2, x;
  int cnt, i, k;

  if ( argc > 1 ) {
    fninp = argv[1];
  }

  if ( argc > 2 ) {
    nerr = atoi( argv[2] );
  }

  rtarr0 = fitpmf(&cnt);

  if ( nerr > 0 ) {
    xnew(x2, cnt);
    for ( k = 0; k < nerr; k++ ) {
      rtarr = varpmf(rtarr0, cnt);
      for ( i = 0; i < cnt; i++ ) {
        x = rtarr[i].pmf - rtarr0[i].pmf;
        x2[i] += x * x;
      }
      free(rtarr);
    }
    for ( i = 0; i < cnt; i++ ) {
      x2[i] /= k;
      rtarr0[i].pmferr = sqrt( x2[i] );
    }
    free(x2);
  }
  savepmf(rtarr0, cnt, fnout);
  return 0;
}
