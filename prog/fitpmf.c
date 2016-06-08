#include "util.h"
#include "lu.h"



char *fninp = "mfrt.log";
char *fnout = "mfrt.out";
double angsp = M_PI / 3;


typedef struct {
  double dis, ang;
  double fr, ft;
  double errfr, errft;
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
    rt[i].errfr = frstd * frstd / (frcnt - 1);
    rt[i].errft = ftstd * ftstd / (ftcnt - 1);
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
  double *invmat, *frc, *ftc, *mm, *mr, *mt;

  xnew(mat, n * n);
  xnew(invmat, n * n);
  xnew(frc, n * n);
  xnew(ftc, n * n);
  xnew(mm, n * n);
  xnew(mr, n * n);
  xnew(mt, n * n);
  xnew(x, n);

  for ( i = 0; i < n * n; i++ ) mat[i] = 0;
  for ( i = 0; i < n * n; i++ ) {
    invmat[i] = frc[i] = ftc[i] = 0;
    mm[i] = mr[i] = mt[i] = 0;
  }
  for ( i = 0; i < n; i++ ) x[i] = 0;

  for ( i = 0; i < n; i++ ) {
    for ( k = 0; k < rt[i].nbcnt; k++ ) {
      j = rt[i].nb[k];
      dx = rt[i].dis - rt[j].dis;
      dy = fmod(rt[i].ang - rt[j].ang + 5 * M_PI, 2 * M_PI) - M_PI;
      wt = 1 / ( (rt[i].errfr + rt[j].errfr) * dx * dx * 0.25
               + (rt[i].errft + rt[j].errft) * dy * dy * 0.25 );
      mf = 0.5 * (rt[i].fr + rt[j].fr) * dx
         + 0.5 * (rt[i].ft + rt[j].ft) * dy;
      mat[i*n + i] -= wt;
      mat[i*n + j] += wt;
      x[i] += mf * wt;
      frc[i*n + i] += wt * 0.5 * dx;
      frc[i*n + j] += wt * 0.5 * dx;
      ftc[i*n + i] += wt * 0.5 * dy;
      ftc[i*n + j] += wt * 0.5 * dy;
    }
  }

  /* replace the last equation by Sum V_i = 0 */
  i = n - 1;
  for ( j = 0; j < n; j++ ) {
    mat[i * n + j] = 1;
  }
  for ( j = 0; j < n; j++ ) {
    frc[i * n + j] = 0;
    ftc[i * n + j] = 0;
  }
  x[i] = 0;

  /* invmat = mat^{-1} */
  luinv(mat, invmat, n, 1e-14);

  //lusolve(mat, x, n, 1e-14);
  for ( i = 0; i < n; i++ ) {
    rt[i].pmf = 0;
    for ( j = 0; j < n; j++ ) {
      rt[i].pmf += invmat[i*n + j] * x[j];
    }
    //rt[i].pmf = x[i];
  }

  /* computing the variance */
  /* mr = invmat * frc
   * mt = invmat * ftc */
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      double yr = 0, yt = 0;
      for ( k = 0; k < n; k++ ) {
        yr += invmat[i*n + k] * frc[k*n + j];
        yt += invmat[i*n + k] * ftc[k*n + j];
      }
      mr[i*n + j] = yr;
      mt[i*n + j] = yt;
    }
  }

  /* mm = mr <fr^2> mr^T + mt <ft^2> mt^T */
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      double y = 0;
      for ( k = 0; k < n; k++ ) {
        y += mr[i*n + k] * rt[k].errfr * mr[j*n + k]
           + mt[i*n + k] * rt[k].errft * mt[j*n + k];
      }
      mm[i*n + j] = y;
    }
  }

  for ( i = 0; i < n; i++ ) {
    rt[i].pmferr = sqrt( mm[i*n+i] );
  }

  free(mat);
  free(invmat);
  free(frc);
  free(ftc);
  free(mm);
  free(mr);
  free(mt);
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

  fprintf(fp, "# %d\n", cnt);
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
  if ( rtarr == NULL ) exit(1);
  findnb(rtarr, *cnt);
  solvepmf(rtarr, *cnt);
  return rtarr;
}

int main(int argc, char **argv)
{
  rtpoint_t *rtarr;
  int cnt;

  if ( argc > 1 ) {
    fninp = argv[1];
  }

  rtarr = fitpmf(&cnt);
  savepmf(rtarr, cnt, fnout);
  return 0;
}
