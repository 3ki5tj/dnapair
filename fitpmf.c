#include "util.h"
#include "lu.h"


char *fninp = "mfrt.log";
char *fnout = "mfrt.out";
double angsp = M_PI / 3;


typedef struct {
  double dis, ang;
  double fr, ft;
  double varfr, varft;
  double pmf;
  int nbcnt;
  int nb[4];
  int nbtype[4];
} rtpoint_t;


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
      rt[i].nbtype[ rt[i].nbcnt ] = 0;
      rt[i].nbcnt += 1;
    }
    j2 = locate(rt, cnt, dis + 0.5, ang);
    if ( j2 >= 0 ) {
      rt[i].nb[ rt[i].nbcnt ] = j2;
      rt[i].nbtype[ rt[i].nbcnt ] = 0;
      rt[i].nbcnt += 1;
    }
    j3 = locate(rt, cnt, dis, ang - angsp);
    if ( j3 >= 0 ) {
      rt[i].nb[ rt[i].nbcnt ] = j3;
      rt[i].nbtype[ rt[i].nbcnt ] = 1;
      rt[i].nbcnt += 1;
    }
    j4 = locate(rt, cnt, dis, ang + angsp);
    if ( j4 >= 0 ) {
      rt[i].nb[ rt[i].nbcnt ] = j4;
      rt[i].nbtype[ rt[i].nbcnt ] = 1;
      rt[i].nbcnt += 1;
    }
    //printf("point %d %g, %g has %d neighbors, %d, %d, %d, %d\n",
    //    i, dis, ang*180/M_PI, rt[i].nbcnt, j1, j2, j3, j4);
    //getchar();
  }
}



static void solvepmf(rtpoint_t *rt, int n)
{
  int i, j, k;
  double *mat, *x, wt, dx, mf;


  xnew(mat, n * n);
  xnew(x, n);

  for ( i = 0; i < n * n; i++ ) mat[i] = 0;
  for ( i = 0; i < n; i++ ) x[i] = 0;

  for ( i = 0; i < n; i++ ) {
    for ( k = 0; k < rt[i].nbcnt; k++ ) {
      j = rt[i].nb[k];
      if ( rt[i].nbtype[k] == 0 ) {
        wt = 1.0 / (rt[i].varfr + rt[j].varfr);
        dx = rt[i].dis - rt[j].dis;
        mf = 0.5 * (rt[i].fr + rt[j].fr);
      } else {
        wt = 1.0 / (rt[i].varft + rt[j].varft);
        dx = fmod(rt[i].ang - rt[j].ang + 5 * M_PI, 2 * M_PI) - M_PI;
        printf("i %2d, j %2d (nb %d), dx %8.3f, ang %8.3f vs %8.3f, dis %8.3f vs %8.3f\n",
            i, j, k, dx*180/M_PI, rt[i].ang*180/M_PI, rt[j].ang*180/M_PI,
            rt[i].dis, rt[j].dis);
        mf = 0.5 * (rt[i].ft + rt[j].ft);
      }
      mat[i*n + i] += wt;
      mat[i*n + j] -= wt;
      x[i] -= mf * dx * wt;
    }
  }

  /* replace the last equation */
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
  int i;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }

  fprintf(fp, "# %d\n", cnt);
  for ( i = 0; i < cnt; i++ ) {
    fprintf(fp, "%g %g %g %g %g\n",
        rt[i].dis, rt[i].ang, rt[i].pmf,
        rt[i].fr, rt[i].ft);
    if ( i < cnt - 1 && fabs(rt[i].dis - rt[i+1].dis) > 0.05 ) {
      fprintf(fp, "\n");
    }
  }
  fclose(fp);

  return 0;
}



static void fitpmf(void)
{
  rtpoint_t *rtarr;
  int cnt;

  rtarr = loaddata(fninp, &cnt);
  findnb(rtarr, cnt);
  solvepmf(rtarr, cnt);
  savepmf(rtarr, cnt, fnout);
}


int main(int argc, char **argv)
{
  if ( argc > 1 ) {
    fninp = argv[1];
  }

  fitpmf();
  return 0;
}
