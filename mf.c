#include "util.h"
#include "mat.h"
#include "corr.h"
#include "param.h"
#include "xf.h"
#include "mf.h"



/* rotate and translate the first helix to the second */
static double rottrans(double (*x)[3], const double *mass,
    int ns, double *ang, double *prmsd, int verbose)
{
  double xc[2][3], rot[3][3], trans[3], rmsd;

  calccom(x,      mass, ns, xc[0]);
  calccom(x + ns, mass, ns, xc[1]);

  rmsd = vrmsd(x, NULL, x + ns, mass, ns, 0, rot, trans);

  if ( verbose ) {
    printf("rmsd %g\n", rmsd);
    printf("trans : %10.5f %10.5f %10.5f\n\n",
          trans[0], trans[1], trans[2]);
    printf("rot   : %10.5f %10.5f %10.5f\n"
           "        %10.5f %10.5f %10.5f\n"
           "        %10.5f %10.5f %10.5f\n\n",
          rot[0][0], rot[0][1], rot[0][2],
          rot[1][0], rot[1][1], rot[1][2],
          rot[2][0], rot[2][1], rot[2][2]);
  }

  if ( ang != NULL ) {
    *ang = atan2(rot[1][0] - rot[0][1], rot[0][0] + rot[1][1]);
  }

  if ( prmsd != NULL ) {
    *prmsd = rmsd;
  }

  return trans[0];
}


/* compute the mean force for the list */
static void mf_dolist(xf_t *xf, char **fns, int cnt,
    const double *mass)
{
  double dis = 0, ang = 0, rmsd = 0;
  double sums[2][3] = {{0, 0, 0}, {0, 0, 0}}, ave[2], std[2];
  int i, j, once = 0, np = xf->np;

  /* load all files in the commandline argument */
  for ( i = 0; i < cnt; i++ ) {
    /* compute the mean force as we scan the data */
    calcmf_inplace(xf, fns[i], mass, sums);

    /* convert the sums to averages and standard deviations */
    for ( j = 0; j < 2; j++ ) {
      ave[j] = sums[j][1] / sums[j][0];
      std[j] = sqrt(sums[j][2] / sums[j][0] - ave[j] * ave[j]);
    }

    if ( !once ) {
      /* compare the geometry of the two helices */
      dis = rottrans(xf->x, mass, np / 2, &ang, &rmsd, 0);
      once = 1;
    }
    printf("dis %g, ang %g/%g, rmsd %g | f %g %g %g | torq %g %g %g\n",
        dis, ang, ang * 180 / M_PI, rmsd,
        ave[0], std[0], sums[0][0],
        ave[1], std[1], sums[1][0]);
  }
}


/* get the file pattern of the directory
 * return the maximal number of blocks */
static int getpat(char *dir, char *head, char *tail)
{
  char fnout[FILENAME_MAX], s[FILENAME_MAX];
  char cmd[FILENAME_MAX*2], *p, *q;
  FILE *fp;
  int i, block = 0;

  tmpnam(fnout);
  sprintf(cmd, "/bin/ls %s/*%s > %s", dir, tail, fnout);
  system(cmd);
  if ( (fp = fopen(fnout, "r")) == NULL ) {
    fprintf(stderr, "cannot open %s for the ls result\n", fnout);
    return -1;
  }
  if ( fgets(s, sizeof s, fp) == NULL ) {
    fclose(fp);
    return -1;
  }
  p = strstr(s, "block.");
  p[0] = '\0'; /* terminate the string */
  strcpy(head, s);

  /* count the number of files */
  rewind(fp);
  while ( fgets(s, sizeof s, fp) ) {
    p = strstr(s, "block.");
    if ( p == NULL ) continue;
    p += 6; /* skip "block." */
    q = strchr(p, '.');
    if ( q == NULL ) continue;
    *q = '\0';
    i = atoi(p);
    if ( i > block ) {
      block = i;
    }
  }

  printf("block %d, head %s\n", block, head);
  fclose(fp);
  remove(fnout);
  return block;
}



/* get a list of force files under the directory */
static char **getlist(param_t *par, int *cnt)
{
  char **fns, *fn;
  char head[FILENAME_MAX] = "sys.160.rot.60.";
  char tail[FILENAME_MAX] = ".fout.dat";
  int i, blkmax;
  FILE *fp;

  /* 1. refine the data directory */
  if ( par->dir[0] == '\0' ) {
    strcpy(par->dir, ".");
  } else {
    i = strlen(par->dir) - 1;
    /* remove the trailing slash, if any */
    if ( par->dir[i] == '/' && i != 0 ) {
      par->dir[i] = '\0';
    }
  }

  /* 2. try to get the pattern of data file */
  blkmax = getpat(par->dir, head, tail);

  /* 3. construct a list of file names */
  xnew(fns, blkmax);
  *cnt = 0;
  for ( i = 0; i < blkmax; i++ ) {
    xnew(fn, FILENAME_MAX);
    sprintf(fn, "%sblock.%d%s", head, i + 1, tail);
    /* try to check if the file exists */
    if ( (fp = fopen(fn, "r")) != NULL ) {
      fclose(fp);
      fns[ (*cnt)++ ] = fn;
      //printf("i %d, cnt %d, %s\n", i + 1, *cnt, fn);
    } else {
      continue;
    }
  }

  return fns;
}


/* scan all force files under the directory */
static int mfscan(param_t *par, const double *mass)
{
  int np = par->np, cnt;
  xf_t *xf;
  char **fns;

  fns = getlist(par, &cnt);
  if ( fns == NULL ) {
    return -1;
  }
  xf = xf_open(np, 1);
  mf_dolist(xf, fns, cnt, mass);
  xf_close(xf);

  return 0;
}




static int do_mf(param_t *par, int argc, char **argv)
{
  xf_t *xf;
  double *mass = NULL, mf[2], var[2];
  double dis, ang, rmsd;
  int i, np = par->np;

  if ( par->usemass ) { /* use the mass as the weight */
    xnew(mass, np);
    loadmass(par->fnpsf, mass, np);
    checkmass(mass, np);
  }

  if ( par->scanf )
  {
    fprintf(stderr, "scanning directory [%s]\n", par->dir);
    mfscan(par, mass);
  }
  else
  {
    xf = xf_open(np, 500);

    if ( par->nargs == 0 ) {
      /* no argument use the default input file
       * this "if" branch for testing */
      xf_load(xf, par->fninp, 0);

      if ( xf->nfr > 0 ) {
        /* compare the geometry of the two helices */
        dis = rottrans(xf->x, mass, np / 2, &ang, &rmsd, 1);
        printf("dis %g, ang %g/%g, rmsd %g\n", dis, ang, ang * 180 / M_PI, rmsd);

        /* compute the mean force */
        calcmf(xf, mass, mf, var, par->docorr);
      }

    } else {
      char **fns = NULL;
      int cnt = 0;

      xnew(fns, argc);
      for ( i = 1; i < argc; i++ ) {
        if ( argv[i][0] != '-' ) { /* not an option */
          fns[cnt++] = argv[i];
        }
      }
      mf_dolist(xf, fns, cnt, mass);
      free(fns);
    }

    xf_close(xf);
  }

  if ( par->usemass ) {
    free(mass);
  }
  return 0;
}


int main(int argc, char **argv)
{
  param_t par[1];

  param_init(par);
  param_doargs(par, argc, argv);
  do_mf(par, argc, argv);
  return 0;
}
