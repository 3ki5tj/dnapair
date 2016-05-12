
const int xf_blksz = 500;
const int xf_nfrx = 1; /* number of x frames to load */


typedef struct {
  double (*x)[3]; /* position, size is np */
  float (*f)[3]; /* force, size is np * nfr */
  int nfr_max; /* capacity */
  int nfr; /* number of frames */
  int np; /* number of atoms in each frame */
} xf_t;



static xf_t *xf_open(int np)
{
  xf_t *xf;

  xnew(xf, 1);
  xf->np = np;
  xf->nfr = 0;
  xf->nfr_max = xf_blksz;
  xnew(xf->x, np * xf_nfrx);
  xnew(xf->f, np * xf_blksz);
  return xf;
}



static void xf_expand(xf_t *xf)
{
  xf->nfr_max += xf_blksz;
  xrenew(xf->f, xf->np * xf->nfr_max);
}



static void xf_close(xf_t *xf)
{
  free(xf->x);
  free(xf->f);
  free(xf);
}



/* load coordinates
 * fn: file name (xxx.fout.dat)
 * nfr_max: maximal number frames to be loaded, 0 to load all
 * */
static int xf_load(xf_t *xf, const char *fn,
    int nfr_max)
{
  FILE *fp;
  char s[256];
  int i, np = xf->np, nfr0 = xf->nfr;
  int id = np * nfr0; /* atom id */
  clock_t starttime = clock();

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }

  for ( ; ; ) {
    if ( fgets(s, sizeof s, fp) == NULL
      || strncmp(s, "timestep", 8) != 0 ) {
      break;
    }

    /* allocate more memory if needed */
    if ( xf->nfr >= xf->nfr_max ) {
      xf_expand(xf);
    }

    /* read in the coordinates */
    for ( i = 0; i < np; i++ ) {
      /* double x[3]; */
      char tok[4][128];

      if ( fgets(s, sizeof s, fp) == NULL ) {
        fprintf(stderr, "cannot read frame %d from %s\n",
            xf->nfr, fn);
        break;
      }
      /*
       sscanf(s, "%lf%lf%lf%lf%lf%lf", &x[0], &x[1], &x[2],
          &xf->f[id][0], &xf->f[id][1], &xf->f[id][2]);
      */
      /* for the coordinates, we only scan them as strings
       * without converting them to real numbers
       * This helps saving time */
      sscanf(s, "%s%s%s%s%f%f%f", tok[0], tok[1], tok[2], tok[3],
          &xf->f[id][0], &xf->f[id][1], &xf->f[id][2]);

      /* save the coordinates only for the first frame
       * because those of a later frame are the same */
      if ( xf->nfr < xf_nfrx ) {
        /*
        xf->x[id][0] = x[0];
        xf->x[id][1] = x[1];
        xf->x[id][2] = x[2];
        */
        xf->x[id][0] = atof(tok[1]);
        xf->x[id][1] = atof(tok[2]);
        xf->x[id][2] = atof(tok[3]);
      }

      id += 1;
    }
    /* something wrong has happened */
    if ( i < np ) break;

    xf->nfr += 1;
    if ( nfr_max > 0 && xf->nfr >= nfr_max )
      break;
  }
  fprintf(stderr, "loaded %s in %.3f seconds, %d -> %d frames\n",
      fn, (double)(clock() - starttime) / CLOCKS_PER_SEC,
      nfr0, xf->nfr);

  fclose(fp);
  return 0;
}



