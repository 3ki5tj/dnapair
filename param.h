/* parameters */
#include "util.h"



typedef struct {
  int nargs;
  int np; /* total number of atoms in the two DNA's */
  char fninp[FILENAME_MAX]; /* force file */
  int usemass; /* use mass as weight */
  char fnpsf[FILENAME_MAX]; /* .psf file */
  int docorr; /* compute correlations */
  int verbose; /* verbose level */
  const char *prog; /* name of the program */
} param_t;


#define MAX_OPT_ALIASES 16

/* initialize the default parameters */
static void param_init(param_t *m)
{
  m->nargs = 0;
  m->np = 3798;
  strcpy(m->fninp, "/Bossman/cllai/share/Cheng/sys.160.rot.120.block.60.fout.dat");
  m->usemass = 0;
  strcpy(m->fnpsf, "/Bossman/cllai/share/Cheng/psfpdb/npara.psf");
  m->docorr = 0;
  m->prog = "mf";
  m->verbose = 0;
}



/* clean up */
/* print help message and die */
static void param_help(const param_t *m)
{
  fprintf(stderr, "Parameters\n");
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  %s [Options]\n\n", m->prog);
  fprintf(stderr, "Description:\n");
  fprintf(stderr, "  Computing mean-force.\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  --np=:         set the number of particles, default %d.\n", m->np);
  fprintf(stderr, "  --inp=:        set the input file, default %s.\n", m->fninp);
  fprintf(stderr, "  --usemass:     use mass as weight, default %d.\n", m->usemass);
  fprintf(stderr, "  --psf=:        set the .psf file, default %s.\n", m->fnpsf);
  fprintf(stderr, "  --corr:        compute correlation functions, default %d.\n", m->docorr);
  fprintf(stderr, "  -v:            be verbose, -vv to be more verbose, etc.\n");
  fprintf(stderr, "  -h, --help:    display this message\n");
  exit(1);
}



/* get integer */
static int param_getint(param_t *m,
    const char *key, const char *val)
{
  if ( val == NULL ) {
    fprintf(stderr, "no value for %s\n", key);
    param_help(m);
  }

  /* if 'e' exists in the string, scan the floating point
   * number, and then convert it to an integer */
  if ( strchr(val, 'e') != NULL
    || strchr(val, 'E') != NULL ) {
    return (int) ( atof(val) + 0.5 );
  } else {
    return atoi(val);
  }
}



/* match string key and value pairs */
static int param_keymatch(param_t *m,
    const char *key, const char *val)
{
  if ( strcmp(key, "np") == 0
    || strcmp(key, "n") == 0 )
  {
    m->np = param_getint(m, key, val);
  }
  else if ( strstartswith(key, "inp")
         || strcmp(key, "fninp") == 0 )
  {
    strcpy(m->fninp, val);
  }
  else if ( strcmp(key, "usemass") == 0
    || strcmp(key, "mass") == 0 )
  {
    m->usemass = 1;
  }
  else if ( strcmp(key, "psf") == 0
         || strcmp(key, "fnpsf") == 0 )
  {
    strcpy(m->fnpsf, val);
  }
  else if ( strcmp(key, "corr") == 0 )
  {
    m->docorr = 1;
  }
  else
  {
    return -1;
  }

  return 0;
}



/* scan parameters from command-line arguments */
static int param_doargs(param_t *m, int argc, char **argv)
{
  int i, j, ch;
  char *p, *q;

  m->nargs = 0;
  for ( i = 1; i < argc; i++ ) {
    /* test if it is an argument */
    if ( argv[i][0] != '-' ) {
      m->nargs += 1;
      continue;
    }

    /* long argument, like --help */
    if ( argv[i][1] == '-' ) {
      /* try to parse the argment
         e.g., `--prog=aaa' is parsed to `--prog' and `aaa' */
      p = argv[i] + 2;
      /* let q point to the argument of the option */
      if ( (q = strchr(p, '=')) != NULL ) {
        *q++ = '\0';
      } else {
        q = NULL;
      }

      if ( param_keymatch(m, p, q) != 0 ) {
        if ( strcmpfuzzy(p, "help") != 0 ) {
          fprintf(stderr, "Unknown option %s, key [%s], val [%s]\n",
              argv[i], p, (q != NULL ? q : "NULL") );
        }
        param_help(m);
      }
      continue;
    }

    /* it is a short option
     * loop over characters in the options
     * in this way, `-vo' is understood as `-v -o' */
    for ( j = 1; (ch = argv[i][j]) != '\0'; j++ ) {
      if ( ch == 'v' ) {
        m->verbose++;
      } else if ( ch == 'h' ) {
        param_help(m);
      } else {
        fprintf(stderr, "unknown option %s, j %d, ch %c\n", argv[i], j, ch);
        param_help(m);
      }
    }
  }

  return 0;
}

