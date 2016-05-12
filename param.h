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



/* compute dependent parameters
 * call this function only once
 * the second call will miss supposedly default parameters */
static void param_compute(param_t *m)
{
}



/* clean up */
static void param_finish(param_t *m)
{
}



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



/* get double */
static double param_getdouble(param_t *m,
    const char *key, const char *val)
{
  if ( val == NULL ) {
    fprintf(stderr, "no value for %s\n", key);
    param_help(m);
  }
  return atof(val);
}



/* select an option */
static int param_selectoption(param_t *m,
    const char *key, const char *val,
    const char *names[][MAX_OPT_ALIASES], int cnt)
{
  int i, j;

  if ( val == NULL ) {
    fprintf(stderr, "no value for %s\n", key);
    param_help(m);
  }

  /* loop over options */
  for ( i = 0; i < cnt; i++ ) {
    /* loop over aliases */
    for ( j = 0; j < MAX_OPT_ALIASES; j++ ) {
      if ( names[i][j] == NULL
        || names[i][j][0] == '\0' ) {
        break;
      }

      if ( strcmpfuzzy(names[i][j], val) == 0 ) {
        return i;
      }
    }
  }

  /* try to treat `val` is a number */
  if ( isdigit(val[0]) ) {
    i = atoi(val);
    if ( i < 0 || i >= cnt ) i = 0;
  } else {
    i = 0;
  }

  return i;
}



/* get a boolean/integer value */
static int param_getbool(param_t *m,
    const char *key, const char *val)
{
  if ( val == NULL || val[0] == '\0' ) {
    return 1;
  }

  if ( isdigit(val[0]) ) {
    return atoi(val);
  }

  if ( strcmpfuzzy(val, "no") == 0
    || strcmpfuzzy(val, "false") == 0
    || strcmpfuzzy(val, "n") == 0
    || strcmpfuzzy(val, "f") == 0 ) {
    return 0;
  } else if ( strcmpfuzzy(val, "yes") == 0
           || strcmpfuzzy(val, "true") == 0
           || strcmpfuzzy(val, "y") == 0
           || strcmpfuzzy(val, "t") == 0 ) {
    return 1;
  } else {
    fprintf(stderr, "unknown value [%s] for %s\n", val, key);
    param_help(m);
  }

  return 1;
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
    m->usemass = param_getbool(m, key, val);
  }
  else if ( strcmp(key, "psf") == 0
         || strcmp(key, "fnpsf") == 0 )
  {
    strcpy(m->fnpsf, val);
  }
  else if ( strcmp(key, "corr") == 0 )
  {
    m->docorr = param_getbool(m, key, val);
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

  param_compute(m);

  return 0;
}

