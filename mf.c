#include "util.h"
#include "mat.h"
#include "corr.h"
#include "param.h"
#include "xf.h"
#include "mf.h"



/* rotate and translate the first helix to the second */
static void rottrans(xf_t *xf, double *mass, int ns)
{
  double xc[2][3], rot[3][3], trans[3], rmsd;

  calccom(xf->x,      mass, ns, xc[0]);
  calccom(xf->x + ns, mass, ns, xc[1]);

  rmsd = vrmsd(xf->x, NULL, xf->x + ns, mass, ns, 0, rot, trans);

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



static int do_mf(param_t *par, int argc, char **argv)
{
  xf_t *xf;
  double *mass = NULL, mf[3];
  int np = par->np;

  xf = xf_open(np);
  if ( par->nargs == 0 ) {
    xf_load(xf, par->fninp, 0);
  } else {
    int i;
    /* load all files in the commandline argument */
    for ( i = 1; i < argc; i++ ) {
      if ( argv[i][0] == '-' ) break;
      xf_load(xf, argv[i], 0);
    }
  }

  if ( xf->nfr <= 0 ) {
    xf_close(xf);
    return -1;
  }

  if ( par->usemass ) { /* use the mass as the weight */
    xnew(mass, np);
    loadmass(par->fnpsf, mass, np);
    checkmass(mass, np);
  }
  /* compare the geometry of the two helices */
  rottrans(xf, mass, np / 2);

  /* compute the mean force */
  calcmf(xf, mass, mf, par->docorr);

  xf_close(xf);
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
  param_finish(par);
  return 0;
}
