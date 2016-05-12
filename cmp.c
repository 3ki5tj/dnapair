/* compare the coordinats from two files */
#include "util.h"
#include "xf.h"
#include "mf.h"
#include "mat.h"


const char *fnin1 = "/Bossman/cllai/share/Cheng/sys.160.rot.60.block.60.fout.dat";
const char *fnin2 = "/Bossman/cllai/share/Cheng/sys.160.rot.120.block.60.fout.dat";
const char *fnpsf = "/Bossman/cllai/share/Cheng/psfpdb/npara.psf";
int ns = 1899; /* number of atoms in single DNA */
int usemass = 0; /* use mass as the weight */


int main(int argc, char **argv)
{
  xf_t *xf1, *xf2;
  double *mass = NULL;
  double xc1[2][3], xc2[2][3];
  double rot[2][3][3], trans[2][3];
  double rmsd1, rmsd2;
  int np = ns * 2, k;

  if ( argc >= 3 ) {
    fnin1 = argv[1];
    fnin2 = argv[2];
  }

  xf1 = xf_open(np);
  xf_load(xf1, fnin1, 1); /* only load one frame */

  xf2 = xf_open(np);
  xf_load(xf2, fnin2, 1);

  if ( usemass ) {
    xnew(mass, np);
    loadmass(fnpsf, mass, np);
    checkmass(mass, np);
  }

  calccom(xf1->x,      mass, ns, xc1[0]);
  calccom(xf1->x + ns, mass, ns, xc1[1]);

  calccom(xf2->x,      mass, ns, xc2[0]);
  calccom(xf2->x + ns, mass, ns, xc2[1]);

  rmsd1 = vrmsd(xf2->x,      NULL, xf1->x,      mass, ns, 0, rot[0], trans[0]);
  rmsd2 = vrmsd(xf2->x + ns, NULL, xf1->x + ns, mass, ns, 0, rot[1], trans[1]);
  printf("rmsd %g, %g\n", rmsd1, rmsd2);

  for ( k = 0; k < 2; k++ ) {
    printf("trans [%d]: %10.5f %10.5f %10.5f\n\n",
        k, trans[k][0], trans[k][1], trans[k][2]);
    printf("rot   [%d]: %10.5f %10.5f %10.5f\n"
            "           %10.5f %10.5f %10.5f\n"
            "           %10.5f %10.5f %10.5f\n\n",
        k,
        rot[k][0][0], rot[k][0][1], rot[k][0][2],
        rot[k][1][0], rot[k][1][1], rot[k][1][2],
        rot[k][2][0], rot[k][2][1], rot[k][2][2]);
  }

  xf_close(xf1);
  xf_close(xf2);
  if ( usemass ) {
    free(mass);
  }
  return 0;
}
