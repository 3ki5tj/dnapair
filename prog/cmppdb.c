/* compare the coordinates of the two DNAs in a single PDB file */
#include "util.h"
#include "mat.h"
#include "com.h"


const char *fnin = "../test/fixed.pdb";
const char *fnpsf = "../test/ionized.psf";
int ns = 1899; /* number of atoms in single DNA */
int usemass = 1;


static vct *loadpdb(const char *fn, int np, char (*anames)[8])
{
  double (*x)[3], r[3];
  FILE *fp;
  char s[1024];
  int id;

  xnew(x, np);
  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot file %s\n", fn);
    return NULL;
  }

  id = 0;
  while ( fgets(s, sizeof s, fp) ) {
    if ( strncmp(s, "ATOM  ", 6) != 0 ) {
      continue;
    }
    /* copy the atom name */
    strncpy(anames[id], s + 12, 4);
    strstrip(anames[id]);

    sscanf(s + 30, "%lf%lf%lf", &r[0], &r[1], &r[2]);
    x[id][0] = r[0];
    x[id][1] = r[1];
    x[id][2] = r[2];
    //printf("%d %s %g %g %g\n", id, anames[id], x[id][0], x[id][1], x[id][2]); getchar();
    id += 1;
    if ( id >= np ) break;
  }

  fclose(fp);

  return x;
}


int main(int argc, char **argv)
{
  double (*x)[3], *mass = NULL;
  double rot[3][3], trans[3], xc1[3], xc2[3];
  double rmsd;
  int np = ns * 2;
  char (*anames)[8];

  if ( argc >= 2 ) {
    fnin = argv[1];
  }

  xnew(anames, np);
  x = loadpdb(fnin, np, anames);

  if ( usemass ) {
    xnew(mass, np);
    loadmass(fnpsf, mass, np);
    checkmass(mass, np);
  }

  calccom(x, mass, ns, xc1);
  calccom(x + ns, mass + ns, ns, xc2);

  rmsd = vrmsd(x + ns, NULL, x, mass, ns, 0, rot, trans);
  printf("RMSD is %g, COM: %g, %g, %g; %g, %g, %g\n",
      rmsd, xc1[0], xc1[1], xc1[2], xc2[0], xc2[1], xc2[2]);

  printf("trans : %10.5f %10.5f %10.5f\n\n",
      trans[0], trans[1], trans[2]);
  printf("rot   : %10.5f %10.5f %10.5f\n"
          "       %10.5f %10.5f %10.5f\n"
          "       %10.5f %10.5f %10.5f\n\n",
      rot[0][0], rot[0][1], rot[0][2],
      rot[1][0], rot[1][1], rot[1][2],
      rot[2][0], rot[2][1], rot[2][2]);

  if ( usemass ) {
    free(mass);
  }
  return 0;
}
