/* compare the coordinats from two files */
#include "util.h"
#include "mat.h"


const char *fnin = "../test/fixed.pdb";
int ns = 1899; /* number of atoms in single DNA */


static vct *loadpdb(const char *fn, int np)
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

    sscanf(s + 30, "%lf%lf%lf", &r[0], &r[1], &r[2]);
    x[id][0] = r[0];
    x[id][1] = r[1];
    x[id][2] = r[2];
    id += 1;
    if ( id >= np ) break;
  }

  fclose(fp);

  return x;
}

int main(int argc, char **argv)
{
  double (*x)[3], *mass = NULL;
  double rot[3][3], trans[3];
  double rmsd;
  int k;

  if ( argc >= 2 ) {
    fnin = argv[1];
  }

  x = loadpdb(fnin, ns * 2);
  rmsd = vrmsd(x + ns, NULL, x, mass, ns, 0, rot, trans);
  printf("RMSD is %g\n", rmsd);

  printf("trans : %10.5f %10.5f %10.5f\n\n",
      trans[0], trans[1], trans[2]);
  printf("rot   : %10.5f %10.5f %10.5f\n"
          "       %10.5f %10.5f %10.5f\n"
          "       %10.5f %10.5f %10.5f\n\n",
      rot[0][0], rot[0][1], rot[0][2],
      rot[1][0], rot[1][1], rot[1][2],
      rot[2][0], rot[2][1], rot[2][2]);

  return 0;
}
