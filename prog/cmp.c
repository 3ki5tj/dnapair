/* compare the coordinates from two files */
#include "util.h"
#include "xf.h"
#include "mf.h"
#include "mat.h"


const char *fnin1 = "/Bossman/cllai/share/Cheng/sys.160.rot.60.block.60.fout.dat";
const char *fnin2 = "/Bossman/cllai/share/Cheng/sys.160.rot.120.block.60.fout.dat";
const char *fnpsf = "/Bossman/cllai/share/Cheng/psfpdb/npara.psf";
int ns = 1899; /* number of atoms in single DNA */
int usemass = 0; /* use mass as the weight */
int dobb = 1; /* compare backbone atoms only */

/* load atom names from .psf file */
static int loadanames(const char *fnpsf,
    char (*anames)[8], int n)
{
  FILE *fp;
  char s[512];
  int i;

  if ( (fp = fopen(fnpsf, "r")) == NULL ) {
    fprintf(stderr, "cannot read %s\n", fnpsf);
    return -1;
  }

  /* search the atom information
   * by looking for the "!NATOM" */
  while ( fgets(s, sizeof s, fp) ) {
    if ( strstr(s, "!NATOM") != NULL )
      break;
  }

  for ( i = 0; i < n; i++ ) {
    char tok[7][128];

    if ( fgets(s, sizeof s, fp) == NULL ) {
      fprintf(stderr, "%s: corrupted in scanning atom %d\n", fnpsf, i);
      fclose(fp);
      return -1;
    }

    sscanf(s, "%s%s%s%s%s",
        tok[0], tok[1], tok[2], tok[3], anames[i]);
  }

  fclose(fp);
  return 0;
}



static int countbackbone(char (*anames)[8], int np, int *isbb)
{
  int i, j, nbb = 0;
  char *bbanames[] = {
    "P",
    "O1P",
    "O2P",
    "O5\'",
    "C5\'",
    "C4\'",
    "C3\'",
    "O3\'",
    ""};


  for ( i = 0; i < np; i++ ) {
    for ( j = 0; bbanames[j][0] != '\0'; j++ )
      if ( strcmp(anames[i], bbanames[j]) == 0 )
        break;
    isbb[i] = ( bbanames[j][0] != '\0' );
    nbb += isbb[i];
    //if(isbb[i]) printf("i %d, atom %s is on the backbone\n", i, anames[i]);
  }
  return nbb;
}


vct *filterbb(vct *x, int np, const int *isbb, int nbb,
    const double *mass, double *mbb)
{
  int i, ibb;
  double (*xbb)[3];

  xnew(xbb, nbb);
  for ( i = ibb = 0; i < np; i++ ) {
    if ( isbb[i] ) {
      vcopy(xbb[ibb], x[i]);
      if ( mbb != NULL && mass != NULL )
        mbb[ibb] = mass[i];
      ibb++;
    }
  }
  return xbb;
}

int main(int argc, char **argv)
{
  xf_t *xf1, *xf2;
  double *mass = NULL;
  double xc1[2][3], xc2[2][3];
  double rot[4][3][3], trans[4][3];
  double rmsd1, rmsd2, rmsd3, rmsd4;
  int np = ns * 2, k;
  char (*anames)[8];
  int *isbb, nbb = 0, nbbs;
  double (*xbb1)[3], (*xbb2)[3], *mbb = NULL;

  if ( argc >= 3 ) {
    fnin1 = argv[1];
    fnin2 = argv[2];
  }
  if ( argc > 3 ) {
    fnpsf = argv[3];
  }

  xf1 = xf_open(np, 1);
  xf_load(xf1, fnin1, 1); /* only load one frame */

  xf2 = xf_open(np, 1);
  xf_load(xf2, fnin2, 1);

  if ( usemass ) {
    xnew(mass, np);
    loadmass(fnpsf, mass, np);
    checkmass(mass, np);
  }

  xnew(anames, np);
  loadanames(fnpsf, anames, np);
  /* count the number of backbone atoms */
  xnew(isbb, np);
  nbb = countbackbone(anames, np, isbb);
  nbbs = nbb / 2;
  printf("%d backbone atoms out of %d; %d out of %d\n", nbb, np, nbbs, ns);

  xbb1 = filterbb(xf1->x, np, isbb, nbb, mass, mbb);
  xbb2 = filterbb(xf2->x, np, isbb, nbb, NULL, NULL);

  if ( dobb ) {
    rmsd1 = vrmsd(xbb2,        NULL, xbb1,        mbb, nbbs, 0, rot[0], trans[0]);
    rmsd2 = vrmsd(xbb2 + nbbs, NULL, xbb1 + nbbs, mbb, nbbs, 0, rot[1], trans[1]);
    rmsd3 = vrmsd(xbb1,        NULL, xbb1 + nbbs, mbb, nbbs, 0, rot[2], trans[2]);
    rmsd4 = vrmsd(xbb2,        NULL, xbb2 + nbbs, mbb, nbbs, 0, rot[3], trans[3]);
  } else {
    calccom(xf1->x,      mass, ns, xc1[0]);
    calccom(xf1->x + ns, mass, ns, xc1[1]);

    calccom(xf2->x,      mass, ns, xc2[0]);
    calccom(xf2->x + ns, mass, ns, xc2[1]);

    rmsd1 = vrmsd(xf2->x,      NULL, xf1->x,      mass, ns, 0, rot[0], trans[0]);
    rmsd2 = vrmsd(xf2->x + ns, NULL, xf1->x + ns, mass, ns, 0, rot[1], trans[1]);
    rmsd3 = vrmsd(xf1->x,      NULL, xf1->x + ns, mass, ns, 0, rot[2], trans[2]);
    rmsd4 = vrmsd(xf2->x,      NULL, xf2->x + ns, mass, ns, 0, rot[3], trans[3]);
  }

  printf("rmsd %g, %g; intra %g, %g\n", rmsd1, rmsd2, rmsd3, rmsd4);
  for ( k = 0; k < 4; k++ ) {
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
  free(anames);
  return 0;
}
