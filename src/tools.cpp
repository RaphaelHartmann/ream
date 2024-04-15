
#include "tools.h"


double unif_L() {

  double u;
  do {
    u = unif_rand();
  } while (u < 0 || u >= 1);
  return u;

}
