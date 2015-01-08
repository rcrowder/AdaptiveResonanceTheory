#ifndef _ART_distance_ALG_H
#define _ART_distance_ALG_H


#include "art_common.h"

namespace ART
{
Clust art_distance(const a_instance& sample, in_param& par);


double euclidean(const prototype& prot, const instance& inst, unsigned long dim);
double euclidean_m(const prototype& prot, const instance& inst, unsigned long dim);
double correl(const prototype& prot, const instance& inst, unsigned long dim);
double manhattan(const prototype& prot, const instance& inst, unsigned long dim);
double minkowski(const prototype& prot, const instance& inst, int p, unsigned long dim);
};

#endif

