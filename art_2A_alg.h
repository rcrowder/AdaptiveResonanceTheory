#ifndef _ART_2A_ALG_H
#define _ART_2A_ALG_H


#include "art_common.h"

void vec_normalize(a_lines& lines, const unsigned int norm_index);

unsigned long denoising(a_lines& lines, const float theta);


Clusters art_2A(a_instance& sample,in_param& par);


#endif

