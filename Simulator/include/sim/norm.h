#ifndef NORM_H
#define NORM_H

extern int ds_seed, ds_b;
extern double ds_m;

//Seeds the RNG using a long as the seed
void ds_srand32(long s);
//Returns a random positive integer between 0 and 2^31
int ds_rand32();
//Returns a random double (uniform distribution in the range [0.0, 1.0)
double ds_uniform();
//Samples from a normal distrubution using a polar box-muller transformation.
double ds_norm();

#endif
