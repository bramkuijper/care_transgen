
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <cmath>
#include <random>
#include "individual.hpp"


// C++ random number generation unsigned int seed = get_nanoseconds();
std::random_device rd;
unsigned seed = rd();
std::mt19937 rng_r{seed};
std::uniform_real_distribution<> uniform(0.0,1.0);

// file base name
std::string base_name;


// number of individuals in population
const int NPatches = 4000;
const int NBreederF = 1;

int number_generations = 20;

const int max_clutch_size = 50;

double init_u = 0.5;
double init_m = 0.5;

double ph_init = 0.5;

// condition dependent costs of care
double k[] = {0.0,0.0};

// patch
struct Patch {
   
    // female and male breeders
    Individual breedersF[NBreederF];
    Individual breedersM[NBreederF];

    int nf;
    int nm;
    
    Individual philopatric_juvs[NBreederF * max_clutch_size];
};

// start a meta population of Npatches
Patch MetaPop[NPatches];

// initialize population
void init_population()
{
    for (int patch_idx = 0; patch_idx < NPatches; ++patch_idx)
    {
        // first initialize females
        for (int breederF_idx = 0; breederF_idx < NBreederF; ++breederF_idx)
        {
            for (int allele_idx = 0; allele_idx < 2; ++allele_idx)
            {
                // parental care loci
                for (int hilo_idx = 0; hilo_idx < 2; ++hilo_idx)
                {
                    MetaPop[patch_idx].breedersF[breederF_idx].u[hilo_idx][allele_idx] = init_u;
                }

                // offspring size loci
                for (int hx_idx = 0; hx_idx < 3; ++hx_idx)
                {
                    MetaPop[patch_idx].breedersF[breederF_idx].m[hx_idx][allele_idx] = init_m;
                }
            }

            // quality
            MetaPop[patch_idx].breedersF[breederF_idx].quality_high = uniform(rng_r) < ph_init;    

        } // end breederF_idx
                
        for (int breederM_idx = 0; breederM_idx < NBreederF; ++breederM_idx)
        {
            for (int allele_idx = 0; allele_idx < 2; ++allele_idx)
            {
                for (int hilo_idx = 0; hilo_idx < 2; ++hilo_idx)
                {
                    MetaPop[patch_idx].breedersM[breederM_idx].u[hilo_idx][allele_idx] = init_u;
                }
                
                for (int hx_idx = 0; hx_idx < 3; ++hx_idx)
                {
                    MetaPop[patch_idx].breedersM[breederM_idx].m[hx_idx][allele_idx] = init_m;
                }
            }
        }

        MetaPop[patch_idx].nf = NBreederF;
        MetaPop[patch_idx].nm = NBreederF;
    }// end for (int patch_idx = 0; patch_idx < NPatches; ++patch_idx)

}

// get parameters from the 
// command line
void get_parameters(int argc, char **argv)
{
    base_name = argv[1];
    init_u = atof(argv[2]);
    init_m = atof(argv[3]);
    k[0] = atof(argv[4]);
    k[1] = atof(argv[5]);
}

// write parameters
void write_parameters(std::ofstream &data_file)
{
    data_file << 
        std::endl
        << std::endl
        << "init_u;" << init_u << std::endl
        << "init_m;" << init_m << std::endl
        << "ph_init;" << ph_init << std::endl
        << "klo;" << k[0] << std::endl
        << "khi;" << k[1] << std::endl
        << "seed;" << seed << std::endl;
}

double mu_mort(bool const quality_hi, double const u)
{
    return(k[quality_hi] + (1.0 - k[quality_hi]) * u * u);
}

void survive()
{
    bool quality_hi;

    double u;

    for (int patch_idx = 0; patch_idx < NPatches; ++patch_idx)
    {
        for (int breederF_idx = 0; breederF_idx < MetaPop[patch_idx].nf; ++breederF_idx)
        {
            quality_hi = MetaPop[patch_idx].breedersF[breederF_idx].quality_high;
            u = MetaPop[patch_idx].breedersF[breederF_idx].u[quality_hi][0]
                    + MetaPop[patch_idx].breedersF[breederF_idx].u[quality_hi][1];

            // individual dies
            if (uniform(rng_r) < mu_mort(quality_hi,u))
            {
                // overwrite individual with last one from the stack
                MetaPop[patch_idx].breedersF[breederF_idx] = 
                    MetaPop[patch_idx].breedersF[MetaPop[patch_idx].nf - 1];

                --breederF_idx;
                --MetaPop[patch_idx].nf;
            }
        }

        for (int breederM_idx = 0; breederM_idx < MetaPop[patch_idx].nm; ++breederM_idx)
        {
            quality_hi = MetaPop[patch_idx].breedersM[breederM_idx].quality_high;
            u = MetaPop[patch_idx].breedersM[breederM_idx].u[quality_hi][0]
                + MetaPop[patch_idx].breedersM[breederM_idx].u[quality_hi][1];

            // individual dies
            if (uniform(rng_r)  < mu_mort(quality_hi, u))
            {
                MetaPop[patch_idx].breedersM[breederM_idx] = 
                    MetaPop[patch_idx].breedersM[MetaPop[patch_idx].nm - 1];
                
                --breederM_idx;
                --MetaPop[patch_idx].nm;
            }
        }
    }
} // end survive()


// replace dead adults by reproducing
void reproduce()
{
    int father;


    for (int patch_idx = 0; patch_idx < NPatches; ++patch_idx)
    {
        if (MetaPop[patch_idx].nf < 1)
        {
            // alas no reproduction, as patch is extinct.
            continue;
        }

        for (int breederF_idx = 0; breederF_idx < MetaPop[patch_idx].nf; ++breederF_idx)
        {
            // randomly choose a male and make offspring
            if 
        }
    }
}


int main(int argc, char **argv)
{
    get_parameters(argc, argv);
    init_population();

    std::ofstream data_file(base_name);

    int generation;

    for (generation = 0; generation < number_generations; ++generation)
    {
        survive();

        reproduce();
    }
}
