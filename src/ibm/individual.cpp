#include "individual.hpp"

// default constructor
Individual::Individual():
    u{{0.0,0.0},{0.0,0.0}},
    m{{0.0,0.0},{0.0,0.0},{0.0,0.0}},
    quality_high(false)
{
}

Individual::Individual(Individual const &other):
    u{{other.u[0][0],other.u[0][1]},{other.u[0][0],other.u[0][1]}},
    m{{other.m[0][0],other.m[0][1]},{other.m[1][0],other.m[1][1]},{other.m[2][0],other.m[2][1]}},
    quality_high(other.quality_high)
{
}

void Individual::operator=(Individual const &other)
{
    for (int allele_idx = 0; allele_idx < 2; ++allele_idx)
    {
        for (int hilo_idx = 0; hilo_idx < 2; ++hilo_idx)
        {
            u[hilo_idx][allele_idx] = other.u[hilo_idx][allele_idx];
        }
    
        for (int hl_idx = 0; hl_idx < 3; ++hl_idx)
        {
            m[hl_idx][allele_idx] = other.m[hl_idx][allele_idx];
        }
    }
    
    quality_high = other.quality_high;
}
