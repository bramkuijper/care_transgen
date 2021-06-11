#ifndef __INDIVIDUAL_INCLUDED__
#define __INDIVIDUAL_INCLUDED__

class Individual
{
    public:

        // parental investment
        // first subscript: high vs low
        // second subscript: allele
        double u[2][2];

        // offspring size 
        //
        double m[3][2];

        bool quality_high;

        Individual();

        Individual(Individual const &other);

        void operator=(Individual const &other);
};

#endif
