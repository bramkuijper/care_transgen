#ifndef CARE_TRANSGEN_HPP_
#define CARE_TRANSGEN_HPP_

#include <vector>

enum Resource
{
    Hi = 0,
    Lo = 1
};

// all combinations of high and low quality parents
// when producing offspring of a certain size
enum M
{
    mHH = 0, 
    mHL = 1,
    mLL = 2
};

class CareTransgen
{
    private:
        double c; // increase in mortality with care
        double a; // power of how kin competition affects survival
        double maxs; // 
        double td;
        double se;
        double r_mother_brood;
        double r_mother_off;
        double mmin;
        std::string base;

        // vector to store stable class frequencies
        double u[2];

        // vector to store reproductive values
        double v[2];

        // vector with the reproductive effort traits (hi, lo)
        double b[2];
        
        // vector with the offspring size (hi, lo)
        double m[3];

        // eigenvalue
        double lambda;

        // max time to iterate
        static const long int max_iter_time = 1e08;
        static constexpr double eul = 0.01;

        static constexpr double vanish_bound = 1e-07;

    // write parameters to the stream str
    void write_parameters(std::ofstream &data_file);

    void init_arguments(int argc, char **argv);

    // calculate left right eigenvectors
    void calc_eigenvectors();

    // probability of becoming a well-fed parent
    double q(double m);
    double dqdm(double m);

    // survival due to size
    double f(double m);
    double dfdm(double m);

    // survival due to kin competition
    double s(double n);
    double dsdx(double n);

    // biparental reproductive effort
    double B(double bx1, double bx2);
    double dBdb(double bx1, double bx2);

    // probability individual dies
    double mort(double bx);
    double dmort(double bx);

    double selgrad_bLo();
    double selgrad_bHi();
    double selgrad_mHH();
    double selgrad_mHL();
    double selgrad_mLL();

    public:
        CareTransgen(int argc, char **argv);

}; //end class CareTransgen

#endif
