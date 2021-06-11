#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cassert>
#include <algorithm>
#include "care_transgen.hpp"

// set a boundary for a value within (0,1)
double bound001(double val)
{
    return(val < 1e-07 ? 1e-07 : val > 1-(1e-07) ? 1-(1e-07) : val);
}
// set a boundary for a value as [0,1]
double bound01(double val)
{
    return(val < 0 ? 0 : val > 1 ? 1 : val);
}



// constructor
CareTransgen::CareTransgen(int argc, char **argv):
    c{0.0}
    ,a{0.0}
    ,maxs{0.0}
    ,td{0.0}
    ,se{0.0}
    ,r_mother_brood{0.0}
    ,r_mother_off{0.0}
    ,mmin{0.0}
    ,base{}
    ,u{0.0,0.0}
    ,v{0.0,0.0}
    ,b{0.0,0.0}
    ,m{0.0, 0.0, 0.0}
    ,lambda{0.0}
{
    // initialize starting values etc
    init_arguments(argc, argv);
    
    assert(base.length() > 2);

    std::ofstream data_file(base);

    write_parameters(data_file);

    double btplus1[2] = {0.0,0.0};
    double mtplus1[3] = {0.0,0.0,0.0};

    for (long int iter = 0; iter < max_iter_time; ++iter)
    {
        // first update eigenvectors
        calc_eigenvectors();

        // then selection gradient
        btplus1[Lo] = b[Lo] + eul * selgrad_bLo();
        btplus1[Hi] = b[Hi] + eul * selgrad_bHi();
        mtplus1[mHH] = m[mHH] + eul * selgrad_mHH();
        mtplus1[mHL] = m[mHL] + eul * selgrad_mHL();
        mtplus1[mLL] = m[mLL] + eul * selgrad_mLL();

        assert(!std::isnan(btplus1[Lo]));

        std::cout << u[Lo] << " " << u[Hi] << " " << btplus1[Lo] << std::endl;
        std::cout << u[Lo] << " " << u[Hi] << " " << btplus1[Hi] << std::endl;

        break;
    }



}

void CareTransgen::write_parameters(std::ofstream &data_file)
{
    data_file 
        << "c;" << c << std::endl
        << "a;" << a << std::endl
        << "r_mother_brood;" << r_mother_brood << std::endl
        << "r_mother_off;" << r_mother_off << std::endl
        << "maxs;" << maxs << std::endl
        << "td;" << td << std::endl
        << "se;" << se << std::endl
        << "mmin;" << mmin << std::endl
        << "uLo;" << u[Lo] << std::endl
        << "uHi;" << u[Hi] << std::endl
        << "vLo;" << v[Lo] << std::endl
        << "vHi;" << v[Hi] << std::endl
        << "bLo;" << b[Lo] << std::endl
        << "bHi;" << b[Hi] << std::endl
        << "mHH;" << m[mHH] << std::endl
        << "mHL;" << m[mHL] << std::endl 
        << "mLL;" << m[mLL] << std::endl 
        << std::endl;
}

// initialize the parameters of the simulation
void CareTransgen::init_arguments(int argc, char **argv)
{
    c = atof(argv[1]);
    a = atof(argv[2]);
    maxs = atof(argv[3]);
    td = atof(argv[4]);
    se = atof(argv[5]);
    mmin = atof(argv[6]);
    u[0] = atof(argv[7]);
    u[1] = atof(argv[8]);
    v[0] = atof(argv[9]);
    v[1] = atof(argv[10]);
    lambda = atof(argv[11]);
    b[Hi] = atof(argv[12]);
    b[Lo] = atof(argv[13]);
    m[mHH] = atof(argv[14]);
    m[mHL] = atof(argv[15]);
    m[mLL] = atof(argv[16]);
    r_mother_brood = atof(argv[17]);
    r_mother_off = atof(argv[18]);
    base = argv[19];
}

double CareTransgen::f(double m)
{
    return(1.0 - exp(-(m - mmin)));
}

double CareTransgen::dfdm(double m)
{
    return(exp(-(m - mmin)));
}

double CareTransgen::dsdx(double n)
{
    if (s(n) == 0.0)
    {
        return(0.0);
    }
    return(- a * 1.0/maxs * pow(n/maxs,a-1));

}

double CareTransgen::s(double n)
{
    return(std::max(1.0 - pow(n/maxs,a),0.0));
}

// derivative wrt to first argument
double CareTransgen::dBdb(double bx1, double bx2)
{
    double retval = 2 - 2 * (bx1 + bx2);

    if (retval < 0)
    {
        return(0.0);
    }

    return(retval);
}

// pair reproductive effort
double CareTransgen::B(double bx1, double bx2)
{
    return(std::max(2 * (bx1 + bx2) - (bx1 + bx2)*(bx1 + bx2),0.0));
}

double CareTransgen::dmort(double bx)
{

    if (mort(bx) < 1.0)
    {
        return(2 * (1.0 - c) * bx);
    }

    return(0.0);    
}

double CareTransgen::mort(double bx)
{
    return(std::min(c + (1.0 - c) * bx * bx, 1.0));
}

double CareTransgen::q(double m)
{
    return(0.5 + 0.5 * erf((m - td)/(sqrt(2.0) * se)));
}

// TODO check this
double CareTransgen::dqdm(double m)
{
    return(1.0/(sqrt(2.0) * se) * 0.5 * 
            exp((m - td)/(sqrt(2.0) * se)*(m - td)/(sqrt(2.0) * se)));
}

// selection gradient on the trait mHH
double CareTransgen::selgrad_mHH()
{
    double dWhh_brood, dWhl_brood;
    double dWhh_off, dWhl_off;

    dWhh_brood = 1.0/lambda * v[Hi] * u[Hi] * 0.5 * (
            u[Hi] * ( 
                // m in last two terms, f(m) * q(m), are functions of offspring m, not brood m
                (-B(b[Hi],b[Hi])/(m[mHH]*m[mHH]) * s(B(b[Hi],b[Hi])/m[mHH]) + 
                + B(b[Hi],b[Hi])/m[mHH] * dsdx(B(b[Hi],b[Hi])/m[mHH]) * -B(b[Hi],b[Hi])/(m[mHH]*m[mHH]) ) * f(m[mHH]) * q(m[mHH])
                )); 

    dWhh_off = 1.0/lambda * v[Hi] * u[Hi] * 0.5 * 
            u[Hi] * ( 
                // m in last two terms, f(m) * q(m), are functions of offspring m, not brood m
                B(b[Hi],b[Hi])/m[mHH] * s(B(b[Hi],b[Hi])/m[mHH]) * (
                    dfdm(m[mHH]) * q(m[mHH])
                    + f(m[mHH]) * dqdm(m[mHH])
                )); 

    dWhl_brood = 1.0/lambda * v[Lo] * u[Hi] * 0.5 * 
                u[Hi] * (
                    (-B(b[Hi],b[Hi])/(m[mHH]*m[mHH]) * s(B(b[Hi],b[Hi])/m[mHH]) +
                        B(b[Hi],b[Hi])/m[mHH] * dsdx(B(b[Hi],b[Hi])/m[mHH]) * -B(b[Hi],b[Hi])/(m[mHH]*m[mHH]) ) * f(m[mHH]) * (1.0 - q(m[mHH]))
                    );

    dWhl_off = 1.0/lambda * v[Lo] * u[Hi] * 0.5 * 
            u[Hi] * B(b[Hi],b[Hi])/m[mHH] * s(B(b[Hi],b[Hi])/m[mHH]) * (
                    dfdm(m[mHH]) * (1.0 - q(m[mHH]))
                    + f(m[mHH]) * -dqdm(m[mHH]));

    return(r_mother_brood * (dWhh_brood + dWhl_brood) + r_mother_off * (dWhh_off + dWhl_off));
} // end double CareTransgen::selgrad_mHH()

// selection gradient on the trait mHL
double CareTransgen::selgrad_mHL()
{
    double lambda_inv = 1.0/lambda;
    double dWhh_brood = lambda_inv * v[Hi] * u[Hi] * 0.5 * 
            u[Lo] *  
                (-B(b[Hi],b[Lo])/(m[mHL]*m[mHL]) * s(B(b[Hi],b[Lo])/m[mHL]) +
                    B(b[Hi],b[Lo])/m[mHL] * dsdx(B(b[Hi],b[Lo])/m[mHL]) * -B(b[Hi],b[Lo])/(m[mHL]*m[mHL])) * f(m[mHL]) * q(m[mHL]);

    double dWhh_off = lambda_inv * v[Hi] * u[Hi] * 0.5 *
        u[Lo] * B(b[Hi],b[Lo])/m[mHL] * s(B(b[Hi],b[Lo])/m[mHL]) * (dfdm(m[mHL]) * q(m[mHL]) + f(m[mHL]) * dqdm(m[mHL]));

    double dWhl_brood = lambda_inv * v[Lo] * u[Hi] * 0.5 *
                u[Lo] * (-B(b[Hi],b[Lo])/(m[mHL]*m[mHL]) * s(B(b[Hi],b[Lo])/m[mHL]) +
                            B(b[Hi],b[Lo])/m[mHL] * dsdx(B(b[Hi],b[Lo])/m[mHL]) * -B(b[Hi],b[Lo])/(m[mHL]*m[mHL])) * f(m[mHL]) * (1.0 - q(m[mHL]));

    double dWhl_off = lambda_inv * v[Lo] * u[Hi] * 0.5 *
        u[Lo] * B(b[Hi],b[Lo])/m[mHL] * s(B(b[Hi],b[Lo])/m[mHL]) * (dfdm(m[mHL]) * (1.0 - q(m[mHL])) + f(m[mHL]) * -dqdm(m[mHL]));

    double dWlh_brood = lambda_inv * v[Hi] * u[Lo] * 0.5 *
                u[Hi] * (
                        -B(b[Hi],b[Lo])/(m[mHL]*m[mHL]) * s(B(b[Hi],b[Lo])/m[mHL]) 
                        + B(b[Hi],b[Lo])/m[mHL] * dsdx(B(b[Hi],b[Lo])/m[mHL]) * -B(b[Hi],b[Lo])/(m[mHL]*m[mHL])
                        )
                        * f(m[mHL]) * q(m[mHL]);

    double dWlh_off = lambda_inv * v[Hi] * u[Lo] * 0.5 * u[Hi] * 
        B(b[Hi],b[Lo])/m[mHL] * s(B(b[Hi],b[Lo])/m[mHL]) * 
        (dfdm(m[mHL]) * q(m[mHL]) + f(m[mHL]) * dqdm(m[mHL]));

    double dWll_brood = lambda_inv * v[Lo] * u[Lo] * 0.5 * u[Hi] *
        (-B(b[Hi],b[Lo])/(m[mHL]*m[mHL]) * s(B(b[Hi],b[Lo])/m[mHL]) 
         + B(b[Hi],b[Lo])/m[mHL] * dsdx(B(b[Hi],b[Lo])/m[mHL]) * -B(b[Hi],b[Lo])/(m[mHL]*m[mHL])) *
                f(m[mHL]) * (1.0 - q(m[mHL]));

    double dWll_off = lambda_inv * v[Lo] * u[Lo] * 0.5 * u[Hi] *
            B(b[Hi],b[Lo])/m[mHL] * s(B(b[Hi],b[Lo])/m[mHL]) * (
                dfdm(m[mHL]) * (1.0 - q(m[mHL])) + f(m[mHL]) * -dqdm(m[mHL])
                );

    return(r_mother_brood * (dWhh_brood + dWhl_brood + dWlh_brood + dWll_brood) 
            + r_mother_off * (dWhh_off + dWhl_off + dWlh_off + dWll_off));
} // end double CareTransgen::selgrad_mHL()

double CareTransgen::selgrad_mLL()
{
    double lambda_inv = 1.0/lambda;

    double dWlh_brood = lambda_inv * v[Hi] * u[Lo] * 0.5 * 
                u[Lo] * (-B(b[Lo],b[Lo])/(m[mLL]*m[mLL]) * s(B(b[Lo],b[Lo])/m[mLL]) 
                        + B(b[Lo],b[Lo])/m[mLL] * dsdx(B(b[Lo],b[Lo])/m[mLL]) * (-B(b[Lo],b[Lo])/(m[mLL]*m[mLL]))
                        )* f(m[mLL]) * q(m[mLL]);

    double dWlh_off = lambda_inv * v[Hi] * u[Lo] * 0.5 * 
                u[Lo] * B(b[Lo],b[Lo])/m[mLL] * s(B(b[Lo],b[Lo])/m[mLL]) * (
                        dfdm(m[mLL]) * q(m[mLL]) + f(m[mLL]) * dqdm(m[mLL]));

    double dWll_brood = lambda_inv * v[Lo] * u[Lo] * 0.5 *
                u[Lo] * (-B(b[Lo],b[Lo])/(m[mLL]*m[mLL]) * s(B(b[Lo],b[Lo])/m[mLL]) +
                            B(b[Lo],b[Lo])/m[mLL] * dsdx(B(b[Lo],b[Lo])/m[mLL]) * -B(b[Lo],b[Lo])/(m[mLL]*m[mLL])
                        )* f(m[mLL]) * (1.0 - q(m[mLL]));

    double dWll_off = lambda_inv * v[Lo] * u[Lo] * 0.5 *
                u[Lo] * B(b[Lo],b[Lo])/m[mLL] * s(B(b[Lo],b[Lo])/m[mLL]) * (
                dfdm(m[mLL]) * (1.0 - q(m[mLL])) + f(m[mLL]) * -q(m[mLL])
                );

    return(r_mother_brood * (dWlh_brood + dWll_brood) 
            + r_mother_off * (dWlh_off + dWll_off));
} // end double CareTransgen::selgrad_mLL()

// selection gradient on the trait bLo
double CareTransgen::selgrad_bLo()
{
    double dWlh, dWll;

//    std::cout << B(b[Lo],b[Hi]) << " "<< B(b[Lo],b[Hi])/m[mHL] << " " << " " << -dBdb(b[Lo],b[Hi])/m[mHL] << std::endl;
//    std::cout << s(B(b[Lo],b[Hi])/m[mHL]) << " " << dsdx(B(b[Lo],b[Hi])/m[mHL]) << std::endl;
//    std::cout << f(m[mHL]) << " " << q(m[mHL]) << std::endl;
//    std::cout << -dmort(b[Lo]) << std::endl;
//    std::cout << u[Hi] << std::endl;
//    std::cout << dBdb(b[Lo],b[Hi])/m[mHL] * s(B(b[Lo],b[Hi])/m[mHL]) << std::endl;

    dWlh = 1.0/lambda * v[Hi] * u[Lo] *  0.5 * ( 
           // parent mates with Hi parent 
                    u[Hi] * (
            dBdb(b[Lo],b[Hi])/m[mHL] * s(B(b[Lo],b[Hi])/m[mHL])
            +
                B(b[Lo],b[Hi])/m[mHL] * dsdx(B(b[Lo],b[Hi])/m[mHL]) * dBdb(b[Lo],b[Hi]) / m[mHL]                ) * f(m[mHL]) * q(m[mHL])

                    // parent mates with other Lo parent
                    + u[Lo] * (
                        dBdb(b[Lo],b[Lo])/m[mLL] * s(B(b[Lo],b[Lo])/m[mLL])
                        +
                        B(b[Lo],b[Lo])/m[mLL] * dsdx(B(b[Lo],b[Lo])/m[mLL]) * dBdb(b[Lo],b[Hi]) / m[mLL]) * f(m[mLL]) * q(m[mLL])
                    );

    dWll = 1.0/lambda * v[Lo] * u[Lo] * (   -dmort(b[Lo]) + 0.5 * (
                u[Hi] * (
           // parent mates with Hi parent 
            dBdb(b[Lo],b[Hi])/m[mHL] * s(B(b[Lo],b[Hi])/m[mHL])
            +
                B(b[Lo],b[Hi])/m[mHL] * dsdx(B(b[Lo],b[Hi])/m[mHL]) * dBdb(b[Lo],b[Hi]) / m[mHL]                ) * f(m[mHL]) * (1.0 - q(m[mHL]))

                    // parent mates with other Lo parent
                    + u[Lo] * (
                        dBdb(b[Lo],b[Lo])/m[mLL] * s(B(b[Lo],b[Lo])/m[mLL])
                        +
                        B(b[Lo],b[Lo])/m[mLL] * dsdx(B(b[Lo],b[Lo])/m[mLL]) * dBdb(b[Lo],b[Hi]) / m[mLL]) * f(m[mLL]) * (1.0 - q(m[mLL]))
                    )
            );

    return(dWlh + dWll);
}

double CareTransgen::selgrad_bHi()
{
    double dWhh, dWhl;

//    std::cout << B(b[Lo],b[Hi]) << " "<< B(b[Lo],b[Hi])/m[mHL] << " " << " " << -dBdb(b[Lo],b[Hi])/m[mHL] << std::endl;
//    std::cout << s(B(b[Lo],b[Hi])/m[mHL]) << " " << dsdx(B(b[Lo],b[Hi])/m[mHL]) << std::endl;
//    std::cout << f(m[mHL]) << " " << q(m[mHL]) << std::endl;
//    std::cout << -dmort(b[Lo]) << std::endl;
//    std::cout << u[Hi] << std::endl;
//    std::cout << dBdb(b[Lo],b[Hi])/m[mHL] * s(B(b[Lo],b[Hi])/m[mHL]) << std::endl;

    dWhl = 1.0/lambda * v[Lo] * u[Hi] *  0.5 * ( 
           // parent mates with Hi parent 
                    u[Hi] * (
            dBdb(b[Hi],b[Hi])/m[mHH] * s(B(b[Hi],b[Hi])/m[mHH])
            +
                B(b[Hi],b[Hi])/m[mHH] * dsdx(B(b[Hi],b[Hi])/m[mHH]) * dBdb(b[Hi],b[Hi]) / m[mHH]                ) * f(m[mHH]) * (1.0 - q(m[mHH]))

                    // parent mates with other Lo parent
                    + u[Lo] * (
                        dBdb(b[Hi],b[Lo])/m[mHL] * s(B(b[Hi],b[Lo])/m[mHL])
                        +
                        B(b[Hi],b[Lo])/m[mHL] * dsdx(B(b[Hi],b[Lo])/m[mHL]) * dBdb(b[Hi],b[Hi]) / m[mHL]) * f(m[mHL]) * (1.0 - q(m[mHL]))
                    );

    dWhh = 1.0/lambda * v[Hi] * u[Hi] * (-dmort(b[Hi]) + 0.5 * (
                u[Hi] * (
           // parent mates with Hi parent 
            dBdb(b[Hi],b[Hi])/m[mHH] * s(B(b[Hi],b[Hi])/m[mHH])
            +
                B(b[Hi],b[Hi])/m[mHH] * dsdx(B(b[Hi],b[Hi])/m[mHH]) * dBdb(b[Hi],b[Hi]) / m[mHH]                ) * f(m[mHH]) * q(m[mHH])

                    // parent mates with other Lo parent
                    + u[Lo] * (
                        dBdb(b[Hi],b[Lo])/m[mHL] * s(B(b[Hi],b[Lo])/m[mHL])
                        +
                        B(b[Hi],b[Lo])/m[mHL] * dsdx(B(b[Hi],b[Lo])/m[mHL]) * dBdb(b[Hi],b[Hi]) / m[mHL]) * f(m[mHL]) * q(m[mHL])
                    )
            );

    std::cout << "dWhh: " << v[Hi] << " " << lambda << " " << dWhh 
        << " " << (dBdb(b[Hi],b[Hi])/m[mHH] * s(B(b[Hi],b[Hi])/m[mHH])
            +
                B(b[Hi],b[Hi])/m[mHH] * dsdx(B(b[Hi],b[Hi])/m[mHH]) * dBdb(b[Hi],b[Hi]) / m[mHH]) 
        << " " <<
                        (dBdb(b[Hi],b[Lo])/m[mHL] * s(B(b[Hi],b[Lo])/m[mHL])
                        +
                        B(b[Hi],b[Lo])/m[mHL] * dsdx(B(b[Hi],b[Lo])/m[mHL]) * dBdb(b[Hi],b[Hi]) / m[mHL]) 
       << " " <<  (-dmort(b[Hi])) << std::endl; 

    return(dWhl + dWhh);
}

// calculate eigenvectors
void CareTransgen::calc_eigenvectors()
{
    double Whh, Whl, Wlh, Wll, trA, detA;

    double uHi_tplus1;
    double vHi_tplus1;
    double vLo_tplus1;
    double vdenom;

    //matrix structure is
    // Whh, Wlh
    // Whl, Wll

    for (int iter = 0; iter < max_iter_time; ++iter)
    {

        // high quality focal parent to high quality offspring
        Whh = 1.0 - mort(b[Hi]) + 0.5 * (
                u[Hi] * B(b[Hi],b[Hi])/m[mHH] * s(B(b[Hi],b[Hi])/m[mHH]) * f(m[mHH]) * q(m[mHH])
                + u[Lo] * B(b[Hi],b[Lo])/m[mHL] * s(B(b[Hi],b[Lo])/m[mHL]) * f(m[mHL]) * q(m[mHL])
                );

        // high quality focal parent to low quality offspring
        Whl = 0.5 * (
                u[Hi] * B(b[Hi],b[Hi])/m[mHH] * s(B(b[Hi],b[Hi])/m[mHH]) * f(m[mHH]) * (1.0 - q(m[mHH]))
                + u[Lo] * B(b[Hi],b[Lo])/m[mHL] * s(B(b[Hi],b[Lo])/m[mHL]) * f(m[mHL]) * (1.0 - q(m[mHL])));

        // low quality parent to high quality offspring
        Wlh = 0.5 * (
                u[Hi] * B(b[Hi],b[Lo])/m[mHL] * s(B(b[Hi],b[Lo])/m[mHL]) * f(m[mHL]) * q(m[mHL])
                + u[Lo] * B(b[Lo],b[Lo])/m[mLL] * s(B(b[Lo],b[Lo])/m[mLL]) * f(m[mLL]) * q(m[mLL])
                );

        // low quality parent to high quality offspring
        Wll = 1.0 - mort(b[Lo]) + 0.5 * (
                u[Hi] * B(b[Hi],b[Lo])/m[mHL] * s(B(b[Hi],b[Lo])/m[mHL]) * f(m[mHL]) * (1.0 - q(m[mHL]))
                + u[Lo] * B(b[Lo],b[Lo])/m[mLL] * s(B(b[Lo],b[Lo])/m[mLL]) * f(m[mLL]) * (1.0 - q(m[mLL]))
                );

//        std::cout << B(b[Hi],b[Hi]) << " " << b[Hi] << " " << b[Lo] << " " << m[mHH] << " " << B(b[Hi],b[Hi])/m[mHH] <<  " " << Whh << " " << Wlh << " " << Whl << " " << Wll << std::endl;

        trA = Whh + Wll;
        detA = Whh * Wll - Whl * Wlh;

        lambda = 0.5 * (trA + sqrt(trA*trA - 4 * detA));

        uHi_tplus1 = Wlh / (lambda  - Whh + Wlh);


        vdenom = (Whl * u[Hi] - Whh * u[Lo] + u[Lo] * lambda);

        // let's set vLo to 1
        vHi_tplus1 = Whl/vdenom;
        vLo_tplus1 = (-Whh + lambda)/vdenom;
        
        std::cout << "ev:" << iter << " " << u[Hi] << " " << uHi_tplus1 << " " << lambda << std::endl;
        std::cout << "ev:" << iter << " " << v[Hi] << " " << vHi_tplus1 << " " << lambda << std::endl;

        if (fabs(uHi_tplus1 - u[Hi]) < vanish_bound
                && fabs(vHi_tplus1 - v[Hi]) < vanish_bound
                && fabs(vLo_tplus1 - v[Lo]) < vanish_bound
                )
        {
            break;
        }

        u[Hi] = uHi_tplus1;
        u[Lo] = 1.0 - uHi_tplus1;
        v[Hi] = vHi_tplus1;
        v[Lo] = vLo_tplus1;
    } // end for
} // end CareTransgen::calc_eigenvectors()
