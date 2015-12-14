#include "individual.hpp"
#include <cstddef>

Individual::Individual(
        double af[2],
        double am[2],
        bool const is_female
) : 
    af{af[0],af[1]}, 
    am{am[0],am[1]}, 
    is_female(is_female)
{
    phen = is_female ? af[0] + af[1] : (is_diploid ? am[0] + am[1] : am[0]);
}

Individual::Individual() 
   : 
        af{0,0},
        am{0,0},
        phen(0),
        is_female(true)
{
}

Individual::Individual(Individual const &Individual) :
    af{Individual.af[0],Individual.af[1]},
    am{Individual.am[0],Individual.am[1]},
    phen(Individual.phen),
    is_female(Individual.is_female)
{
}

void Individual::operator=(Individual const &other)
{
    is_female = other.is_female;
    af[0] = other.af[0];
    af[1] = other.af[1];
    am[0] = other.am[0];
    am[1] = other.am[1];
    phen = other.phen;
}

bool Individual::is_diploid;


