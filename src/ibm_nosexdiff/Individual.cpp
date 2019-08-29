#include "individual.hpp"
#include <cstddef>

Individual::Individual(
        bool const is_female
) : 
    is_female(is_female),
    phen(0.0)
{
}

Individual::Individual() 
   : 
        phen(0),
        is_female(true)
{
}

Individual::Individual(Individual const &Individual) :
    a{Individual.a[0],Individual.a[1]},
    phen(Individual.phen),
    is_female(Individual.is_female)
{
}

void Individual::operator=(Individual const &other)
{
    is_female = other.is_female;
    a[0] = other.a[0];
    a[1] = other.a[1];
    phen = other.phen;
}
