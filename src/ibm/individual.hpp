
// individual class

#ifndef __INDIVIDUAL_INCLUDED__
#define __INDIVIDUAL_INCLUDED__

#include <cstddef>

class Individual
{
    public:
        static bool is_diploid;
        double af[2];
        double am[2];
        double phen;
        bool is_female;

    Individual(
        double af[2],
        double am[2],
        bool const is_female
    );
    
    Individual();

    Individual(Individual const &Individual);

    void operator=(Individual const &other);
};

#endif
