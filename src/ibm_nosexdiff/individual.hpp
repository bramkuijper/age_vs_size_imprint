
// individual class

#ifndef __INDIVIDUAL_INCLUDED__
#define __INDIVIDUAL_INCLUDED__

#include <vector>

class Individual
{
    public:
        static bool is_diploid; // whether individual is diploid or haplodiploid
        static int nloci_f; // female-only loci
        static int nloci_male_only; // male-only loci
        static int nloci_shared; // shared loci

        // two genome copies
        std::vector <double> a[2];

        // resulting phenotype
        double phen;

        // sex of the individual
        bool is_female;

    // constructor, which assumes that
    // af and am are both initialized as
    // empty vectors
    Individual(
        bool const is_female
    );
    
    Individual();

    Individual(Individual const &Individual);

    void operator=(Individual const &other);
};

#endif
