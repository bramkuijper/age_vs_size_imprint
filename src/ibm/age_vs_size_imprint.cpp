/* age versus size at maturity: imprinting */

#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <cmath>
#include <cassert>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "bramauxiliary.h"
#include "individual.hpp"


//#define NDEBUG

using namespace std;

// random number generator 
// see http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html#Random-Number-Generation 
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *r; // gnu scientific rng 

const size_t Npatches = 1000; 
const size_t Nfp = 5; // females per patch
const size_t Nmp = 5; // males per patch
const size_t numgen = 20000;
const size_t Clutch = 200;

// mutation rates
double mu = 0;
double sdmu = 0;

double df = 0;
double dm = 0;
double q = 0;
double k = 0;
size_t type = 0;
bool diploid = false;


// tally of dispersers
size_t NdispF = 0;
size_t NdispM = 0;
size_t NsurvM = 0;
size_t NsurvF = 0;

// runtime for stats
time_t total_time; 

size_t generation = 0;

int seed = -1;

// skip the number of generations in output
// to prevent output files from becoming too large
size_t skip = 1;

// 
struct Patch
{
    Individual localsF[Nfp]; // all the local female breeders
    Individual localsM[Nmp]; // all the local male breeders

    // philopatric female offspring
    vector <Individual> philsF;     
    vector <Individual> philsM;     

    // (note that dispersing offspring
    // ends up in global pool, see below)

    // total number of kids 
    size_t NkidsF; 
    size_t NkidsM; 
    size_t NhelpersF; 
    size_t NhelpersM; 

    // variable that allows for correct
    // rounding of the number of immigrants
    // per patch (see below)
    size_t immigrant_bonusF;
    size_t immigrant_bonusM;

};

// generate the complete population
Patch MetaPop[Npatches];
vector<Individual> DispersersF;
vector<Individual> DispersersM;

// give the outputfile a unique name
// by using the create_filename function (see 
// "bramauxiliary.h")
string filename("sim_sa_sexbias");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

// initialize the command line arguments to vary 
// the parameters
void init_arguments(int argc, char *argv[])
{
    mu = atof(argv[1]); // mutation prob
    sdmu = atof(argv[2]); // mutational distribution
    df = atof(argv[3]);
    dm = atof(argv[4]);
    k = atof(argv[5]);
    q = atof(argv[6]);
    diploid = atoi(argv[7]);
    type = atoi(argv[8]);

    Individual::is_diploid = diploid;
}

void init_pop()
{
    // start the time
    total_time = time(NULL);

    // obtain a seed from current nanosecond count
	seed = get_nanoseconds();
    // set the seed to the random number generator
    // stupidly enough, this can only be done by setting
    // a shell environment parameter
    stringstream s;
    s << "GSL_RNG_SEED=" << setprecision(10) << seed;
    putenv(const_cast<char *>(s.str().c_str()));

    // set up the random number generators
    // (from the gnu gsl library)
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // go through all patches
    for (size_t i = 0; i < Npatches; ++i)
    {
        for (size_t j = 0; j < Nfp; ++j)
        {
            // initialize phenotypes
            MetaPop[i].localsF[j].af[0] = 3.0;
            MetaPop[i].localsF[j].af[1] = 3.0;
            MetaPop[i].localsF[j].am[0] = 3.0;
            MetaPop[i].localsF[j].am[1] = 3.0;
            MetaPop[i].localsF[j].phen = 3.0;
            MetaPop[i].localsF[j].is_female = true;
        }
        
        for (size_t j = 0; j < Nmp; ++j)
        {
            // initialize phenotypes
            MetaPop[i].localsM[j].af[0] = 3.0;
            MetaPop[i].localsM[j].af[1] = 3.0;
            MetaPop[i].localsM[j].am[0] = 3.0;
            MetaPop[i].localsM[j].am[1] = 3.0;
            MetaPop[i].localsM[j].phen = 3.0;
            MetaPop[i].localsM[j].is_female = false;
        }
    }
}

// mutate alleles
double Mut(double allele)
{
    if (gsl_rng_uniform(r) < mu)
    {
        allele += gsl_ran_gaussian(r, sdmu);

        if (allele < 0)
        {
            allele = 0;
        }
    }

    return(allele);
}


// allocate a kid and give it genes
void create_kid(Individual &mother, Individual &father, Individual &Kid, bool is_female)
{
    Kid.phen = -1;
    Kid.is_female = is_female;

    // first inherit af
    
    // maternal allele
    size_t maternal_allele_f = gsl_rng_uniform_int(r,2);
    Kid.af[0] = Mut(mother.af[maternal_allele_f]); 

    // paternal allele
    // if in haplodiploids, it is the first allele of the father
    size_t paternal_allele_f = diploid ? gsl_rng_uniform_int(r,2) : 0;

    // if in a haplodiploid male 
    // assign no value
    Kid.af[1] = !is_female && !diploid ? -1000 : Mut(father.af[paternal_allele_f]);

    // if female, determine phenotype
    if (is_female)
    {
        switch (type)
        {
            case 1: // madumnal allele
                Kid.phen = mother.af[maternal_allele_f];
                break;
            case 2: // padumnal allele
                Kid.phen = father.af[paternal_allele_f];
                break;
            default: // normal allele
                Kid.phen = .5 * (Kid.af[0] + Kid.af[1]);
                break;
        }
    }

    // inherit am
    size_t maternal_allele_m = gsl_rng_uniform_int(r,2);
    Kid.am[0] = Mut(mother.am[maternal_allele_m]); 

    // paternal allele am
    size_t paternal_allele_m = diploid ? gsl_rng_uniform_int(r, 2) : 0;
    Kid.am[1] = !is_female && !diploid ? -1000 : Mut(father.am[paternal_allele_m]);

    if (!is_female)
    {
        switch(type)
        {
            case 1: // madumnal allele
                Kid.phen = mother.am[maternal_allele_m];
                break;
            case 2: // padumnal allele
                Kid.phen = diploid ? father.am[paternal_allele_m] : Kid.am[0];
                break;
            default: // normal allele
                Kid.phen = diploid ? .5 * (Kid.am[0] + Kid.am[1]) : Kid.am[0];
                break;
        }
    }

    assert(Kid.phen >= 0);
}

// mate and create kids across all patches...
void make_juveniles()
{
    // clear vectors
    DispersersF.clear();
    DispersersM.clear();


    // reset surviving individual counters
    NsurvM = 0;
    NsurvF = 0;

    // variable to calculate clutch size;
    double exact_clutch;
    size_t clutch_size;

    // store current father
    size_t father;

    // store current sex of offspring
    bool is_female;

    // cumulative distribution of father's effort
    double cumul_male_effort[Nmp];

    double sum_male_effort = 0;

    double cumul_deviate;
    
    // generate offspring on each patch
    for (size_t i = 0; i < Npatches; ++i)
    {
        // clear anything from the previous generation
        MetaPop[i].philsF.clear();
        MetaPop[i].philsM.clear();

        // variable that takes into account
        // any beyond-the-minimum number of immigrants
        // (see later)
        MetaPop[i].immigrant_bonusF = 0;
        MetaPop[i].immigrant_bonusM = 0;

        // make cumulative distribution of male efforts
        sum_male_effort = 0;

        for (size_t j = 0; j < Nmp; ++j)
        {
            cumul_male_effort[j] = sum_male_effort + pow(k * MetaPop[i].localsM[j].phen,3);
            sum_male_effort = cumul_male_effort[j];
        }

        clutch_size = 0;

        // reproduce
        for (size_t j = 0; j < Nfp; ++j)
        {
            // calculate female clutch size
            exact_clutch = pow(k * MetaPop[i].localsF[j].phen,3) * Clutch; 
            clutch_size = floor(exact_clutch);
           
            // stochastic rounding
            if (gsl_rng_uniform(r) < exact_clutch - clutch_size)
            {
                ++clutch_size;
            }


            // create kids and reduce parental resources
            for (size_t egg_i = 0; egg_i < clutch_size; ++egg_i)
            {
                // select male based on his reproductive effort
                cumul_deviate = gsl_rng_uniform(r) * sum_male_effort;

                father = 1000;

                for (size_t male_i = 0; male_i < Nmp; ++male_i)
                {
                    if (cumul_deviate <= cumul_male_effort[male_i])
                    {
                        father = male_i;
                        break;
                    }
                }

                assert(father >= 0 && father < Nmp);

                Individual Kid; 

                is_female = gsl_rng_uniform(r) < 0.5;

                create_kid(MetaPop[i].localsF[j],MetaPop[i].localsM[father],Kid, is_female);


                // female or male
                if (is_female)
                {
                    ++NsurvF;
                    // dispersal or not
                    if (gsl_rng_uniform(r) < df)
                    {                            
                        // DispersersF[NdispF++] = Kid;
                        DispersersF.push_back(Kid);
                    }
                    else
                    {
                        // MetaPop[i].philsF[MetaPop[i].NkidsF++] = Kid;
                        MetaPop[i].philsF.push_back(Kid);
                    }
                }
                else
                {
                    ++NsurvM;
                    // dispersal or not
                    if (gsl_rng_uniform(r) < dm)
                    {                            
                        // DispersersM[NdispM++] = Kid;
                        DispersersM.push_back(Kid);
                    }
                    else
                    {
                        // MetaPop[i].philsM[MetaPop[i].NkidsM++] = Kid;
                        MetaPop[i].philsM.push_back(Kid);
                    }
                }
            } // end clutch

//            if (MetaPop[i].philsM.size() < Nmp)
//            {
//                cout << MetaPop[i].philsM.size() << " " << MetaPop[i].philsF.size() << " " << total_offspring_produced << " " << MetaPop[i].localsF[Nfp-1].phen << " " << clutch_size << endl;
//            }
        } // end Nfp
    }//Npatches

}

// replacement of adults with juveniles
void replace_adults()
{
    // okay, we have Ndisp/Npatches dispersers per patch
    size_t dispersers_per_patchF = floor((double) DispersersF.size() / Npatches);
    size_t dispersers_per_patchM = floor((double) DispersersM.size() / Npatches);

    // however, we need to correctly round this rational number
    // to the actual number of immigrants for 
    // a given patch. To this end, 
    // we just randomly distribute the dispersers that remain after the
    // previous rounding over the different patches
    for (size_t i = 0; i < DispersersF.size() - Npatches * dispersers_per_patchF; ++i)
    {
        // randomly picked patch receives additional immigrant
        MetaPop[gsl_rng_uniform_int(r,Npatches)].immigrant_bonusF++;
    }
    for (size_t i = 0; i < DispersersM.size() - Npatches * dispersers_per_patchM; ++i)
    {
        // randomly picked patch receives additional immigrant
        MetaPop[gsl_rng_uniform_int(r,Npatches)].immigrant_bonusM++;
    }

    size_t random_ind;

    double rand_deviate;

    bool new_breeder_chosen = false;
    
    // now replace local breeders on each patch
    for (size_t i = 0; i < Npatches; ++i)
    {
        if (dispersers_per_patchF + MetaPop[i].immigrant_bonusF > 0)
        {
            assert(DispersersF.size() > 0);
        }

        if (dispersers_per_patchM + MetaPop[i].immigrant_bonusM > 0)
        {
            assert(DispersersM.size() > 0);
        }

        size_t arriving_immigrantsF = dispersers_per_patchF + MetaPop[i].immigrant_bonusF;
        size_t arriving_immigrantsM = dispersers_per_patchM + MetaPop[i].immigrant_bonusM;

        for (size_t j = 0; j < Nfp; ++j)
        {
            assert(arriving_immigrantsF + MetaPop[i].philsF.size() > 0);
            // make cumulative distribution of female immigrants and philopatric individuals
            double cumul_dist[arriving_immigrantsF + MetaPop[i].philsF.size()];
            size_t immigrant_id[arriving_immigrantsF];
            double sum_juv = 0;
            size_t juv_i = 0;
            size_t phil_i = 0;

            // first work through the immigrant offspring
            for (; juv_i < arriving_immigrantsF; ++juv_i)
            {
                // draw a random immigrant
                random_ind = gsl_rng_uniform_int(r,DispersersF.size());
                immigrant_id[juv_i] = random_ind;
                assert(DispersersF[random_ind].phen >= 0 && DispersersF[random_ind].phen <= 100);
//                assert(DispersersF[random_ind].phen == .5 * (DispersersF[random_ind].af[0] + DispersersF[random_ind].af[1]));
                cumul_dist[juv_i] = sum_juv + exp(-q * DispersersF[random_ind].phen);
                sum_juv = cumul_dist[juv_i];
            }

            // then work through the philopatric offspring
            for (; juv_i < arriving_immigrantsF + MetaPop[i].philsF.size(); ++juv_i)
            {
                phil_i = juv_i - arriving_immigrantsF;

                assert(phil_i >= 0 && phil_i < MetaPop[i].philsF.size());
                assert(MetaPop[i].philsF[phil_i].phen >= 0 && MetaPop[i].philsF[phil_i].phen <= 100);
//                assert(MetaPop[i].philsF[phil_i].phen == .5 * (MetaPop[i].philsF[phil_i].af[0] + MetaPop[i].philsF[phil_i].af[1]));
                cumul_dist[juv_i] = sum_juv + exp(-q * MetaPop[i].philsF[phil_i].phen);
                sum_juv = cumul_dist[juv_i];
            }

            // finally, choose the individual
            rand_deviate = sum_juv * gsl_rng_uniform(r);

            new_breeder_chosen = false;

            for (juv_i = 0; juv_i < arriving_immigrantsF + MetaPop[i].philsF.size(); ++juv_i)
            {
                if (rand_deviate <= cumul_dist[juv_i])
                {
                    if (juv_i < arriving_immigrantsF)
                    {
                        assert(arriving_immigrantsF > 0);
                        assert(immigrant_id[juv_i] >= 0 && immigrant_id[juv_i] < DispersersF.size());
                        MetaPop[i].localsF[j] = DispersersF[immigrant_id[juv_i]];
                        DispersersF[immigrant_id[juv_i]] = DispersersF.back();
                        DispersersF.pop_back();
                        --arriving_immigrantsF;
                    }
                    else
                    {
                        assert(MetaPop[i].philsF.size() > 0);
                        phil_i = juv_i - arriving_immigrantsF;
                        assert(phil_i >= 0 && phil_i < MetaPop[i].philsF.size());
                        MetaPop[i].localsF[j] = MetaPop[i].philsF[phil_i];
                        MetaPop[i].philsF[phil_i] = MetaPop[i].philsF.back();
                        MetaPop[i].philsF.pop_back();
                    }
                    new_breeder_chosen = true;
                    break;
                }
            }
            assert(new_breeder_chosen);
        }
        
        for (size_t j = 0; j < Nmp; ++j)
        {
            assert(arriving_immigrantsM + MetaPop[i].philsM.size() > 0);
            // make cumulative distribution of female immigrants and philopatric individuals
            double cumul_dist[arriving_immigrantsM + MetaPop[i].philsM.size()];
            size_t immigrant_id[arriving_immigrantsM];
            double sum_juv = 0;
            size_t juv_i = 0;
            size_t phil_i = 0;

            // first work through the immigrant offspring
            for (; juv_i < arriving_immigrantsM; ++juv_i)
            {
                // draw a random immigrant
                random_ind = gsl_rng_uniform_int(r,DispersersM.size());
                immigrant_id[juv_i] = random_ind;
                assert(DispersersM[random_ind].phen >= 0 && DispersersM[random_ind].phen <= 100);
//                assert(DispersersM[random_ind].phen == DispersersM[random_ind].am[0]);
                cumul_dist[juv_i] = sum_juv + exp(-q * DispersersM[random_ind].phen);
                sum_juv = cumul_dist[juv_i];
            }

            // then work through the philopatric offspring
            for (; juv_i < arriving_immigrantsM + MetaPop[i].philsM.size(); ++juv_i)
            {
                phil_i = juv_i - arriving_immigrantsM;

                assert(phil_i >= 0 && phil_i < MetaPop[i].philsM.size());

                assert(MetaPop[i].philsM[phil_i].phen >= 0 && MetaPop[i].philsM[phil_i].phen <= 100);
//                assert(MetaPop[i].philsM[phil_i].phen == MetaPop[i].philsM[phil_i].am[0]);
                cumul_dist[juv_i] = sum_juv + exp(-q * MetaPop[i].philsM[phil_i].phen);
                sum_juv = cumul_dist[juv_i];
            }

            // finally, choose the individual
            rand_deviate = sum_juv * gsl_rng_uniform(r);

            new_breeder_chosen = false;


            for (juv_i = 0; juv_i < arriving_immigrantsM + MetaPop[i].philsM.size(); ++juv_i)
            {
                if (rand_deviate <= cumul_dist[juv_i])
                {
                    if (juv_i < arriving_immigrantsM)
                    {
                        assert(arriving_immigrantsM > 0);
                        assert(immigrant_id[juv_i] >= 0 && immigrant_id[juv_i] < DispersersM.size());
                        MetaPop[i].localsM[j] = DispersersM[immigrant_id[juv_i]];
                        DispersersM[immigrant_id[juv_i]] = DispersersM.back();
                        DispersersM.pop_back();
                        --arriving_immigrantsM;
                    }
                    else
                    {
                        assert(MetaPop[i].philsM.size() > 0);
                        phil_i = juv_i - arriving_immigrantsM;
                        assert(phil_i >= 0 && phil_i < MetaPop[i].philsM.size());
                        MetaPop[i].localsM[j] = MetaPop[i].philsM[phil_i];
                        MetaPop[i].philsM[phil_i] = MetaPop[i].philsM.back();
                        MetaPop[i].philsM.pop_back();
                    }
                    new_breeder_chosen = true;
                    break;
                }
            }

            assert(new_breeder_chosen);
        }
        
        for (size_t j = 0; j < arriving_immigrantsF; ++j)
        {
            random_ind = gsl_rng_uniform_int(r, DispersersF.size());
            DispersersF[random_ind] = DispersersF.back();
            DispersersF.pop_back();
        }
        // remove all remaining male immigrants from the global dispersal pool
        for (size_t j = 0; j < arriving_immigrantsM; ++j)
        {
            random_ind = gsl_rng_uniform_int(r, DispersersM.size());
            DispersersM[random_ind] = DispersersM.back();
            DispersersM.pop_back();
        }
    }
    assert(DispersersF.size()==0);
    assert(DispersersM.size()==0);
}

void write_data_headers()
{
    DataFile << "generation;meanaf;varaf;meanam;varam;nsurvf;nsurvm;" << endl;
}

void write_data()
{
    double meanaf = 0;
    double ssaf = 0;
    double meanam = 0;
    double ssam = 0;

    for (size_t i = 0; i < Npatches; ++i)
    {
        for (size_t j = 0; j < Nfp; ++j)
        {
            meanaf += .5 * (MetaPop[i].localsF[j].af[0]+ MetaPop[i].localsF[j].af[1]);
            ssaf += .5 * (MetaPop[i].localsF[j].af[0]+ MetaPop[i].localsF[j].af[1]) * 
                .5 * (MetaPop[i].localsF[j].af[0]+ MetaPop[i].localsF[j].af[1]);

            meanam += .5 * (MetaPop[i].localsF[j].am[0]+ MetaPop[i].localsF[j].am[1]);
            ssam += .5 * (MetaPop[i].localsF[j].am[0]+ MetaPop[i].localsF[j].am[1]) * 
                .5 * (MetaPop[i].localsF[j].am[0]+ MetaPop[i].localsF[j].am[1]);
        }
        
        for (size_t j = 0; j < Nmp; ++j)
        {
            if (diploid)
            {
                meanaf += .5 * (MetaPop[i].localsM[j].af[0]+ MetaPop[i].localsM[j].af[1]);
                ssaf += .5 * (MetaPop[i].localsM[j].af[0]+ MetaPop[i].localsM[j].af[1]) * 
                    .5 * (MetaPop[i].localsM[j].af[0]+ MetaPop[i].localsM[j].af[1]);

                meanam += .5 *(MetaPop[i].localsM[j].am[0]+ MetaPop[i].localsM[j].am[1]);
                ssam += .5 * (MetaPop[i].localsM[j].am[0]+ MetaPop[i].localsM[j].am[1]) * 
                    .5 * (MetaPop[i].localsM[j].am[0]+ MetaPop[i].localsM[j].am[1]);
            } else
            {
                meanaf += MetaPop[i].localsM[j].af[0];
                ssaf += MetaPop[i].localsM[j].af[0] * MetaPop[i].localsM[j].af[0];

                meanam += MetaPop[i].localsM[j].am[0];
                ssam += MetaPop[i].localsM[j].am[0] * MetaPop[i].localsM[j].am[0];
            }
        }
    }

    double varaf = 0;
    double varam = 0;
    
    meanaf /= Npatches * (Nfp + Nmp);
    meanam /= Npatches * (Nfp + Nmp);
    varaf = ssaf / (Npatches * (Nfp + Nmp)) - meanaf * meanaf;
    varam = ssam / (Npatches * (Nfp + Nmp)) - meanam * meanam;

        
    DataFile << generation 
                        << ";" << meanaf 
                        << ";" << varaf 
                        << ";" << meanam 
                        << ";" << varam  
                        << ";" << NsurvF
                        << ";" << NsurvM
                        << ";" << endl;
}

void write_parameters()
{
    DataFile << endl << endl << "system;" << (diploid ? "diploid" : "haplodiploid") << endl
                << "imprint;" << (type == 0 ? "offspring" : (type == 1 ? "madumnal" : "padumnal")) << endl
                << "patch;" << Npatches << endl
                << "nfp;" << Nfp << endl
                << "nmp;" << Nmp << endl
                << "seed;" << seed << endl
                << "df;" << df << endl
                << "dm;" << dm << endl
                << "numgen;" << numgen << endl
                << "mu;" << mu << endl
                << "sdmu;" << sdmu << endl
                << "q;" << q << endl
                << "k;" << k << endl
                << "clutch;" << Clutch << endl
                << "runtime;" << total_time << endl;
}


int main(int argc, char * argv[])
{
    init_arguments(argc,argv);
    init_pop();

    write_data_headers();


    for (generation = 0; generation < numgen; ++generation)
    {
        make_juveniles();

        if (generation % skip == 0)
        {
            write_data();
        }

        replace_adults();
    }

    write_data();
    write_parameters();
}
