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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "bramauxiliary.h"


//#define NDEBUG

using namespace std;

// random number generator 
// see http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html#Random-Number-Generation 
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *r; // gnu scientific rng 

const size_t Npatches = 1000; 
const size_t Nfp = 5; // females per patch
const size_t Nmp = 5; // males per patch
const size_t numgen = 30000;
const size_t Clutch = 500;

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
size_t skip = 10;

// haploid individual
struct Individual
{
    double af[2];
    double am[2];
    double phen;
};

struct Patch
{
    Individual localsF[Nfp]; // all the local female breeders
    Individual localsM[Nmp]; // all the local male breeders

    // philopatric female offspring
    Individual philsF[Nfp * Clutch];     

    // philopatric male offspring
    Individual philsM[Nfp * Clutch];     

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
Individual DispersersF[Npatches * Nfp * Clutch];
Individual DispersersM[Npatches * Nfp * Clutch];

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
            MetaPop[i].localsF[j].af[0] = 1.5;
            MetaPop[i].localsF[j].af[1] = 1.5;
            MetaPop[i].localsF[j].am[0] = 0.75;
            MetaPop[i].localsF[j].am[1] = 0.75;
            MetaPop[i].localsF[j].phen = 3.0;
        }
        
        for (size_t j = 0; j < Nmp; ++j)
        {
            // initialize phenotypes
            MetaPop[i].localsM[j].af[0] = 1.5;
            MetaPop[i].localsM[j].af[1] = 1.5;
            MetaPop[i].localsM[j].am[0] = 0.75;
            MetaPop[i].localsM[j].am[1] = 0.75;
            MetaPop[i].localsM[j].phen = 1.5;
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
void Create_Kid(Individual &mother, Individual &father, Individual &Kid, bool is_female)
{
    Kid.phen = -1;
    // first inherit af
    
    // maternal allele
    size_t maternal_allele_f = gsl_rng_uniform_int(r,2);
    Kid.af[0] = Mut(mother.af[maternal_allele_f]); 

    // paternal allele
    size_t paternal_allele_f = 0;

    if (diploid)
    {
        paternal_allele_f = gsl_rng_uniform_int(r,2);
        Kid.af[1] = Mut(father.af[paternal_allele_f]);

        // if female, determine phenotype
        if (is_female)
        {
            switch(type)
            {
                case 1: // madumnal allele
                    Kid.phen = mother.af[maternal_allele_f];
                    break;
                case 2: // padumnal allele
                    Kid.phen = father.af[paternal_allele_f];
                    break;
                default: // normal allele
                    Kid.phen = Kid.af[0] + Kid.af[1];
                    break;
            }
        }
    }
    else if (is_female) // haplodiploid and female
    {
        Kid.af[1] = Mut(father.af[paternal_allele_f]);
        
        switch(type)
        {
            case 1: // madumnal allele
                Kid.phen = mother.af[maternal_allele_f];
                break;
            case 2: // padumnal allele
                Kid.phen = father.af[paternal_allele_f];
                break;
            default: // normal allele
                Kid.phen = Kid.af[0] + Kid.af[1];
                break;
        }
    } // we're done here; male's phenotypes come later

    // inherit am
    size_t maternal_allele_m = gsl_rng_uniform_int(r,2);
    Kid.am[0] = Mut(mother.am[maternal_allele_m]); 

    // paternal allele am
    size_t paternal_allele_m = 0;

    if (diploid)
    {
        paternal_allele_m = gsl_rng_uniform_int(r, 2);
        Kid.am[1] = Mut(father.am[paternal_allele_m]);

        if (!is_female)
        {
            switch(type)
            {
                case 1: // madumnal allele
                    Kid.phen = mother.am[maternal_allele_m];
                    break;
                case 2: // padumnal allele
                    Kid.phen = father.am[paternal_allele_m];
                    break;
                default: // normal allele
                    Kid.phen = Kid.am[0] + Kid.am[1];
                    break;
            }
        }
    }
    else if (is_female) // haplodiploid female
    {
        Kid.am[1] = Mut(father.am[paternal_allele_m]);
    }
    else // haplodiploid male
    {
        switch(type)
        {
            case 1: // madumnal allele
                Kid.phen = mother.am[maternal_allele_m];
                break;
            case 2: // padumnal allele
                Kid.phen = Kid.am[0];
                break;
            default: // normal allele
                Kid.phen = Kid.am[0] + Kid.am[1];
                break;
        }
    }

    assert(Kid.phen >= 0);
}

// mate and create kids across all patches...
void make_juveniles()
{
    // reset counters of number of
    // dispersing individuals
    NdispF = 0;
    NdispM = 0;

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
        // reset local offspring counter
        MetaPop[i].NkidsF = 0;
        MetaPop[i].NkidsM = 0;

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

        // reproduce
        for (size_t j = 0; j < Nfp; ++j)
        {
            // select male based on his reproductive effort
            cumul_deviate = gsl_rng_uniform(r) * sum_male_effort;

            father = gsl_rng_uniform_int(r, Nmp);

            for (size_t male_i = 0; male_i < Nmp; ++male_i)
            {
                if (cumul_deviate <= cumul_male_effort[male_i])
                {
                    father = male_i;
                    break;
                }
            }

            // calculate female clutch size
            exact_clutch = pow(k * MetaPop[i].localsF[j].phen,3) * Clutch; 
            clutch_size = floor(exact_clutch);
           
            // stochastic rounding
            if (gsl_rng_uniform(r) < exact_clutch - clutch_size)
            {
                ++clutch_size;
            }


            // create kids and reduce parental resources
            for (size_t k = 0; k < clutch_size; ++k)
            {
                Individual Kid; 

                is_female = gsl_rng_uniform(r) < 0.5;

                Create_Kid(MetaPop[i].localsF[j],MetaPop[i].localsM[father],Kid, is_female);

    
                
                // juvenile survival
                if (gsl_rng_uniform(r) < exp(-q * Kid.phen))
                {
                    // female or male
                    if (is_female)
                    {
                        ++NsurvF;
                        // dispersal or not
                        if (gsl_rng_uniform(r) < df)
                        {                            
                            DispersersF[NdispF++] = Kid;
                        }
                        else
                        {
                            MetaPop[i].philsF[MetaPop[i].NkidsF++] = Kid;
                        }
                    }
                    else
                    {
                        ++NsurvM;
                        // dispersal or not
                        if (gsl_rng_uniform(r) < dm)
                        {                            
                            DispersersM[NdispM++] = Kid;
                        }
                        else
                        {
                            MetaPop[i].philsM[MetaPop[i].NkidsM++] = Kid;
                        }
                    }
                }// end if (survival)
            } // end clutch
        } // end Nfp
    }//Npatches

    assert(NdispF < Npatches * Nfp * Clutch);
    assert(NdispM < Npatches * Nfp * Clutch);
}

// replacement of adults with juveniles
void replace_adults()
{
    // okay, we have Ndisp/Npatches dispersers per patch
    size_t dispersers_per_patchF = floor((double) NdispF / Npatches);
    size_t dispersers_per_patchM = floor((double) NdispM / Npatches);

    // however, we need to correctly round this rational number
    // to the actual number of immigrants for 
    // a given patch. To this end, 
    // we just randomly distribute the dispersers that remain after the
    // previous rounding over the different patches
    for (size_t i = 0; i < NdispF - Npatches * dispersers_per_patchF; ++i)
    {
        // randomly picked patch receives additional immigrant
        MetaPop[gsl_rng_uniform_int(r,Npatches)].immigrant_bonusF++;
    }
    for (size_t i = 0; i < NdispM - Npatches * dispersers_per_patchM; ++i)
    {
        // randomly picked patch receives additional immigrant
        MetaPop[gsl_rng_uniform_int(r,Npatches)].immigrant_bonusM++;
    }
    
    // now replace local breeders on each patch
    for (size_t i = 0; i < Npatches; ++i)
    {
        if (dispersers_per_patchF + MetaPop[i].immigrant_bonusF > 0)
        {
            assert(NdispF > 0);
        }

        if (dispersers_per_patchM + MetaPop[i].immigrant_bonusM > 0)
        {
            assert(NdispM > 0);
        }

        size_t arriving_immigrantsF = dispersers_per_patchF + MetaPop[i].immigrant_bonusF;
        size_t arriving_immigrantsM = dispersers_per_patchM + MetaPop[i].immigrant_bonusM;

        for (size_t j = 0; j < Nfp; ++j)
        {
            //cout << arriving_immigrantsF << " " << MetaPop[i].NkidsF << endl;
            assert((double) arriving_immigrantsF / (arriving_immigrantsF + MetaPop[i].NkidsF) >= 0.0 
                    && (double) arriving_immigrantsF / (arriving_immigrantsF + MetaPop[i].NkidsF) <= 1.0);
            
            if (gsl_rng_uniform(r) < (double) arriving_immigrantsF / (arriving_immigrantsF + MetaPop[i].NkidsF))
            {
                size_t rand_disp = gsl_rng_uniform_int(r,NdispF);
                MetaPop[i].localsF[j] = DispersersF[rand_disp];
                DispersersF[rand_disp] = DispersersF[NdispF-1];
                --NdispF;
                --arriving_immigrantsF;
            }
            else
            {
                size_t rand_phil = gsl_rng_uniform_int(r,MetaPop[i].NkidsF);
                MetaPop[i].localsF[j] = MetaPop[i].philsF[rand_phil];
                MetaPop[i].philsF[rand_phil] = MetaPop[i].philsF[MetaPop[i].NkidsF-1];
                --MetaPop[i].NkidsF;
            }
        }
        
        for (size_t j = 0; j < Nmp; ++j)
        {
            //cout << seed << " gen: " << generation << " patch: " << i << " nmp: " << j << " nkidsM: " << MetaPop[i].NkidsM << " ndispM: " << NdispM << " nsurvM: " << NsurvM <<  " " << NsurvF << endl;
            assert((double) arriving_immigrantsM / (arriving_immigrantsM + MetaPop[i].NkidsM) >= 0.0 
                    && (double) arriving_immigrantsM / (arriving_immigrantsM + MetaPop[i].NkidsM) <= 1.0);

            if (gsl_rng_uniform(r) < (double) arriving_immigrantsM / (arriving_immigrantsM + MetaPop[i].NkidsM))
            {
                size_t rand_disp = gsl_rng_uniform_int(r,NdispM);
                MetaPop[i].localsM[j] = DispersersM[rand_disp];
                DispersersM[rand_disp] = DispersersM[NdispM-1];
                --NdispM;
                --arriving_immigrantsM;
            }
            else
            {
                size_t rand_phil = gsl_rng_uniform_int(r,MetaPop[i].NkidsM);
                MetaPop[i].localsM[j] = MetaPop[i].philsM[rand_phil];
                MetaPop[i].philsM[rand_phil] = MetaPop[i].philsM[MetaPop[i].NkidsM-1];
                --MetaPop[i].NkidsM;
            }
        }

        // remove all remaining female immigrants from the global dispersal pool
        for (size_t j = 0; j < arriving_immigrantsF; ++j)
        {
            size_t rand_disp = gsl_rng_uniform_int(r,NdispF);
            DispersersF[rand_disp] = DispersersF[NdispF-1];
            --NdispF;
        }
        // remove all remaining male immigrants from the global dispersal pool
        for (size_t j = 0; j < arriving_immigrantsM; ++j)
        {
            size_t rand_disp = gsl_rng_uniform_int(r,NdispM);
            DispersersM[rand_disp] = DispersersM[NdispM-1];
            --NdispM;
        }
    }
    assert(NdispF==0);
    assert(NdispM==0);
}

void write_data_headers()
{
    DataFile << "generation;meanaf;varaf;meanam;varam;" << endl;
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
            meanaf += MetaPop[i].localsF[j].af[0]+ MetaPop[i].localsF[j].af[1];
            ssaf += ( MetaPop[i].localsF[j].af[0]+ MetaPop[i].localsF[j].af[1]) * 
                ( MetaPop[i].localsF[j].af[0]+ MetaPop[i].localsF[j].af[1]);

            meanam += MetaPop[i].localsF[j].am[0]+ MetaPop[i].localsF[j].am[1];
            ssam += ( MetaPop[i].localsF[j].am[0]+ MetaPop[i].localsF[j].am[1]) * 
                ( MetaPop[i].localsF[j].am[0]+ MetaPop[i].localsF[j].am[1]);
        }
        
        for (size_t j = 0; j < Nmp; ++j)
        {
            if (diploid)
            {
                meanaf += MetaPop[i].localsM[j].af[0]+ MetaPop[i].localsM[j].af[1];
                ssaf += ( MetaPop[i].localsM[j].af[0]+ MetaPop[i].localsM[j].af[1]) * 
                    ( MetaPop[i].localsM[j].af[0]+ MetaPop[i].localsM[j].af[1]);

                meanam += MetaPop[i].localsM[j].am[0]+ MetaPop[i].localsM[j].am[1];
                ssam += ( MetaPop[i].localsM[j].am[0]+ MetaPop[i].localsM[j].am[1]) * 
                    ( MetaPop[i].localsM[j].am[0]+ MetaPop[i].localsM[j].am[1]);
            } else
            {
                meanaf += MetaPop[i].localsM[j].af[0];
                ssaf += MetaPop[i].localsM[j].af[0] * MetaPop[i].localsM[j].af[0];

                meanam += MetaPop[i].localsM[j].am[0];
                ssam += MetaPop[i].localsM[j].am[0] * MetaPop[i].localsM[j].am[0];
            }
        }
    }

    double varaf, varam;
    
    if (diploid)
    {
        meanaf /= Npatches * (Nfp + Nmp);
        meanam /= Npatches * (Nfp + Nmp);
        varaf = ssaf / (Npatches * (Nfp + Nmp)) - meanaf * meanaf;
        varam = ssam / (Npatches * (Nfp + Nmp)) - meanam * meanam;
    }
    else
    {
        meanaf /= Npatches * (2*Nfp + Nmp);
        meanam /= Npatches * (2*Nfp + Nmp);
        varaf = ssaf / (Npatches * (2*Nfp + Nmp)) - meanaf * meanaf;
        varam = ssam / (Npatches * (2*Nfp + Nmp)) - meanam * meanam;
    }

        
    DataFile << generation 
                        << ";" << meanaf 
                        << ";" << varaf 
                        << ";" << meanam 
                        << ";" << varam  << ";" << endl;
}

void write_parameters()
{
    DataFile << endl << endl << "system;" << (diploid ? "diploid" : "haplodiploid") << endl
                << "imprint;" << (type == 0 ? "offspring" : (type == 1 ? "madumnal" : "padumnal")) << endl
                << "patch;" << Npatches << endl
                << "nfp;" << Nfp << endl
                << "nmp;" << Nmp << endl
                << "df;" << df << endl
                << "dm;" << dm << endl
                << "numgen;" << numgen << endl
                << "mu;" << mu << endl
                << "sdmu;" << sdmu << endl
                << "q;" << q << endl
                << "k;" << k << endl
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
