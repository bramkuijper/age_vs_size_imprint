#ifndef AGEMATSOL_HPP_ 
#define AGEMATSOL_HPP_

#include <vector>

class AgeMatSol
{
    private:
        double nm;
        double nf;
        double df;
        double dm;
        double l;

        // current values of af and am
        double af;
        double am;

        // vector with the equilibrium values of af and am
        // for different genetic systems and expression levels
        std::vector <double> af_vec;
        std::vector <double> am_vec;


        // the solver which takes two functions as their value and solves them
        void solveAfAm(
                double (*af_selgrad)() // function specifying the selection gradient on af
                ,double (*am_selgrad)() // function specifying the selection gradient on am
                );

        selgradfunc 

    public:
        AgeMatSol(int argc, char **argv);



#endif


