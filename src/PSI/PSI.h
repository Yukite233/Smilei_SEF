/*

PSI class - Frederic Perez - 03/2015

This is based on the work described here
http://dx.doi.org/10.1063/1.4742167

Binary collisions, between macro-particles, are treated according to a scheme
first described by Nanbu (http://dx.doi.org/10.1103/PhysRevE.55.4642).

To include collisions in the simulations, add a block in the input file,
similar to the following:

# PSI
# species1    = list of strings, the names of the first species that collide
# species2    = list of strings, the names of the second species that collide
#               (can be the same as species1)
# coulomb_log = float, Coulomb logarithm. If negative or zero, then automatically computed.
PSI(
	species1 = ["ion1"],
	species2 = ["electron1"],
	coulomb_log = 2.0
)

Several collision types can be defined. For each type, add a group "PSI()".

*/

#ifndef PSI_H
#define PSI_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "H5.h"

class PSI
{

public:
    //! Constructor for PSI between two species
    PSI(){};
    ~PSI(){};

    //! Group of the species numbers that are associated for PSI.
    //> for PSI_Injection, only species1 is used;
    //> for sputtering and secondary electron emission, species1 is the incident particle. 
    std::vector<unsigned int> species1, species2;


    //! Method called in the main smilei loop to apply collisions at each timestep
    virtual void performPSI(PicParams&,std::vector<Species*>&,int){};





private:



};


#endif
