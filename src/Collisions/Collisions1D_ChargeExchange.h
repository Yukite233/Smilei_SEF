/*
Collisions1D class
*/

#ifndef COLLISIONS1D_CHARGEEXCHANGE_H
#define COLLISIONS1D_CHARGEEXCHANGE_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "Collisions1D.h"
#include "H5.h"

class Collisions1D_ChargeExchange : public Collisions1D
{

public:
    //! Constructor for Collisions1D between two species
    Collisions1D_ChargeExchange(PicParams&,std::vector<Species*>&,SmileiMPI*,unsigned int,std::vector<unsigned int>,std::vector<unsigned int>,double,bool,int);
    ~Collisions1D_ChargeExchange();


    //! Coulomb logarithm (zero or negative means automatic)
    double coulomb_log;

    //! Method to calculate the Debye length in each cluster
    void calculate_debye_length(PicParams&,std::vector<Species*>&);

    //! Method called in the main smilei loop to apply collisions at each timestep
    void collide(PicParams&,std::vector<Species*>&,int);

    virtual double cos_chi(double);
private:
    //>the ionization threshold energy
    double energy_ion;

};


#endif
