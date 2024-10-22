#ifndef ELECTROMAGN_H
#define ELECTROMAGN_H

#include <vector>
#include <string>
#include <map>

#include "Tools.h"
#include "Profile.h"


class PicParams;
class Species;
class Projector;
class Field;
class SmileiMPI;
class Solver;

//! class ElectroMagn: generic class containing all information on the electromagnetic fields and currents

class ElectroMagn
{

public:
    //! Constructor for Electromagn
    ElectroMagn( PicParams &params, InputData &input_data, SmileiMPI* smpi );

    //! Destructor for Electromagn
    virtual ~ElectroMagn();

    std::vector<unsigned int> dimPrim;
    std::vector<unsigned int> dimDual;

    std::vector<unsigned int> dim_global;



    std::vector<unsigned int> index_bc_min;
    std::vector<unsigned int> index_bc_max;

    //! time-step (from picparams)
    const double timestep;

    //! cell length (from picparams)
    const std::vector<double> cell_length;

    //! \todo Generalise this to none-cartersian geometry (e.g rz, MG & JD)

    //! x-component of the electric field
    Field* Ex_;

    //! y-component of the electric field
    Field* Ey_;

    //! z-component of the electric field
    Field* Ez_;

    //! x-component of the magnetic field
    Field* Bx_;

    //! y-component of the magnetic field
    Field* By_;

    //! z-component of the magnetic field
    Field* Bz_;

    //! x-component of the time-centered magnetic field
    Field* Bx_m;

    //! y-component of the time-centered magnetic field
    Field* By_m;

    //! z-component of the time-centered magnetic field
    Field* Bz_m;

    //! x-component of the total charge current
    Field* Jx_;

    //! y-component of the total charge current
    Field* Jy_;

    //! z-component of the total charge current
    Field* Jz_;

    //! Total charge density
    Field* rho_;

    //! time-average x-component of the electric field
    Field* Ex_avg;

    //! time-average y-component of the electric field
    Field* Ey_avg;

    //! time-average z-component of the electric field
    Field* Ez_avg;

    //! time-average x-component of the magnetic field
    Field* Bx_avg;

    //! time-average y-component of the magnetic field
    Field* By_avg;

    //! time-average z-component of the magnetic field
    Field* Bz_avg;


    //! the global electric potential for solving poisson eqution using SuperLU
    Field* phi_global;

    //! the global Ex for solving poisson eqution using SuperLU
    Field* Ex_global;

    //! the global Ey for solving poisson eqution using SuperLU
    Field* Ey_global;

    //! the global Ez for solving poisson eqution using SuperLU
    Field* Ez_global;

    //! the global rho for solving poisson eqution using SuperLU
    Field* rho_global;

    //! all Fields in electromagn (filled in ElectromagnFactory.h)
    std::vector<Field*> allFields;

    //! all Fields in electromagn (filled in ElectromagnFactory.h)
    std::vector<Field*> allFields_avg;

    //! Vector of charge density and currents for each species
    const unsigned int n_species;
    std::vector<Field*> Jx_s;
    std::vector<Field*> Jy_s;
    std::vector<Field*> Jz_s;
    std::vector<Field*> rho_s;

    //! nDim_field (from params)
    const unsigned int nDim_field;

    //! Volume of the single cell (from params)
    const double cell_volume;

    //! n_space (from params) always 3D
    const std::vector<unsigned int> n_space;
    const std::vector<unsigned int> n_space_global;
    //! Index of starting elements in arrays without duplicated borders
    //! By constuction 1 element is shared in primal field, 2 in dual
    //! 3 : Number of direction (=1, if dim not defined)
    //! 2 : isPrim/isDual
    unsigned int istart[3][2];
    //! Number of elements in arrays without duplicated borders
    unsigned int bufsize[3][2];

    //!\todo should this be just an integer???
    //! Oversize domain to exchange less particles (from params)
    const std::vector<unsigned int> oversize;

    //! Method used to dump data contained in ElectroMagn
    void dump();

    //! Method used to initialize the total charge currents and densities
    virtual void restartRhoJ() = 0;
    //! Method used to initialize the total charge currents and densities of species
    virtual void restartRhoJs(int ispec, bool currents) = 0;

    //! Method used to initialize the total charge density
    void initRhoJ(std::vector<Species*>& vecSpecies, Projector* Proj);

    //! Method used to sum all species densities and currents to compute the total charge density and currents
    virtual void computeTotalRhoJ() = 0;

    //! \todo check time_dual or time_prim (MG)
    //! method used to solve Maxwell's equation (takes current time and time-step as input parameter)
    virtual void solveMaxwellAmpere() = 0;
    //! Maxwell Faraday Solver
    Solver* MaxwellFaradaySolver_;
    virtual void saveMagneticFields() = 0;
    virtual void centerMagneticFields() = 0;

    virtual void incrementAvgFields(unsigned int time_step, unsigned int ntime_step_avg) = 0;

    //! compute Poynting on borders
    virtual void computePoynting() = 0;

    //! pointing vector on borders
    //! 1D: poynting[0][0]=left , poynting[1][0]=right
    //! 2D: poynting[0][0]=west , poynting[1][0]=east
    //!     poynting[1][0]=south, poynting[1][0]=north
    std::vector<double> poynting[2];

    //same as above but instantaneous
    std::vector<double> poynting_inst[2];

    //! Check if norm of charge denisty is not null
    bool isRhoNull(SmileiMPI* smpi);

    double computeNRJ(unsigned int shift, SmileiMPI *smpi);
    double getLostNrjMW() const {return nrj_mw_lost;}

    double getNewFieldsNRJ() const {return nrj_new_fields;}

    void reinitDiags() {
        nrj_mw_lost = 0.;
        nrj_new_fields = 0.;
    }

    inline int getMemFootPrint() {

	// Size of temporary arrays in Species::createParticles
	/*
	int N1(1), N2(1);
	if (nDim_field>1) {
	    N1 = dimPrim[1];
	    if (nDim_field>2) N2 = dimPrim[2];
	}
	int tmpArrayInit = (dimPrim[0]*N1*N2)*sizeof(double) //
	    + dimPrim[0]*sizeof(double**)
	    + dimPrim[0] * N1 * sizeof(double*);
	tmpArrayInit *= 9;
	std::cout << tmpArrayInit  << std::endl;
	*/

	int emSize = 9+4; // 3 x (E, B, Bm) + 3 x J, rho
	if (true) // For now, no test to compute or not per species
	    emSize += n_species * 4; // 3 x J, rho
	if (true) // For now, no test to compute or not average
	    emSize += 6; // 3 x (E, B)

	for (size_t i=0 ; i<nDim_field ; i++)
	    emSize *= dimPrim[i];

	emSize *= sizeof(double);
	return emSize;
    }


protected:
    //! Vector of boundary-condition per side for the fields
    //std::vector<ElectroMagnBC*> emBoundCond;

private:

    //! Accumulate nrj lost with moving window
    double nrj_mw_lost;

    //! Accumulate nrj added with new fields
    double nrj_new_fields;

};

#endif
