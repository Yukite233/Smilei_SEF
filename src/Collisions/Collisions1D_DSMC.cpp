#include "Collisions1D_DSMC.h"
#include "SmileiMPI.h"
#include "Field2D.h"
#include "H5.h"


#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <sstream>

using namespace std;


// Constructor
Collisions1D_DSMC::Collisions1D_DSMC(PicParams& param, vector<Species*>& vecSpecies, SmileiMPI* smpi,
                       unsigned int n_collisions,
                       vector<unsigned int> species_group1,
                       vector<unsigned int> species_group2,
                       double coulomb_log,
                       bool intra_collisions,
                       int debug_every)
{

    n_collisions    = (n_collisions    );
    species_group1  = (species_group1  );
    species_group2  = (species_group2  );
    coulomb_log     = (coulomb_log     );
    intra_collisions= (intra_collisions);
    debug_every     = (debug_every     );
    start           = (0               );



    // Calculate total number of bins
    int nbins = vecSpecies[0]->bmin.size();
    totbins = nbins;
    //MPI_Allreduce( smpi->isMaster()?MPI_IN_PLACE:&totbins, &totbins, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
    MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&totbins, &totbins, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);

    // if debug requested, prepare hdf5 file
    fileId = 0;
    if( debug_every>0 ) {
        ostringstream mystream;
        mystream.str("");
        mystream << "Collisions1D" << n_collisions << ".h5";
        // Create the HDF5 file
        hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
        fileId = H5Fcreate(mystream.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, pid);
        H5Pclose(pid);
        // write all parameters as HDF5 attributes
        string ver(__VERSION);
        H5::attr(fileId, "Version", ver);
        mystream.str("");
        mystream << species_group1[0];
        for(unsigned int i=1; i<species_group1.size(); i++) mystream << "," << species_group1[i];
        H5::attr(fileId, "species1" , mystream.str());
        mystream.str("");
        mystream << species_group2[0];
        for(unsigned int i=1; i<species_group2.size(); i++) mystream << "," << species_group2[i];
        H5::attr(fileId, "species2" , mystream.str());
        H5::attr(fileId, "coulomb_log" , coulomb_log);
        H5::attr(fileId, "debug_every"  , debug_every);

        // Find out where the proc will start writing in the overall array
        MPI_Status status;
        // Receive the location where to start from the previous node
        if (smpi->getRank()>0) MPI_Recv( &(start), 1, MPI_INTEGER, smpi->getRank()-1, 0, MPI_COMM_WORLD, &status );
        // Send the location where to end to the next node
        int end = start+nbins;
        if (smpi->getRank()!=smpi->getSize()-1) MPI_Send( &end, 1, MPI_INTEGER, smpi->getRank()+1, 0, MPI_COMM_WORLD );
    }

}

Collisions1D_DSMC::~Collisions1D_DSMC()
{
    if (fileId != 0) H5Fclose(fileId);
}


// Declare other static variables here
bool               Collisions1D::debye_length_required;
vector<double>     Collisions1D::debye_length_squared;



// Calculates the debye length squared in each cluster
// The formula for the inverse debye length squared is sumOverSpecies(density*charge^2/temperature)
void Collisions1D_DSMC::calculate_debye_length(PicParams& params, vector<Species*>& vecSpecies)
{

    // get info on particle binning
    unsigned int nbins = vecSpecies[0]->bmin.size(); // number of bins
    unsigned int nspec = vecSpecies.size(); // number of species
    unsigned int bmin, bmax;
    double p2, density, density_max, charge, temperature, rmin2;
    Species   * s;
    Particles * p;
    double coeff = params.wavelength_SI/(6.*M_PI*2.8179403267e-15); // normLength/(3*electronRadius) = wavelength/(6*pi*electronRadius)

    debye_length_squared.resize(nbins);

    // Loop on bins
    //! \todo Make OpenMP parallelization (MF & JD)
    for (unsigned int ibin=0 ; ibin<nbins ; ibin++) {

        density_max = 0.;
        debye_length_squared[ibin] = 0.;
        for (unsigned int ispec=0 ; ispec<nspec ; ispec++) { // loop all species
            // Calculation of particles density, mean charge, and temperature
            // Density is the sum of weights
            // Temperature basic definition is the average <v*p> divided by 3
            //    (instead of v*p, we use p^2/gamma)
            s  = vecSpecies[ispec];
            p  = &(s->particles);
            bmin = s->bmin[ibin];
            bmax = s->bmax[ibin];
            density     = 0.;
            charge      = 0.;
            temperature = 0.;
            // loop particles to calculate average quantities
            for (unsigned int iPart=bmin ; iPart<bmax ; iPart++ ) {
                p2 = pow(p->momentum(0,iPart),2)+pow(p->momentum(1,iPart),2)+pow(p->momentum(2,iPart),2);
                density     += p->weight(iPart);
                charge      += p->weight(iPart) * p->charge(iPart);
                temperature += p->weight(iPart) * p2/sqrt(1.+p2);
            }
            if (density <= 0.) continue;
            charge /= density; // average charge
            temperature *= (s->species_param.mass) / (3.*density); // Te in units of me*c^2
            density /= params.n_cell_per_cluster; // density in units of critical density
            // compute inverse debye length squared
            if (temperature>0.) debye_length_squared[ibin] += density*charge*charge/temperature;
            // compute maximum density of species
            if (density>density_max) density_max = density;
        }

        // if there were particles,
        if (debye_length_squared[ibin] > 0.) {
            // compute debye length squared in code units
            debye_length_squared[ibin] = 1./(debye_length_squared[ibin]);
            // apply lower limit to the debye length (minimum interatomic distance)
            rmin2 = pow(coeff*density_max, -2./3.);
            if (debye_length_squared[ibin] < rmin2) debye_length_squared[ibin] = rmin2;
        }

    }

#ifdef  __DEBUG
    // calculate and print average debye length
    double mean_debye_length = 0.;
    for (unsigned int ibin=0 ; ibin<nbins ; ibin++)
        mean_debye_length += sqrt(debye_length_squared[ibin]);
    mean_debye_length /= (double)nbins;
    //DEBUG("Mean Debye length in code length units = " << scientific << setprecision(3) << mean_debye_length);
    mean_debye_length *= params.wavelength_SI/(2.*M_PI); // switch to SI
    DEBUG("Mean Debye length in meters = " << scientific << setprecision(3) << mean_debye_length );
#endif

}


// Calculates the collisions for a given Collisions1D object
void Collisions1D_DSMC::collide(PicParams& params, vector<Species*>& vecSpecies, int itime)
{

    unsigned int nbins = vecSpecies[0]->bmin.size(); // number of bins
    vector<unsigned int> *sg1, *sg2, *sgtmp, bmin1, bmax1, bmin2, bmax2, index1, index2;


    vector<vector<<vector int> > > index1, index2;
    vector<int> index1_tot, index2_tot;
    vector<int> n1, n2;
    vector<double> momentum_unit(3, 0.0), momentum_temp(3, 0.0);
    int idNew;
    int totNCollision = 0;
    vector<int> bmin1, bmax1, bmin2, bmax2, bmin3, bmax3;





    unsigned int nspec1, nspec2; // numbers of species in each group
    unsigned int npart1, npart2; // numbers of macro-particles in each group
    unsigned int npairs; // number of pairs of macro-particles
    vector<unsigned int> np1, np2; // numbers of macro-particles in each species, in each group
    double n1, n2, n12, n123, n223; // densities of particles
    unsigned int i1, i2, ispec1, ispec2;
    Species   *s1, *s2;
    Particles *p1, *p2;
    double m1, m2, m12, W1, W2, qqm, qqm2, gamma1, gamma2, gamma12, gamma12_inv,
           COM_vx, COM_vy, COM_vz, COM_vsquare, COM_gamma,
           term1, term2, term3, term4, term5, term6, coeff1, coeff2, coeff3, coeff4, twoPi,
           vcv1, vcv2, px_COM, py_COM, pz_COM, p2_COM, p_COM, gamma1_COM, gamma2_COM,
           logL, bmin, s, vrel, smax,
           cosX, sinX, phi, sinXcosPhi, sinXsinPhi, p_perp, inv_p_perp,
           newpx_COM, newpy_COM, newpz_COM, U, vcp;
    Field2D *smean, *logLmean, *ncol;//, *temperature
    ostringstream name;
    hid_t did;

    sg1 = &species_group1;
    sg2 = &species_group2;


    // Loop on bins
    for (unsigned int ibin=0 ; ibin<nbins ; ibin++) {


        for (ispec1=0 ; ispec1<nspec1 ; ispec1++) {
            s1 = vecSpecies[(*sg1)[ispec1]];
            bmin1[ispec1] = s1->bmin[ibin];
            bmax1[ispec1] = s1->bmax[ibin];
            np1[ispec1] = bmax1[ispec1] - bmin1[ispec1];
            npart1 += np1[ispec1];
        }



        //>calculate the particle number of species1 in each cell, and the indexs of particles in the cell
        index1.resize(params.n_space[1]);
        np1.resize(nspec1);
        for(int i = 0; i < index1.size(); i++)
        {
            index1[i].resize(nspec1);
            for(int j = 0; j < index1[i].size(); j++)
            {
                index1[i][j].resize(0);
            }
        }


        //> calculate the particle number of each species in each cell, and the indexs of particles in the cell
        for (ispec1=0 ; ispec1<nspec1 ; ispec1++)
        {
            s1 = vecSpecies[(*sg1)[ispec1]];
            p1 = &(s1->particles);
            m1 = s1->species_param.mass;
            W1 = p1->weight(i1);
            bmin1 = s1->bmin;
            bmax1 = s1->bmax;

            for(int iPart = bmin1[ibin]; iPart < bmax1[ibin]; iPart++)
            {
                double ypn = p1->position(1, ipart)*dy_inv_;
                jp_ = floor(ypn);
                jp_ = jp_ - j_domain_begin;
                index1[jp_][ispec1].push_back(iPart);
            }
        }



        // Now start the real loop
        //> the code is mainly from the website: https://www.particleincell.com/2012/html5-dsmc/
        //> the original code language is javascript
        // ----------------------------------------------------
        for (int iCell = 0; iCell < params.n_space[1]; iCell++)
        {
            sigma_cr_max_temp = 0.0;
            //> calculate the total number of all neutral particles (all species) in the cell
            npart1 = 0;
            for(ispec1=0 ; ispec1<nspec1 ; ispec1++)
            {
                np1[ispec1] = index1[iCell][ispec1].size();
                npart1 += np1[ispec1];
            }
            //> for now, the weights of all neutral species are the same (W1)
            npairs = int (0.5 * npart1 * npart1 * W1 * sigma_cr_max * params.timestep  / params.cell_volume);

            //> for updating sigma_cr_max
            sigma_cr_max_temp = 0;

            for(int i =0; i < npairs; i++)
            {
                //> select the first particle
                i1 = int( ((double)rand() / RAND_MAX) * npart1 );
                //> select the second one, making sure the two particles are unique
                do{
                    i2 = int( ((double)rand() / RAND_MAX) * npart1 );
                }
                while(i1 == i2);


                // find species and index i1 of particle "1"
                for (ispec1=0 ; ispec1<nspec1 ; ispec1++) {
                    if (i1 < np1[ispec1]) break;
                    i1 -= np1[ispec1];
                }
                i1 += bmin1[ispec1];
                // find species and index i2 of particle "2"
                i2 = index2[i];
                for (ispec2=0 ; ispec2<nspec1 ; ispec2++) {
                    if (i2 < np1[ispec2]) break;
                    i2 -= np1[ispec2];
                }
                i2 += bmin2[ispec2];

                //> only one species group, but one or two particles (the particle is a c++ class)
                s1 = vecSpecies[(*sg1)[ispec1]]; s2 = vecSpecies[(*sg1)[ispec2]];
                p1 = &(s1->particles);           p2 = &(s2->particles);
                m1 = s1->species_param.mass;     m2 = s2->species_param.mass;
                W1 = p1->weight(i1);             W2 = p2->weight(i2);

                //> relative velocity
                cr = relative_velocity(p1, i1, p2, i2);
                sigma = evalSigma(cr);
                sigma_cr = sigma * cr;
                if(sigma_cr > sigma_cr_max_temp){
                    sigma_cr_max_temp = sigma_cr
                }
                P_collision = sigma_cr / sigma_cr_max;
                double ran_p = (double)rand() / RAND_MAX;

                if(ran_P < P_collision){
                    scatter_particles(p1, i1, p2, i2);
                }

            }


        } //> end loop on cells




    } // end loop on bins

}


inline double relative_velocity(Particles* particle1, int iPart1, Particles* particle2, iPart2)
{
    double rv;
    rv = sqrt( pow((particle1->momentum(0, iPart1) - particle2->momentum(0, iPart2)), 2)
            +  pow((particle1->momentum(1, iPart1) - particle2->momentum(1, iPart2)), 2)
            +  pow((particle1->momentum(2, iPart1) - particle2->momentum(2, iPart2)), 2)
            );
    return rv;
}



inline double scatter_particles(Particles* particle1, int iPart1, Particles* particle2, iPart2)
{
    double rv;
    rv = sqrt( pow((particle1->momentum(0, iPart1) - particle2->momentum(0, iPart2)), 2)
            +  pow((particle1->momentum(1, iPart1) - particle2->momentum(1, iPart2)), 2)
            +  pow((particle1->momentum(2, iPart1) - particle2->momentum(2, iPart2)), 2)
            );
    return rv;
}



// Technique given by Nanbu in http://dx.doi.org/10.1103/PhysRevE.55.4642
//   to pick randomly the deflection angle cosine, in the center-of-mass frame.
// It involves the "s" parameter (~ collision frequency * deflection expectation)
//   and a random number "U".
// Technique slightly modified in http://dx.doi.org/10.1063/1.4742167
inline double Collisions1D_DSMC::cos_chi(double s)
{

    double A, invA;
    //!\todo make a faster rand by preallocating ??
    double U = (double)rand() / RAND_MAX;

    if( s < 0.1 ) {
        if ( U<0.0001 ) U=0.0001; // ensures cos_chi > 0
        return 1. + s*log(U);
    }
    if( s < 3.  ) {
        // the polynomial has been modified from the article in order to have a better form
        invA = 0.00569578 +(0.95602 + (-0.508139 + (0.479139 + ( -0.12789 + 0.0238957*s )*s )*s )*s )*s;
        A = 1./invA;
        return  invA  * log( exp(-A) + 2.*U*sinh(A) );
    }
    if( s < 6.  ) {
        A = 3.*exp(-s);
        return (1./A) * log( exp(-A) + 2.*U*sinh(A) );
    }
    return 2.*U - 1.;

}
