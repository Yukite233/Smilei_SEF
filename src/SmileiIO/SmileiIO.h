/*
 * SmileiIO.h
 *
 *  Created on: 3 juil. 2013
 */

#ifndef SMILEIIO_H
#define SMILEIIO_H

#include <string>
#include <vector>

#include <hdf5.h>
#include <Tools.h>

class PicParams;
class InputData;
class SmileiMPI;
class ElectroMagn;
class Field;
class Species;

#include <csignal>

//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiIO
//  --------------------------------------------------------------------------------------------------------------------
class SmileiIO {
public:
    //! Create // HDF5 environment
    //! @see global_file_id_
    //! @see global_file_id_avg
    SmileiIO( PicParams& params, SmileiMPI* smpi );
    //! Destructor for SmileiIO
    virtual ~SmileiIO();

    //! Basic write field on its own file (debug)
    virtual void write( ElectroMagn* fields, SmileiMPI* smpi ){};


    //! Id of "Fields.h5", contains all fields per timestep
    hid_t global_file_id_;

    hsize_t     dims_global[4];
    double* data_;

    //! Id of "Fields_avg.h5", contains time-averaged fields per timestep
    hid_t global_file_id_avg;

    //! Property list for collective dataset write, set for // IO.
    hid_t write_plist;

    //! Id of "particles-mpirank.h5", contains particles of current mpirank
    //! Disabled for now
    hid_t  partFile_id;

    int ndims_t;
    int ndims_[4];
    //! Space dimension of a particle
    unsigned int nDim_particle;



private:
    //! incremental number of times we've done a dump
    unsigned int dump_times;

    //! name of the fields to dump
    std::vector<std::string> fieldsToDump;


};

#endif /* SMILEI_OUTPUT_H_ */
