/*
 * SmileiIO.cpp
 *
 *  Created on: 3 juil. 2013
 */

#include "SmileiIO.h"

#include <sstream>
#include <iomanip>

#include <mpi.h>

#include "PicParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Species.h"

using namespace std;

// static varable must be defined and initialized here

SmileiIO::SmileiIO( PicParams& params, SmileiMPI* smpi )
{
    // Fields_global.h5
    global_file_id_  = H5Fcreate( "Fields_global.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    ndims_t = 0;

}

SmileiIO::~SmileiIO()
{
    // Management of global IO file
    H5Fclose( global_file_id_ );
}
