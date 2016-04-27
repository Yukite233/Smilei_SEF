/*
 * SmileiIO_Cart1D.cpp
 *
 *  Created on: 3 juil. 2013
 */

#include "SmileiIO_Cart1D.h"

#include <sstream>

#include "PicParams.h"
#include "SmileiMPI_Cart1D.h"
#include "Field1D.h"
#include "ElectroMagn.h"

using namespace std;

SmileiIO_Cart1D::SmileiIO_Cart1D( PicParams& params, SmileiMPI* smpi )
: SmileiIO( params, smpi )
{
    createPattern(params,smpi);
    offset[0] = 0;
    offset[1] = 0;
    offset[2] = 0;
    offset[3] = 0;

    stride[0] = 1;
    stride[1] = 1;
    stride[2] = 1;
    stride[3] = 1;

    block[0] = 1;
    block[1] = 1;
    block[2] = 1;
    block[3] = 1;

}

SmileiIO_Cart1D::~SmileiIO_Cart1D()
{
}

//> create hdf5 data hierarchical structure: datespace, dateset and so on
void SmileiIO_Cart1D::createPattern( PicParams& params, SmileiMPI* smpi )
{
    hsize_t     dims;

    dims_global[3] = params.n_space_global[0] + 1;
    dims_global[2] = 1;
    dims_global[1] = 1;
    dims_global[0] = params.n_time / params.dump_step;

    ndims_[0] = dims_global[0];
    ndims_[1] = dims_global[1];
    ndims_[2] = dims_global[2];
    ndims_[3] = dims_global[3];

    data_ =  (double*)malloc(dims_global[3] * dims_global[2] * dims_global[1] * dims_global[0] * sizeof(double));
    for( int i = 0; i < dims_global[3] * dims_global[2] * dims_global[1] * dims_global[0]; i++)
    {
      data_[i] = 20.0;
    }


    group_id = H5Gcreate(global_file_id_, "/1d_global", H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
    H5Gcreate(global_file_id_, "/2d_global", H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);

    /* Create a datagroup attribute. */
    dims = 4;
    //> dataspace is to descript the structure of data: the number of data dimension and the size of each dimension
    //> the first parameter rank=1: is the number of dimensions used in the dataspace
    dataspace_id = H5Screate_simple(1, &dims, NULL);
    attribute_id = H5Acreate2 (group_id, "dims_global", H5T_STD_I32BE, dataspace_id,
                                 H5P_DEFAULT, H5P_DEFAULT);
    /* Write the attribute data. */
    status = H5Awrite(attribute_id, H5T_NATIVE_INT, ndims_);
    /* Close the attribute. */
    status = H5Aclose(attribute_id);



    //> output initial 1d_global data===========================================
    //> if without below process, the method write() will go wrong, no ideas now!!!

    /* Create the data space for the dataset. */
    dataspace_id = H5Screate_simple(4, dims_global, NULL);

    /* Create the second dataset in group "Group_A". */
    dataset_id = H5Dcreate2(group_id, "potential", H5T_NATIVE_DOUBLE, dataspace_id,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    data_);
    status = H5Dclose(dataset_id);


    dataset_id = H5Dcreate2(group_id, "rho", H5T_NATIVE_DOUBLE, dataspace_id,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      data_);
    status = H5Dclose(dataset_id);

    free(data_);


} // END createPattern




//! write potential, rho and so on into hdf5 file every some timesteps
void SmileiIO_Cart1D::write(  ElectroMagn* fields, SmileiMPI* smpi )
{
    SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);
    // Static-cast of the fields
    Field1D* Ex1D = static_cast<Field1D*>(fields->Ex_);


    Field1D* rho1D           = static_cast<Field1D*>(fields->rho_);
    Field1D* rho1D_global    = static_cast<Field1D*>(fields->rho_global);
    Field1D* phi1D_global    = static_cast<Field1D*>(fields->phi_global);
    Field1D* Ex1D_global    = static_cast<Field1D*>(fields->Ex_global);

    offset[0] = ndims_t;


    count[0]  = 1;
    count[1]  = dims_global[1];
    count[2]  = dims_global[2];
    count[3]  = dims_global[3];


    //>write the ===global potential=== in time ndims_t to the file
    dataset_id = H5Dopen2 (group_id, "potential", H5P_DEFAULT);

    memspace_id = H5Screate_simple (4, count, NULL);
    dataspace_id = H5Dget_space (dataset_id);
    status = H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset,
                                      stride, count, block);
    status = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, memspace_id,
                         dataspace_id, H5P_DEFAULT, phi1D_global->data_);

    status = H5Sclose (memspace_id);
    status = H5Sclose (dataspace_id);
    status = H5Dclose (dataset_id);


    //>write the ===global rho=== in time ndims_t to the file
    dataset_id = H5Dopen2 (group_id, "rho", H5P_DEFAULT);

    memspace_id = H5Screate_simple (4, count, NULL);
    dataspace_id = H5Dget_space (dataset_id);
    status = H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset,
                                      stride, count, block);

    status = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, memspace_id,
                         dataspace_id, H5P_DEFAULT, rho1D_global->data_);

    status = H5Sclose (memspace_id);
    status = H5Sclose (dataspace_id);
    status = H5Dclose (dataset_id);


    ndims_t++;

} // END write
