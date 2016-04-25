/*
 * SmileIO_Cart1D.h
 *
 *  Created on: 3 juil. 2013
 */

#ifndef SMILEIO_CART1D_H
#define SMILEIO_CART1D_H

#include <string>
#include <vector>

#include "SmileiIO.h"
#include "Tools.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiIO_Cart1D
//  --------------------------------------------------------------------------------------------------------------------
class SmileiIO_Cart1D : public SmileiIO {
public:
    //! Create // HDF5 environment
    SmileiIO_Cart1D( PicParams& params, SmileiMPI* smpi );
    //! Destructor for SmileiIO
    ~SmileiIO_Cart1D();

    //! Build memory and file space for // HDF5 write/read
    void createPattern( PicParams& params, SmileiMPI* smpi ){};

    //! Basic write field on its own file (debug)
    void write( ElectroMagn* fields, SmileiMPI* smpi ){};

private:


};

#endif /* SMILEIO_CART1D_H */
