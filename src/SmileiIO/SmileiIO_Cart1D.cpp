/*
 * SmileiIO_Cart1D.cpp
 *
 *  Created on: 3 juil. 2013
 */

#include "SmileiIO_Cart1D.h"

#include "PicParams.h"
#include "SmileiMPI_Cart1D.h"
#include "Field1D.h"

using namespace std;

SmileiIO_Cart1D::SmileiIO_Cart1D( PicParams& params, SmileiMPI* smpi )
    : SmileiIO( params, smpi )
{
}

SmileiIO_Cart1D::~SmileiIO_Cart1D()
{
}
