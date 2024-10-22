#ifndef SOLVERFACTORY_H
#define SOLVERFACTORY_H

#include "EF_Solver1D_TDMA.h"
#include "MF_Solver1D_Yee.h"
#include "MF_Solver2D_Yee.h"
#include "MF_Solver2D_Cowan.h"
#include "EF_Solver2D_SLU.h"

#include "PicParams.h"

#include "Tools.h"

class SolverFactory {
public:
    static Solver* create(PicParams& params, Grid* grid, SmileiMPI* smpi) {
        Solver* solver = NULL;
        if ( params.geometry == "1d3v" ) {
            solver = new EF_Solver1D_TDMA(params, smpi);
        }
        else if ( params.geometry == "2d3v" ) {
	    //if ()
            //solver = new MF_Solver2D_Yee(params);
            solver = new EF_Solver2D_SLU(params, grid, smpi);
	    //elseif()
	    //solver = new MF_Solver1D_Cowan(params);
        }
        else {}

        return solver;
    }

};


#endif
