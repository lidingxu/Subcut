/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_sub.cpp
 * @brief  submodular problem reader (to implement)
 * @author Liding XU
 * This file implements the reader/parser used to read the cbp input data. For more details see \ref NETWORKROUTING_READER.
 *
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <assert.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <utility>
#include <tuple>
//#include <armadillo>

#include "objscip/objscip.h"

#include "probdata.h"
#include "reader_sub.h"


using namespace scip;
using namespace std;

/** destructor of file reader to free user data (called when SCIP is exiting) */
SCIP_DECL_READERFREE(ReaderSubmodular::scip_free)
{
	return SCIP_OKAY;
} /*lint !e715*/


/** problem writing method of reader; NOTE: if the parameter "genericnames" is TRUE, then
 *  SCIP already set all variable and constraint names to generic names; therefore, this
 *  method should always use SCIPvarGetName() and SCIPconsGetName();
 *
 *  possible return values for *result:
 *  - SCIP_SUCCESS    : the reader read the file correctly and created an appropritate problem
 *  - SCIP_DIDNOTRUN  : the reader is not responsible for given input file
 *
 *  If the reader detected an error in the writing to the file stream, it should return
 *  with RETCODE SCIP_WRITEERROR.
 */
SCIP_DECL_READERWRITE(ReaderSubmodular::scip_write)
{
	*result = SCIP_DIDNOTRUN;

	return SCIP_OKAY;
} /*lint !e715*/

/**@} */

/** problem reading method of reader
 *
 *  possible return values for *result:
 *  - SCIP_SUCCESS    : the reader read the file correctly and created an appropritate problem
 *  - SCIP_DIDNOTRUN  : the reader is not responsible for given input file
 *
 *  If the reader detected an error in the input file, it should return with RETCODE SCIP_READERR or SCIP_NOFILE.
 */
SCIP_DECL_READERREAD(ReaderSubmodular::scip_read) {
	*result = SCIP_DIDNOTRUN;

   	SCIPdebugMessage("Start read!\n");
	// open the file
	ifstream filedata(filename);
	if (!filedata) {
		return SCIP_READERROR;
	}
	filedata.clear();

	// read parameters
	int numvars, numedges;
	filedata >> numvars >> numedges; 

	vector<tuple<int, int, SCIP_Real>> weights(numedges);
	// read graph weights
	for(int i = 0; i < numedges; i++){
		int u, v;
		SCIP_Real weight;
		filedata >> u >> v >> weight;
		weights[i] = make_tuple(u-1, v-1, weight);
	}

	SCIPdebugMessage("numvars:%d numedges:%d \n", numvars, numedges);
	// create the problem's data structure
	
	ProbData * probdata = NULL;
	probdata = new ProbData(numvars, numedges, weights);
	SCIP_CALL(SCIPgetBoolParam(scip,  "submodular/is_nature", &probdata->is_nature ));
	SCIP_CALL(SCIPgetBoolParam(scip,  "submodular/gradient_cut", &probdata->gradient_cut ));
	assert(probdata != NULL);
	SCIPdebugMessage("--problem data completed!\n");
	SCIP_CALL(SCIPcreateObjProb(scip, filename, probdata, FALSE));


	SCIP_CALL(probdata->createInitial(scip));
   
   	*result = SCIP_SUCCESS;

	SCIPdebugMessage("--reader read completed!\n");
	return SCIP_OKAY;
}


/**@} */
