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

/**@file   sepa_interlattice.h
 * @brief  lattice separator with intersection cuts
 * @author Liding Xu
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_INTERLATTICE_H__
#define __SCIP_SEPA_INTERLATTICE_H__

#define DEBUG_TIME_

#include "objscip/objscip.h"
#include <unordered_map>
#include <list>
#include <set>
#include <tuple>
#include <vector>
#include <algorithm>
#include <ctime>
#include "utility.h"


using namespace scip;
using namespace std;





/** C++ constraint handler for  constraints */
class SepaInterLattice : public ObjSepa
{
public:
	/** default constructor */
	SepaInterLattice(
		SCIP* scip /**< SCIP data structure */
	)
		: ObjSepa(scip, /**< SCIP data structure */
			"interlattice", /**< 	name of cut separator  */
			"stores the lattice constraints", /**< desc	description of cut separator */
			100000,	/**< priority	priority of the cut separator */
   			1,	/** freq	frequency for calling separator  */
			1.0,	/**< maxbounddist	maximal relative distance from current node's dual bound to primal bound compared to best node's dual bound for applying separation */
			FALSE,	/**< usessubscip	does the separator use a secondary SCIP instance? */
			TRUE	/**< delay	should separator be delayed, if other separators found cuts? */)
	{
	}

	SCIP_SepaData sepadata;

	virtual SCIP_DECL_SEPAEXECLP(scip_execlp);

	virtual SCIP_DECL_SEPAEXITSOL(scip_exitsol);

	virtual SCIP_DECL_SEPAINITSOL(scip_initsol);

	virtual SCIP_DECL_SEPAINIT(scip_init);

	virtual SCIP_DECL_SEPAEXIT(scip_exit);

	virtual SCIP_DECL_SEPACOPY(scip_copy);

	virtual SCIP_DECL_SEPAFREE(scip_free);
}; /*lint !e1712*/


#endif
