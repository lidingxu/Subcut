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

/**@file   probdata.cpp
 * @brief  Problem data for max cut problem
 * @author Liding XU
 * /
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <string.h>
#include <utility>  
#include <math.h> 
#include <tuple>

#include "probdata.h"
#include "objscip/objscip.h"
#include "scip/struct_cons.h"
#include "scip/cons_linear.h"
#include "scip/cons_nonlinear.h"
#include "scip/var.h"
#include "scip/scip.h"
#include <algorithm>

using namespace scip;
using namespace std;



/** ProbData destructor */
ProbData::~ProbData()
{

}


/** release scip reference in probelme data*/
SCIP_RETCODE ProbData::releaseAll(
	SCIP*              scip                /**< SCIP data structure */
) {
	// release
	for (int i = 0; i < bin_vars.size(); i++) {
		SCIP_CALL(SCIPreleaseVar(scip, &bin_vars[i]));
	}

	SCIP_CALL(SCIPreleaseVar(scip, &obj_var));

	for (int i = 0; i < conss.size(); i++){
		SCIP_CALL(SCIPreleaseCons(scip, &conss[i]));
	}

	return SCIP_OKAY;
}


/** destructor of user problem data to free original user data (called when original problem is freed)
 *
 *  If the "deleteobject" flag in the SCIPcreateObjProb() method was set to TRUE, this method is not needed,
 *  because all the work to delete the user problem data can be done in the destructor of the user problem
 *  data object. If the "deleteobject" flag was set to FALSE, and the user problem data object stays alive
 *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
 *  longer needed.
 */
SCIP_RETCODE ProbData::scip_delorig(
   SCIP*              scip                /**< SCIP data structure */
   )
{
	SCIP_CALL(releaseAll(scip));
	return SCIP_OKAY;
}


/**< evaluate the ray */
SCIP_Real ProbData::evalRay(
	SCIP * scip,
	const vector<SCIP_Real> & ray,
	const vector<SCIP_Real> & sol,
	const SCIP_Real t
){
	for(int i = 0; i < numvars; i++){
		w_val[i] = sol[i] + t * ray[i];
	}
	SCIP_Real underobj = evalUnderSub(scip, w_val);
	SCIP_Real obj =  sol[numvars] + t * ray[numvars];
	return obj - underobj;
}

/**< evaluate the ray with the gradient information*/
SCIP_Real ProbData::evalRayWith(
	SCIP * scip,
	const vector<SCIP_Real> & ray,
	const vector<SCIP_Real> & sol,
	const SCIP_Real  t,
	SCIP_Real & gradient_
){
	for(int i = 0; i < numvars; i++){
		w_val[i] = sol[i] + t * ray[i];
	}
	gradient = ray[numvars];
	SCIP_Real underobj = evalUnderSubWith(scip, w_val, ray);
	SCIP_Real obj =  sol[numvars] + t * ray[numvars];
	gradient_ = gradient;
	return obj - underobj;	
}

void ProbData::useBranch(
	SCIP * scip
){
	startvalue = emptyvalue;
	is_branch_opt = true;
	for(int i = 0; i < numvars; i++){
		if(SCIPvarGetUbLocal(bin_vars[i]) < 0.5){
			branch_indicator[i] = 0;
		}
		else if(SCIPvarGetLbLocal(bin_vars[i]) > 0.5){
			branch_indicator[i] = 1;
			SCIP_Real delta_val = 0;
			auto & incident_edges = incident_list[i];
			for(auto & p: incident_edges){
				SCIP_Real wt = p.second;
				delta_val += branch_indicator[p.first] == 1 ? -wt : wt;
			}
			startvalue += delta_val;
		}
		else{
			branch_indicator[i] = 2;
		}
	}
}

// evaluate the underestimator of submodular function
SCIP_Real ProbData::evalUnderSub(
   SCIP*              scip,    /**< SCIP data structure */
   const vector<SCIP_Real> & varval  /**< vector value */
){
	if(is_branch_opt){
		int nactive = 0;
		for(int i = 0; i < numvars; i++){
			if(branch_indicator[i] == 2){
				w_valinds[i].first  = varval[i]; 
				w_valinds[i].second = i;
				nactive++;
			}
		}
		sort(w_valinds.begin(), w_valinds.begin() + nactive, [](const pair<SCIP_Real, int> & p1, const  pair<SCIP_Real, int> & p2){ return p1.first > p2.first || (p1.first == p2.first && p1.second > p2.second);});
		
		std::fill(w_indicator.begin(), w_indicator.begin() + nactive, false);
		SCIP_Real underobj = startvalue;
		for(int i = 0; i < nactive; i++){
			int id  = w_valinds[i].second;
			w_indicator[id] = true;
			SCIP_Real delta_val = 0;
			auto & incident_edges = incident_list[id];
			for(auto & p: incident_edges){
				SCIP_Real wt = p.second;
				delta_val += (branch_indicator[p.first] == 1 || (branch_indicator[p.first] == 2 && w_indicator[p.first]))  ? -wt : wt;
			}
			underobj += delta_val * w_valinds[i].first;
		}
		return underobj;
	}{
		for(int i = 0; i < numvars; i++){
			w_valinds[i].first  = varval[i]; 
			w_valinds[i].second = i;
		}
		sort(w_valinds.begin(), w_valinds.end(), [](const pair<SCIP_Real, int> & p1, const  pair<SCIP_Real, int> & p2){ return p1.first > p2.first || (p1.first == p2.first && p1.second > p2.second);});
		
		std::fill(w_indicator.begin(), w_indicator.end(), false);
		SCIP_Real underobj = emptyvalue;
		for(int i = 0; i < numvars; i++){
			int id  = w_valinds[i].second;
			w_indicator[id] = true;
			SCIP_Real delta_val = 0;
			auto & incident_edges = incident_list[id];
			for(auto & p: incident_edges){
				SCIP_Real wt = p.second;
				delta_val += w_indicator[p.first] ? -wt : wt;
			}
			underobj += delta_val * w_valinds[i].first;
		}
		return underobj;
	}
}

// evaluate the underestimator of submodular function with the linear underestimator
SCIP_Real ProbData::evalUnderSubWith(
   	SCIP*              		scip,    	/**< SCIP data structure */
   	const vector<SCIP_Real> & varval,    /**< variable value*/
		const vector<SCIP_Real> & ray
){
	if(is_branch_opt){
		for(int i = 0; i < numvars; i++){
			w_valinds[i].first  = varval[i]; 
			w_valinds[i].second = i;
		}
		sort(w_valinds.begin(), w_valinds.end(), [](const pair<SCIP_Real, int> & p1, const  pair<SCIP_Real, int> & p2){ return p1.first > p2.first || (p1.first == p2.first && p1.second > p2.second);});
		
		std::fill(w_indicator.begin(), w_indicator.end(), false);
		SCIP_Real underobj = emptyvalue;
		for(int i = 0; i < numvars; i++){
			int id  = w_valinds[i].second;
			w_indicator[id] = true;
			SCIP_Real delta_val = 0;
			auto & incident_edges = incident_list[id];
			for(auto & p: incident_edges){
				SCIP_Real wt = p.second;
				delta_val += w_indicator[p.first] ? -wt : wt;
			}
			gradient -= delta_val * ray[id];
			underobj += delta_val * w_valinds[i].first;
		}
		return underobj;
	}
	else{
		int nactive = 0;
		for(int i = 0; i < numvars; i++){
			if(branch_indicator[i] == 2){
				w_valinds[i].first  = varval[i]; 
				w_valinds[i].second = i;
				nactive++;
			}
		}
		sort(w_valinds.begin(), w_valinds.begin() + nactive, [](const pair<SCIP_Real, int> & p1, const  pair<SCIP_Real, int> & p2){ return p1.first > p2.first || (p1.first == p2.first && p1.second > p2.second);});
		
		std::fill(w_indicator.begin(), w_indicator.begin() + nactive, false);
		SCIP_Real underobj = startvalue;
		for(int i = 0; i < nactive; i++){
			int id  = w_valinds[i].second;
			w_indicator[id] = true;
			SCIP_Real delta_val = 0;
			auto & incident_edges = incident_list[id];
			for(auto & p: incident_edges){
				SCIP_Real wt = p.second;
				delta_val += (branch_indicator[p.first] == 1 || (branch_indicator[p.first] == 2 && w_indicator[p.first])) ? -wt : wt;
			}
			gradient -= delta_val * ray[id];
			underobj += delta_val * w_valinds[i].first;
		}
		return underobj;
	}
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed)
 *
 *  The user has two possibilities to implement this method:
 *   1. Return the pointer to the original problem data object (this) as pointer to the transformed problem data
 *      object. The user may modify some internal attributes, but he has to make sure, that these modifications are
 *      reversed in the scip_deltrans() method, such that the original problem data is restored. In this case,
 *      he should set *deleteobject to FALSE, because the problem data must not be destructed by SCIP after the
 *      solving process is terminated.
 *   2. Call the copy constructor of the problem data object and return the created copy as transformed problem
 *      data object. In this case, he probably wants to set *deleteobject to TRUE, thus letting SCIP call the
 *      destructor of the object if the transformed problem data is no longer needed.
 */
SCIP_RETCODE ProbData::scip_trans(
   SCIP*              scip,               /**< SCIP data structure */
   ObjProbData**      objprobdata,        /**< pointer to store the transformed problem data object */
   SCIP_Bool*         deleteobject        /**< pointer to store whether SCIP should delete the object after solving */
   )
{  /*lint --e{715}*/
   	assert( objprobdata != NULL );
   	assert( deleteobject != NULL );


   	// create and cpature transformed path varibles 
   	SCIPdebugMessage("start transform !!!!!!!!!!!\n");

    // allocate memory for target prob data
	ProbData * transprobdata = new ProbData(numvars, numedges, weights);
	transprobdata->fullvalue = fullvalue;
	transprobdata->emptyvalue = emptyvalue;
	transprobdata->card = card;
	transprobdata->knapweights = knapweights;
	transprobdata->is_nature  = is_nature;

	for (int i = 0; i < bin_vars.size(); i++) {
		SCIP_VAR * var;
		SCIP_CALL(SCIPtransformVar(scip, bin_vars[i], &var));
		transprobdata->bin_vars.push_back(var);
	}

	SCIP_VAR * var;
	SCIP_CALL(SCIPtransformVar(scip, obj_var, &var));
	transprobdata->obj_var = var;

	for (int i = 0; i < conss.size(); i++){
		SCIP_CONS * cons;
		SCIP_CALL(SCIPtransformCons(scip, conss[i], &cons ));
		transprobdata->conss.push_back(cons);
	}

	
	SCIPdebugMessage("transformed data check!");
  	// transform and capture transformed set partition constraints

   SCIPdebugMessage("end transform \n");
   assert( transprobdata != NULL );
   *objprobdata = transprobdata;           
   
   *deleteobject = FALSE;

   return SCIP_OKAY;
}      

/** destructor of user problem data to free original user data (called when original problem is freed)
 *
 *  If the "deleteobject" flag in the SCIPcreateObjProb() method was set to TRUE, this method is not needed,
 *  because all the work to delete the user problem data can be done in the destructor of the user problem
 *  data object. If the "deleteobject" flag was set to FALSE, and the user problem data object stays alive
 *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
 *  longer needed.
 */
SCIP_RETCODE ProbData::scip_deltrans(
   SCIP*              scip                /**< SCIP data structure */
   )
{ 
   SCIP_CALL(releaseAll(scip));
   return SCIP_OKAY;
}





/** create variables and initial  constraints */
SCIP_RETCODE ProbData::createInitial(
	SCIP*                 scip               /**< SCIP data structure */
) {   
	// add binary variables
	for(int i = 0; i < numvars; i++){
		SCIP_VAR* bin_var;
		string tmp = "b"+ std::to_string(i);
		SCIP_CALL(SCIPcreateVar(
			scip, /**<	SCIP data structure*/
			&bin_var, /**< 	pointer to variable object*/
			tmp.c_str(), /**< name of variable, or NULL for automatic name creation*/
			0.0, /**<	lower bound of variable*/
			1.0, /**< 	upper bound of variable */
			0, /**<	objective function value */
			SCIP_VARTYPE_BINARY, /**< type of variable */
			TRUE, /**<	should var's column be present in the initial root LP?*/
			FALSE, /**<	is var's column removable from the LP (due to aging or cleanup)?*/
			NULL, NULL, NULL, NULL, NULL
		));
		SCIP_CALL(SCIPaddVar(scip, bin_var));
		SCIP_CALL(SCIPcaptureVar(scip, bin_var));
		bin_vars.push_back(bin_var);
		SCIP_CALL(SCIPreleaseVar(scip, &bin_var));
	}
	// add objective variable (in minimization form, max to min - )
	SCIP_VAR * obj_var_;
	SCIP_CALL(SCIPcreateVar(
		scip, /**<	SCIP data structure*/
		&obj_var_, /**< 	pointer to variable object*/
		"obj_var", /**< name of variable, or NULL for automatic name creation*/
		-SCIPinfinity(scip), /**<	lower bound of variable*/
		SCIPinfinity(scip), /**< 	upper bound of variable */
		-1, /**<	objective function value */
		SCIP_VARTYPE_CONTINUOUS, /**< type of variable */
		TRUE, /**<	should var's column be present in the initial root LP?*/
		FALSE, /**<	is var's column removable from the LP (due to aging or cleanup)?*/
		NULL, NULL, NULL, NULL, NULL
	));
	SCIP_CALL(SCIPaddVar(scip, obj_var_));
	obj_var = obj_var_;
	SCIP_CALL(SCIPcaptureVar(scip, obj_var_));	
	SCIP_CALL(SCIPreleaseVar(scip, &obj_var_));


	// add initial submodular constraints
	emptyvalue = 0;

	//SCIPdebugMessage("variable created\n");
	SCIP_Cons * cons = NULL;
	if(!is_nature){
		// S = e_j and e_{-j}
		/*
		for(int j = 0; j < numvars; j++){
			std::fill(v.begin(), v.end(), 0);
			v[j] = 1;
			//SCIPdebugMessage("con 1 to created\n");
			//cout<< v << endl;
			SCIP_CALL(addSubCons(scip, v));
			//SCIPdebugMessage("con 1 created\n");

			std::fill(v.begin(), v.end(), 1);
			v[j] = 0;
			SCIP_CALL(addSubCons(scip, v));
		}*/
	}
	else{
		int nlinvars = numvars + 1;
		vector<SCIP_VAR *> linvars(bin_vars);
		linvars.push_back(obj_var);
		vector<SCIP_Real> lincoefs(nlinvars, 0);
		vector<SCIP_VAR *> quadvars1(numedges, NULL);
		vector<SCIP_VAR *> quadvars2(numedges, NULL);
		vector<SCIP_Real> quadcoeffs(numedges, 0);
		int quad_id = 0;
		for(int i  = 0; i < numedges; i++){
			int u = get<0>(weights[i]);
			int v = get<1>(weights[i]);
			assert(u < numvars && v < numvars);
			SCIP_Real weight = get<2>(weights[i]);
			lincoefs[u] += weight;
			lincoefs[v] += weight;
			quadvars1[quad_id] = linvars[u];
			quadvars2[quad_id] = linvars[v];
			assert(quadvars1[quad_id] != NULL && quadvars2[quad_id] != NULL);
			quadcoeffs[quad_id] = -2*weight;
			quad_id++;
		}
		lincoefs[numvars] = -1;
		//SCIPdebugMessage("1\n");
		// 0 <= lin + quad - t <= inf 
		SCIP_CONS * quad_cons = NULL;
		SCIP_CALL(SCIPcreateConsQuadraticNonlinear( 
			scip, 
			&quad_cons,
			"objcons",
			nlinvars,
			linvars.data(),
			lincoefs.data(),
			numedges,
			quadvars1.data(),
			quadvars2.data(),
			quadcoeffs.data(),
			0,
			SCIPinfinity(scip),
			TRUE,               /**< should the LP relaxation of constraint be in the initial LP?
														*   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
			TRUE,                /**< should the constraint be separated during LP processing?
														*   Usually set to TRUE. */
			TRUE,               /**< should the constraint be enforced during node processing?
														*   TRUE for model constraints, FALSE for additional, redundant constraints. */
			TRUE,               /**< should the constraint be checked for feasibility?
														*   TRUE for model constraints, FALSE for additional, redundant constraints. */
			TRUE,               /**< should the constraint be propagated during node processing?
														*   Usually set to TRUE. */
			FALSE,              /**< is constraint only valid locally?
														*   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
			FALSE,              /**< is constraint subject to aging?
														*   Usually set to FALSE. Set to TRUE for own cuts which
														*   are separated as constraints. */
			FALSE,              /**< should the relaxation be removed from the LP due to aging or cleanup?
														*   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
			FALSE               /**< should the constraint always be kept at the node where it was added, even
														*   if it may be moved to a more global node?
														*   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
		)); 	
		SCIP_CALL(SCIPaddCons(scip, quad_cons));
		SCIP_CALL(SCIPcaptureCons(scip, quad_cons));
		conss.push_back(quad_cons);
		SCIP_CALL(SCIPreleaseCons(scip, &quad_cons));
	}	
	return SCIP_OKAY;
}

/**@} */
