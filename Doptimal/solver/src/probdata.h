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

/**@file   probdata.h
 * @brief  Problem data for D-optimal design problem
 * @author Liding XU
 * /

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROBDATA_SUB__
#define __SCIP_PROBDATA_SUB__

#include "objscip/objscip.h"
#include "scip/cons_linear.h"
#include <map>
#include <list>
#include <vector>
#include <utility>
#include <armadillo>

using namespace scip;
using namespace std;
using namespace arma;


/** SCIP user problem data for a Basic Submodular Maximization Problem */
class ProbData : public ObjProbData
{
public:
   /** default constructor */
   ProbData(
		const int numvars_,  /**< the number of items */
      const SCIP_Real dim_, /**< the problem dimension */
      const mat A_, /**<  A: dim_ * numvars_ */
      const int card_, /**< cardinality */
      const SCIP_Real epsilon_ /**<  epsilon: it is already sqrt, so the real epsilon in consideration is epsilon^2*/
   ): numvars(numvars_), dim(dim_), A(A_), card(card_), epsilon(epsilon_){
      E = epsilon * epsilon * eye(dim_, dim_);
      w_valinds = vector<pair<SCIP_Real, int>> (numvars_);
      w_indicator = vec(numvars_,fill::randu);
      w_val = vector<SCIP_Real> (numvars_);
      w_coeffs = vector<SCIP_Real> (numvars_);
      has_cardcons = card_ >= 0;
      has_knapcons = card < 0;
   };

   /**< destructor */
   ~ProbData();

   /**< utility functions*/
   SCIP_Real orcaleValue(
      const vec & values,
	   bool iszero = false
   );


   /**< evaluate the ray */
   SCIP_Real evalRay(
      SCIP * scip,
      vector<SCIP_Real> & ray,
      vector<SCIP_Real> & sol,
      const SCIP_Real & t
   );

   /**< evaluate the ray with gradient information*/
   SCIP_Real  evalRayWith(
      SCIP * scip,
      vector<SCIP_Real> & ray,
      vector<SCIP_Real> & sol,
      const SCIP_Real & t,
      SCIP_Real & gradient
   );

   /**< evaluate the underestimator of the submodular function */
   SCIP_Real evalUnderSub(
      SCIP * scip,                         /**< SCIP data structure */
      const vector<SCIP_Real> & varval    /**< variable value*/
   );

   /**< evaluate the underestimator of the submodular function with linear underestimator */
   SCIP_Real evalUnderSubWith(
      SCIP*              scip,    			/**< SCIP data structure */
      const vector<SCIP_Real> & varval    /**< variable value*/
   );

   /** add submodular constraint associated with S*/
   SCIP_RETCODE addSubCons(
      SCIP * scip, 
      const vec & indicator
   );

   /** create variable and initial constraints */
   SCIP_RETCODE createInitial(
	   SCIP*                 scip               /**< SCIP data structure */
   );

   /** release all */
   SCIP_RETCODE releaseAll(
	   SCIP*                 scip               /**< SCIP data structure */
   );


   /** destructor of user problem data to free original user data (called when original problem is freed)
    *
    *  If the "deleteobject" flag in the SCIPcreateObjProb() method was set to TRUE, this method is not needed,
    *  because all the work to delete the user problem data can be done in the destructor of the user problem
    *  data object. If the "deleteobject" flag was set to FALSE, and the user problem data object stays alive
    *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
    *  longer needed.
    */
   virtual SCIP_RETCODE scip_delorig(
      SCIP*              scip                /**< SCIP data structure */
      );
   
   /** destructor of user problem data to free transformed user data (called when transformed problem is freed)
    *
    *  If the "*deleteobject" flag in the scip_trans() method was set to TRUE, this method is not needed,
    *  because all the work to delete the user problem data can be done in the destructor of the user problem
    *  data object. If the "*deleteobject" flag was set to FALSE, and the user problem data object stays alive
    *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
    *  longer needed.
    */
   virtual SCIP_RETCODE scip_deltrans(
      SCIP*              scip                /**< SCIP data structure */
      );

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
   virtual SCIP_RETCODE scip_trans(
      SCIP*              scip,               /**< SCIP data structure */
      ObjProbData**      objprobdata,        /**< pointer to store the transformed problem data object */
      SCIP_Bool*         deleteobject        /**< pointer to store whether SCIP should delete the object after solving */
      );


   // problem relevant data
   int dim; // the dimension
   mat A, E; // the data matrix : dim * numvars, dim * dim
   vector<vector<SCIP_VAR*>> J; // lower-trigianle: dim * dim 
   vector<vector<SCIP_VAR*>> Z; // numvars * dim 
   vector<vector<SCIP_VAR*>> epsZ; // dim * dim 
   vector<vector<SCIP_VAR*>> epsZ2; // dim * dim 
   vector<vector<SCIP_VAR*>> t; // (numvars + 1)* dim 
   vector<SCIP_VAR *> w; // numvars
   SCIP_Real epsilon = 1e-3; // epsilon

   // probelem irrelevant data
   int numvars;  // the number of 0-1 variables
   SCIP_Bool has_cardcons; // has a cardinality constraint, assume to have the variable card
   SCIP_Real card;
   SCIP_Bool has_knapcons; // has a knapsack constraint, assume capacity = 1
   vector<SCIP_Real> knapweights;
   SCIP_Real fullvalue; 
   SCIP_Real emptyvalue;
   vector<SCIP_VAR*> bin_vars;
   SCIP_VAR* obj_var;
   vector<SCIP_CONS*> conss; // submodular constraints, model constraints.

   // settings
   SCIP_Bool is_nature;
   SCIP_Bool gradient_cut;


   // working vectors
   vector<pair<SCIP_Real, int>> w_valinds;
   vec w_indicator;
   vector<SCIP_Real> w_val;
   vector<SCIP_Real> w_coeffs;
};/*lint !e1712*/


#endif

