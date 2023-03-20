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
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa_intersub.cpp
 * @brief  submodular mamization separator with intersection cuts
 * @author Liding Xu
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <assert.h>
#include <string.h>

#include "sepa_intersub.h"




/** constructs map between lp position of a basic variable and its row in the tableau,
 * A x + I s = 0, and A_b x_b + I_b s_b + A_n x_n + I_n s_n = 0.
*/
static
SCIP_RETCODE constructBasicVars2TableauRowMap(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  map                 /**< buffer to store the map */
   )
{
   int* basisind;
   int nrows;
   int i;

   nrows = SCIPgetNLPRows(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &basisind, nrows) );

   SCIP_CALL( SCIPgetLPBasisInd(scip, basisind) );
   for( i = 0; i < nrows; ++i )
   {
      if( basisind[i] >= 0 )
         map[basisind[i]] = i;
   }

   SCIPfreeBufferArray(scip, &basisind);

   return SCIP_OKAY;
}


/** get the tableau rows of the variables in vars,
 * A x + I s = 0, and A_b x_b + I_b s_b + A_n x_n + I_n s_n = 0. B = [A_b, I_b]
 * B^{-1} (Ax + Is) = 0,  x_b +  B^{-1}A_n x_n + s_b + B^{-1} I_n s_n = 0.
 */
static
SCIP_RETCODE getTableauRows(
   SCIP*       scip,               /**< SCIP data structure */
   const vector<SCIP_VAR *>  &    vars,       /**< variables */
   int * basicvarpos2tableaurow,  /**< map between basic var and its tableau row */
   vector< SCIP_Real*> & tableau,     /**< map between all var and its tableau row */
   SCIP_Bool*     success             /**< set to TRUE if no variable had basisstat = ZERO */
   )
{

   *success = TRUE;

   int nrows = SCIPgetNLPRows(scip);
   int ncols = SCIPgetNLPCols(scip);

   /* check if we have the tableau row of the variable and if not compute it */
   for(int i = 0; i < vars.size(); i++){
      SCIP_VAR * var = vars[i];
      if( tableau[i] == NULL)
      {
         /* get column of variable */
         SCIP_COL* col = SCIPvarGetCol(var);
         /* if variable is basic, then get its tableau row and insert it in the hashmap */
         if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_BASIC )
         {
            int lppos = SCIPcolGetLPPos(col);
            SCIP_Real* densetableaurow;
            SCIP_CALL( SCIPallocBufferArray(scip, &densetableaurow, ncols + nrows));
            SCIP_CALL( SCIPgetLPBInvRow(scip, basicvarpos2tableaurow[lppos], &densetableaurow[ncols], NULL, NULL) );
            SCIP_CALL( SCIPgetLPBInvARow(scip, basicvarpos2tableaurow[lppos], &densetableaurow[ncols], densetableaurow, NULL, NULL) );

            /* insert tableau row in hashmap*/
            tableau[i] = densetableaurow;
         }
         else if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_ZERO)
         {
            *success = FALSE;
            return SCIP_OKAY; /* don't even bother */
         }
         else
         {
            tableau[i] = (SCIP_Real *) NULL;
         }
      }
   }
   return SCIP_OKAY;
}

/** adds cutcoef * (col - col*) to rowprep */
static
SCIP_RETCODE addColToCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to store intersection cut */
   SCIP_Real             cutcoef,            /**< cut coefficient */
   SCIP_COL*             col                 /**< column to add to rowprep */
   )
{
   assert(col != NULL);

#ifdef DEBUG_INTERS_NUMERICS
   SCIPinfoMessage(scip, NULL, "adding col %s to cut. %g <= col <= %g\n", SCIPvarGetName(SCIPcolGetVar(col)),
      SCIPvarGetLbLocal(SCIPcolGetVar(col)), SCIPvarGetUbLocal(SCIPcolGetVar(col)));
   SCIPinfoMessage(scip, NULL, "col is active at %s. Value %.15f\n", SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_LOWER ? "lower bound" :
      "upper bound" , SCIPcolGetPrimsol(col));
#endif

   SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, SCIPcolGetVar(col), cutcoef) );
   SCIProwprepAddConstant(rowprep, -cutcoef * SCIPcolGetPrimsol(col) );

   return SCIP_OKAY;
}

/** The maximal S-free set is gamma(z) <= 0; we find the intersection point of the ray `ray` starting from zlp with the
 * boundary of the S-free set.
 * That is, we find t >= 0 such that gamma(zlp + t * ray) = 0 and zlp+t*ray >= 0.
 */
static
SCIP_Real computeIntersectionPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SepaData * sepadata,
   bool iszero,                 // is all binaries in ray zero?
   SCIP_Bool*            success             /**< pointer to store whether the generation of cutcoefs was successful */
   )
{
   // to do binary search;
   //SCIPdebugMessage("bin search\n");
   auto & mlf = sepadata->mlf;
   //SCIPdebugMessage("%f %f\n", (mlf.val), mlf.ray[mlf.ct_ind]);
   if(iszero){
      SCIP_Real t = (mlf.val /  mlf.ray[mlf.ct_ind]) * (1 - 1e-9);
      if(t < 0){
         t = SCIPinfinity(scip);
         //SCIPdebugMessage("%f %f %f\n",  t, mlf.val, mlf.linearray);
      }
      assert( t >= 0);
      return t; 
   }


   SCIP_Real lb = 0.0;
   SCIP_Real ub = SCIPinfinity(scip) / 2;
   SCIP_Real gradient = 0;
   SCIP_Real gamma_t;
   SCIP_Real t = 0.2;

   // speed up for the case that ray has no intersection with S-free set
   gamma_t = mlf.evalRayWith(ub, gradient, false); 
   if(gamma_t > 0){
      t = SCIPinfinity(scip);
      return t;
   }


   for(int i = 0; i < sepadata->BINSEARCH_MAXITERS; ++i )
   {
      gamma_t = mlf.evalRayWith(t, gradient, false); 
      //SCIPdebugMessage("bin search_ %d %f %f %f\n", i, t, gamma_t, gradient);

      if( SCIPisFeasZero(scip, gamma_t))
      {  
         *success = TRUE;
         //SCIPdebugMessage("%f %d %f %f\n", gradient, i, t, gamma_t);
         break;
      }
      else if(SCIPisNegative(scip, gradient)){
         //SCIPdebugMessage("%f %f %f %f\n", gamma_t, ( mlf.val), gradient, mlf.ray[cont_ind]);
         t -= gamma_t / gradient;
      }
      else{
         t *= 2;
         if( t >= ub -1 ){
            t = SCIPinfinity(scip);
            //assert(false);
            break;
         }
      }


   }
   //assert(SCIPisFeasZero(scip, gamma_t));
   //abort();
   //SCIPdebugMessage("%f\n", t);
   assert(t >= 0);
   return t * (1 - 1e-9);
}



/** The maximal S-free set is gamma(z) <= 0; we find the intersection point of the ray `ray` starting from zlp with the
 * boundary of the S-free set.
 * That is, we find t >= 0 such that gamma(zlp + t * ray) = 0 and zlp+t*ray >= 0.
 */
static
SCIP_Real computeIntersectionPointf(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SepaData * sepadata,
   bool iszero,                              // is all binaries in ray zero?
   SCIP_Bool*            success             /**< pointer to store whether the generation of cutcoefs was successful */
   )
{
   auto & mlf = sepadata->mlf;
   if(iszero){
      SCIP_Real t = (mlf.val / mlf.linearray) * (1 - 1e-9); // (0 - mlf.val) / (- mlf.linearray)
      if(t < 0){
         t = SCIPinfinity(scip);
         //SCIPdebugMessage("%f %f %f\n",  t, mlf.val, mlf.linearray);
      }
      assert(t >= 0);
      return t; 
   }
   SCIP_Real lb = 0.0;
   SCIP_Real ub = SCIPinfinity(scip) / 2;
   SCIP_Real gradient = 0;
   SCIP_Real gamma_t;
   SCIP_Real t = 0.2;

   // speed up for the case that ray has no intersection with S-free set
   gamma_t = mlf.evalRayfWith(ub, gradient); 
   if(gamma_t < 0){
      t = SCIPinfinity(scip);
      return t;
   }

   for(int i = 0; i < sepadata->BINSEARCH_MAXITERS; ++i )
   {
      gamma_t = mlf.evalRayfWith(t, gradient); 

      //SCIPdebugMessage("bin search_ %d %f %f %f\n", i, t, gamma_t, gradient);

      if( SCIPisFeasZero(scip, gamma_t))
      {  
         //SCIPdebugMessage("%f %f\n", gamma_t, t);
         *success = TRUE;
         break;
      }
      else if(SCIPisPositive(scip, gradient)){
         t -= gamma_t / gradient;
      }
      else{
         t *= 2;
         if( t >= ub -1 ){
            t = SCIPinfinity(scip);
            //assert(false);
            break;
         }
      }
   }
   //SCIPdebugMessage("%f %f\n",  t, gamma_t);
   return t * (1 - 1e-9);
}

/** computes the cut coefs of the  non-basic (non-slack) variables (correspond to cols) and adds them to the
 * A x + I s = 0, and A_b x_b + I_b s_b + A_n x_n + I_n s_n = 0. B = [A_b, I_b]
 * B^{-1} (Ax + Is) = 0,  x_b +  B^{-1}A_n x_n + s_b + B^{-1} I_n s_n = 0.
 */
static
SCIP_RETCODE addCols(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA * sepadata,
   vector<SCIP_Real*> & tableau, /**< tableau rows corresponding to the variables in vars */
   SCIP_ROWPREP*         rowprep,            /**< store cut */
   SCIP_Bool*            success             /**< pointer to store whether the generation of cutcoefs was successful */
   )
{
   //SCIPdebugMessage("add col %d\n", int(tableaurows.size()));
   *success = TRUE;

   int nrays = 0;

   auto & mlf = sepadata->mlf;
   /* loop over non-basic (non-slack) variables */
   SCIP_COL** cols = SCIPgetLPCols(scip);
   int ncols = SCIPgetNLPCols(scip);
   SCIPdebugMessage("ncols: %d\n", ncols);
   for(int i = 0; i < ncols; ++i )
   {

      SCIP_COL* col = cols[i];

      SCIP_Real factor;
      /* set factor to store entries of ray as = [-BinvL, BinvU] */
      if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_LOWER )
         factor = -1.0;
      else if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_UPPER )
         factor = 1.0;
      else if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_ZERO )
      {
         *success = FALSE;
         return SCIP_OKAY;
      }
      else
         continue;

      SCIP_Real interpoint = 0;
      //SCIPdebugMessage("step %d %d\n", i, ncols);
      /* build the ray and compute the intersection point */
      if(sepadata->disaggregate){
         bool israyzero = true;
         bool iszero = true;
         auto & ray = mlf.ray;
         auto & extendedfvars = *mlf.extendedfvars;
         int  varsize = mlf.size;
         int ftype =  mlf.typef ? 1 : 2;
         mlf.linearray = 0;
         //unordered_map<SCIP_VAR*,SCIP_Real> raysbycol;
         for(int j = 0; j < extendedfvars.size(); j++)
         {
            SCIP_VAR * var = extendedfvars[j];
            SCIP_Real ray_val = 0;
            //SCIPdebugMessage("add col 2.1 %d %d\n",  int(tableaurows.find(var) != tableaurows.end()), int(tableaurows[var] != NULL));
            if( tableau[j] != NULL )
            {
               //SCIPdebugMessage("%f %f add col 2.11\n", factor,  tableau[j][i]);
               ray_val = factor * ( SCIPisZero(scip, tableau[j][i]) ? 0.0 : tableau[j][i] );
               //SCIPdebugMessage("add col 2.12\n");
            }
            else
            {
               //SCIPdebugMessage("add col 2.2\n");
               if( col == SCIPvarGetCol(var) )
                  ray_val= (-factor);
               else
                  ray_val =(0.0);
            }
            if(j < varsize){
               ray[j] = ray_val;
               bool affected  = (mlf.probvartypes[j] == ftype || mlf.probvartypes[j] == 3 );
               bool rayjzero =  (  affected  ? SCIPisZero(scip, ray[j]) : true);
               israyzero = israyzero && rayjzero;
               iszero = iszero && rayjzero;
            }
            else{
               //SCIPdebugMessage("%f %f\n", ray_val, mlf.mlf_term[j - varsize].second);
               mlf.linearray += ray_val * (*mlf.mlf_term)[j - varsize].second;
               //israyzero = israyzero && SCIPisZero(scip, ray_val);
               //SCIPdebugMessage("%f %d\n", ray_val, SCIPisZero(scip, ray_val));
            }
            //SCIPdebugMessage("%f\n", ray_val);
         }
         israyzero = israyzero && SCIPisZero(scip, mlf.linearray);
         //SCIPdebugMessage("add col 3 %d %d\n", int(israyzero), int(israyzero));
         /* do nothing if ray is 0 */
         if( israyzero )
            continue;


         interpoint = computeIntersectionPointf(scip, sepadata, iszero, success);
      }
      else{
         bool israyzero = true;
         bool iszero = true;
         auto & ray = mlf.ray;
         auto & probvars = mlf.probvars;

         #ifdef DEBUG_TIME
         time(&t5);
         #endif
         for(int j = 0; j < probvars.size(); j++)
         {
            SCIP_VAR * var = probvars[j];
            //SCIPdebugMessage("add col 2.1 %d %d\n",  int(tableaurows.find(var) != tableaurows.end()), int(tableaurows[var] != NULL));
            if( tableau[j] != NULL )
            {
               //SCIPdebugMessage("add col 2.11\n");
               ray[j] = factor * ( SCIPisZero(scip, tableau[j][i]) ? 0.0 : tableau[j][i] );
               //SCIPdebugMessage("add col 2.12\n");
            }
            else
            {
               //SCIPdebugMessage("add col 2.2\n");
               if( col == SCIPvarGetCol(var) )
                  ray[j]= (-factor);
               else
                  ray[j]=(0.0);
            }
            //SCIPdebugMessage("%f\n", raysbycol[var]);
            israyzero = israyzero &&  SCIPisZero(scip, ray[j]);
            iszero = iszero && ( (mlf.probvartypes[j] == 2 ||  mlf.probvartypes[j] == 3)  ? SCIPisZero(scip, ray[j]) : true);
         }
         #ifdef DEBUG_TIME
         time(&t6);
         memory_time += difftime(t6,t5);
         #endif
         
         /* do nothing if ray is 0 */
         if(israyzero)
            continue;

         mlf.linearf1ray  = 0;
         for(int ind: mlf.set1)
         {
				mlf.linearf1ray += ray[ind] * mlf.linearize1[ind];
         }         
         mlf.linearray = 0;
         for(auto p: mlf.linear_term)
         {
				mlf.linearray += ray[p.first] * p.second;
         }

         //nonzeros++;
         //nonrayzeros+= !israyzero;
         #ifdef DEBUG_TIME
         time(&t5);
         #endif
         interpoint = computeIntersectionPoint(scip, sepadata, iszero, success);
         #ifdef DEBUG_TIME
         time(&t6);
         compute_time += difftime(t6,t5);
         #endif
         //SCIPdebugMessage("%d/%d %d/%d %f\n", iszero, nonzeros, israyzero, nonrayzeros, interpoint);
      }


      if(SCIPisZero(scip, interpoint)){
         *success = FALSE;
         break;
      }

      //SCIPdebugMessage("add %d/%d4\n", i, ncols);

      if( *success == FALSE )
         return SCIP_OKAY;

      /* count nonzero rays */
      nrays += 1;

      /* compute cut coef */
      SCIP_Real cutcoef = SCIPisInfinity(scip, interpoint) ? 0.0 : 1.0 / interpoint;

      /* add var to cut: if variable is nonbasic at upper we have to flip sign of cutcoef */
      assert(SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_UPPER || SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_LOWER);
      SCIP_CALL( addColToCut(scip, rowprep, SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_UPPER ? -cutcoef :
            cutcoef, col) );
   }
   //SCIPdebugMessage("add col__ nray:%d\n", nrays);
   //abort();
   return SCIP_OKAY;
}



/** adds cutcoef * (slack - slack*) to rowprep
  *
  * row is lhs <= <coefs, vars> + constant <= rhs, thus slack is defined by
  * slack + <coefs, vars> + constant = side
  * If row (slack) is at upper, it means that <coefs,vars*> + constant = rhs, and so
  * slack* = side - rhs --> slack - slack* = rhs - <coefs, vars> - constant.
  * If row (slack) is at lower, then <coefs,vars*> + constant = lhs, and so
  * slack* = side - lhs --> slack - slack* = lhs - <coefs, vars> - constant.
  */
static
SCIP_RETCODE addRowToCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to store intersection cut */
   SCIP_Real             cutcoef,            /**< cut coefficient */
   SCIP_ROW*             row,                /**< row, whose slack we are adding to rowprep */
   SCIP_Bool*            success             /**< buffer to store whether the row is nonbasic enough */
   )
{
   int i;
   SCIP_COL** rowcols;
   SCIP_Real* rowcoefs;
   int nnonz;

   assert(row != NULL);

   rowcols = SCIProwGetCols(row);
   rowcoefs = SCIProwGetVals(row);
   nnonz = SCIProwGetNLPNonz(row);

#ifdef DEBUG_INTERS_NUMERICS
   SCIPinfoMessage(scip, NULL, "adding slack var row_%d to cut. %g <= row <= %g\n", SCIProwGetLPPos(row), SCIProwGetLhs(row), SCIProwGetRhs(row));
   SCIPinfoMessage(scip, NULL, "row is active at %s = %.15f Activity %.15f\n", SCIProwGetBasisStatus(row) == SCIP_BASESTAT_LOWER ? "lhs" :
   "rhs" , SCIProwGetBasisStatus(row) == SCIP_BASESTAT_LOWER ? SCIProwGetLhs(row) : SCIProwGetRhs(row),
   SCIPgetRowActivity(scip, row));
#endif

   if( SCIProwGetBasisStatus(row) == SCIP_BASESTAT_LOWER )
   {
      assert(!SCIPisInfinity(scip, -SCIProwGetLhs(row)));
      if( ! SCIPisFeasEQ(scip, SCIProwGetLhs(row), SCIPgetRowActivity(scip, row)) )
      {
         *success = FALSE;
         return SCIP_OKAY;
      }

      SCIProwprepAddConstant(rowprep, SCIProwGetLhs(row) * cutcoef);
   }
   else
   {
      assert(!SCIPisInfinity(scip, SCIProwGetRhs(row)));
      if( ! SCIPisFeasEQ(scip, SCIProwGetRhs(row), SCIPgetRowActivity(scip, row)) )
      {
         *success = FALSE;
         return SCIP_OKAY;
      }

      SCIProwprepAddConstant(rowprep, SCIProwGetRhs(row) * cutcoef);
   }

   for( i = 0; i < nnonz; i++ )
   {
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, SCIPcolGetVar(rowcols[i]), -rowcoefs[i] * cutcoef) );
   }

   SCIProwprepAddConstant(rowprep, -SCIProwGetConstant(row) * cutcoef);

   return SCIP_OKAY;
}

/** computes the cut coefs of the non-basic slack variables (correspond to rows) and adds them to the
 * A x + I s = 0, and A_b x_b + I_b s_b + A_n x_n + I_n s_n = 0. B = [A_b, I_b]
 * B^{-1} (Ax + Is) = 0,  x_b +  B^{-1}A_n x_n + s_b + B^{-1} I_n s_n = 0. 
 */
static
SCIP_RETCODE addRows(
   SCIP*                    scip,               /**< SCIP data structure */
   SCIP_SEPADATA *      sepadata,
   vector< SCIP_Real*> & tableau, /**< tableau rows corresponding to the variables in vars */
   SCIP_ROWPREP*         rowprep,            /**< store cut */
   SCIP_Bool*            success             /**< pointer to store whether the generation of cutcoefs was successful */
   )
{
   //SCIPdebugMessage("add row\n");
   *success = TRUE;

   int nrays = 0;
   auto & mlf = sepadata->mlf;
   /* loop over non-basic (slack) variables */
   SCIP_ROW** rows = SCIPgetLPRows(scip);
   int nrows = SCIPgetNLPRows(scip);
   int ncols = SCIPgetNLPCols(scip);
   SCIPdebugMessage("nrows: %d\n", nrows);
   for(int i = 0; i < nrows; ++i )
   {

      //SCIPdebugMessage("add row 1\n");
      SCIP_ROW* row = rows[i];

      SCIP_Real factor;
      /* set factor to store entries of ray as = [-BinvL, BinvU] */
      if( SCIProwGetBasisStatus(row) == SCIP_BASESTAT_LOWER )
         factor = 1.0;
      else if( SCIProwGetBasisStatus(row) == SCIP_BASESTAT_UPPER )
         factor = -1.0;
      else if( SCIProwGetBasisStatus(row) == SCIP_BASESTAT_ZERO )
      {
         *success = FALSE;
         return SCIP_OKAY;
      }
      else
         continue;


      //SCIPdebugMessage("add row 2\n");
      /* build the ray  and compute intersection point*/
      SCIP_Real interpoint = 0;
      if(sepadata->disaggregate){
         bool israyzero = true;
         bool iszero = true;
         auto & ray = mlf.ray;
         auto & extendedfvars = *mlf.extendedfvars;
         int ftype = mlf.typef ? 1 : 2;
         int  varsize = mlf.size;
         mlf.linearray = 0;
         //unordered_map<SCIP_VAR*,SCIP_Real> raysbycol;
         for(int j = 0; j < extendedfvars.size(); j++)
         {
            SCIP_VAR * var = extendedfvars[j];
            SCIP_Real ray_val = 0;
            if( tableau[j] != NULL )
               ray_val = factor * ( SCIPisZero(scip, tableau[j][ncols + i]) ? 0.0 : tableau[j][ncols + i] );
            else
            {
               ray_val = 0.0;
            }
            //israynonzero = israynonzero || (j != mlf && ray_val != 0.0);
            if(j < varsize){
               ray[j] = ray_val;
               bool affected  = ( mlf.probvartypes[j] == ftype || mlf.probvartypes[j] == 3 );
               bool rayjzero =  (  affected  ? SCIPisZero(scip, ray[j]) : true);
               israyzero = israyzero && rayjzero;
               iszero = iszero && rayjzero;
            }
            else{
               mlf.linearray += ray_val * (*mlf.mlf_term)[j - varsize].second;
            }
         }
         israyzero = israyzero && SCIPisZero(scip, mlf.linearray);
         //SCIPdebugMessage("add col 3 %d\n", int(israynonzero));
         /* do nothing if ray is 0 */
         if( israyzero )
            continue;


         interpoint = computeIntersectionPointf(scip, sepadata, iszero, success);
      }
      else {
         bool iszero = true;
         bool israyzero = true;
         auto & ray = mlf.ray;
         auto & probvars = mlf.probvars;
         #ifdef DEBUG_TIME
         time(&t5);
         #endif
         for(int j = 0; j < probvars.size(); j++)
         {
            SCIP_VAR * var = probvars[j];
            if( tableau[j] != NULL )
               ray[j] = factor * ( SCIPisZero(scip, tableau[j][ncols + i]) ? 0.0 : tableau[j][ncols + i] );
            else
            {
               ray[j]= 0.0;
            }
            israyzero = israyzero &&  SCIPisZero(scip, ray[j]);
            iszero = iszero && ( (mlf.probvartypes[j] == 2 ||  mlf.probvartypes[j] == 3)  ? SCIPisZero(scip, ray[j]) : true);
         }
         #ifdef DEBUG_TIME
         time(&t6);
         memory_time += difftime(t6,t5);
         #endif

         //SCIPdebugMessage("add col 3 %d\n", int(israynonzero));
         /* do nothing if ray is 0 */
         if(israyzero)
            continue;

         mlf.linearf1ray  = 0;
         for(int ind: mlf.set1)
         {
				mlf.linearf1ray += ray[ind] * mlf.linearize1[ind];
         }         
         mlf.linearray = 0;
         for(auto p: mlf.linear_term)
         {
				mlf.linearray += ray[p.first] * p.second;
         }
         
         /* compute intersection point */
         #ifdef DEBUG_TIME
         time(&t5);
         #endif
         interpoint = computeIntersectionPoint(scip, sepadata, iszero, success);
         #ifdef DEBUG_TIME
         time(&t6);
         compute_time += difftime(t6,t5);
         #endif
      }

      if(SCIPisZero(scip, interpoint)){
         *success = FALSE;
         break;
      }

  
      if( *success == FALSE )
         return SCIP_OKAY;
         

      //SCIPdebugMessage("add row 4\n");

      /* count nonzero rays */
      nrays += 1;

      /* compute cut coef */
      SCIP_Real cutcoef = SCIPisInfinity(scip, interpoint) ? 0.0 : 1.0 / interpoint;



      /* add var to cut: if variable is nonbasic at upper we have to flip sign of cutcoef */
      assert(SCIProwGetBasisStatus(row) == SCIP_BASESTAT_UPPER || SCIProwGetBasisStatus(row) == SCIP_BASESTAT_LOWER);
      SCIP_CALL( addRowToCut(scip, rowprep, SCIProwGetBasisStatus(row) == SCIP_BASESTAT_UPPER ? cutcoef :
            -cutcoef, row, success) ); /* rows have flipper base status! */
   }

   //SCIPdebugMessage("add row__ nay:%d\n", nrays);
   return SCIP_OKAY;
}

/** separates cuts for intersection cuts */
SCIP_RETCODE separateSubInter(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA *       sepadata,
   SCIP_SEPA*            sepa,               /**< separator */
   int *                 basicvarpos2tableaurow,/**< unordered_map from basic var to its tableau row */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   //SCIPdebugMessage("start sepinter!\n");
   SCIP_ROWPREP* rowprep;

   int ncols = SCIPgetNLPCols(scip);
   int nrows = SCIPgetNLPRows(scip);

   SCIP_Bool success = TRUE;

   /* cut (in the nonbasic space) is of the form alpha^T x >= 1 */
   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, SCIP_SIDETYPE_LEFT, TRUE) );
   SCIProwprepAddSide(rowprep, 1.0);

   vector<SCIP_Real*> tableau;
   if(sepadata->disaggregate){
      auto & mlf = sepadata->mlf;
      tableau = vector<SCIP_Real*>((*mlf.extendedfvars).size(), NULL);
      SCIP_CALL( getTableauRows(scip, *mlf.extendedfvars, basicvarpos2tableaurow, tableau,  &success) );
   }
   else{
      auto & mlf = sepadata->mlf;
      tableau = vector<SCIP_Real*>(mlf.size, NULL);
      SCIP_CALL( getTableauRows(scip, mlf.probvars, basicvarpos2tableaurow, tableau,  &success) );
   }
   assert(tableau.size() != 0 );
   //SCIPdebugMessage("get tableau finish!\n");

   if(! success )
      goto CLEANUP;



   //SCIPdebugMessage("col to add!\n");

   /* loop over each non-basic var; get the ray; compute cut coefficient */
   SCIP_CALL( addCols(scip, sepadata, tableau, rowprep, &success) );

   if( ! success )
      goto CLEANUP;
   //SCIPdebugMessage("row to add!\n");

   /* loop over non-basic slack variables */
   SCIP_CALL( addRows(scip, sepadata, tableau, rowprep, &success) );

   //SCIPdebugMessage("row added!\n");
   if( ! success )
      goto CLEANUP;


   /* merge coefficients that belong to same variable */

   //SCIPprintRowprep( scip,   rowprep,  NULL);
   
   //SCIPdebugMessage("row added0.5!\n");

   SCIPmergeRowprepTerms(scip, rowprep);

   SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, NULL, sepadata->mincutviol, NULL, &success) );

   //SCIPprintRowprep( scip,   rowprep,  NULL);

   //SCIPdebugMessage("row added1!\n");
   /* if cleanup was successfull, create row out of rowprep and add it */
   if( success )
   {
      SCIP_ROW* row;
      SCIP_Bool infeasible;

      /* create row */
      SCIP_CALL( SCIPgetRowprepRowSepa(scip, &row, rowprep, sepa) );

      assert(SCIPgetCutEfficacy(scip, NULL, row) > 0.0);

      /* add row */
      SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );


      SCIPdebugMessage("%f %f<=%f<=%f %d \n", SCIPgetCutEfficacy(scip, NULL, row), SCIProwGetLhs(row), SCIPgetRowActivity(scip, row), SCIProwGetRhs(row), infeasible);

      if( infeasible )
         *result = SCIP_CUTOFF;
      else
         *result = SCIP_SEPARATED;

      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

   //SCIPdebugMessage("row added2!\n");
CLEANUP:
   SCIPfreeRowprep(scip, &rowprep);
   //SCIPfreeBuffer(scip, &rayslppos);
   //SCIPfreeBuffer(scip, &rays);
   //SCIPfreeBuffer(scip, &interpoints);

   for(int i = 0; i < tableau.size(); i++){
      if(tableau[i] != NULL)
         SCIPfreeBufferArrayNull(scip, &tableau[i]);
   }
   return SCIP_OKAY;
}

/** separates cuts for stored signomial terms */
static
SCIP_RETCODE separateSubOuter(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_RESULT*          result,              /**< pointer to store the result of the separation call */
   SCIP_Bool &           success 
   )
{
   //SCIPdebugMessage("start sep!\n");
   SCIP_ROWPREP* rowprep;
   success = TRUE;

   /* cut (in the nonbasic space) is of the form alpha^T x <= -constant */
   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, SCIP_SIDETYPE_RIGHT, TRUE));

   //unordered_map<SCIP_VAR *, SCIP_Real> coefficients;
   /* check if we have the tableau row of the variable and if not compute it */
   //SCIPdebugMessage("get tableau!\n");

   //SCIP_Real cut = 0;
   auto & mlf = sepadata->mlf;

   if(! success )
      goto CLEANUP;

   if(sepadata->disaggregate){
      // f1(x) \le w*lift(x)
      for(int varind: *mlf.fset){
         SCIP_CALL(SCIPaddRowprepTerm(scip, rowprep, mlf.probvars[varind], (*mlf.linearize)[varind]));
      }
      for(int i = 0; i < mlf.fsize; i++){
         SCIP_CALL(SCIPaddRowprepTerm(scip, rowprep, (*mlf.liftfvars)[i], -(*mlf.mlf_term)[i].second));
      }
      SCIProwprepAddSide(rowprep, 0);
   }
   else{
      //SCIPdebugMessage("S-free construct!\n");
      for(int varind: mlf.set1){
         SCIP_CALL(SCIPaddRowprepTerm(scip, rowprep, mlf.probvars[varind], mlf.linearize1[varind]));
      }

      for(auto &p: mlf.linear_term){
         SCIP_CALL(SCIPaddRowprepTerm(scip, rowprep, mlf.probvars[p.first], p.second));
      }

      SCIProwprepAddSide(rowprep, -mlf.constant);
   }

   //SCIPdebugMessage("row added!\n");

   /* merge coefficients that belong to same variable */
   SCIPmergeRowprepTerms(scip, rowprep);

   SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, NULL, sepadata->mincutviol, NULL, &success) );

   /* if cleanup was successfull, create row out of rowprep and add it */
   if( success )
   {
      SCIP_ROW* row;
      SCIP_Bool infeasible;

      /* create row */
      SCIP_CALL( SCIPgetRowprepRowSepa(scip, &row, rowprep, sepa) );

      assert(SCIPgetCutEfficacy(scip, NULL, row) > 0.0);

      /* add row */
      SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );

      //SCIPdebugMessage("mis\n");
      //SCIPprintRow ( scip,row, NULL ) ;
      //SCIPdebugMessage("mis\n");

      SCIPdebugMessage("%f %d \n", SCIPgetCutEfficacy(scip, NULL, row), infeasible);
      if( infeasible )
         *result = SCIP_CUTOFF;
      else
         *result = SCIP_SEPARATED;

      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

CLEANUP:
   SCIPfreeRowprep(scip, &rowprep);
   //SCIPfreeBuffer(scip, &rayslppos);
   //SCIPfreeBuffer(scip, &rays);
   //SCIPfreeBuffer(scip, &interpoints);

   return SCIP_OKAY;
}

/** separate cuts for stored signomial */
static
SCIP_RETCODE separatePoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SepaData *        sepadata,
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   assert(sepa != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

	//ProbData * probdata = NULL;
   //SCIPdebugMessage("%p\n", scip);
   //return SCIP_OKAY;
	//probdata = dynamic_cast<ProbData *>(SCIPgetObjProbData(scip));
	//assert(probdata != NULL);

   int *  basicvarpos2tableaurow = NULL; /* map between basic var and its tableau row */

   // check whether to separate
   auto & mlf = sepadata->mlf;

   //vector<SCIP_Real> vsol(mlf.probvars.size());
   //SCIPdebugMessage(">>>>\n");
   int id = 0;
   SCIP_Real gradient;
   for(auto var: mlf.probvars){
      mlf.sol[id] = SCIPvarGetLPSol(var);
      id++;
   }


   #ifdef DEBUG_TIME
   time_t start,finish;
   time(&start);
   memory_time = 0;
   compute_time = 0;
   ct_zero = 0;
   ct_full = 0;
   sort_time = 0;
   coef_time = 0;
   loop_time = 0;
   inner_loop_time = 0;
   #endif

   //SCIPdebugMessage("%d >>>>\n", sepadata->disaggregate);

   if(sepadata->disaggregate)
   {
      //SCIP_Real f1val = mlf.evalRayf1WithZero();
      SCIP_Bool success;
      // separate outer approximation cuts for f1

      if(mlf.mlf1_term.size() != 0){
         id = 0;
         for(auto var: mlf.liftf1vars){
            mlf.solf1[id] = SCIPvarGetLPSol(var);
            id++;
         }
         SCIP_Real f1val = mlf.evalRayfWithZero(TRUE);
         SCIPdebugMessage("f1:%f\n", f1val);
         //SCIPdebugMessage("%f\n", f1val);
         if(SCIPisGT(scip, f1val, 0)){
            SCIPdebugMessage("outer\n");
            SCIP_CALL(separateSubOuter(scip, sepa, sepadata, result, success));
         }
         else if(SCIPisLT(scip, f1val, 0)){
            SCIP_CALL( SCIPallocBufferArray(scip, &basicvarpos2tableaurow, SCIPgetNLPCols(scip)) );

            /* construct basicvar to tableau row map */
            SCIP_CALL( constructBasicVars2TableauRowMap(scip, basicvarpos2tableaurow));
            
            //SCIPdebugMessage("1\n");
            SCIP_CALL(separateSubInter(scip, sepadata, sepa,  basicvarpos2tableaurow, result));
            //SCIPdebugMessage("2\n");

            /* all feasible, so no memory to free */
            if( basicvarpos2tableaurow == NULL )
               return SCIP_OKAY;
            SCIPfreeBufferArray(scip, &basicvarpos2tableaurow);
           
         }
      }


      if(mlf.mlf2_term.size() != 0){
         id = 0;
         for(auto var: mlf.liftf2vars){
            mlf.solf2[id] = SCIPvarGetLPSol(var);
            id++;
         }
         SCIP_Real f2val = mlf.evalRayfWithZero(FALSE);
         SCIPdebugMessage("f2:%f\n", f2val);
         //SCIPdebugMessage("%f\n", f2val);
         if(SCIPisGT(scip, f2val, 0)){
            SCIPdebugMessage("outer\n");
            SCIP_CALL(separateSubOuter(scip, sepa, sepadata, result, success));
         }
         else if(SCIPisLT(scip, f2val, 0)){
            SCIP_CALL( SCIPallocBufferArray(scip, &basicvarpos2tableaurow, SCIPgetNLPCols(scip)) );

            /* construct basicvar to tableau row map */
            SCIP_CALL( constructBasicVars2TableauRowMap(scip, basicvarpos2tableaurow));
            
            //SCIPdebugMessage("1\n");
            SCIP_CALL(separateSubInter(scip, sepadata, sepa, basicvarpos2tableaurow, result));
            //SCIPdebugMessage("2\n");

            /* all feasible, so no memory to free */
            if( basicvarpos2tableaurow == NULL )
               return SCIP_OKAY;
            SCIPfreeBufferArray(scip, &basicvarpos2tableaurow);
           
         }
      }
      //abort();
      
   }
   else{
      SCIP_Real fval = mlf.evalRayWithZero();

      //SCIPdebugMessage("----%f %f %f\n", fval, mlf.evalRayWith(0, gradient, false), mlf.val);

      // separate if both forms are violated
      if(SCIPisLT(scip, fval, 0)){
         *result = SCIP_DIDNOTFIND;
         return SCIP_OKAY;
      }
      // f_2 null
      SCIP_Bool success;
      if(mlf.mlf2_term.size() == 0){
         //SCIPdebugMessage("start sep 2\n");
         SCIP_CALL(separateSubOuter(scip, sepa, sepadata, result, success));
         return SCIP_OKAY;
      }      

      // separate
      /* allocate memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &basicvarpos2tableaurow, SCIPgetNLPCols(scip)) );

      /* construct basicvar to tableau row map */
      SCIP_CALL( constructBasicVars2TableauRowMap(scip, basicvarpos2tableaurow));
      
      //SCIPdebugMessage("1\n");
      SCIP_CALL(separateSubInter(scip, sepadata, sepa, basicvarpos2tableaurow, result));
      //SCIPdebugMessage("2\n");

      /* all feasible, so no memory to free */
      if( basicvarpos2tableaurow == NULL )
         return SCIP_OKAY;
      SCIPfreeBufferArray(scip, &basicvarpos2tableaurow);
   }

   #ifdef DEBUG_TIME
   time(&finish);
   SCIP_Real t =  difftime(finish,start);
   SCIPdebugMessage("%d,%d mem:%f, compute:%f, total:%f  sort:%f, coef:%f, inner:%f, loop:%f\n", ct_zero, ct_full, memory_time, compute_time, t, sort_time, coef_time, inner_loop_time, loop_time);
   #endif
   return SCIP_OKAY;
}



/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
SCIP_DECL_SEPACOPY(SepaInterSub::scip_copy)
{  /*lint --e{715}*/
   //SCIPdebugMessage("copy %p\n", scip);
   return SCIP_OKAY;
}


/** destructor of separator to free user data (called when SCIP is exiting) */
SCIP_DECL_SEPAFREE(SepaInterSub::scip_free)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/** initialization method of separator (called after problem was transformed) */
SCIP_DECL_SEPAINIT(SepaInterSub::scip_init)
{  /*lint --e{715}*/
   /* get separator data */
	SCIP_CALL(SCIPgetRealParam(scip,  "separating/intersub/mincutviol", &sepadata.mincutviol ));
	SCIP_CALL(SCIPgetIntParam(scip,  "separating/intersub/maxroundsroot", &sepadata.maxroundsroot ));
   SCIP_CALL(SCIPgetIntParam(scip,  "separating/intersub/maxrounds", &sepadata.maxrounds ));
   SCIP_CALL(SCIPgetBoolParam(scip,  "separating/intersub/disaggregate", &sepadata.disaggregate ));
   SCIPdebugMessage("scip_init out\n");
   return SCIP_OKAY;
}


/** deinitialization method of separator (called before transformed problem is freed) */
SCIP_DECL_SEPAEXIT(SepaInterSub::scip_exit)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
SCIP_DECL_SEPAINITSOL(SepaInterSub::scip_initsol)
{  /*lint --e{715}*/
   // collect terms
   SCIPdebugMessage("scip_init sol\n");
   SCIP_CALL(detectMLFCons(scip, sepadata.mlf));
   return SCIP_OKAY;
}


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
SCIP_DECL_SEPAEXITSOL(SepaInterSub::scip_exitsol)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** LP solution separation method of separator */
SCIP_DECL_SEPAEXECLP(SepaInterSub::scip_execlp)
{  /*lint --e{715}*/
   int ncalls;
   int currentdepth;

   currentdepth = SCIPgetDepth(scip);
   ncalls = SCIPsepaGetNCallsAtNode(sepa);

   *result = SCIP_DIDNOTRUN;

   /* only call the separator a given number of times at each node */
   if( (currentdepth == 0 && sepadata.maxroundsroot >= 0 && ncalls >= sepadata.maxroundsroot)
      || (currentdepth > 0 && sepadata.maxrounds >= 0 && ncalls >= sepadata.maxrounds) )
   {
      SCIPdebugMsg(scip, "reached round limit for node\n");
      return SCIP_OKAY;
   }

   /* call separation method */
   //SCIPdebugMessage("execlp\n");
   SCIP_CALL( separatePoint(scip, &sepadata, sepa, result) );

   return SCIP_OKAY;
}


