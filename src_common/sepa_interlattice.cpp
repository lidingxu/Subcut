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

/**@file   sepa_interlattice.cpp
 * @brief  lattice separator with intersection cuts
 * @author Liding Xu
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <assert.h>
#include <string.h>

#include "sepa_interlattice.h"
#include "scip/cons_nonlinear.h"
#include "scip/struct_misc.h"


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
   list<SCIP_VAR *>  &    vars,       /**< signomial term */
   int * basicvarpos2tableaurow,  /**< map between basic var and its tableau row */
   unordered_map<SCIP_VAR*, SCIP_Real*> & tableau,     /**< map between all var and its tableau row */
   unordered_map<SCIP_VAR*, SCIP_Real*> & tableaurows,    /**< buffer to store tableau row associated with signomial vars */
   SCIP_Bool*     success             /**< set to TRUE if no variable had basisstat = ZERO */
   )
{

   *success = TRUE;

   int nrows = SCIPgetNLPRows(scip);
   int ncols = SCIPgetNLPCols(scip);

   /* check if we have t   abort();he tableau row of the variable and if not compute it */
   for(auto var: vars)
   {
      if( tableau.find(var) == tableau.end())
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
            tableau[var] = densetableaurow;
         }
         else if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_ZERO)
         {
            *success = FALSE;
            return SCIP_OKAY; /* don't even bother */
         }
         else
         {
            tableau[var] = (SCIP_Real *)NULL;
         }
      }
      /* get tableau row of var */
      tableaurows[var] = tableau[var];
   }
   return SCIP_OKAY;
}

/** adds cutcoef * (col    abort();- col*) to rowprep */
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
   int BINSEARCH_MAXITERS,
   unordered_map<SCIP_VAR*, SCIP_Real> & ray,          /**< iterator of ray */
   ProbData *       probdata,                  /**< problem data structure */
   int var_id,
   SCIP_Bool*            success             /**< pointer to store whether the generation of cutcoefs was successful */
   )
{
   SCIP_Real vray = ray[probdata->bin_vars[var_id]];
   bool iszero =  SCIPisZero(scip,  vray);
   SCIP_Real vsol  = SCIPvarGetLPSol(probdata->bin_vars[var_id]);
   SCIP_Real t;
   if(iszero){
      t = SCIPinfinity(scip);
   }
   else{
      if(vray > 0){
         t = (1 - vsol) / vray;
      }
      else{
         t = -vsol / vray;
      }
   }
   return t;
}

/** computes the cut coefs of the  non-basic (non-slack) variables (correspond to cols) and adds them to the
 * A x + I s = 0, and A_b x_b + I_b s_b + A_n x_n + I_n s_n = 0. B = [A_b, I_b]
 * B^{-1} (Ax + Is) = 0,  x_b +  B^{-1}A_n x_n + s_b + B^{-1} I_n s_n = 0.
 */
static
SCIP_RETCODE addCols(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SepaData_Lattice * sepadata,
   list<SCIP_VAR*> & probvars,             /**< problem variables */
   ProbData * probdata,                  /**< problem data structure */
   int var_id,
   unordered_map<SCIP_VAR*, SCIP_Real*> & tableaurows, /**< tableau rows corresponding to the variables in vars */
   SCIP_ROWPREP*         rowprep,            /**< store cut */
   SCIP_Bool*            success             /**< pointer to store whether the generation of cutcoefs was successful */
   )
{
   //SCIPdebugMessage("add col %d\n", int(tableaurows.size()));
   *success = TRUE;

   int nrays = 0;

   /* loop over non-basic (non-slack) variables */
   SCIP_COL** cols = SCIPgetLPCols(scip);
   int ncols = SCIPgetNLPCols(scip);
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

      /* build the ray */
      SCIP_Bool israynonzero = FALSE;
      unordered_map<SCIP_VAR*,SCIP_Real> raysbycol;
      for(auto var: probvars)
      {
         //SCIPdebugMessage("add col 2.1 %d %d\n",  int(tableaurows.find(var) != tableaurows.end()), int(tableaurows[var] != NULL));
         if( tableaurows.find(var) != tableaurows.end() && tableaurows[var] != NULL )
         {
            //SCIPdebugMessage("add col 2.11\n");
            raysbycol[var] = factor * ( SCIPisZero(scip, tableaurows[var][i]) ? 0.0 : tableaurows[var][i] );
            //SCIPdebugMessage("add col 2.12\n");
         }
         else
         {
            //SCIPdebugMessage("add col 2.2\n");
            if( col == SCIPvarGetCol(var) )
               raysbycol[var]= (-factor);
            else
               raysbycol[var]=(0.0);
         }
         //SCIPdebugMessage("%f\n", raysbycol[var]);
         israynonzero = israynonzero || (raysbycol[var] != 0.0);
      }

      //SCIPdebugMessage("add col 3 %d\n", int(israynonzero));
      /* do nothing if ray is 0 */
      if( ! israynonzero )
         continue;

      //SCIPdebugMessage("add col 3\n");
      

      /* compute intersection point */
      SCIP_Real interpoint = computeIntersectionPoint(scip, sepadata->BINSEARCH_MAXITERS, raysbycol, probdata, var_id, success);


      if(SCIPisZero(scip, interpoint)){
         *success = FALSE;
         break;
      }

      //SCIPdebugMessage("add col 4\n");

      if( *success == FALSE )
         return SCIP_OKAY;

      /* count nonzero rays */
      nrays += 1;

      /* compute cut coef */
      SCIP_Real cutcoef = SCIPisInfinity(scip, interpoint) ? 0.0 : 1.0 / interpoint;


      //SCIPdebugMessage("%f\n", cutcoef); cutcoef *= 0.01;SCIPdebugMessage("%f\n", cutcoef); 
      //cutcoef = 0;

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
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SepaData_Lattice * sepadata,
   list<SCIP_VAR*> & probvars,             /**< problem variables */
   ProbData * probdata,                  /**< problem data structure */
   int var_id,
   unordered_map<SCIP_VAR*, SCIP_Real*> & tableaurows, /**< tableau rows corresponding to the variables in vars */
   SCIP_ROWPREP*         rowprep,            /**< store cut */
   SCIP_Bool*            success             /**< pointer to store whether the generation of cutcoefs was successful */
   )
{
   //SCIPdebugMessage("add row\n");
   *success = TRUE;

   int nrays = 0;
   /* loop over non-basic (slack) variables */
   SCIP_ROW** rows = SCIPgetLPRows(scip);
   int nrows = SCIPgetNLPRows(scip);
   int ncols = SCIPgetNLPCols(scip);
   for(int i = 0; i < nrows; ++i )
   {
      
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
      /* build the ray */
      SCIP_Bool israynonzero = FALSE;
      unordered_map<SCIP_VAR*,SCIP_Real> raysbyrow;
      for(auto var: probvars)
      {
         if( tableaurows.find(var) != tableaurows.end() && tableaurows[var] != NULL )
            raysbyrow[var] = factor * ( SCIPisZero(scip, tableaurows[var][ncols + i]) ? 0.0 : tableaurows[var][ncols + i] );
         else
         {
            raysbyrow[var]= 0.0;
         }
         israynonzero = israynonzero || (raysbyrow[var] != 0.0);
      }

      //SCIPdebugMessage("add row 3\n");
      /* do nothing if ray is 0 */
      if( ! israynonzero )
         continue;


      /* compute intersection point */
      
      SCIP_Real interpoint = computeIntersectionPoint(scip, sepadata->BINSEARCH_MAXITERS, raysbyrow, probdata, var_id, success);

      //SCIPdebugMessage("row: %f %d\n", interpoint, SCIPisZero(scip, interpoint));


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

      //SCIPdebugMessage("%f\n", cutcoef); cutcoef *= 0.01;SCIPdebugMessage("%f\n", cutcoef); 

      /* add var to cut: if variable is nonbasic at upper we have to flip sign of cutcoef */
      assert(SCIProwGetBasisStatus(row) == SCIP_BASESTAT_UPPER || SCIProwGetBasisStatus(row) == SCIP_BASESTAT_LOWER);
      SCIP_CALL( addRowToCut(scip, rowprep, SCIProwGetBasisStatus(row) == SCIP_BASESTAT_UPPER ? cutcoef :
            -cutcoef, row, success) ); /* rows have flipper base status! */
   }

   //abort();
   //run++;
   //SCIPdebugMessage("add row__ nay:%d\n", nrays);
   return SCIP_OKAY;
}



static SCIP_RETCODE rowprepCleanupSortTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep             /**< rowprep to be sorted */
   )
{
   int i;

   assert(scip != NULL);
   assert(rowprep != NULL);

   /* special treatment for cuts with few variables */
   switch( rowprep->nvars )
   {
      case 0:
      case 1:
         break;

      case 2:
      {
         if( REALABS(rowprep->coefs[0]) < REALABS(rowprep->coefs[1]) )
         {
            SCIP_Real tmp1;
            SCIP_VAR* tmp2;

            tmp1 = rowprep->coefs[0];
            rowprep->coefs[0] = rowprep->coefs[1];
            rowprep->coefs[1] = tmp1;

            tmp2 = rowprep->vars[0];
            rowprep->vars[0] = rowprep->vars[1];
            rowprep->vars[1] = tmp2;
         }
         break;
      }

      default :
      {
         SCIP_Real* abscoefs;

         SCIP_CALL( SCIPallocBufferArray(scip, &abscoefs, rowprep->nvars) );
         for( i = 0; i < rowprep->nvars; ++i )
            abscoefs[i] = REALABS(rowprep->coefs[i]);
         SCIPsortDownRealRealPtr(abscoefs, rowprep->coefs, (void**)rowprep->vars, rowprep->nvars);
         SCIPfreeBufferArray(scip, &abscoefs);
      }
   }

   /* forget about coefs that are exactly zero (unlikely to have some) */
   while( rowprep->nvars > 0 && rowprep->coefs[rowprep->nvars-1] == 0.0 )
      --rowprep->nvars;

   return SCIP_OKAY;
}

/** separates cuts for intersection cuts */
SCIP_RETCODE separateLatticeInter(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SepaData_Lattice *       sepadata,
   SCIP_SEPA*            sepa,               /**< separator */
   ProbData *        probdata,           /**< problem data */
   int               var_id,
   int *                 basicvarpos2tableaurow,/**< unordered_map from basic var to its tableau row */
   unordered_map<SCIP_VAR *, SCIP_Real*> &  tableau,            /**< unordered_map from var to its tableau row */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   //SCIPdebugMessage("start sep!\n");
   SCIP_ROWPREP* rowprep;

   int ncols = SCIPgetNLPCols(scip);
   int nrows = SCIPgetNLPRows(scip);

   SCIP_Bool success = TRUE;

   SCIP_VAR** vars;
   SCIP_Real* coeffs;

   /* cut (in the nonbasic space) is of the form alpha^T x >= 1 */
   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, SCIP_SIDETYPE_LEFT, TRUE) );
   SCIProwprepAddSide(rowprep, 1.0);

   unordered_map<SCIP_VAR *, SCIP_Real*> tableaurows;
   /* check if we have the tableau row of the variable and if not compute it */
   //SCIPdebugMessage("get tableau!\n");

   //SCIPdebugMessage("31\n");

   list<SCIP_VAR * > probvars(probdata->bin_vars.begin(), probdata->bin_vars.end()); 

   probvars.push_back(probdata->obj_var);
   SCIP_CALL( getTableauRows(scip, probvars, basicvarpos2tableaurow, tableau, tableaurows, &success) );
   //SCIPdebugMessage("get tableau finish!\n");

   if(! success )
      goto CLEANUP;



   //SCIPdebugMessage("col to add!\n");

   /* loop over each non-basic var; get the ray; compute cut coefficient */
   SCIP_CALL( addCols(scip, sepadata, probvars, probdata, var_id, tableaurows, rowprep, &success) );

   if( ! success )
      goto CLEANUP;
   //SCIPdebugMessage("col added!\n");

   /* loop over non-basic slack variables */
   SCIP_CALL( addRows(scip, sepadata, probvars, probdata, var_id, tableaurows, rowprep, &success) );

   //SCIPdebugMessage("row added %d!\n", success);
   if( ! success )
      goto CLEANUP;


   /* merge coefficients that belong to same variable */

   
   //SCIPdebugMessage("row added0.5!\n");
   //f = SCIPgetRowprepViolation ( scip, 	rowprep, NULL,   &	reliable   ); 
   //SCIPdebugMessage("%f %d\n", f, reliable);	

   SCIPmergeRowprepTerms(scip, rowprep);

   //f = SCIPgetRowprepViolation ( scip, 	rowprep, NULL,   &	reliable   ); 
   //SCIPdebugMessage("%f %d\n", f, reliable);	
   

   //SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, NULL, 0, NULL, &success) );
   //SCIP_CALL( SCIPcleanupRowprep2(scip, rowprep, NULL, 0, NULL, &success) );

   //f = SCIPgetRowprepViolation ( scip, 	rowprep, NULL,   &	reliable   ); 
   //SCIPdebugMessage("%f %d\n", f, reliable);	

   //SCIPdebugMessage("%d\n", success);
   //SCIPprintRowprep(scip, rowprep, NULL);
   if(sepadata->cleanrowprep)
   {
      SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, NULL, 0, NULL, &success) );
   }
   else{
      vars = SCIProwprepGetVars(rowprep);
      coeffs = SCIProwprepGetCoefs(rowprep);
      for(int i = 0; i < SCIProwprepGetNVars(rowprep); i++){
         if(fabs(coeffs[i]) < 1e-6 ){
            coeffs[i] = 0;// SCIPdebugMessage("%s %f\n", SCIPvarGetName(vars[i]), coeffs[i]);
         }
      }
      SCIP_CALL(rowprepCleanupSortTerms(scip, rowprep) );
   }
  

   

   //SCIPprintRowprep(scip, rowprep, NULL);
   //success = true;

   //SCIPdebugMessage("row added1!\n");
   /* if cleanup was successfull, create row out of rowprep and add it */
   if( success )
   {
      SCIP_ROW* row;
      SCIP_Bool infeasible;



      /* create row */
      SCIP_CALL( SCIPgetRowprepRowSepa(scip, &row, rowprep, sepa) );

      //SCIPprintRow(scip, row, NULL);

      assert(SCIPgetCutEfficacy(scip, NULL, row) > 0.0);

      /* add row */
      SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );


      SCIPdebugMessage("%f %f<=%f<=%f %d \n", SCIPgetCutEfficacy(scip, NULL, row), SCIProwGetLhs(row), SCIPgetRowActivity(scip, row), SCIProwGetRhs(row), infeasible);

      //abort();

      if( infeasible )
         *result = SCIP_CUTOFF;
      else
         *result = SCIP_SEPARATED;

      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }
   //abort();

   //SCIPdebugMessage("row added2!\n");
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
   SCIP_SepaData_Lattice *        sepadata,
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{

   assert(sepa != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

	ProbData * probdata = NULL;
   //SCIPdebugMessage("%p\n", scip);
   //return SCIP_OKAY;
	probdata = dynamic_cast<ProbData *>(SCIPgetObjProbData(scip));
	assert(probdata != NULL);

   *result = SCIP_DIDNOTFIND;



   unordered_map<SCIP_VAR *, SCIP_Real *> tableau;
   int *  basicvarpos2tableaurow = NULL; /* map between basic var and its tableau row */

   // check whether to separate
   int numvars = probdata->numvars;

   int id = 0;
   SCIP_Real fractional = 0;
   for(int i = 0; i < numvars; i++){
      SCIP_Real vsol = SCIPvarGetLPSol(probdata->bin_vars[i]);
      SCIP_Real fractional_ = min(vsol, 1 - vsol);
      if(fractional_ > fractional){
         fractional = fractional_;
         id = i;
      }
   }


   bool display = false; 
   SCIPdebugMessage("%f %d\n", fractional, id);
   //return SCIP_OKAY;
   if(SCIPisFeasZero(scip, fractional)){
      return SCIP_OKAY;
   }

   // separate
   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &basicvarpos2tableaurow, SCIPgetNLPCols(scip)) );

   /* construct basicvar to tableau row map */
   SCIP_CALL( constructBasicVars2TableauRowMap(scip, basicvarpos2tableaurow) );
   
   //SCIPdebugMessage("sep inter\n");
   SCIP_CALL(separateLatticeInter(scip, sepadata, sepa, probdata, id, basicvarpos2tableaurow, tableau, result));
   //SCIPdebugMessage("2\n");

   for(auto & table: tableau){
      SCIPfreeBufferArrayNull(scip, &table.second);
   }

   /* all minors were feasible, so no memory to free */
   if( basicvarpos2tableaurow == NULL )
      return SCIP_OKAY;
   SCIPfreeBufferArray(scip, &basicvarpos2tableaurow);

   return SCIP_OKAY;
}

/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
SCIP_DECL_SEPACOPY(SepaInterLattice::scip_copy)
{  /*lint --e{715}*/
   //SCIPdebugMessage("copy %p\n", scip);
   return SCIP_OKAY;
}


/** destructor of separator to free user data (called when SCIP is exiting) */
SCIP_DECL_SEPAFREE(SepaInterLattice::scip_free)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** initialization method of separator (called after problem was transformed) */
SCIP_DECL_SEPAINIT(SepaInterLattice::scip_init)
{  /*lint --e{715}*/
   //SCIPdebugMessage("init %p\n", scip);
   /* get separator data */
	SCIP_CALL(SCIPgetRealParam(scip,  "separating/interlattice/mincutviol", &sepadata.mincutviol ));
	SCIP_CALL(SCIPgetIntParam(scip,  "separating/interlattice/maxroundsroot", &sepadata.maxroundsroot ));
   SCIP_CALL(SCIPgetIntParam(scip,  "separating/interlattice/maxrounds", &sepadata.maxrounds ));
   SCIP_CALL(SCIPgetBoolParam(scip,  "separating/interlattice/cleanrowprep", &sepadata.cleanrowprep));
   return SCIP_OKAY;
}


/** deinitialization method of separator (called before transformed problem is freed) */
SCIP_DECL_SEPAEXIT(SepaInterLattice::scip_exit)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
SCIP_DECL_SEPAINITSOL(SepaInterLattice::scip_initsol)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
SCIP_DECL_SEPAEXITSOL(SepaInterLattice::scip_exitsol)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** LP solution separation method of separator */
SCIP_DECL_SEPAEXECLP(SepaInterLattice::scip_execlp)
{  /*lint --e{715}*/
   int ncalls;
   int currentdepth;

   //SCIPdebugMessage("execl sep\n");
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
   SCIP_CALL( separatePoint(scip, &sepadata, sepa, result) );

   return SCIP_OKAY;
}


