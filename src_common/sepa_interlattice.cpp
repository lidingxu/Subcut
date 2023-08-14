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

#ifdef INTERCUTS_VERBOSE
#define INTER_LOG
#endif

#ifdef INTER_LOG
#define INTERLOG(x) if( SCIPgetSubscipDepth(scip) == 0 && SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_NORMAL ) { x }
#else
#define INTERLOG(x)
#endif

/* fundamental nonlinear handler properties */
#define SEPA_NAME                      "interlattice"

/* some default values */
#define INTERCUTS_MINVIOL              1e-4
#define BINSEARCH_MAXITERS             500
#define DEFAULT_NCUTSROOT              80
#define DEFAULT_NCUTS                  2

/* structure to store rays. note that for a given ray, the entries in raysidx are sorted. */
struct Rays
{
   SCIP_Real*            rays;               /**< coefficients of rays */
   int*                  raysidx;            /**< to which var the coef belongs; vars are [quadvars, linvars, auxvar] */
   int*                  raysbegin;          /**< positions of rays: coefs of i-th ray [raybegin[i], raybegin[i+1]) */
   int*                  lpposray;           /**< lp pos of var associated with the ray;
                                               >= 0 -> lppos of var; < 0 -> var is row -lpposray -1 */
   int                   rayssize;           /**< size of rays and rays idx */
   int                   nrays;              /**< size of raysbegin is nrays + 1; size of lpposray */
};
typedef struct Rays RAYS;




 /*
  * Callback methods of separator
  */


/*
 * static methods
 */

/** adds cutcoef * (col - col*) to rowprep */
static
SCIP_RETCODE addColToCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to store intersection cut */
   SCIP_SOL*             sol,                /**< solution to separate */
   SCIP_Real             cutcoef,            /**< cut coefficient */
   SCIP_COL*             col                 /**< column to add to rowprep */
   )
{
   assert(col != NULL);

#ifdef DEBUG_INTERCUTS_NUMERICS
   SCIPinfoMessage(scip, NULL, "adding col %s to cut. %g <= col <= %g\n", SCIPvarGetName(SCIPcolGetVar(col)),
      SCIPvarGetLbLocal(SCIPcolGetVar(col)), SCIPvarGetUbLocal(SCIPcolGetVar(col)));
   SCIPinfoMessage(scip, NULL, "col is active at %s. Value %.15f\n", SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_LOWER ? "lower bound" :
      "upper bound" , SCIPcolGetPrimsol(col));
#endif

   SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, SCIPcolGetVar(col), cutcoef) );
   SCIProwprepAddConstant(rowprep, -cutcoef * SCIPgetSolVal(scip, sol, SCIPcolGetVar(col)) );

   return SCIP_OKAY;
}

/** adds cutcoef * (slack - slack*) to rowprep
 *
  * row is lhs &le; <coefs, vars> + constant &le; rhs, thus slack is defined by
  * slack + <coefs.vars> + constant = side
  *
  * If row (slack) is at upper, it means that <coefs,vars*> + constant = rhs, and so
  * slack* = side - rhs --> slack - slack* = rhs - <coefs, vars> - constant.
  *
  * If row (slack) is at lower, then <coefs,vars*> + constant = lhs, and so
  * slack* = side - lhs --> slack - slack* = lhs - <coefs, vars> - constant.
  */
static
SCIP_RETCODE addRowToCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to store intersection cut */
   SCIP_Real             cutcoef,            /**< cut coefficient */
   SCIP_ROW*             row,                /**< row, whose slack we are ading to rowprep */
   SCIP_Bool*            success             /**< if the row is nonbasic enough */
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

#ifdef DEBUG_INTERCUTS_NUMERICS
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

/** constructs map between lp position of a basic variable and its row in the tableau */
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

/** counts the number of basic variables in the signomial expr */
static
int countBasicVars(
   ProbData *      probdata,
   SCIP_VAR*             auxvar,             /**< aux var of expr or NULL if not needed (e.g. separating real cons) */
   SCIP_Bool*            nozerostat          /**< whether there is no variable with basis status zero */
   )
{
 
   SCIP_COL* col;
   int i;
   int nbasicvars = 0;

   *nozerostat = TRUE;

   /* loop over signomial vars */
   for( i = 0; i < probdata->numvars; ++i )
   {
      col = SCIPvarGetCol(probdata->bin_vars[i]);
      if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_BASIC )
         nbasicvars += 1;
      else if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_ZERO )
      {
         *nozerostat = FALSE;
         return 0;
      }
   }

   /* finally consider the aux var (if it exists) */
   if( auxvar != NULL )
   {
      col = SCIPvarGetCol(auxvar);
      if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_BASIC )
         nbasicvars += 1;
      else if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_ZERO )
      {
         *nozerostat = FALSE;
         return 0;
      }
   }

   return nbasicvars;
}

/** stores the row of the tableau where `col` is basic
 *
 *  In general, we will have
 *
 *      basicvar1 = tableaurow var1
 *      basicvar2 = tableaurow var2
 *      ...
 *      basicvarn = tableaurow varn
 *
 *  However, we want to store the the tableau row by columns. Thus, we need to know which of the basic vars `col` is.
 *
 *  Note we only store the entries of the nonbasic variables
 */
static
SCIP_RETCODE storeDenseTableauRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COL*             col,                /**< basic columns to store its tableau row */
   int*                  basicvarpos2tableaurow,/**< map between basic var and its tableau row */
   int                   nbasiccol,          /**< which basic var this is */
   int                   raylength,          /**< the length of a ray (the total number of basic vars) */
   SCIP_Real*            binvrow,            /**< buffer to store row of Binv */
   SCIP_Real*            binvarow,           /**< buffer to store row of Binv A */
   SCIP_Real*            tableaurows         /**< pointer to store the tableau rows */
   )
{
   int nrows;
   int ncols;
   int lppos;
   int i;
   SCIP_COL** cols;
   SCIP_ROW** rows;
   int nray;

   assert(nbasiccol < raylength);
   assert(col != NULL);
   assert(binvrow != NULL);
   assert(binvarow != NULL);
   assert(tableaurows != NULL);
   assert(basicvarpos2tableaurow != NULL);
   assert(SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_BASIC);

   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );

   lppos = SCIPcolGetLPPos(col);

   assert(basicvarpos2tableaurow[lppos] >= 0);

   SCIP_CALL( SCIPgetLPBInvRow(scip, basicvarpos2tableaurow[lppos], binvrow, NULL, NULL) );
   SCIP_CALL( SCIPgetLPBInvARow(scip, basicvarpos2tableaurow[lppos], binvrow, binvarow, NULL, NULL) );

   nray = 0;
   for( i = 0; i < ncols; ++i )
      if( SCIPcolGetBasisStatus(cols[i]) != SCIP_BASESTAT_BASIC )
      {
         tableaurows[nbasiccol + nray * raylength] = binvarow[i];
         nray++;
      }
   for( ; i < ncols + nrows; ++i )
      if( SCIProwGetBasisStatus(rows[i - ncols]) != SCIP_BASESTAT_BASIC )
      {
         tableaurows[nbasiccol + nray * raylength] = binvrow[i - ncols];
         nray++;
      }

   return SCIP_OKAY;
}

/** stores the rows of the tableau corresponding to the basic variables in the quadratic expression
 *
 * Also return a map storing to which var the entry of a ray corresponds, i.e., if the tableau is
 *
 *     basicvar_1 = ray1_1 nonbasicvar_1 + ...
 *     basicvar_2 = ray1_2 nonbasicvar_1 + ...
 *     ...
 *     basicvar_n = ray1_n nonbasicvar_1 + ...
 *
 * The map maps k to the position of basicvar_k in the variables of the constraint assuming the variables are sorted as
 * [vars, auxvar].
 */
static
SCIP_RETCODE storeDenseTableauRowsByColumns(
   SCIP*                 scip,               /**< SCIP data structure */
   ProbData *      probdata,
   int                   raylength,          /**< length of a ray of the tableau */
   SCIP_VAR*             auxvar,             /**< aux var of expr or NULL if not needed (e.g. separating real cons) */
   SCIP_Real*            tableaurows,        /**< buffer to store the tableau rows */
   int*                  rayentry2conspos    /**< buffer to store the map */
   )
{
   SCIP_Real* binvarow;
   SCIP_Real* binvrow;
   SCIP_COL* col;
   int* basicvarpos2tableaurow; /* map between basic var and its tableau row */
   int nrayentries;
   int nrows;
   int ncols;
   int i;

   nrows = SCIPgetNLPRows(scip);
   ncols = SCIPgetNLPCols(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &basicvarpos2tableaurow, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &binvrow, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &binvarow, ncols) );

   for( i = 0; i < ncols; ++i )
      basicvarpos2tableaurow[i] = -1;
   SCIP_CALL( constructBasicVars2TableauRowMap(scip, basicvarpos2tableaurow) );


   /* entries of vars */
   nrayentries = 0;
   for( i = 0; i < probdata->numvars; ++i )
   {
      col = SCIPvarGetCol(probdata->bin_vars[i]);
      if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_BASIC )
      {
         SCIP_CALL( storeDenseTableauRow(scip, col, basicvarpos2tableaurow, nrayentries, raylength, binvrow, binvarow,
                  tableaurows) );
         rayentry2conspos[nrayentries] = i;
         nrayentries++;
      }
   }


   /* entry of aux var (if it exists) */
   if( auxvar != NULL )
   {
      col = SCIPvarGetCol(auxvar);
      if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_BASIC )
      {
         SCIP_CALL( storeDenseTableauRow(scip, col, basicvarpos2tableaurow, nrayentries, raylength, binvrow, binvarow,
                  tableaurows) );

         rayentry2conspos[nrayentries] = probdata->numvars;
         nrayentries++;
      }
   }
   assert(nrayentries == raylength);

#ifdef  DEBUG_INTERSECTIONCUT
   for( i = 0; i < ncols; ++i )
   {
      SCIPinfoMessage(scip, NULL, "%d column of tableau is: ",i);
      for( int j = 0; j < raylength; ++j )
         SCIPinfoMessage(scip, NULL, "%g ", tableaurows[i * raylength + j]);
      SCIPinfoMessage(scip, NULL, "\n");
   }
#endif

   SCIPfreeBufferArray(scip, &binvarow);
   SCIPfreeBufferArray(scip, &binvrow);
   SCIPfreeBufferArray(scip, &basicvarpos2tableaurow);

   return SCIP_OKAY;
}

/** initializes rays data structure */
static
SCIP_RETCODE createRays(
   SCIP*                 scip,               /**< SCIP data structure */
   RAYS**                rays                /**< rays data structure */
   )
{
   SCIP_CALL( SCIPallocBuffer(scip, rays) );
   BMSclearMemory(*rays);

   SCIP_CALL( SCIPallocBufferArray(scip, &(*rays)->rays, SCIPgetNLPCols(scip)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*rays)->raysidx, SCIPgetNLPCols(scip)) );

   /* overestimate raysbegin and lpposray */
   SCIP_CALL( SCIPallocBufferArray(scip, &(*rays)->raysbegin, SCIPgetNLPCols(scip) + SCIPgetNLPRows(scip) + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*rays)->lpposray, SCIPgetNLPCols(scip) + SCIPgetNLPRows(scip)) );
   (*rays)->raysbegin[0] = 0;

   (*rays)->rayssize = SCIPgetNLPCols(scip);

   return SCIP_OKAY;
}



/** frees rays data structure */
static
void freeRays(
   SCIP*                 scip,               /**< SCIP data structure */
   RAYS**                rays                /**< rays data structure */
   )
{
   if( *rays == NULL )
      return;

   SCIPfreeBufferArray(scip, &(*rays)->lpposray);
   SCIPfreeBufferArray(scip, &(*rays)->raysbegin);
   SCIPfreeBufferArray(scip, &(*rays)->raysidx);
   SCIPfreeBufferArray(scip, &(*rays)->rays);

   SCIPfreeBuffer(scip, rays);
}

/** inserts entry to rays data structure; checks and resizes if more space is needed */
static
SCIP_RETCODE insertRayEntry(
   SCIP*                 scip,               /**< SCIP data structure */
   RAYS*                 rays,               /**< rays data structure */
   SCIP_Real             coef,               /**< coefficient to insert */
   int                   coefidx,            /**< index of coefficient (conspos of var it corresponds to) */
   int                   coefpos             /**< where to insert coefficient */
   )
{
   /* check for size */
   if( rays->rayssize <= coefpos + 1 )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, coefpos + 1);
      SCIP_CALL( SCIPreallocBufferArray(scip, &(rays->rays), newsize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &(rays->raysidx), newsize) );
      rays->rayssize = newsize;
   }

   /* insert entry */
   rays->rays[coefpos] = coef;
   rays->raysidx[coefpos] = coefidx;

   return SCIP_OKAY;
}

/** constructs map between the lppos of a variables and its position in the constraint assuming the constraint variables
 * are sorted as [quad vars, lin vars, aux var (if it exists)]
 *
 * If a variable doesn't appear in the constraint, then its position is -1.
 */
static
void constructLPPos2ConsPosMap(
   ProbData *      probdata,
   SCIP_VAR*             auxvar,             /**< aux var of the expr */
   int*                  map                 /**< buffer to store the mapping */
   )
{

   int lppos;
   int i;

   for( i = 0; i < probdata->numvars; ++i )
   {
      lppos = SCIPcolGetLPPos(SCIPvarGetCol(probdata->bin_vars[i]));
      map[lppos] = i;
   }

   /* set pos of aux var (if it exists) */
   if( auxvar != NULL )
   {
      lppos = SCIPcolGetLPPos(SCIPvarGetCol(auxvar));
      map[lppos] = probdata->numvars;
   }

   return;
}

/** inserts entries of factor * nray-th column of densetableaucols into rays data structure */
static
SCIP_RETCODE insertRayEntries(
   SCIP*                 scip,               /**< SCIP data structure */
   RAYS*                 rays,               /**< rays data structure */
   SCIP_Real*            densetableaucols,   /**< column of the tableau in dense format */
   int*                  rayentry2conspos,   /**< map between rayentry and conspos of associated var */
   int                   raylength,          /**< length of a tableau column */
   int                   nray,               /**< which tableau column to insert */
   int                   conspos,            /**< conspos of ray's nonbasic var in the cons; -1 if not in the cons */
   SCIP_Real             factor,             /**< factor to multiply each tableau col */
   int*                  nnonz,              /**< position to start adding the ray in rays and buffer to store nnonz */
   SCIP_Bool*            success             /**< we can't separate if there is a nonzero ray with basis status ZERO */
   )
{
   int i;

   *success = TRUE;

   for( i = 0; i < raylength; ++i )
   {
      SCIP_Real coef;

      /* we have a nonzero ray with base stat zero -> can't generate cut */
      if( factor == 0.0 && ! SCIPisZero(scip, densetableaucols[nray * raylength + i]) )
      {
         *success = FALSE;
         return SCIP_OKAY;
      }

      coef = factor * densetableaucols[nray * raylength + i];

      /* this might be a source of numerical issues
       * TODO: in case of problems, an idea would be to scale the ray entries; compute the cut coef and scale it back down
       * another idea would be to check against a smaller epsilion.
       * The problem is that if the cut coefficient is 1/t where lpsol + t*ray intersects the S-free set.
       * Now if t is super big, then a super small coefficient would have had an impact...
       */
      if( SCIPisZero(scip, coef) )
         continue;

      /* check if nonbasic var entry should come before this one */
      if( conspos > -1 && conspos < rayentry2conspos[i] )
      {
         /* add nonbasic entry */
         assert(factor != 0.0);

#ifdef  DEBUG_INTERSECTIONCUT
         SCIPinfoMessage(scip, NULL, "ray belongs to nonbasic var pos %d in cons\n", conspos);
#endif

         SCIP_CALL( insertRayEntry(scip, rays, -factor, conspos, *nnonz) );
         (*nnonz)++;

         /* we are done with nonbasic entry */
         conspos = -1;
      }

      SCIP_CALL( insertRayEntry(scip, rays, coef, rayentry2conspos[i], *nnonz) );
      (*nnonz)++;
   }

   /* if nonbasic entry was not added and should still be added, then it should go at the end */
   if( conspos > -1 )
   {
      /* add nonbasic entry */
      assert(factor != 0.0);

#ifdef  DEBUG_INTERSECTIONCUT
      SCIPinfoMessage(scip, NULL, "ray belongs to nonbasic var pos %d in cons\n", conspos);
#endif

      SCIP_CALL( insertRayEntry(scip, rays, -factor, conspos, *nnonz) );
      (*nnonz)++;
   }

   /* finished ray entry; store its end */
   rays->raysbegin[rays->nrays + 1] = *nnonz;

#ifdef  DEBUG_INTERSECTIONCUT
   SCIPinfoMessage(scip, NULL, "entries of ray %d are between [%d, %d):\n", rays->nrays, rays->raysbegin[rays->nrays], *nnonz);
   for( i = rays->raysbegin[rays->nrays]; i < *nnonz; ++i )
      SCIPinfoMessage(scip, NULL, "(%d, %g), ", rays->raysidx[i], rays->rays[i]);
   SCIPinfoMessage(scip, NULL, "\n");
#endif

   return SCIP_OKAY;
}


static
SCIP_RETCODE createAndStoreSparseRays(
   SCIP*                 scip,               /**< SCIP data structure */
   ProbData *      probdata,
   SCIP_SepaData_Lattice * sepadata,
   SCIP_VAR*             auxvar,             /**< aux var of expr or NULL if not needed (e.g. separating real cons) */
   RAYS**                raysptr,            /**< buffer to store rays datastructure */
   SCIP_Bool*            success             /**< we can't separate if there is a var with basis status ZERO */
   )
{
   SCIP_Real* densetableaucols;
   SCIP_COL** cols;
   SCIP_ROW** rows;
   RAYS* rays;
   int* rayentry2conspos;
   int* lppos2conspos;
   int nnonbasic;
   int nrows;
   int ncols;
   int nnonz;
   int raylength;
   int i;

   nrows = SCIPgetNLPRows(scip);
   ncols = SCIPgetNLPCols(scip);

   *success = TRUE;

   raylength = countBasicVars(probdata, auxvar, success);
   if( ! *success )
   {
      SCIPdebugMsg(scip, "failed to store sparse rays: there is a var with base status zero\n");
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &densetableaucols, raylength * (ncols + nrows)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rayentry2conspos, raylength) );

   /* construct dense tableau and map between ray entries and position of corresponding var  */
   SCIP_CALL( storeDenseTableauRowsByColumns(scip, probdata, raylength, auxvar,
            densetableaucols, rayentry2conspos) );

   /* build rays sparsely now */
   SCIP_CALL( SCIPallocBufferArray(scip, &lppos2conspos, ncols) );
   for( i = 0; i < ncols; ++i )
      lppos2conspos[i] = -1;

   constructLPPos2ConsPosMap(probdata, auxvar, lppos2conspos);

   /* store sparse rays */
   SCIP_CALL( createRays(scip, raysptr) );
   rays = *raysptr;

   nnonz = 0;
   nnonbasic = 0;

   /* go through nonbasic variables */
   cols = SCIPgetLPCols(scip);
   for( i = 0; i < ncols; ++i )
   {
      int oldnnonz = nnonz;
      SCIP_COL* col;
      SCIP_Real factor;

      col = cols[i];

      /* set factor to store basic entries of ray as = [-BinvL, BinvU] */
      if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_LOWER )
         factor = -1.0;
      else if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_UPPER )
         factor = 1.0;
      else if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_ZERO )
         factor = 0.0;
      else
         continue;

      SCIP_CALL( insertRayEntries(scip, rays, densetableaucols, rayentry2conspos, raylength, nnonbasic,
               lppos2conspos[SCIPcolGetLPPos(col)], factor, &nnonz, success) );


      if( ! (*success) )
      {
         goto CLEANUP;
      }

      /* if ray is non zero remember who it belongs to */
      assert(oldnnonz <= nnonz);
      if( oldnnonz < nnonz )
      {
         rays->lpposray[rays->nrays] = SCIPcolGetLPPos(col);
         (rays->nrays)++;
      }
      nnonbasic++;
   }

   /* go through nonbasic slack variables */
   rows = SCIPgetLPRows(scip);
   for( i = 0; i < nrows; ++i )
   {
      int oldnnonz = nnonz;
      SCIP_ROW* row;
      SCIP_Real factor;

      row = rows[i];

      /* set factor to store basic entries of ray as = [-BinvL, BinvU]; basic status of rows are flipped! See lpi.h! */
      assert(SCIProwGetBasisStatus(row) != SCIP_BASESTAT_ZERO);
      if( SCIProwGetBasisStatus(row) == SCIP_BASESTAT_LOWER )
         factor = 1.0;
      else if( SCIProwGetBasisStatus(row) == SCIP_BASESTAT_UPPER )
         factor =-1.0;
      else
         continue;

      SCIP_CALL( insertRayEntries(scip, rays, densetableaucols, rayentry2conspos, raylength, nnonbasic, -1, factor,
               &nnonz, success) );
      assert(*success);

      /* if ray is non zero remember who it belongs to */
#ifdef  DEBUG_INTERSECTIONCUT
      SCIPinfoMessage(scip, NULL, "looked at ray of row %d, it has %d nonzeros\n-----------------\n", i, nnonz - oldnnonz);
#endif
      assert(oldnnonz <= nnonz);
      if( oldnnonz < nnonz )
      {
         rays->lpposray[rays->nrays] = -SCIProwGetLPPos(row) - 1;
         (rays->nrays)++;
      }
      nnonbasic++;
   }

CLEANUP:
   SCIPfreeBufferArray(scip, &lppos2conspos);
   SCIPfreeBufferArray(scip, &rayentry2conspos);
   SCIPfreeBufferArray(scip, &densetableaucols);

   if( ! *success )
   {
      freeRays(scip, &rays);
   }

   return SCIP_OKAY;
}



static
SCIP_Real computeIntersectionPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   ProbData *      probdata,
   SCIP_SepaData_Lattice * sepadata,
   SCIP_Real*            ray,        
   int*                  rayidx,           
   int                   raylen,
   SCIP_SOL*             sol                /**< solution we want to separate */    
   )
{
   int  i = sepadata->var_id;
   SCIP_Real vsol = SCIPgetSolVal(scip, sol, probdata->bin_vars[i]);
   int j = 0;
   bool isfind = false;
   for(j = 0; j < raylen; j++){
      if(rayidx[j] == i){
         isfind = true;
         break;
      }
   }

   SCIP_Real t = 0;
   if(isfind){
      SCIP_Real vray = ray[j];
      if(vray > 0){
         t = (1 - vsol) / vray;
      }
      else{
         t = -vsol / vray;
      }
   }
   else{
      t = SCIPinfinity(scip);
   }
   return t;
}

/** computes intersection cut cuts off sol (because solution sol violates the signomial expr)
 * and stores it in rowprep. Here, we don't use any strengthening.
 */
static
SCIP_RETCODE computeIntercut(
   SCIP*                 scip,               /**< SCIP data structure */
   ProbData *      probdata,
   SCIP_SepaData_Lattice * sepadata,
   RAYS*                 rays,               /**< rays */
   SCIP_ROWPREP*         rowprep,            /**< rowprep for the generated cut */
   SCIP_Real*            interpoints,        /**< array to store intersection points for all rays or NULL if nothing
                                                  needs to be stored */
   SCIP_SOL*             sol,                /**< solution we want to separate */
   SCIP_Bool*            success             /**< if a cut candidate could be computed */
   )
{
   SCIP_COL** cols;
   SCIP_ROW** rows;
   int i;

   cols = SCIPgetLPCols(scip);
   rows = SCIPgetLPRows(scip);

   /* for every ray: compute cut coefficient and add var associated to ray into cut */
   for( i = 0; i < rays->nrays; ++i )
   {
      SCIP_Real interpoint;
      SCIP_Real cutcoef;
      int lppos;

      SCIP_Real * ray = &rays->rays[rays->raysbegin[i]];
      int * raysidx = &rays->raysidx[rays->raysbegin[i]];
      int raylen =  rays->raysbegin[i + 1] - rays->raysbegin[i];

      /* compute intersection point */
      interpoint = computeIntersectionPoint(scip, probdata, sepadata, ray, raysidx, raylen, NULL);
#ifdef  DEBUG_INTERSECTIONCUT
      SCIPinfoMessage(scip, NULL, "interpoint for ray %d is %g\n", i, interpoint);
#endif

      /* store intersection point */
      if( interpoints != NULL )
         interpoints[i] = interpoint;

      /* compute cut coef */
      cutcoef = SCIPisInfinity(scip, interpoint) ? 0.0 : 1.0 / interpoint;

      /* add var to cut: if variable is nonbasic at upper we have to flip sign of cutcoef */
      lppos = rays->lpposray[i];
      if( lppos < 0 )
      {
         lppos = -lppos - 1;

         assert(SCIProwGetBasisStatus(rows[lppos]) == SCIP_BASESTAT_LOWER || SCIProwGetBasisStatus(rows[lppos]) ==
               SCIP_BASESTAT_UPPER);

         SCIP_CALL( addRowToCut(scip, rowprep, SCIProwGetBasisStatus(rows[lppos]) == SCIP_BASESTAT_UPPER ? cutcoef :
                  -cutcoef, rows[lppos], success) ); /* rows have flipper base status! */

         if( ! *success )
         {
            INTERLOG(printf("Bad numeric: now not nonbasic enough\n");)
            return SCIP_OKAY;
         }
      }
      else
      {
            assert(SCIPcolGetBasisStatus(cols[lppos]) == SCIP_BASESTAT_UPPER || SCIPcolGetBasisStatus(cols[lppos]) ==
                  SCIP_BASESTAT_LOWER);
            SCIP_CALL( addColToCut(scip, rowprep, sol, SCIPcolGetBasisStatus(cols[lppos]) == SCIP_BASESTAT_UPPER ? -cutcoef :
                  cutcoef, cols[lppos]) );
      }
   }

   return SCIP_OKAY;
}


/** generates intersection cut that cuts off sol (which violates the signomial expr) */
static
SCIP_RETCODE generateIntercut(
   SCIP*                 scip,               /**< SCIP data structure */
   ProbData *      probdata,
   SCIP_SepaData_Lattice * sepadata,
   SCIP_ROWPREP*         rowprep,            /**< rowprep for the generated cut */
   SCIP_Bool*            success             /**< whether separation was successfull or not */
   )
{
   RAYS* rays;
   SCIP_VAR* auxvar;
   SCIP_SOL* soltoseparate;


   *success = TRUE;

   /* in nonbasic space cut is >= 1 */
   assert(SCIProwprepGetSide(rowprep) == 0.0);
   SCIProwprepAddSide(rowprep, 1.0);
   SCIProwprepSetSidetype(rowprep, SCIP_SIDETYPE_LEFT);
   assert(SCIProwprepGetSide(rowprep) == 1.0);

   rays = NULL;
   auxvar = probdata->obj_var;

   /* check if we use tableau or bounds as rays */
   SCIP_CALL( createAndStoreSparseRays(scip, probdata, sepadata, auxvar, &rays, success) );

   if( ! *success )
   {
      return SCIP_OKAY;
   }

   soltoseparate = NULL;
   
   SCIP_CALL( computeIntercut(scip, probdata, sepadata, rays, rowprep, NULL, soltoseparate, success) );

   freeRays(scip, &rays);

   return SCIP_OKAY;
}

/** separate cuts for stored signomial */
static
SCIP_RETCODE separatePoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SepaData_Lattice * sepadata,
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{

   assert(sepa != NULL);
   assert(result != NULL);

	ProbData * probdata = NULL;
	probdata = dynamic_cast<ProbData *>(SCIPgetObjProbData(scip));
	assert(probdata != NULL);


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

   if(SCIPisFeasLT(scip, fractional, 0.005)){
      return SCIP_OKAY;
   }

   sepadata->var_id = id;
   SCIP_ROWPREP* rowprep;
   /* cut (in the nonbasic space) is of the form alpha^T x >= 1 */
   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, SCIP_SIDETYPE_LEFT, TRUE) );

   SCIP_Bool success = FALSE;
   SCIP_CALL( generateIntercut(scip, probdata, sepadata, rowprep,  &success) );
   SCIP_Real violation = 0;

   /* we generated something, let us see if it survives the clean up */
   if( success )
   {
      sepadata->ncutsgenerated += 1;
      sepadata->ncutsadded += 1;

      /* merge coefficients that belong to same variable */
      SCIPmergeRowprepTerms(scip, rowprep);

      SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, NULL, sepadata->mincutviolation, &violation, &success) );
   }

   /* if cut looks good (numerics ok and cutting off solution), then turn into row and add to sepastore */
   if( success )
   {
      SCIP_ROW* row;
      SCIP_Bool infeasible;

      SCIP_CALL( SCIPgetRowprepRowSepa(scip, &row, rowprep, sepa) );

      
      printf("## New cut\n");
      printf(" -> found maxquad-free cut <%s>: act=%f, lhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n\n",
            SCIProwGetName(row), SCIPgetRowLPActivity(scip, row), SCIProwGetLhs(row), SCIProwGetNorm(row),
            SCIPgetCutEfficacy(scip, NULL, row),
            SCIPgetRowMinCoef(scip, row), SCIPgetRowMaxCoef(scip, row),
            SCIPgetRowMaxCoef(scip, row)/SCIPgetRowMinCoef(scip, row)); 
      //SCIP_CALL(SCIPprintRow(scip, row, NULL));

      printf("SCIP DEPTH %d got a cut with violation %g, efficacy %g and r/e %g\n", SCIPgetSubscipDepth(scip),
        violation, SCIPgetCutEfficacy(scip, NULL, row), SCIPgetRowMaxCoef(scip, row) / SCIPgetRowMinCoef(scip, row) /
        SCIPgetCutEfficacy(scip, NULL, row));
      
      
      assert(SCIPgetCutEfficacy(scip, NULL, row) > 0.0);
      if( ! sepadata->ignorehighre || SCIPgetRowMaxCoef(scip, row) / SCIPgetRowMinCoef(scip, row) / SCIPgetCutEfficacy(scip, NULL, row) < 1e9 )
      {

         SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
      
         if( infeasible )
         {
            *result = SCIP_CUTOFF;
         }
         else
         {
            *result = SCIP_SEPARATED;
            sepadata->ncutsadded += 1;
         }
      }
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

   SCIPfreeRowprep(scip, &rowprep);


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
   SCIP_CALL( SCIPgetIntParam(scip, "separating/" SEPA_NAME "/ncutslimit", &sepadata.ncutslimit) );
   SCIP_CALL( SCIPgetIntParam(scip, "separating/" SEPA_NAME "/ncutslimitroot", &sepadata.ncutslimitroot) );
   SCIP_CALL( SCIPgetRealParam(scip, "separating/" SEPA_NAME "/mincutviolation", &sepadata.mincutviolation) );
   SCIP_CALL( SCIPgetRealParam(scip, "separating/" SEPA_NAME "/minviolation", &sepadata.mincutviolation) );
   SCIP_CALL( SCIPgetIntParam(scip, "separating/" SEPA_NAME "/atwhichnodes", &sepadata.atwhichnodes) );
   SCIP_CALL( SCIPgetBoolParam(scip, "separating/" SEPA_NAME "/ignorebadrayrestriction", &sepadata.ignorebadrayrestriction) );
   SCIP_CALL( SCIPgetBoolParam(scip, "separating/" SEPA_NAME "/ignorenhighre", &sepadata.ignorehighre) );
   
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

   *result = SCIP_DIDNOTRUN;
   SCIP_Longint nodenumber;

   /* only separate at selected nodes */
   SCIP_NODE* node = SCIPgetCurrentNode(scip);

   if( (sepadata.atwhichnodes == -1 && depth != 0) || (sepadata.atwhichnodes != -1 && depth % sepadata.atwhichnodes != 0) )
   {
      return SCIP_OKAY;
   }

   /* do not add more than ncutslimitroot cuts in root node and ncutslimit cuts in the non-root nodes */
   nodenumber = SCIPnodeGetNumber(node);
   if( sepadata.lastnodenumber != nodenumber )
   {
      sepadata.lastnodenumber = nodenumber;
      sepadata.lastncuts = sepadata.ncutsadded;
      sepadata.ncutsadded = 0;
   }

   /*else if( (depth > 0 && sepadata.ncutsadded - sepadata.lastncuts >= sepadata.ncutslimit) || (depth == 0 &&
            sepadata.ncutsadded - sepadata.lastncuts >= sepadata.ncutslimitroot)) */
   /* allow the addition of a certain number of cuts per signomials */
   long long nodelimit;
   SCIP_CALL(SCIPgetLongintParam(scip,  "limits/nodes", &nodelimit));
   SCIP_Bool rootmode = nodelimit == 1;
   if(!rootmode)
   {
      //ProbData * probdata = NULL;
	   //probdata = dynamic_cast<ProbData *>(SCIPgetObjProbData(scip));
	   //assert(probdata != NULL);
      //int ncutslimit = sepadata.ncutslimitroot / 128.0 * probdata->numvars;
      //ncutslimit = int(ncutslimit / (1<< (depth + 1)));
      int ncutslimit = sepadata.ncutslimit;
      if( (depth > 0 && sepadata.ncutsadded >= ncutslimit) || (depth == 0 &&
         sepadata.ncutsadded >= sepadata.ncutslimitroot) )
      {
         return SCIP_OKAY;
      }
   }

   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL || !SCIPisLPSolBasic(scip) )
   {
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   /* call separation method */
   SCIP_CALL( separatePoint(scip, &sepadata, sepa, result) );

   return SCIP_OKAY;
}


