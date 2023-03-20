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

/**@file   main.c
 * @brief  Main file for Doptimal solver
 * @author Liding Xu
 */

#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/scipdefplugins.h"


#include "sepa_intersub.h"
#include "sepa_interlattice.h"

/** creates a SCIP instance with default plugins, evaluates command line parameters, runs SCIP appropriately,
 *  and frees the SCIP instance
 */
static
SCIP_RETCODE runShell(
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array with shell parameters */
   const char*           defaultsetname      /**< name of default settings file */
   )
{
   SCIP* scip = NULL;

   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );
   
   /* we explicitly enable the use of a debug solution for this main SCIP instance */
   SCIPenableDebugSol(scip);

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );


   /* include submodular problem plugins */
   SCIP_CALL( SCIPincludeObjSepa(scip, new SepaInterSub(scip), TRUE) );
   SCIP_CALL( SCIPincludeObjSepa(scip, new SepaInterLattice(scip), TRUE) );

   /* parameter setting */


   SCIP_CALL( SCIPaddRealParam(scip,
      "separating/intersub/mincutviol",
      "minimum required violation of a cut",
      NULL, FALSE, 1e-4, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/intersub/maxrounds",
         "maximal number of separation rounds per node (-1: unlimited)",
         NULL, FALSE, 10, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/intersub/maxroundsroot",
         "maximal number of separation rounds in the root node (-1: unlimited)",
         NULL, FALSE, -1, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/intersub/disaggregate",
         "disaggregate the cut separation be  (TRUE: disaggregate)",
         NULL, FALSE, TRUE,  NULL, NULL) );


   SCIP_CALL( SCIPaddRealParam(scip,
      "separating/interlattice/mincutviol",
      "minimum required violation of a cut",
      NULL, FALSE, 1e-4, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/interlattice/maxrounds",
         "maximal number of separation rounds per node (-1: unlimited)",
         NULL, FALSE, 10, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/interlattice/maxroundsroot",
         "maximal number of separation rounds in the root node (-1: unlimited)",
         NULL, FALSE, -1, -1, INT_MAX, NULL, NULL) );


   /**********************************
    * Process command line arguments *
    **********************************/
   SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, defaultsetname) );

   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

int
main(
   int                        argc,
   char**                     argv
   )
{
   SCIP_RETCODE retcode;

   retcode = runShell(argc, argv, "scip.set");
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}

