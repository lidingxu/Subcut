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
 * @brief  Main file for submodular maxmization problem solver
 * @author Liding Xu
 */

#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/scipdefplugins.h"

#include "reader_sub.h"
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
   SCIP_CALL( SCIPincludeObjReader(scip, new ReaderSubmodular(scip), TRUE));
   SCIP_CALL( SCIPincludeObjSepa(scip, new SepaInterSub(scip), TRUE) );
   SCIP_CALL( SCIPincludeObjSepa(scip, new SepaInterLattice(scip), TRUE) );

   /* parameter setting */

   /* add parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "submodular/is_nature","use the submodular problem's nature formulation or its overestimating formulation", NULL, FALSE, TRUE, NULL,  NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "submodular/gradient_cut","use gradient cut for the submodular problem's nature formulation", NULL, FALSE, FALSE, NULL,  NULL) );
   
   string SEPA_NAME1= "interlattice";
   string SEPA_NAME2= "intersub";
   for(int i = 0; i < 2; i++){
      string SEPA_NAME;
      if(i == 0)
         SEPA_NAME = SEPA_NAME1;
      else
         SEPA_NAME = SEPA_NAME2;

      SCIP_CALL( SCIPaddIntParam(scip, ("separating/" + SEPA_NAME + "/ncutslimit").c_str(),
            "limit for number of cuts generated consecutively",
            NULL, FALSE, 2, 0, INT_MAX, NULL, NULL) );

      SCIP_CALL( SCIPaddIntParam(scip, ("separating/" + SEPA_NAME + "/ncutslimitroot").c_str(),
            "limit for number of cuts generated at root node",
            NULL, FALSE, 40, 0, INT_MAX, NULL, NULL) );

      SCIP_CALL( SCIPaddRealParam(scip, ("separating/" + SEPA_NAME + "/mincutviolation").c_str(),
            "minimal cut violation the generated cuts must fulfill to be added to the LP",
            NULL, FALSE, 1e-4, 0.0, SCIPinfinity(scip), NULL, NULL) );

      SCIP_CALL( SCIPaddRealParam(scip, ("separating/" + SEPA_NAME  + "/minviolation").c_str(),
            "minimal violation the constraint must fulfill such that a cut is generated",
            NULL, FALSE, 1e-4, 0.0, SCIPinfinity(scip), NULL, NULL) );

      SCIP_CALL( SCIPaddIntParam(scip, ("separating/" + SEPA_NAME + "/atwhichnodes").c_str(),
            "determines at which nodes cut is used (if it's -1, it's used only at the root node, if it's n >= 0, it's used at every multiple of n",
            NULL, FALSE, 1, -1, INT_MAX, NULL, NULL) );

      SCIP_CALL( SCIPaddBoolParam(scip, ("separating/" + SEPA_NAME + "/ignorebadrayrestriction").c_str(),
            "should cut be generated even with bad numerics when restricting to ray?",
            NULL, FALSE, FALSE, NULL, NULL) );

      SCIP_CALL( SCIPaddBoolParam(scip, ("separating/" + SEPA_NAME + "/ignorenhighre").c_str(),
            "should cut be added even when range / efficacy is large?",
            NULL, FALSE, FALSE, NULL, NULL) );
   }


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

