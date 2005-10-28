/*****************************************************************************\
 *  Function: ll_fetch
 *
 *  Description: This function will return the data value from the specified 
 *  data object. The field to be retrieved within the object is referenced by
 *  the dataField.
 *
 *  Arguments:
 *    IN perfObject : Pointer to the object to be read.
 *    IN dataField : Enum which references the data field within the
 *                  object to be retrieved.
 *    OUT data : Pointer to the data to be retrieved.
 *    RET Success: 0
 *        Failure: -1 : Invalid perfObject.
 *                 -2 : Version check.
 *                 -3 : Invalid dataField.
 *****************************************************************************
 *  Copyright (C) 2004 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory.
 *  Written by Morris Jette <jette1@llnl.gov>
 * 
 *  This file is part of slurm_ll_api, a collection of LoadLeveler-compatable
 *  interfaces to Simple Linux Utility for Resource Managment (SLURM).  These 
 *  interfaces are used by POE (IBM's Parallel Operating Environment) to 
 *  initiated SLURM jobs. For details, see <http://www.llnl.gov/linux/slurm/>.
 *
 *  This notice is required to be provided under our contract with the U.S.
 *  Department of Energy (DOE).  This work was produced at the University
 *  of California, Lawrence Livermore National Laboratory under Contract
 *  No. W-7405-ENG-48 with the DOE.
 * 
 *  Neither the United States Government nor the University of California
 *  nor any of their employees, makes any warranty, express or implied, or
 *  assumes any liability or responsibility for the accuracy, completeness,
 *  or usefulness of any information, apparatus, product, or process
 *  disclosed, or represents that its use would not infringe
 *  privately-owned rights.
 *
 *  Also, reference herein to any specific commercial products, process, or
 *  services by trade name, trademark, manufacturer or otherwise does not
 *  necessarily constitute or imply its endorsement, recommendation, or
 *  favoring by the United States Government or the University of
 *  California.  The views and opinions of authors expressed herein do not
 *  necessarily state or reflect those of the United States Government or
 *  the University of California, and shall not be used for advertising or
 *  product endorsement purposes.
 * 
 *  The precise terms and conditions for copying, distribution and
 *  modification are specified in the file "COPYING".
\*****************************************************************************/

#include "common.h"
#include "config.h"
#include "llapi.h"

extern int ll_fetch(LL_element *perfObject,enum LLAPI_Specification dataField, 
		void *data)
{
	int rc = -2;

	VERBOSE("++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	VERBOSE("ll_fetch\n");
	VERBOSE("Returns: %d\n", rc);
	VERBOSE("--------------------------------------------------\n");
	return rc;
}
