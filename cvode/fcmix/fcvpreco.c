/**********************************************************************
 * File          : fcvpreco.c                                         *
 * Programmers   : Alan C. Hindmarsh and Radu Serban @ LLNL           *
 * Version of    : 30 March 2003                                      *
 *--------------------------------------------------------------------*
 * This C function CVPreco is to interface between the CVSPGMR module *
 * and the user-supplied preconditioner setup routine CVPRECO.        *
 * Note the use of the generic name FCV_PRECO below.                  *
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h" /* definitions of types realtype and integertype   */
#include "nvector.h"       /* definitions of type N_Vector and vector macros  */
#include "fcvode.h"        /* actual function names, prototypes, global vars. */
#include "cvspgmr.h"       /* CVSpgmr prototype                               */

/*********************************************************************/

/* Prototype of the Fortran routine */
void FCV_PRECO(realtype*, realtype*, realtype*, booleantype*, 
               booleantype*, realtype*, realtype*, realtype*, realtype*, 
               long int *, realtype*, realtype*, realtype*, int*);

/***************************************************************************/

/* C function CVPreco to interface between CVODE and a Fortran subroutine
   CVPRECO for setup of a Krylov preconditioner.
   Addresses of Nlocal, t, jok, gamma, h, uround, y, fy, ewt, vtemp1, vtemp2, 
   vtemp3, and the addresses jcurPtr and nfePtr, are passed to CVPRECO, using
   the routine N_VGetData from NVECTOR.  A return flag ier from CVPRECO
   is returned by CVPreco.
   Auxiliary data is assumed to be communicated by Common. */

int CVPreco(realtype t, N_Vector y, N_Vector fy, booleantype jok,
            booleantype *jcurPtr, realtype gamma, N_Vector ewt, realtype h,
            realtype uround, long int *nfePtr, void *P_data,
            N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  realtype *ydata, *fydata, *ewtdata, *v1data, *v2data, *v3data;
  int ier = 0;

  ydata   = N_VGetData(y);
  fydata  = N_VGetData(fy);
  ewtdata = N_VGetData(ewt);
  v1data  = N_VGetData(vtemp1);
  v2data  = N_VGetData(vtemp2);
  v3data  = N_VGetData(vtemp3);

  FCV_PRECO (&t, ydata, fydata, &jok, jcurPtr, &gamma, ewtdata,
           &h, &uround, nfePtr, v1data, v2data, v3data, &ier);

  /* Note: there is no need to use N_VSetData since we are not getting back any
    information that should go into an N_Vector */

  return(ier);
}
