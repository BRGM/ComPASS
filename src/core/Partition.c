#include <stdio.h>
#include "metis.h"

void Metis_C( int nMaille, int *ptMaillebyMaille, int *numMaillebyMaille, int npart, int objval, int *Part )
{
  int i, pt;
  int Metis_Message;
  // idx_t options[METIS_NOPTIONS];
  idx_t nvtxs =  nMaille;
  idx_t ncon  = 1;
  idx_t nparts = npart;

  //Recursive, partition par defaut

  //diagnostic
  if ( 0 )
  {
    printf( "appel metis (c) avec\n" );
    for ( i = 0; i < nMaille; ++i )
    {
      printf( "maille = %d\n", i );
      for ( pt = ptMaillebyMaille[i]; pt < ptMaillebyMaille[i + 1]; ++pt )
        printf( "maille voisine %d \n", numMaillebyMaille[pt] );
      printf( "\n" );
    }
  }

  Metis_Message = METIS_PartGraphRecursive( &nvtxs, &ncon, ptMaillebyMaille, numMaillebyMaille, \
                                            NULL, NULL, NULL, &nparts, NULL, \
                                            NULL, NULL, &objval, Part );

  //Kway, autre partition possible
  if ( 0 )
  {
    Metis_Message = METIS_PartGraphKway( &nvtxs, &ncon, ptMaillebyMaille, numMaillebyMaille, \
                                         NULL, NULL, NULL, &nparts, NULL, \
                                         NULL, NULL, &objval, Part );
  }



  //diagnostic
  if ( 0 )
  {
    printf( "NbCell = %d\n", nMaille );
    printf( "npart  = %d\n", npart );
    printf( "%d\n", numMaillebyMaille[0] );
    printf( "Part:\n" );
    for ( i = 0; i < nMaille; ++i )
      printf( "%d ", Part[i] );
    printf( "\n" );
    printf( "objval = %d\n", objval );
  }
}
