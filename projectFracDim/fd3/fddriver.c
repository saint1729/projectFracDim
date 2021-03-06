/* BEGIN NOTICE

Copyright (c) 1992 by John Sarraille and Peter DiFalco
(john@ishi.csustan.edu)

Permission to use, copy, modify, and distribute this software
and its documentation for any purpose and without fee is hereby
granted, provided that the above copyright notice appear in all
copies and that both that copyright notice and this permission
notice appear in supporting documentation.

The algorithm used in this program was inspired by the paper
entitled "A Fast Algorithm To Determine Fractal Dimensions By
Box Counting", which was written by Liebovitch and Toth, and
which appeared in the journal "Physics Letters A", volume 141,
pp 386-390, (1989).

This program is not warranteed: use at your own risk.

END NOTICE */

#include <stdio.h>
#include "fd.h"
#include <math.h>
#include <stdlib.h>


main(argc,argv)
        int     argc;
        char    *argv[];
{
  extern FILE  *infile;
  extern int   embed_dim,         /* How many coordinates points have. */
               *marker;           /*  To mark queues during radix sort */
  extern ulong **data,    /* initial array of data points */
               dataLines,        /* number of data points in input file */
              *diff_test;        /* XOR's of pairs of coordinates */
  extern long int  *next_after,    /*  array for pointers to data. */
                   avail;             /*  avail list pointer */
  extern QHeader Q[2];       /* array of two queues for the radix sort */


  int    i,j,m;
  ulong  maxdiam=0,              /* 2^numbits - 1: determines the scaling. */
         boxCount[numbits+1],  /* Array to keep track of the box
                                         counts of each size.  Each index
                                         is the base-2 log of the
                                         corresponding box size. */
         pointCount[numbits+1],  /*  Array to keep track of how many points
                                     of the data set are contained in each
                                     box */
	 mark;                  /* Marks smallest usable box count */
  double    max, min,   /* Maximum and minimum values in input file-- we
                           need to know these in order to scale the
                           input points. */
            buf,        /* A buffer for inputs to stay in until they are
                           scaled and converted to integers.  */
	    negLogBoxCount[numbits+1],
            logSumSqrFreq[numbits+1],
            information[numbits+1],
	    capDim, infDim, corrDim;
  long int  pointer;    /* holds indexes used as pointers to records */

printf("\n\n") ;
printf("******************************************************************\n");
printf("  FRACTAL DIMENSION REPORT -- by fd software (DiFalco/Sarraille)\n");
printf("******************************************************************\n");
printf("\n") ;

     /*  IDENTIFY THE INPUT FILE */
  printf ("\nReporting on file named: %s.\n\n", argv[1]) ;

     /* FIND OUT HOW MANY COORDINATES POINTS HAVE -- 1?, 2?, 3?, MORE? */
  printf("Getting the embedding dimension ...\n");
  embed_dim = get_e_dim(argv[1]);
  printf("Embedding Dimension taken to be: %d\n\n",embed_dim);

/*
   find max and min in input file.  This function should open the input file
   for reading, scan it to get the maximum and minimum data values it
   contains, and then seek to the beginning of the input file so that the
   code below can proceed to read the data sequentially from the beginning.
   We are having trouble getting that to work, so we have temporarily taken
   the expedient of closing the file in the function max_min, and opening it
   again below.
*/
  printf("Finding max and min values in data ...");
  fflush(stdout);
  max_min (argv[1], &max, &min) ;
  printf ("\n\nMinimum value in input file is: %lf ...\n", min) ;
  printf ("Maximum value in input file is: %lf ...\n\n", max) ;
     /*
        set maxdiam=2^numbits-1 -- this is the largest value that an
        unsigned integer can have with the host compiler.
        (I checked the value obtained in this manner, and it
         is correct -- j.s.)
     */
  if (debugging || checking) printf("numbits is: %d\n",numbits);
  printf("%d different cell sizes will be used ...\n",numbits);
  for (m=0;m<numbits;m++)       maxdiam = maxdiam + (1<<m);
  printf("Data will be shifted and re-scaled so that minimum coordinate\n");
  printf("value is ZERO and maximum coordinate value is: %lu ...\n\n",maxdiam);

     /* open input file */
  if((infile = fopen(argv[1],"r")) == NULL)
    {
      printf("Cannot open input file, %s\n",argv[1]);
      exit(-1);
    }

     /* get number of data points from input file */
  fscanf(infile,"%lu",&dataLines);
  if (debugging || checking) printf("dataLines is: %ld\n",dataLines);

/*
    In the next few statements, we allocate "parallel arrays" which are
    being used to implement records containing fields for as many
    coordinates as the embedding dimension requires, for a "next_after"
    pointer field to be used to link the records into lists, and for a
    "marker" field to indicate whether the element is being used to mark the
    end of a pass in the radix sort.
*/
/*  allocate memory to point to a coordinate array for each coordinate --
    embedding dimension determines how many coordinates.
*/

  printf("Allocating storage ...");
  fflush(stdout);

  if ( (data = (ulong **) malloc (embed_dim*sizeof(ulong *))) == NULL)
    {
      printf("Memory allocation failure: \"data\".\n");
      exit(-1);
    }

/*  Here we allocate the storage for the coordinates of data points, plus
    extra room for 2 queue markers needed in the radix sort.
*/
  for (i=0;i<embed_dim;i++)
     if ((data[i] = (ulong *) malloc ((dataLines+2)*sizul)) == NULL)
        {
          printf("Memory allocation failure: \"data[%d]\").\n",i);
          exit(-1);
        }

/*  Here we allocate the storage for the pointers.  These are implemented as
    long integers.  The values of these pointers will be used as indexes into
    the arrays being declared here.  Note that here we also allocate room for
    two markers.
*/
  if((next_after=(long *)malloc((dataLines+2)*sizeof(long))) == NULL)
     {
       printf("Memory allocation failure: \"next_after\" array).\n");
       exit(-1);
     }

/* Here we allocate storage for a tag field that tells whether a record is a
   marker.
*/

   if ( (marker = (int *) malloc( (dataLines+2)*sizeof(int) ) ) == NULL)
       {
         printf("Memory allocation failure: \"marker\" array).\n");
         exit(-1);
       }

/* Here we allocate storage for a word that describes the bit changes
   between two different unsigned long int's.  We have one of these words
   for each of the embedding dimensions.  These will come in handy when the
   sweep is done that gets the box counts.
*/

   if ( (diff_test = (ulong *) malloc ((embed_dim)*sizeof(ulong)) ) == NULL)
       {
         printf("Memory allocation failure: \"diff_test\" array).\n");
         exit(-1);
       }

   printf(" Done ...\n");

/* Initialize the queues to an empty state.  We load the data directly into
   these queues, and use them afterwards to do the radix sort and the sweep.
*/
   printf("Initializing queues ... ");
   fflush(stdout);

   create_Q(&Q[0]); create_Q(&Q[1]);
   if (debugging || checking) printf ("Just created queues.\n");
   if (debugging) print_Q(&Q[0]);
   if (debugging) print_Q(&Q[1]);

/* Initialize the avail list to hold all the records implemented by the
   parallel arrays.
*/
   create_Av(&avail)   ;
   if (debugging || checking) printf ("Just created avail.\n");
   if (debugging) print_Av(avail);

/* Place a marker in each queue. Markers will be in front of each queue
   after all the data has been loaded.
*/
   for (i=0;i<2;i++)
      {
        pop_Av(&avail, &pointer);
        marker[pointer] = 1;
        en_Q(&(Q[i]),pointer);
      }

   printf("Done ...\n");

/* Place all the data points in the appropriate queue. */

   printf ("Loading data ... ");
   fflush(stdout);

   for(i=0;i<dataLines;i++)
      {
        pop_Av(&avail, &pointer);   /* get next record */
        marker[pointer] = 0;        /* it is not a marker */
        for (j=0; j<embed_dim;j++)
          { fscanf(infile,"%lf",&buf);  /* put scaled data in array */
            if (debugging) printf ("Just read %lf\n",buf);
            data[j][pointer]=(ulong)((buf-min)*(maxdiam/(max-min)));
          }
           /* Get started on the radix sort by putting the data in the
              queue corresponding to the least significant bit of the
              last coordinate. */
        if (IS_A_ONE(data[embed_dim-1][pointer],0))
        en_Q(&(Q[1]),pointer);       /* last coord ends in "1" */
        else en_Q(&(Q[0]),pointer);   /* last coord ends in "0" */
      }

  /* close the input file */
  fclose(infile);
  printf ("Done ...\n");
  if (debugging) print_Av(avail);
  if (debugging) print_Q(&Q[0]);
  if (debugging) print_Q(&Q[1]);

     /* radix sort queues */

  printf("Sorting the data ... ");
  fflush(stdout);

  radixsort();

  printf("Done ... \n");

    /* sweep data */

  printf("Doing sweep to get counts ... ");
  fflush(stdout);

  sweep(boxCount, negLogBoxCount, logSumSqrFreq, information);

  printf("Done ...\n\n");

  findMark(&mark, boxCount) ;

  /* GET RID OF THIS LINE WHEN DONE WITH TEST!!! */
  /* mark = 24; */

     /* WRITE COUNTS AND DERIVED DATA TO STANDARD OUTPUT. */
  printf("\n");
  printf("[log(epsl)]"); printf("  [CellCount]");
  printf("  [log(CellCount)]");  printf("  [informtn]");
  printf("  [-log(SumSqrFreqs)]\n");

  printf(" = [log(e)]"); printf("   = [N(e)] ");
  printf("     = [logN(e)]  ");  printf("    = [I(e)] ");
  printf("    = [-logSSF(e)]\n\n");

  for(i=0;i<=numbits;i++)
  {
   if ( (mark < numbits-2) && (i == mark) )
      { printf ("**************************************");
        printf ("**************************************\n");	}
    printf("%9d", i );
    printf("%13lu", boxCount[i]);
    printf("%18.5f", -negLogBoxCount[i]);
    printf("%12.5f", -information[i]);
    printf("%21.5f\n", -logSumSqrFreq[i]);
    if ( (mark < numbits-2) && (i == numbits-2) )
      { printf ("**************************************");
        printf ("**************************************\n");	}
  }

  printf ("\n\n\nTwo-Point Estimates of FD's:\n\n");
  printf("[log(e)]");      printf("  [logN(e)-logN(2e)]");
  printf("  [I(e)-I(2e)]");  printf("  [logSSF(2e)-logSSF(e)]\n\n");

  for(i=0; i<=numbits-1; i++)
    {
      if ( (mark < numbits-2) && (i == mark) )
        { printf ("**************************************");
          printf ("**************************************\n");	}
      printf ("%8d",i);
      printf ("%20.5f",negLogBoxCount[i+1]-negLogBoxCount[i]);
      printf ("%14.5f",information[i+1]-information[i]);
      printf ("%24.5f\n",logSumSqrFreq[i+1]-logSumSqrFreq[i]);
      if ( (mark < numbits-2) && (i==numbits-3) )
        { printf ("**************************************");
          printf ("**************************************\n");	}
    }

  if (mark >= numbits-2)
  {
    printf("\nInsufficient Data.  More points are needed.\n") ;
    printf("Cannot assign a Fractal Dimension to this set.\n");
    printf("Examine the data above to make your own conjecture.\n");
    exit(-1) ;
  }

  printf ("\n\n\n%ludd is the smallest cell size used in ", mark);
  printf ("the overall dimension estimates\n");
  printf ("below.  The largest cell size is ");
  printf ("%d.  Data above corresponding to\n", numbits-2);
  printf ("this range is between rows of asterisks.\n\n");

  GetDims(negLogBoxCount, logSumSqrFreq, information, mark,
                &capDim, &infDim, &corrDim) ;

  printf("\n\nLeast-Square Estimates based on Indicated Cell Range:\n\n");
  printf("Fractal Dimension  (Capacity)   =  %.5f\n", capDim );
  printf("Fractal Dimension (Information) =  %.5f\n", infDim );
  printf("Fractal Dimension (Correlation) =  %.5f\n", corrDim);
  printf("\n\n****************************************\n");
  
  FILE *fp2;
  fp2 = fopen("/home/tsn/Desktop/cse_project/fd3/results/results", "a");
  if (fp2!=NULL)
    fprintf(fp2, "%s\t\t\t\t\t%f\t\t\t\t\t\t%f\t\t\t\t\t\t\t%f\n" ,argv[1], capDim, infDim, corrDim);
}
