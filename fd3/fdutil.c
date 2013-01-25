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

int get_e_dim();
void max_min();
void rad_pass();
void radixsort();
void sweep();
float fracdim();

/* ################################################################## */
/* ################################################################## */
int get_e_dim(userinputfile)
        char *userinputfile;
{
   extern FILE  *infile;
   int c, count ;   double temp ;
   if (debugging || checking) printf("Now inside get_e_dim.\n");
   if (debugging || checking) printf("About to open input file.\n");
   if((infile = fopen(userinputfile,"r")) == NULL)  /* open input file */
     {
       printf("Cannot open input file, %s\n",userinputfile);
       exit(-1);
     }
      /* get past "dataLines", and the first number on the first line of
         actual data. */
   fscanf (infile, "%lf%lf", &temp, &temp);
   count = 1 ;        /* One number on this line, so far. */
   do      /*  Read and count the rest of the numbers on this line. */
    {
      do     c=getc(infile) ;          /* Skip over same-line white space */
      while  ( (c=='\t') || (c==' ')  ) ;
      if   ( (c != '\n') && (c != '\r') && (c !='\f') )
       {                       /* if we don't seem to be on a new line */
        count = count + 1 ;    /* then increment # of numbers on this line */
        do   c=getc(infile) ;  /* and read past that next number */
        while  ( (c != '\n') && (c != '\r') && (c !='\f')
                 && (c !='\t') && (c != ' ') );
       }
    }  while ( (c != '\n') && (c != '\r') && (c !='\f') ) ;
   fclose (infile);
   return (count);
}
/* ################################################################## */
/* ################################################################## */
void max_min(userinputfile,pmax, pmin)
        char *userinputfile;
        double *pmax, *pmin ;
{
        extern int   embed_dim;      /* How many coordinates points have. */
        FILE    *infile;        /* input file */
        double  temp;
        ulong   i, numToRead;           /* control variable */
        ulong   dataLines;      /* number of data points in input file */
        /* open input file */
        if((infile = fopen(userinputfile,"r")) == NULL)
        {
                printf("Cannot open input file, %s\n",userinputfile);
                exit(-1);
        }
        /* get dataLines */
        if ((fscanf(infile,"%lu",&dataLines)) == 0)
          {
             printf ("Format of input file, %s, ", userinputfile);
             printf ("is incorrect -- please check.\n");
             exit(-1);
          }
        if (debugging || checking)
	  { printf ("Now control is in max_min.\n");
            printf ("Number of data lines is %lu ...\n", dataLines);
            printf ("Starting to read data.\n");
	  }
	      /* find maximum and minimum data values in the input file.
		 These are needed so that the data can be scaled to be a set
		 of non-negative integers between 0 and maxdiam, where
		 maxdiam will be the largest integer value expressible as an
		 element of the data type "unsigned long int". */

        fscanf(infile,"%lf",&temp);
        if (debugging || checking) printf ("First datum is %lf.\n", temp);
        *pmin=*pmax=temp;
        if (debugging)
	  printf ("Pmin is now: %lf, Pmax is now: %lf\n", *pmin, *pmax);
        numToRead = ((ulong) (embed_dim)) * dataLines;
        if (debugging)
	  printf ("Total number of coordinates to read is: %lu\n", numToRead);
        for(i=1;i<numToRead;i++)
        {
                fscanf(infile,"%lf",&temp);
                if (debugging) printf ("Next datum read is: %lf\n", temp);
                if(temp > *pmax)
                        *pmax = temp;
                else
                        if(temp < *pmin)
                                *pmin = temp;
 if (debugging) printf ("Pmin is now: %lf, Pmax is now: %lf\n", *pmin, *pmax);
        }
        /* check to see if the maximum equals the minimum -- this is a
           degenerate case -- or it might mean that the input file is
}          faulty. */
        if(*pmax == *pmin)
        {
          printf ("The input file, %s, is confusing!\n\n", userinputfile);
          printf ("Either all the points in %s have the same ", userinputfile);
          printf ("coordinates,\n");
          printf ("(THE FRACTAL DIMENSION IS ZERO IN THIS CASE.)\n\n");
          printf ("or %s is simply of the wrong form for an\n",userinputfile);
          printf ("input file to this program -- please check.\n");
          exit(-1);
        }
        /*  close the input file */
        fclose(infile);
}
/* ################################################################## */
/* ################################################################## */
/*
   This procedure performs one pass of the radix sort.  It assumes that the
   queues each have a marker in front at the time of the call, and this is
   the condition the procedure LEAVES the queues in when it terminates.
*/
void rad_pass(coord, bitpos)
        int coord,   /*  The coordinate we are doing this pass on */
            bitpos;  /*  The bit position we are doing this pass on */
{  extern int  *marker;           /*  To mark queues during radix sort */
   extern QHeader Q[2];       /* array of two queues for the radix sort */
   extern ulong **data;   /* initial array of data points */
   int queue_num, index ;
    /*  Move the markers to the rear of the queues */
 if (debugging)
    printf ("Starting pass on bit %d of coord %d.\n",bitpos,coord);
 if (debugging) printf ("Moving markers to the rear of queues.\n");
 for (queue_num=0;queue_num<2;queue_num++)
   {  de_Q(&(Q[queue_num]),&index); en_Q(&(Q[queue_num]),index);  }
 if (debugging) print_Q(&Q[0]);
 if (debugging) print_Q(&Q[1]);

    /*  Move all non-marker elements to the appropriate queue. */
 if (debugging) printf ("Starting main part of pass.\n");
 for (queue_num=0;queue_num<2;queue_num++)
   {    /* Peek first to see if Qfront is a marker -- this violates
           information hiding, and would be "cleaner" if done with a
           procedure in the queue package.  */
     while ( marker[Q[queue_num].Qfront] == 0 )
       {     /*  Peeking again! */
         if (IS_A_ONE(data[coord][Q[queue_num].Qfront],bitpos))
               /* Directly transfer from the source to target queue */
          transf_Q( &(Q[queue_num]),&(Q[1]) );
         else transf_Q( &(Q[queue_num]),&(Q[0]) ) ;
        if (debugging) printf("A queue transfer is done.  \n");
        if (debugging) print_Q(&Q[0]);
        if (debugging) print_Q(&Q[1]);
       }  }  }

/* ################################################################## */
/* ################################################################## */
/*
THIS SORT IS TO BE USED DIRECTLY AFTER FDDRIVER DOES THE INITIAL LOADING OF
THE DATA INTO THE QUEUES.  IT LEAVES THE DATA ESSENTIALLY SORTED, WITH ALL
THE DATA WHOSE X-COORD STARTS WITH 0 IN Q[0], IN ORDER, AND ALL THE DATA
WHOSE X-COORD STARTS WITH 1 IN Q[1], IN ORDER.  THUS THE SWEEP THAT COMES
NEXT MUST TRAVERSE Q[0], AND THEN Q	[1].
*/
void radixsort()
{
  extern QHeader Q[2];       /* ARRAY OF TWO QUEUES FOR THE RADIX SORT */
  extern int   embed_dim;      /* HOW MANY COORDINATES POINTS HAVE. */
  extern long int  avail;             /*  AVAIL LIST POINTER */
  int  bitpos, coord, queue_num;
  long int index;

    /* FINISH UP ON THE ZERO-BIT -- LOADING DATA TOOK CARE OF FIRST PASS. */
  for (coord=embed_dim-2;coord>=0;coord--) rad_pass(coord,0);

      /* NOW SORT ON THE REST OF THE BITS. */
  for (bitpos=1;bitpos<numbits;bitpos++)
      for (coord=embed_dim-1;coord>=0;coord--) rad_pass(coord,bitpos);

      /*  LASTLY, GET THE MARKERS OUT OF THE QUEUES. */
  for (queue_num=0;queue_num<2;queue_num++)
    {
      de_Q(&(Q[queue_num]),&index);
      push_Av(&avail,index);
    }}
/* ################################################################## */
/* ################################################################## */
/*
THIS PROCEDURE TRAVERSES THE SORTED DATA POINT LIST, AND EXTRACTS THE
INFORMATION NEEDED TO COMPUTE THE CAPACITY, INFORMATION, AND CORRELATION
DIMENSIONS OF THE DATA.  DATA IS CORRECTED ON THE FLY DURING THE TRAVERSAL,
AND THEN "MASSAGED".  AFTER "SWEEP" RUNS, WE NEED ONLY FIT A SLOPE TO THE
DATA TO OBTAIN THE FRACTAL DIMENSIONS.
*/
void sweep(boxCountS, negLogBoxCountS, logSumSqrFreqS, informationS)

   ulong   boxCountS[numbits+1];  /* COUNT BOXES OF EACH SIZE */
   
   double
               /* FOR EACH BOX SIZE #s, negLogBoxCountS[s] WILL EVENTUALLY BE
	          SET TO THE NEGATIVE OF THE LOG (BASE TWO) OF THE NUMBER OF
		  BOXES OF SIZE #s THAT ARE OCCUPIED BY ONE OR MORE DATA
		  POINTS.
		*/

           negLogBoxCountS[numbits+1],
   
               /* FOR EACH BOX SIZE #s, logSumSqrFreqS[s] WILL EVENTUALLY BE
		  SET TO THE LOG (BASE TWO) OF THE PROBABILITY THAT TWO
		  RANDOMLY CHOSEN DATA POINTS ARE IN THE SAME BOX OF SIZE
		  #s.
	       */

           logSumSqrFreqS[numbits+1],

              /* FOR EACH BOX SIZE #s, informationS[s] WILL EVENTUALLY BE
		 SET TO THE INFORMATION (BASE TWO) IN THE DISTRIBUTION OF
		 DATA POINTS IN THE BOXES OF SIZE #s.
	      */
           informationS[numbits+1] ;
	   
{ extern int  embed_dim;         /* HOW MANY COORDINATES POINTS HAVE. */
  extern ulong **data,           /* INITIAL ARRAY OF DATA POINTS */
               dataLines ;
  extern QHeader Q[2];           /* ARRAY OF TWO QUEUES FOR THE RADIX SORT */
  extern ulong  *diff_test;      /* XOR'S OF PAIRS OF COORDINATES */
  extern long int  *next_after;  /*  ARRAY FOR POINTERS TO DATA. */

  int       bitpos, countpos, coord, queue_num, found;
  ulong     pointCount[numbits+1] ;
  long int  current, previous;
  double    sumSqrFreq[numbits+1], freq, log2=log(2.0) ;

    /* GET A POINTER TO THE FIRST DATA VALUE */
  if    (Q[0].Qfront != -1) previous = Q[0].Qfront;
  else  previous = Q[1].Qfront;

     /* INIT boxCountS, pointCountS, sumSqrFreq, AND informationS.*/
  for(countpos=0;countpos<=numbits;countpos++)
    {  boxCountS [countpos]=1;   sumSqrFreq  [countpos]=0.0;
       pointCount[countpos]=1;   informationS[countpos]=0.0; }

  for (queue_num=0;queue_num<2;queue_num++)
  { current = Q[queue_num].Qfront;
    while (current != -1)
    { found = 0 ;
        /* START BY LOOKING AT THE BIGGEST BOX SIZE */
      bitpos=numbits-1;
      for (coord=embed_dim-1;coord>=0;coord--)
        diff_test[coord] = data[coord][previous]^ data[coord][current];
      do {coord = embed_dim - 1;
          do {/* IF THE CURRENT POINT AND PREVIOUS POINTS ARE IN DIFFERENT
	         BOXES OF THIS SIZE, */
	       if ( IS_A_ONE(diff_test[coord],bitpos) )
               {/* THEN THE CURRENT POINT IS IN NEW BOXES OF ALL SMALLER
		   SIZES TOO, AND STILL IN THE SAME BOXES OF LARGER SIZES,
		   SO ... */
	        for (countpos=bitpos;countpos>=0;countpos--)
		  {/* CALCULATE FREQUENCY OF POINTS IN THE BOX, ASSUMING FOR
		      NOW THAT THE NUMBER OF DATA LINES IN THE INPUT FILE IS
		      THE NUMBER OF DISTINCT POINTS IN THE DATA SET.  WE
		      ADJUST THIS AT THE END OF THIS FUNCTION. */
                   if (debugging)
		    printf("pointCount[%d] is %lud...\n", countpos, pointCount[countpos]);
		   freq = pointCount[countpos]/(double)dataLines ;
		     /* WE WILL ENCOUNTER NO MORE OF THE POINTS IN THE BOX
			WE JUST LEFT (THE SPECIAL ORDERING OF THE SORT WE
			USED ABOVE GUARANTEES THIS!), SO WE COMPUTE WHAT
			THIS BOX CONTRIBUTES TO THE RUNNING SUMS. */
                   sumSqrFreq[countpos] += (freq * freq) ;
                   informationS[countpos] += ( freq*log(freq)/log2 ) ;
 		     /* WE HAVE GOTTEN INTO A NEW BOX AT THIS LEVEL, SO WE
			REFLECT THE NEW BOX IN THE COUNT */
                   boxCountS[countpos]++ ;
                     /* SINCE WE HAVE A NEW BOX AT THIS LEVEL, THERE IS ONLY
		        ONE KNOWN POINT IN IT SO FAR -- THE CURRENT POINT */
		   pointCount[countpos]=1;
                  }
                for (countpos=bitpos+1;countpos<=numbits;countpos++)
		    /* THE CURRENT POINT IS IN THE BOXES AT THESE LEVELS, SO
		       JUST INCREMENT THE POINT COUNTER.  */
                  pointCount[countpos]++ ;
                found = 1;
               }
               else coord-- ;
             }
          while ( (found == 0) && (coord > -1) );
          bitpos-- ;
         }
      while ( (found == 0) && (bitpos > -1) );
      previous = current; current = next_after[current];
    } }
    
    /* NOW ADD IN THE CONTRIBUTION DUE TO THE COUNTS REMAINING AFTER THE
       LAST POINT HAS BEEN FOUND, RENORMALIZE WITH BOXCOUNT[0], AND MASSAGE
       THE RAW DATA FROM THE TRAVERSAL SO THAT IS IS READY FOR THE LEAST
       SQUARES FITTING. */
       
  for (countpos=numbits;countpos>=0;countpos--)
  {
    negLogBoxCountS[countpos] = -log((double)boxCountS[countpos])/log(2.0);

    if (debugging)
      printf("pointCount[%d] is %lud...\n", countpos, pointCount[countpos]);
    freq = pointCount[countpos]/(double)dataLines ;

    sumSqrFreq[countpos] += (freq * freq) ;
    sumSqrFreq[countpos] *= (dataLines/(double) boxCountS[0]) ;
    sumSqrFreq[countpos] *= (dataLines/(double) boxCountS[0]) ;

       /* sumSqrFreq[countpos] NOW CONTAINS THE SUM OF THE SQUARES OF THE
	  FREQUENCIES OF POINTS IN ALL OCCUPIED BOXES OF THE SIZE
	  CORRESPONDING TO countpos. */
	  
    logSumSqrFreqS[countpos] = log(sumSqrFreq[countpos])/log(2.0) ;
    informationS[countpos] += ( freq*log(freq)/log(2.0) ) ;
    informationS[countpos] *= (dataLines/(double)boxCountS[0]) ;
    informationS[countpos] +=
       ( log((double)dataLines)-log((double)boxCountS[0]) ) / log(2.0) ;

       /* information[countpos] NOW CONTAINS THE INFORMATION SUM FOR ALL THE
          OCCUPIED BOXES OF THIS SIZE. */

   }
}
/* ################################################################## */
/* ################################################################## */
  /*  FIT LEAST SQUARE LINE TO DATA IN X,Y.  NO PROTECTION AGAINST OVERFLOW
      HERE.  IT IS ASSUMED THAT LAST > FIRST AND THAT THE X'S ARE NOT ALL THE
      SAME -- ELSE DIVISION BY ZERO WILL OCCUR.  */
void fitLSqrLine (first, last, X, Y, slopePtr, interceptPtr)
     long   first, last ;
     double X[], Y[], *slopePtr, *interceptPtr ;
{
  int    index , pointCount ;
  double Xsum=0, Ysum=0, XYsum=0, XXsum=0, Xmean=0, Ymean=0,
         Xtemp, Ytemp;
  for (index=first; index<=last; index++)
    { Xtemp = X[index]; Ytemp = Y[index];
      Xsum += Xtemp; Ysum += Ytemp;
      XYsum += (Xtemp * Ytemp); XXsum += (Xtemp * Xtemp);  }
  pointCount = last - first + 1 ;
  Xmean = Xsum/pointCount;  Ymean = Ysum/pointCount;
  *slopePtr = (XYsum - Xsum * Ymean)/(XXsum - Xsum * Xmean) ;
  *interceptPtr = Ymean - *slopePtr * Xmean ;
}
/* ################################################################## */
/* ################################################################## */
/*
MARK GREATEST INDEX WHERE COUNT > boxCountF[0]/cutOff_factor.

COUNTS AT LESSER INDEXES WILL NOT BE USED IN THE ESTIMATE OF FRACTAL
DIMENSION -- DISTORTION DUE TO SATURATION IS THE CONCERN.

NOTE THAT boxCountF[0] IS THE NUMBER OF BOXES OF SIZE 1 (THE SMALLEST SIZE)
THAT CONTAIN A POINT OF THE SET.  FOR ALL PRACTICAL PURPOSES, boxCountF[0]
WILL EQUAL THE NUMBER OF DISTINCT POINTS IN THE INPUT FILE, BECAUSE THESE
BOXES ARE REALLY SMALL COMPARED TO THE SIZE OF THE BIGGEST BOX (ABOUT 4
BILLION IF AN UNSIGNED LONG INT IS 32 BITS TO THE PLATFORM COMPILER.  THE
POINTS ARE SCALED BY THE PROGRAM SO THAT THE SET IS TOO "LARGE" TO FIT IN
THE NEXT SMALLEST BOX SIZE, SO THAT "1" IS THE SMALLEST DIFFERENCE IN
VALUE THAT CAN BE RESOLVED.) ONE BOX, IN EFFECT, COVERS ONLY A SINGLE POINT
OF THE INPUT SET BECAUSE THE PROGRAM CAN'T RESOLVE POINTS WITH A SMALLER
DIFFERENCE.

WE THINK IT WOULD BE A BAD IDEA TO USE dataLines/cutOff_factor AS THE LIMIT
BECAUSE IN CASES WHERE THERE WERE MANY DUPLICATE POINTS, WE WOULD SERIOUSLY
OVER-ESTIMATE THE NUMBER OF DISTINCT POINTS, AND THUS USE SATURATED DATA TO
BASE THE ESTIMATE OF FRACTAL DIMENSION UPON.  WHEN TESTING THE PROGRAM WITH
RANDOM DATA SAMPLED WITH REPLACEMENT, THIS COULD THROW THE RESULTS WAY OFF.
(THIS HAPPENED TO US, AND IT TOOK US A WHILE TO FIGURE OUT WHY.  AFTERWARDS,
WE STOPPED USING dataLines/cutOff_factor, AND CHANGED TO
boxCountF[0]/cutOff_factor.)
*/
void findMark(markPtr, boxCountM)
   ulong *markPtr, boxCountM[numbits+1] ;
{
    int     i,  cutOff_factor=1;
    
   /* Calculate cutOff_factor = 2^(embed_dim) + 1 */
 for (i=1;i<=embed_dim;i++) cutOff_factor = cutOff_factor * 2;
 cutOff_factor++;

 *markPtr=0;
 for(i=0;i<numbits;i++)
   { if(boxCountM[i] > boxCountM[0]/cutOff_factor) *markPtr=i; }
}

/* ################################################################## */
/* ################################################################## */
float GetDims(negLogBoxCountF, logSumSqrFreqF, informationF, markF,
                capDimPtr, infDimPtr, corrDimPtr)
        ulong   markF ;
        double  negLogBoxCountF[numbits+1], informationF[numbits+1],
                logSumSqrFreqF[numbits+1],
	        *capDimPtr, *infDimPtr, *corrDimPtr;
{
    int     i;
    double  logEps[numbits+1], slope, intercept;
    
    /* GET LOG (BASE 2) OF THE DIAMETER OF THE I'TH SIZE OF BOX. */
  for(i=numbits; i>=0; i--)  logEps[i] = i;

/* fitLSqrLine (markF, numbits-4, logEps, negLogBoxCountF, &slope,&intercept);*/
  fitLSqrLine (markF, numbits-2, logEps, negLogBoxCountF, &slope, &intercept);
  *capDimPtr = slope ;
/*fitLSqrLine(markF, numbits-4, logEps, informationF, &slope, &intercept);*/
  fitLSqrLine(markF, numbits-2, logEps, informationF, &slope, &intercept);
  *infDimPtr = slope ;
/*fitLSqrLine(markF,numbits-4, logEps, logSumSqrFreqF, &slope, &intercept);*/
  fitLSqrLine(markF,numbits-2, logEps, logSumSqrFreqF, &slope, &intercept);
  *corrDimPtr = slope ;
}
/* ################################################################## */
/* ################################################################## */
