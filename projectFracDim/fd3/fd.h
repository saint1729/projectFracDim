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

/*
	fd.h  Header file for Fractal Dimension Programs
*/

#define MASK(n)		(((unsigned long int) 1)<<(n))
#define IS_A_ONE(x,n)	((x) & MASK(n))
#define false		(0)
#define true		(1)
#define sizul		sizeof(unsigned long int)
#define numbits		(8*sizul)
/* Use this (change to true) for a small amount of checking. */
#define checking   (false)  

/* This gives EXCESSIVE output -- even with tiny input sets -- use in a way
   that allows you to abort when you have seen enough. */
#define debugging (false) 

typedef unsigned long int	ulong;

typedef struct QHeader
{
   long int  Qfront,	/* indexes used as pointers to the front */
	     Qrear;     /* and rear of the queue */
} QHeader;

/*  Some things that are global variables. */
 FILE  *infile;
 int   embed_dim,         /* How many coordinates points have. */
             *marker;           /*  To mark queues during radix sort */
 ulong **data,	/* initial array of data points */
             dataLines,        /* number of data points in input file */
            *diff_test; 	 /* XOR's of pairs of coordinates */
 long int  *next_after,    /*  array for pointers to data. */
                 avail;             /*  avail list pointer */
 QHeader Q[2];       /* array of two queues for the radix sort */

