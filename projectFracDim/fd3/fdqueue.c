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

int  empty_Av();
void push_Av();
void pop_Av();
void create_Av();
void print_Av();

int empty_Q();
void en_Q();
void de_Q();
void create_Q();
void transf_Q();
void print_Q();


int empty_Av(avail)
       long int avail;
{
  return (avail == -1);
}

/*
    Puts the info pointed to by "pointer" into the avail list.  Use
    this after dequeuing an element if you want to recycle it.  Note the
    need to access the global array "next_after", which holds all the
    pointers.
*/
void push_Av(p_avail, pointer)
	long int *p_avail, pointer ;
{
  extern long int  *next_after;   /*  array for pointers to data. */
  next_after[pointer] = *p_avail;
  *p_avail = pointer;
}

/*
   Pops the info pointed to by "p_pointer" from the avail list.  Use
   this if you need an element to place on a queue.
*/
void pop_Av(p_avail, p_pointer)
	long int *p_avail, *p_pointer;
{
  extern long int  *next_after;   /*  array for pointers to data. */
  *p_pointer = *p_avail;
  *p_avail = next_after[*p_avail];
}

/*
    Creates an avail list with room for all the data, plus two more
    elements to be used as markers.
*/
void create_Av(p_avail)
	long int *p_avail;
{
  extern long int  *next_after;   /*  array for pointers to data. */
  extern ulong dataLines;        /* number of data points in input file */
  unsigned long int i;
  for (i=0;i<dataLines+1;i++)
    next_after[i] = i+1;
  next_after[dataLines+1] = -1;
  *p_avail = 0;
}

/*
   This may prove to be useful when debugging.
*/
void print_Av(avail)
	long int avail;
{
  extern long int  *next_after;   /*  array for pointers to data. */
  extern int  *marker;    /*  To mark queues during radix sort */
  extern int   embed_dim;         /* How many coordinates points have. */
  extern ulong **data;	/* array of data points */
  long int temp;
  int coord;
  printf("Starting to print the avail list:\n\n");
  printf("Index\tNextField\tMarker?\t\tCoords\n\n");
  temp = avail;
  while (temp != -1)
    {
        printf("%ld\t%ld\t\t%d\t", temp, next_after[temp], marker[temp]);
	for (coord=0;coord<embed_dim;coord++)
	  printf("\t%lu", data[coord][temp]);
        printf("\n");
        temp = next_after[temp];
     }
 printf("End of printing of avail list.\n");
}


int empty_Q(Q)
	QHeader Q; 	/* The queue to be tested. */
{
   return(Q.Qrear == -1);
}

/*
    This function enqueues the info at "pointer" into "*p_QHeader".  It does
    NOT remove the info from the avail list, NOR does it check to see if the
    queue is full.  Use with caution.  BEFORE calling this function, you
    should free the info from the avail list and be sure somehow that
    "*p_QHeader" is not full.  Note also that the function needs to access
    the array called "next_after" globally.
*/

void en_Q(p_QHeader,pointer)
	QHeader *p_QHeader;		/* Queue to enqueue. */
	long int pointer;    /* array index of info to go into the queue. */
{
  extern long int  *next_after;   /*  array for pointers to data. */
  next_after[pointer] = -1;
  if (empty_Q(*p_QHeader))
     p_QHeader->Qfront = pointer;
  else
     next_after[p_QHeader->Qrear] = pointer;
  p_QHeader->Qrear = pointer;
}

/*
    This Dequeuing function does NOT check to see if the queue is
    empty first, nor does it return the dequeued info to the avail
    list.  The caller MUST take care of these things!  Note also
    that the function needs to access the array called "next_after"
    globally.
*/
void de_Q(p_QHeader,p_pointer)
	QHeader *p_QHeader;		/* Queue to be dequeued. */
	long int *p_pointer;	/* Pointer to the data to be returned */
				/* to the calling function. */	
{
  extern long int  *next_after;   /*  array for pointers to data. */
  *p_pointer = p_QHeader->Qfront;
  if ((p_QHeader->Qfront) == (p_QHeader->Qrear))
     p_QHeader->Qrear = -1 ;
  p_QHeader->Qfront = next_after[p_QHeader->Qfront];
}

void create_Q (p_QHeader)
	QHeader *p_QHeader;	/* Queue to be created. */
{
  p_QHeader->Qfront = -1;
  p_QHeader->Qrear = -1;
}

/*
    This is a function to rapidly move an "element" from the source
    to the target queue.  It bypasses the actions of placing the info
    into the avail list, and taking it back out again.  It also
    avoids testing to see if the source queue is being emptied, and
    testing to see if the target queue is empty before the enqueue
    is done.

    Such a function is useful because there will be a great many of
    these operations done in the radix sort of a large set.  For
    further time-savings, one should think about putting this code
    in-line, so that the overhead of function-calling can be saved.

    USE with EXTREME CAUTION.  The effect of the function will be
    INCORRECT if the source has one element left, or the target is
    empty!

    Note the reliance on the global array "next_after".
*/
void transf_Q(p_Q1,p_Q2)
        QHeader *p_Q1, *p_Q2;  /* source and target queues */
{
  extern long int  *next_after;   /*  array for pointers to data. */
  long int temp;
  temp = p_Q1->Qfront;                     /* save front of Q1 */
  p_Q1->Qfront = next_after[p_Q1->Qfront];  /* fast dequeue */
  next_after[temp] = -1;               /* fast enqueue */
  next_after[p_Q2->Qrear] = temp;        /* fast enqueue */
  p_Q2->Qrear = temp;                    /* fast enqueue */
}

/*
    This may come in handy when de-bugging.  Note the reliance on
    the global arrays "next_after" and "data".
*/
void print_Q(p_Q)
     QHeader *p_Q;
{
  extern long int  *next_after;   /*  array for pointers to data. */
  extern int  *marker;    /*  To mark queues during radix sort */
  extern int   embed_dim;         /* How many coordinates points have. */
  extern ulong **data;	/* array of data points */  
  long int temp;
  int coord;
  printf("Starting to print the queue:\n\n");
  printf("Index\tNextField\tMarker?\t\tCoords\n\n");
  temp = p_Q->Qfront;
  while (temp != -1)
    {
        printf("%ld\t%ld\t\t%d\t", temp, next_after[temp], marker[temp]);
	for (coord=0;coord<embed_dim;coord++)
	  printf("\t%lu", data[coord][temp]);
        printf("\n");
        temp = next_after[temp];
     }
 printf("End of printing of queue.\n");
}
