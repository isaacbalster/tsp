#include "hi_pr.hpp"

typedef double excessType;

typedef double cType;

typedef struct arcSt {
  cType           resCap;          /* residual capacity */
  struct nodeSt   *head;           /* arc head */
  struct arcSt    *rev;            /* reverse arc */
} arc;

typedef struct nodeSt {
  arc             *first;           /* first outgoing arc */
  arc             *current;         /* current outgoing arc */
  excessType      excess;           /* excess at the node 
				       change to double if needed */
  long            d;                /* distance label */
  struct nodeSt   *bNext;           /* next node in bucket */
  struct nodeSt   *bPrev;           /* previous node in bucket */
} node;


typedef struct bucketSt {
  node             *firstActive;      /* first node with positive excess */
  node             *firstInactive;    /* first node with zero excess */
} bucket;

#define EPS 0.00001

int parse( double  **x,      /* capa */
	   long    n_nodes, /* nb nodes */
	   long    s,       /* source */
	   long    t,       /* sink */
	   long    *n_ad,
	   long    *m_ad,
	   node    **nodes_ad,
	   arc     **arcs_ad,
	   double  **cap_ad,
	   node    **source_ad,
	   node    **sink_ad,
	   long    *node_min_ad) {

  long    n,                      /* internal number of nodes */
          node_min=0,             /* minimal no of node  */
          node_max=0,             /* maximal no of nodes */
         *arc_first=NULL,         /* internal array for holding
                                     - node degree
                                     - position of the first outgoing arc */
         *arc_tail=NULL,          /* internal array: tails of the arcs */
          source=0,               /* no of the source */
          sink=0,                 /* no of the sink */
          /* temporary variables carrying no of nodes */
          head, tail, i;

  long    m,                      /* internal number of arcs */
          /* temporary variables carrying no of arcs */
          last, arc_num, arc_new_num;

  node    *nodes=NULL,            /* pointer to the node structure */
          *head_p,
          *ndp;

  arc     *arcs=NULL,             /* pointer to the arc structure */
          *arc_current=NULL,
          *arc_new,
          *arc_tmp;

  double  *acap=NULL,             /* array of capasities */
          cap;                    /* capasity of the current arc */

  long    no_alines=0,            /* no of arc-lines */
          pos_current=0;          /* 2*no_alines */

/* The main loop:
        -  reads the line of the input,
        -  analyses its type,
        -  checks correctness of parameters,
        -  puts data to the arrays,
        -  does service functions
*/
  int count = 0;
  int count1 = 0;
  int count2 = 0;
  int count3 = 0;

  for (int jj = 0; jj < n_nodes; jj++){
    for (int kk = 0; kk < n_nodes; kk++){
      if (x[jj][kk] > EPS){
/* 	printf("%d --- %d : %f\n",jj,kk,x[jj][kk]); */
	count++;
/* 	if (jj == s){ */
/* 	  count1++; */
/* 	} */
	if (jj == n_nodes-1){
	  count2++;
	}
	if (jj == 0){
	  count3++;
	}
      }
    }
  }
/*   printf("source : %d    puits : %d  \n\n\n",s,t); */
/*   count = n_arcs; */

/*   if (count1 == 0){ */
/*     count++; */
/*   } */
  if (count2 == 0){
    count++;
  }
  if (count3 == 0){
    count++;
  }

  n=n_nodes;
  m=count;

  /* allocating memory for  'nodes', 'arcs'  and internal arrays */
  nodes     = (node*) calloc ( n+2, sizeof(node) );
  arcs      = (arc*)  calloc ( 2*m+1, sizeof(arc) );
  arc_tail  = (long*) calloc ( 2*m,   sizeof(long) ); 
  arc_first = (long*) calloc ( n+2, sizeof(long) );
  acap      = (double*) calloc ( 2*m, sizeof(double) );
  /* arc_first [ 0 .. n+1 ] = 0 - initialized by calloc */
  
  /* setting pointer to the first arc */
  arc_current = arcs;

  source = s;
  sink = t;
  node_max = 0;
  node_min = n;

  for (int ii = 0; ii < n_nodes; ii++){
    for (int jj = 0; jj < n_nodes; jj++){
/*       if (x[ii][jj] > EPS){ */
      if ((x[ii][jj] > EPS)
/* 	  || ((ii == s) && (jj == t) && (count1 == 0)) */
	  || ((ii == n_nodes-1) && (jj == 0) && (count2 == 0))
	  || ((ii == 0) && (jj == n_nodes-1) && (count3 == 0))){
	tail = ii;
	head = jj;
	cap = x[ii][jj];

	arc_first[tail + 1] ++; 
	arc_first[head + 1] ++;

	/* storing information about the arc */
	arc_tail[pos_current]        = tail;
	arc_tail[pos_current+1]      = head;
	arc_current       -> head    = nodes + head;
	arc_current       -> resCap    = cap;
	arc_current       -> rev  = arc_current + 1;
	( arc_current + 1 ) -> head    = nodes + tail;
	( arc_current + 1 ) -> resCap    = 0.;
	( arc_current + 1 ) -> rev  = arc_current;

	/* searching minimum and maximum node */
	if ( head < node_min ) node_min = head;
	if ( tail < node_min ) node_min = tail;
	if ( head > node_max ) node_max = head;
	if ( tail > node_max ) node_max = tail;

	no_alines   ++;
	arc_current += 2;
	pos_current += 2;
      }
    }
  }

  /********** ordering arcs - linear time algorithm ***********/

  /* first arc from the first node */
  ( nodes + node_min ) -> first = arcs;

  /* before below loop arc_first[i+1] is the number of arcs outgoing from i;
     after this loop arc_first[i] is the position of the first 
     outgoing from node i arcs after they would be ordered;
     this value is transformed to pointer and written to node.first[i]
   */
 
  for ( i = node_min + 1; i <= node_max + 1; i ++ ) {
    arc_first[i]          += arc_first[i-1];
    ( nodes + i ) -> first = arcs + arc_first[i];
  }


  for ( i = node_min; i < node_max; i ++ ) /* scanning all the nodes  
					      exept the last*/
    {

      last = ( ( nodes + i + 1 ) -> first ) - arcs;
      /* arcs outgoing from i must be cited    
	 from position arc_first[i] to the position
	 equal to initial value of arc_first[i+1]-1  */

      for ( arc_num = arc_first[i]; arc_num < last; arc_num ++ )
	{ tail = arc_tail[arc_num];

	while ( tail != i )
          /* the arc no  arc_num  is not in place because arc cited here
             must go out from i;
             we'll put it to its place and continue this process
             until an arc in this position would go out from i */

	  { arc_new_num  = arc_first[tail];
	  arc_current  = arcs + arc_num;
	  arc_new      = arcs + arc_new_num;
	    
	  /* arc_current must be cited in the position arc_new    
	     swapping these arcs:                                 */

	  head_p               = arc_new -> head;
	  arc_new -> head      = arc_current -> head;
	  arc_current -> head  = head_p;

	  cap                 = arc_new -> resCap;
	  arc_new -> resCap     = arc_current -> resCap;
	  arc_current -> resCap = cap;

	  if ( arc_new != arc_current -> rev )
	    {
	      arc_tmp                = arc_new -> rev;
	      arc_new  -> rev     = arc_current -> rev;
	      arc_current -> rev  = arc_tmp;

	      ( arc_current -> rev ) -> rev = arc_current;
	      ( arc_new     -> rev ) -> rev = arc_new;
	    }

	  arc_tail[arc_num] = arc_tail[arc_new_num];
	  arc_tail[arc_new_num] = tail;

	  /* we increase arc_first[tail]  */
	  arc_first[tail] ++ ;

	  tail = arc_tail[arc_num];
	  }
	}
      /* all arcs outgoing from  i  are in place */
    }       

  /* -----------------------  arcs are ordered  ------------------------- */

  /*----------- constructing lists ---------------*/


  for ( ndp = nodes + node_min; ndp <= nodes + node_max;  ndp ++ )
    ndp -> first = (arc*) NULL;

  for ( arc_current = arcs + (2*m-1); arc_current >= arcs; arc_current -- )
    {
      arc_num = arc_current - arcs;
      tail = arc_tail [arc_num];
      ndp = nodes + tail;
      /* avg
	 arc_current -> next = ndp -> first;
      */
      ndp -> first = arc_current;
    }

  for ( ndp = nodes + node_max; ndp >= nodes + node_min; ndp -- )
    {
      if (ndp -> first == (arc*) NULL)
	{
	  ndp -> first = (ndp + 1) -> first;
	}
    }

  /* ----------- assigning output values ------------*/
  *m_ad = m;
  *n_ad = node_max - node_min + 1;
  *source_ad = nodes + source;
  *sink_ad   = nodes + sink;
  *node_min_ad = node_min;
  *nodes_ad = nodes + node_min;
  *arcs_ad = arcs;
  *cap_ad = acap;

  for ( arc_current = arcs, arc_num = 0; 
	arc_num < 2*m;
	arc_current ++, arc_num ++
	)
    acap [ arc_num ] = arc_current -> resCap; 

  /* free internal memory */
  free ( arc_first ); free ( arc_tail );

/*   if (n_nodes > *n_ad){ */
/*     for (int jj = 0; jj < n_nodes; jj++){ */
/*       for (int kk = 0; kk < n_nodes; kk++){ */
/* 	if (x[jj][kk] > EPS){ */
/* 	  printf("%d --- %d : %f\n",jj,kk,x[jj][kk]); */
/* 	} */
/*       } */
/*     } */
/*   } */
/*   printf("avant %d   apres  %d\n",count,*m_ad); */


  /* Thanks God! all is done */
  return (0);

}

#define GLOB_UPDT_FREQ 0.5
#define ALPHA 6
#define BETA 12

#define WHITE 0
#define GREY 1
#define BLACK 2

/* global variables */

long   n;                    /* number of nodes */
long   m;                    /* number of arcs */
long   nm;                   /* n + ALPHA * m */
long   nMin;                 /* smallest node id */
node   *nodes;               /* array of nodes */
arc    *arcs;                /* array of arcs */
bucket *buckets;             /* array of buckets */
cType  *cap;                 /* array of capacities */
node   *source;              /* source node pointer */
node   *sink;                /* sink node pointer */
node   **queue;              /* queue for BFS */
node   **qHead, **qTail, **qLast;     /* queue pointers */
long   dMax;                 /* maximum label */
long   aMax;                 /* maximum actie node label */
long   aMin;                 /* minimum active node label */
double flow;                 /* flow value */
long pushCnt  = 0;           /* number of pushes */
long relabelCnt   = 0;       /* number of relabels */
long updateCnt    = 0;       /* number of updates */
long gapCnt   = 0;           /* number of gaps */
long gNodeCnt = 0;           /* number of nodes after gap */  
node   *sentinelNode;        /* end of the node list marker */
arc *stopA;                  /* used in forAllArcs */
long workSinceUpdate=0;      /* the number of arc scans since last update */
float globUpdtFreq;          /* global update frequency */

/* macros */

#define forAllNodes(i) for ( i = nodes; i != sentinelNode; i++ )
#define forAllArcs(i,a) for (a = i->first, stopA = (i+1)->first; a != stopA; a++)

#define nNode( i ) ( (i) - nodes + nMin )
#define nArc( a )  ( ( a == NULL )? -1 : (a) - arcs )

#define min( a, b ) ( ( (a) < (b - EPS) ) ? a : b )

/* FIFO queue for BFS macros */
#define qInit() \
{\
  qHead = qTail = queue;\
}

#define qEmpty ( qHead == qTail )

#define qEnqueue(i) \
{\
  *qTail = i;\
  if ( qTail == qLast ) qTail = queue;\
  else qTail++;\
}

#define qDequeue(i) \
{\
  i = *qHead;\
  if ( qHead == qLast ) qHead = queue;\
  else qHead++;\
}


/* 
   bucket macros:
   bucket's active node list is singly-linked
     operations aAdd, aRemove (from the front)
   bucket's inactive list is doubly-linked
     operations iAdd, iDelete (from arbitrary position)
*/

long i_dist;

#define aAdd(l,i)\
{\
  i->bNext = l->firstActive;\
  l->firstActive = i;\
  i_dist = i->d;\
  if (i_dist < aMin)\
    aMin = i_dist;\
  if (i_dist > aMax)\
    aMax = i_dist;\
}

/* i must be the first element */
#define aRemove(l,i)\
{\
  l->firstActive = i->bNext;\
}

node *i_next, *i_prev;
#define iAdd(l,i)\
{\
  i_next = l->firstInactive;\
  i->bNext = i_next;\
  i->bPrev = sentinelNode;\
  i_next->bPrev = i;\
  l->firstInactive = i;\
}

#define iDelete(l,i)\
{\
  i_next = i->bNext;\
  if (l->firstInactive == i) {\
    l->firstInactive = i_next;\
    i_next->bPrev = sentinelNode;\
  }\
  else {\
    i_prev = i->bPrev;\
    i_prev->bNext = i_next;\
    i_next->bPrev = i_prev;\
  }\
}

/* allocate datastructures, initialize related variables */

int allocDS( )

{

  nm = ALPHA * n + m;
  queue = (node**) calloc ( n, sizeof (node*) );
  if ( queue == NULL ) return ( 1 );
  qLast = queue + n - 1;
  qInit();

  buckets = (bucket*) calloc ( n+2, sizeof (bucket) );
  if ( buckets == NULL ) return ( 1 );

  sentinelNode = nodes + n;
  sentinelNode->first = arcs + 2*m;

  return ( 0 );

} /* end of allocate */


void init( )

{
  node  *i;        /* current node */
  int overflowDetected;
  double delta;
  bucket *l;
  arc *a;

  // initialize excesses

  forAllNodes(i) {
    i->excess = 0.0;
    i->current = i->first;
/*     printf("%d \n", i->first); */
    forAllArcs(i, a){
/*       printf("%d --- \n",nNode(i)); */
/*       printf("%f \n",a->resCap); */
      a->resCap = cap[a-arcs];
    }
  }

  overflowDetected = 0;
  source->excess = 0.0;
  forAllArcs(source,a) {
    if (a->head != source) {
      pushCnt ++;
      delta = a -> resCap;
      a -> resCap -= delta;
      (a -> rev) -> resCap += delta;
      a->head->excess += delta;
    }
  }

  /*  setup labels and buckets */
  for (l = buckets; l <= buckets + n-1; l++) {
    l -> firstActive   = sentinelNode;
    l -> firstInactive  = sentinelNode;
  }
    
  l = buckets + 1;
    
  aMax = 0;
  aMin = n;
    
  forAllNodes(i) {
    if (i == sink) {
      i->d = 0;
      continue;
    }
    if ((i == source) && (!overflowDetected)) {
      i->d = n;
    }
    else
      i->d = 1;
    if (i->excess > EPS) {
      /* put into active list */
      aAdd(l,i);
    }
    else { /* i -> excess == 0 */
      /* put into inactive list */
      if (i->d < n)
	iAdd(l,i);
    }
  }
  dMax = 1;

  dMax = n-1;
  flow = 0.0;

} /* end of init */



/* global update via backward breadth first search from the sink */

void globalUpdate ()

{

  node  *i, *j;       /* node pointers */
  arc   *a;           /* current arc pointers  */
  bucket *l;          /* bucket */
  long  jD;           /* d of node j */


  updateCnt ++;

  /* initialization */

  forAllNodes(i)
    i -> d = n;
  sink -> d = 0;

  for (l = buckets; l <= buckets + dMax; l++) {
    l -> firstActive   = sentinelNode;
    l -> firstInactive  = sentinelNode;
  }

  dMax = aMax = 0;
  aMin = n;

  /* breadth first search */

  qEnqueue (sink);
  do {
    qDequeue (i);
    jD = (i->d) + 1;

    /* scanning arcs incident to node i */
    forAllArcs(i,a) {
/*       printf("coucou ? %d %d \n",a,sentinelNode->first); */
      if (a->rev->resCap > EPS ) {
	j = a->head;
	if (j->d == n) {
	  j->d = jD;
	  j->current = j->first;
	  l = buckets + jD;
	  if (jD > dMax) dMax = jD;

	  if (j->excess > EPS) {
	    /* put into active list */
	    aAdd(l,j);
	  }
	  else
	    /* put into inactive list */
	    iAdd(l,j);
    	  qEnqueue (j);
    	}
      }
/*   printf("coucou ! \n"); */
    } /* node i is scanned */ 
  } while (!qEmpty);

  assert(source->d == n);

} /* end of global update */

/* second stage -- preflow to flow */
void stageTwo ( )
/*
   do dsf in the reverse flow graph from nodes with excess
   cancel cycles if found
   return excess flow in topological order
*/

/*
   i->d is used for dfs labels 
   i->bNext is used for topological order list
   buckets[i-nodes]->firstActive is used for DSF tree
*/

{
  node *i, *j, *tos, *bos, *restart, *r;
  arc *a;
  cType delta;

  /* deal with self-loops */
  forAllNodes(i) {
    forAllArcs(i,a)
      if ( a -> head == i ) {
	a -> resCap = cap[a - arcs];
      }
  }

  /* initialize */
  tos = bos = NULL;
  forAllNodes(i) {
    i -> d = WHITE;
    buckets[i-nodes].firstActive = NULL;
    i -> current = i -> first;
  }

  /* eliminate flow cycles, topologically order vertices */
  forAllNodes(i)
    if (( i -> d == WHITE ) && ( i -> excess > EPS ) &&
	( i != source ) && ( i != sink )) {
      r = i;
      r -> d = GREY;
      do {
	for ( ; i->current != (i+1)->first; i->current++) {
	  a = i -> current;
	  if (( (-EPS <= cap[a - arcs]) && (cap[a - arcs] <= EPS) )
	      && ( a -> resCap > EPS )) { 
	    j = a -> head;
	    if ( j -> d == WHITE ) {
	      /* start scanning j */
	      j -> d = GREY;
	      buckets[j-nodes].firstActive = i;
	      i = j;
	      break;
	    }
	    else
	      if ( j -> d == GREY ) {
		/* find minimum flow on the cycle */
		delta = a -> resCap;
		while ( 1 ) {
		  delta = min ( delta, j -> current -> resCap );
		  if ( j == i )
		    break;
		  else
		    j = j -> current -> head;
		}

		/* remove delta flow units */
		j = i;
		while ( 1 ) {
		  a = j -> current;
		  a -> resCap -= delta;
		  a -> rev -> resCap += delta;
		  j = a -> head;
		  if ( j == i )
		    break;
		}
	  
		/* backup DFS to the first saturated arc */
		restart = i;
		for ( j = i -> current -> head; j != i; j = a -> head ) {
		  a = j -> current;
		  if (( j -> d == WHITE ) || ( (-EPS <= a -> resCap) && (a -> resCap <= EPS) )) {
		    j -> current -> head -> d = WHITE;
		    if ( j -> d != WHITE )
		      restart = j;
		  }
		}
	  
		if ( restart != i ) {
		  i = restart;
		  i->current++;
		  break;
		}
	      }
	  }
	}

	if (i->current == (i+1)->first) {
	  /* scan of i complete */
	  i -> d = BLACK;
	  if ( i != source ) {
	    if ( bos == NULL ) {
	      bos = i;
	      tos = i;
	    }
	    else {
	      i -> bNext = tos;
	      tos = i;
	    }
	  }

	  if ( i != r ) {
	    i = buckets[i-nodes].firstActive;
	    i->current++;
	  }
	  else
	    break;
	}
      } while ( 1 );
    }


  /* return excesses */
  /* note that sink is not on the stack */
  if ( bos != NULL ) {
    for ( i = tos; i != bos; i = i -> bNext ) {
      a = i -> first;
      while ( i -> excess > EPS ) {
	if (( (-EPS <= cap[a - arcs]) && (cap[a -arcs] <= EPS) ) && ( a -> resCap > EPS )) {
	  if (a->resCap < i->excess - EPS)
	    delta = a->resCap;
	  else
	    delta = i->excess;
	  a -> resCap -= delta;
	  a -> rev -> resCap += delta;
	  i -> excess -= delta;
	  a -> head -> excess += delta;
	}
	a++;
      }
    }
    /* now do the bottom */
    i = bos;
    a = i -> first;
    while ( i -> excess > EPS ) {
      if (( (-EPS <= cap[a - arcs]) && (cap[a-arcs] <= EPS) ) && ( a -> resCap > EPS )) {
	if (a->resCap < i->excess - EPS)
	  delta = a->resCap;
	else
	  delta = i->excess;
	a -> resCap -= delta;
	a -> rev -> resCap += delta;
	i -> excess -= delta;
	a -> head -> excess += delta;
      }
      a++;
    }
  }
}


/* gap relabeling */

int gap (bucket *emptyB){

  bucket *l;
  node  *i; 
  long  r;           /* index of the bucket before l  */
  int   cc;          /* cc = 1 if no nodes with positive excess before
		      the gap */

  gapCnt ++;
  r = ( emptyB - buckets ) - 1;

  /* set labels of nodes beyond the gap to "infinity" */
  for ( l = emptyB + 1; l <= buckets + dMax; l ++ ) {
    /* this does nothing for high level selection 
    for (i = l -> firstActive; i != sentinelNode; i = i -> bNext) {
      i -> d = n;
      gNodeCnt++;
    }
    l -> firstActive = sentinelNode;
    */

    for ( i = l -> firstInactive; i != sentinelNode; i = i -> bNext ) {
      i -> d = n;
      gNodeCnt ++;
    }

    l -> firstInactive = sentinelNode;
  }

  cc = ( aMin > r ) ? 1 : 0;

  dMax = r;
  aMax = r;

  return ( cc );

}

/*--- relabelling node i */

long relabel (node *i){

  node  *j;
  long  minD;     /* minimum d of a node reachable from i */
  arc   *minA;    /* an arc which leads to the node with minimal d */
  arc   *a;

  assert(i->excess > EPS);

  relabelCnt++;
  workSinceUpdate += BETA;

  i->d = minD = n;
  minA = NULL;

  /* find the minimum */
  forAllArcs(i,a) {
    workSinceUpdate++;
    if (a -> resCap > EPS) {
      j = a -> head;
      if (j->d < minD) {
	minD = j->d;
	minA = a;
      }
    }
  }

  minD++;
      
  if (minD < n) {

    i->d = minD;
    i->current = minA;

    if (dMax < minD) dMax = minD;

  } /* end of minD < n */
      
  return ( minD );

} /* end of relabel */


/* discharge: push flow out of i until i becomes inactive */

void discharge (node  *i){

  node  *j;                 /* sucsessor of i */
  long  jD;                 /* d of the next bucket */
  bucket *lj;               /* j's bucket */
  bucket *l;                /* i's bucket */
  arc   *a;                 /* current arc (i,j) */
  cType  delta;
  arc *stopA;

  assert(i->excess > EPS);

  do {

    jD = i->d - 1;
    l = buckets + i->d;

    /* scanning arcs outgoing from  i  */
    for (a = i->current, stopA = (i+1)->first; a != stopA; a++) {
      if (a -> resCap > EPS) {
	j = a -> head;

	if (j->d == jD) {
	  pushCnt ++;
	  if (a->resCap < i->excess - EPS)
	    delta = a->resCap;
	  else
	    delta = i->excess;
	  a->resCap -= delta;
	  a->rev->resCap += delta;

	  if (j != sink) {
	    lj = buckets + jD;

	    if ((-EPS <= j->excess) && (j->excess <= EPS)) {
	      /* remove j from inactive list */
	      iDelete(lj,j);
	      /* add j to active list */
	      aAdd(lj,j);
	    }
	  }
	  
	  j -> excess += delta;
	  i -> excess -= delta;
	  
	  if ((-EPS <= i->excess) && (i->excess <= EPS)) break;

	} /* j belongs to the next bucket */
      } /* a  is not saturated */
    } /* end of scanning arcs from  i */

    if (a == stopA) {
      /* i must be relabeled */
      relabel (i);
      if ((l -> firstActive == sentinelNode) && 
	  (l -> firstInactive == sentinelNode)
	  ) /* gap is found */
	gap (l);
      if (i->d == n) break;
    }
    else {
      /* i no longer active */
      i->current = a;
      /* put i on inactive list */
      iAdd(l,i);
      break;
    }
  } while (1);
}

/* first stage  -- maximum preflow*/

void stageOne ( )

{

  node   *i;
  bucket  *l;             /* current bucket */


  workSinceUpdate = 0;

  /* main loop */
  while ( aMax >= aMin ) {
    l = buckets + aMax;
    i = l->firstActive;

    if (i == sentinelNode)
      aMax--;
    else {
      aRemove(l,i);

      assert(i->excess > EPS);
      discharge (i);

      /* is it time for global update? */
      if (workSinceUpdate * globUpdtFreq > nm) {
	globalUpdate ();
	workSinceUpdate = 0;
      }

    }
    
  } /* end of the main loop */
    
  flow = sink -> excess;

} 


/* Compute the directed minimum cut separating the sink from the source.
 * Parameters :
 * - (in)  double **x : arc capacity matrix
 * - (in)  long n_nodes : number of nodes (size of the graph)
 * - (in)  long source : source node between 0 and n_nodes-1
 * - (in)  long sin : sink node between 0 and n_nodes-1
 * - (out) double & val_flow : value of the resulting min cut / max flow
 * - (out) double *& dist : distance of the nodes to the sink in the residual graph,
 *                          needs to be initialized to n_nodes,
 *                          the minimum cut V\S -> S is given by the set of nodes S
 *                          having a distance less than or equal to n_nodes-1
 */

void directed_min_cut(double **x,
		              long n_nodes,
		              long sourc,
		              long sin,
		              double & val_flow,
		              long *& dist)
{
  node *j;
  int  cc;

  globUpdtFreq = GLOB_UPDT_FREQ;

  parse( x, n_nodes, sourc, sin, &n, &m, &nodes, &arcs, &cap, &source, &sink, &nMin );

  cc = allocDS();

  init();
  stageOne ( );

  val_flow = flow;

  stageTwo ( );

  globalUpdate();
  forAllNodes(j) {
    dist[nNode(j)] = j->d;
  }


  free(nodes-nMin);
  free(arcs);
  free(cap);
  free(queue);
  free(buckets);

}
