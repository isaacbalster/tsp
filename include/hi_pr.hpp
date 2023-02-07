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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
//#include <values.h>

void directed_min_cut(double **x,
		              long n_nodes,
		              long sourc,
		              long sin,
		              double & val_flow,
		              long *& dist);

