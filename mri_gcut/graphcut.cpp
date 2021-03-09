/*
 * Original Author: Vitali Zagorodnov, ZHU Jiaqi (September, 2009)
 *
 * Copyright © 2009 Nanyang Technological University, Singapore
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

#include "graphcut.h"

#define TERMINAL ( (arc *) 1 )  /* to terminal */
#define ORPHAN   ( (arc *) 2 )  /* orphan */
#define INFINITE_D 1000000000  /* infinite distance to the terminal */

#define NEW_VECTOR(X,N,TYPE,MSG) Test_Pointer(X = new TYPE[N], MSG)

#define max_graph_size 20000000

struct edgeW
{
  int edge1;
  int edge2;
  int weight;
  int fsw;
  int bsw;
};

//no include
extern bool matrix_alloc(int ****pointer, int z, int y, int x);
extern bool matrix_free(int ***pointer, int z, int y, int x);



inline void Graph::set_active(node *i)
{
  if (!i->next)
  {
    /* it's not in the list yet */
    if (queue_last[1]) queue_last[1] -> next = i;
    else               queue_first[1]        = i;
    queue_last[1] = i;
    i -> next = i;
  }
}

/*
 Returns the next active node.
 If it is connected to the sink, it stays in the list,
 otherwise it is removed from the list
*/
inline Graph::node * Graph::next_active()
{
  node *i;

  while ( 1 )
  {
    if (!(i=queue_first[0]))
    {
      queue_first[0] = i = queue_first[1];
      queue_last[0]  = queue_last[1];
      queue_first[1] = NULL;
      queue_last[1]  = NULL;
      if (!i) return NULL;
    }

    /* remove it from the active list */
    if (i->next == i) queue_first[0] = queue_last[0] = NULL;
    else              queue_first[0] = i -> next;
    i -> next = NULL;

    /* a node in the list is active iff it has a parent */
    if (i->parent) return i;
  }
}

/***********************************************************************/

void Graph::maxflow_init()
{
  node *i;

  queue_first[0] = queue_last[0] = NULL;
  queue_first[1] = queue_last[1] = NULL;
  orphan_first = NULL;

  for (i=node_block->ScanFirst(); i; i=node_block->ScanNext())
  {
    i -> next = NULL;
    i -> mark_count = 0;
    i -> is_sink = 0;
    if (i->tr_cap > 0)
    {
      /* i is connected to the source */
      i -> is_sink = 0;
      i -> parent = TERMINAL;
      set_active(i);
      i -> mark_count = 1;
      i -> mark_d = 1;
    }
    else if (i->tr_cap < 0)
    {
      /* i is connected to the sink */
      i -> is_sink = 1;
      i -> parent = TERMINAL;
      set_active(i);
      i -> mark_count = 1;
      i -> mark_d = 1;
    }
    else
    {
      i -> parent = NULL;
    }
  }
  mark_count = 2;
}

/***********************************************************************/

void Graph::augment(arc *middle_arc)
{
  node *i;
  arc *a;
  captype bottleneck;
  nodeptr *np;


  /* 1. Finding bottleneck capacity */
  /* 1a - the source tree */
  bottleneck = middle_arc -> r_cap;
  for (i=middle_arc->sister->head; ; i=a->head)
  {
    a = i -> parent;
    if (a == TERMINAL) break;
    if (bottleneck > a->sister->r_cap) bottleneck = a -> sister -> r_cap;
  }
  if (bottleneck > i->tr_cap) bottleneck = i -> tr_cap;
  /* 1b - the sink tree */
  for (i=middle_arc->head; ; i=a->head)
  {
    a = i -> parent;
    if (a == TERMINAL) break;
    if (bottleneck > a->r_cap) bottleneck = a -> r_cap;
  }
  if (bottleneck > - i->tr_cap) bottleneck = - i -> tr_cap;


  /* 2. Augmenting */
  /* 2a - the source tree */
  middle_arc -> sister -> r_cap += bottleneck;
  middle_arc -> r_cap -= bottleneck;
  for (i=middle_arc->sister->head; ; i=a->head)
  {
    a = i -> parent;
    if (a == TERMINAL) break;
    a -> r_cap += bottleneck;
    a -> sister -> r_cap -= bottleneck;
    if (!a->sister->r_cap)
    {
      /* add i to the adoption list */
      i -> parent = ORPHAN;
      np = nodeptr_block -> New();
      np -> ptr = i;
      np -> next = orphan_first;
      orphan_first = np;
    }
  }
  i -> tr_cap -= bottleneck;
  if (!i->tr_cap)
  {
    /* add i to the adoption list */
    i -> parent = ORPHAN;
    np = nodeptr_block -> New();
    np -> ptr = i;
    np -> next = orphan_first;
    orphan_first = np;
  }
  /* 2b - the sink tree */
  for (i=middle_arc->head; ; i=a->head)
  {
    a = i -> parent;
    if (a == TERMINAL) break;
    a -> sister -> r_cap += bottleneck;
    a -> r_cap -= bottleneck;
    if (!a->r_cap)
    {
      /* add i to the adoption list */
      i -> parent = ORPHAN;
      np = nodeptr_block -> New();
      np -> ptr = i;
      np -> next = orphan_first;
      orphan_first = np;
    }
  }
  i -> tr_cap += bottleneck;
  if (!i->tr_cap)
  {
    /* add i to the adoption list */
    i -> parent = ORPHAN;
    np = nodeptr_block -> New();
    np -> ptr = i;
    np -> next = orphan_first;
    orphan_first = np;
  }


  flow += bottleneck;
}

/***********************************************************************/

void Graph::process_source_orphan(node *i)
{
  node *j;
  arc *a0, *a0_min = NULL, *a;
  nodeptr *np;
  int d, d_min = INFINITE_D;

  /* trying to find a new parent */
  for (a0=i->first; a0; a0=a0->next)
    if (a0->sister->r_cap)
    {
      j = a0 -> head;
      if (!j->is_sink && (a=j->parent))
      {
        /* checking the origin of j */
        d = 0;
        while ( 1 )
        {
          if (j->mark_count == mark_count)
          {
            d += j -> mark_d;
            break;
          }
          a = j -> parent;
          d ++;
          if (a==TERMINAL)
          {
            j -> mark_count = mark_count;
            j -> mark_d = 1;
            break;
          }
          if (a==ORPHAN)
          {
            d = INFINITE_D;
            break;
          }
          j = a -> head;
        }
        if (d<INFINITE_D) /* j originates from the source - done */
        {
          if (d<d_min)
          {
            a0_min = a0;
            d_min = d;
          }
          /* set marks along the path */
          for (j=a0->head; j->mark_count!=mark_count; j=j->parent->head)
          {
            j -> mark_count = mark_count;
            j -> mark_d = d --;
          }
        }
      }
    }

  if ((i->parent = a0_min))
  {
    i -> mark_count = mark_count;
    i -> mark_d = d_min + 1;
  }
  else
  {
    /* no parent is found */
    i -> mark_count = 0;

    /* process neighbors */
    for (a0=i->first; a0; a0=a0->next)
    {
      j = a0 -> head;
      if (!j->is_sink && (a=j->parent))
      {
        if (a0->sister->r_cap) set_active(j);
        if (a!=TERMINAL && a!=ORPHAN && a->head==i)
        {
          /* add j to the adoption list */
          j -> parent = ORPHAN;
          np = nodeptr_block -> New();
          np -> ptr = j;
          if (orphan_last) orphan_last -> next = np;
          else             orphan_first        = np;
          orphan_last = np;
          np -> next = NULL;
        }
      }
    }
  }
}

void Graph::process_sink_orphan(node *i)
{
  node *j;
  arc *a0, *a0_min = NULL, *a;
  nodeptr *np;
  int d, d_min = INFINITE_D;

  /* trying to find a new parent */
  for (a0=i->first; a0; a0=a0->next)
    if (a0->r_cap)
    {
      j = a0 -> head;
      if ((j->is_sink) && (a=j->parent))
      {
        /* checking the origin of j */
        d = 0;
        while ( 1 )
        {
          if (j->mark_count == mark_count)
          {
            d += j -> mark_d;
            break;
          }
          a = j -> parent;
          d ++;
          if (a==TERMINAL)
          {
            j -> mark_count = mark_count;
            j -> mark_d = 1;
            break;
          }
          if (a==ORPHAN)
          {
            d = INFINITE_D;
            break;
          }
          j = a -> head;
        }
        if (d<INFINITE_D) /* j originates from the sink - done */
        {
          if (d<d_min)
          {
            a0_min = a0;
            d_min = d;
          }
          /* set marks along the path */
          for (j=a0->head; j->mark_count!=mark_count; j=j->parent->head)
          {
            j -> mark_count = mark_count;
            j -> mark_d = d --;
          }
        }
      }
    }

  if ((i->parent = a0_min))
  {
    i -> mark_count = mark_count;
    i -> mark_d = d_min + 1;
  }
  else
  {
    /* no parent is found */
    i -> mark_count = 0;

    /* process neighbors */
    for (a0=i->first; a0; a0=a0->next)
    {
      j = a0 -> head;
      if (j->is_sink && (a=j->parent))
      {
        if (a0->r_cap) set_active(j);
        if (a!=TERMINAL && a!=ORPHAN && a->head==i)
        {
          /* add j to the adoption list */
          j -> parent = ORPHAN;
          np = nodeptr_block -> New();
          np -> ptr = j;
          if (orphan_last) orphan_last -> next = np;
          else             orphan_first        = np;
          orphan_last = np;
          np -> next = NULL;
        }
      }
    }
  }
}

/***********************************************************************/

Graph::flowtype Graph::maxflow()
{
  node *i = NULL, *j = NULL, *current_node = NULL;
  arc *a = NULL;
  nodeptr *np = NULL, *np_next = NULL;

  maxflow_init();
  nodeptr_block = new DBlock<nodeptr>(NODEPTR_BLOCK_SIZE, error_function);

  while ( 1 )
  {
    if ((i=current_node))
    {
      i -> next = NULL; /* remove active flag */
      if (!i->parent) i = NULL;
    }
    if (!i)
    {
      if (!(i = next_active())) break;
    }

    /* growth */
    if (!i->is_sink)
    {
      /* grow source tree */
      for (a=i->first; a; a=a->next)
        if (a->r_cap)
        {
          j = a -> head;
          if (!j->parent)
          {
            j -> is_sink = 0;
            j -> parent = a -> sister;
            j -> mark_count = i -> mark_count;
            j -> mark_d = i -> mark_d + 1;
            set_active(j);
          }
          else if (j->is_sink) break;
          else if (j->mark_count &&
                   j->mark_count <= i->mark_count &&
                   j->mark_d > i->mark_d)
          {
            /* heuristic - trying to make the distance from 
               j to the source shorter */
            j -> parent = a -> sister;
            j -> mark_count = i -> mark_count;
            j -> mark_d = i -> mark_d + 1;
          }
        }
    }
    else
    {
      /* grow sink tree */
      for (a=i->first; a; a=a->next)
        if (a->sister->r_cap)
        {
          j = a -> head;
          if (!j->parent)
          {
            j -> is_sink = 1;
            j -> parent = a -> sister;
            j -> mark_count = i -> mark_count;
            j -> mark_d = i -> mark_d + 1;
            set_active(j);
          }
          else if (!j->is_sink)
          {
            a = a -> sister;
            break;
          }
          else if (j->mark_count &&
                   j->mark_count <= i->mark_count &&
                   j->mark_d > i->mark_d)
          {
            /* heuristic - trying to make the distance 
               from j to the sink shorter */
            j -> parent = a -> sister;
            j -> mark_count = i -> mark_count;
            j -> mark_d = i -> mark_d + 1;
          }
        }
    }

    if (a)
    {
      i -> next = i; /* set active flag */
      current_node = i;

      /* augmentation */
      augment(a);
      /* augmentation end */

      /* adoption */
      while ((np=orphan_first))
      {
        np_next = np -> next;
        np -> next = NULL;

        while ((np=orphan_first))
        {
          orphan_first = np -> next;
          i = np -> ptr;
          nodeptr_block -> Delete(np);
          if (!orphan_first) orphan_last = NULL;
          if (i->is_sink) process_sink_orphan(i);
          else            process_source_orphan(i);
        }

        orphan_first = np_next;
      }
      mark_count ++;
      /* adoption end */
    }
    else current_node = NULL;
  }

  delete nodeptr_block;

  return flow;
}

/***********************************************************************/

Graph::termtype Graph::what_segment(node_id i)
{
  if (((node*)i)->parent && !((node*)i)->is_sink) return SOURCE;
  return SINK;
}


Graph::Graph(void (*err_function)(char *))
{
  error_function = err_function;
  node_block = new Block<node>(NODE_BLOCK_SIZE, error_function);
  arc_block  = new Block<arc>(NODE_BLOCK_SIZE, error_function);
  flow = 0;
}

Graph::~Graph()
{
  delete node_block;
  delete arc_block;
}

Graph::node_id Graph::add_node()
{
  node *i = node_block -> New();

  i -> first = NULL;
  i -> tr_cap = 0;

  return (node_id) i;
}

void Graph::add_edge(node_id from, node_id to, captype cap, captype rev_cap)
{
  arc *a, *a_rev;

  a = arc_block -> New(2);
  a_rev = a + 1;

  a -> sister = a_rev;
  a_rev -> sister = a;
  a -> next = ((node*)from) -> first;
  ((node*)from) -> first = a;
  a_rev -> next = ((node*)to) -> first;
  ((node*)to) -> first = a_rev;
  a -> head = (node*)to;
  a_rev -> head = (node*)from;
  a -> r_cap = cap;
  a_rev -> r_cap = rev_cap;
}

void Graph::set_tweights(node_id i, captype cap_source, captype cap_sink)
{
  flow += (cap_source < cap_sink) ? cap_source : cap_sink;
  ((node*)i) -> tr_cap = cap_source - cap_sink;
}


void Test_Pointer(void * p, char * _or)
{
  if (p == NULL)
  {
    printf("\n\n\a -> Error: insufficient memory.");
    if (_or != NULL) printf("\nOrigin = %s\n", _or);
    printf("\n Execution terminated!\n");
    exit(1);
  }
}

void map(int index, int *x, int *y, int *z, int xV, int yV, int zV)
{
  *z = index / (yV * xV);
  *y = (index - (*z) * yV * xV) / xV;
  *x = index - (*z) * yV * xV - (*y) * xV;
}

void map(int index, int *x, int *y, int *z, 
         int xV, int yV, int zV, 
         int xoff, int yoff, int zoff)
{
  *z = index / (yV * xV);
  *y = (index - (*z) * yV * xV) / xV;
  *x = index - (*z) * yV * xV - (*y) * xV;
  //plus offset
  *z = *z + zoff;
  *y = *y + yoff;
  *x = *x + xoff;
}

void do_cityblock(int ***pointer, int zVol, int yVol, int xVol)
{
  //set a boundary
  for (int z = 0; z < zVol; z++)
  {
    for (int y = 0; y < yVol; y++)
    {
      for (int x = 0; x < xVol; x++)
      {
        if ( x == 0 || x == xVol-1 || 
             y == 0 || y == yVol-1 || 
             z == 0 || z == zVol-1 )
          pointer[z][y][x] = 0;
      }
    }
  }
  //fill
  for (int time = 0; time <= 19; time++)
  {
    for (int z = 1; z < zVol-1; z++)
    {
      for (int y = 1; y < yVol-1; y++)
      {
        for (int x = 1; x < xVol-1; x++)
        {
          if ( pointer[z][y][x] == -1 )//need to do
          {
            if ( pointer[z][y][x-1]==time || pointer[z][y][x+1]==time
                 || pointer[z][y-1][x]==time || pointer[z][y+1][x]==time
                 || pointer[z-1][y][x]==time || pointer[z+1][y][x]==time )
              pointer[z][y][x] = time + 1;
          }//end of if
        }
      }
    }
  }
  //fill rest
  for (int z = 1; z < zVol-1; z++)
  {
    for (int y = 1; y < yVol-1; y++)
    {
      for (int x = 1; x < xVol-1; x++)
      {
        if ( pointer[z][y][x] == -1 )//need to do
        {
          pointer[z][y][x] = 21;
        }//end of if
      }
    }
  }
}

void mincut(int ***im_gcut, edgeW *hor, edgeW *ver, edgeW *tra, 
            int length_h, int length_v, int length_t, 
            int x_start, int y_start, int z_start, 
            int x_end, int y_end, int z_end)
{
  int i, k, x, y, z, m1, m2, m3;
  int *marked;

  x = x_end - x_start + 1;
  y = y_end - y_start + 1;
  z = z_end - z_start + 1;
  m1 = length_h;
  m2 = length_v;
  m3 = length_t;

  /* if ((m1 >= m2) && (m1 >= m3))
    m = m1;
   if ((m2 >= m1) && (m2 >= m3))
    m = m2;
   if ((m3 >= m1) && (m3 >= m2))
    m = m3;
  */
  // Flag for the indices in the matrix.
  int sizeofMarked = x * y * z;
  marked = new int[sizeofMarked];
  for (i=0; i<sizeofMarked; i++) marked[i]=0;

  // Graph construction
  Graph::node_id *nodes = new Graph::node_id[x * y * z + 1];
  //Graph::node_id nodes[max_graph_size];
  Graph *g = new Graph();

  // Create the first node
  k = 0;
  nodes[k] = g -> add_node();
  g -> set_tweights(nodes[k], 400, 0); //set_tweights(node, source_weight, sink_weight);
  marked[0] = 1;
  k = k + 1;

  //edgewts1 = (float **) calloc(m, sizeof(float *));
  //for (i = 1; i <= m; i++)
  // edgewts1[i] = (float *) calloc(n, sizeof(float));

  /* for (i = 1; i <= m1; i++)
   {
    for (j = 1; j <= n; j++)
    {
     fscanf (fp,"%f",&edgewts1[i][j]);
    }
   }
  */

  // Create nodes and edges (from right to left)
  for (i = 0; i < m1; i++)
  {
    if (marked[hor[i].edge1 - 1] == 0)
    {
      marked[hor[i].edge1 - 1] = 1;
      nodes[k] = g -> add_node();
      g -> set_tweights(nodes[k], hor[i].bsw, hor[i].fsw);
      k = k + 1;
    }

    if (marked[hor[i].edge2 - 1] == 0)
    {
      marked[hor[i].edge2 - 1] = 1;
      nodes[k] = g -> add_node();
      g -> set_tweights(nodes[k], hor[i].bsw, hor[i].fsw);
      k = k + 1;
    }

    g -> add_edge(nodes[hor[i].edge2 - 1], 
                  nodes[hor[i].edge1 - 1], 
                  hor[i].weight, 
                  hor[i].weight);
//  if( i % 10000 == 0 )
//    printf("doing hor:%d; k=%d\n", i, k);
  }


  /* for (i = 1; i <= m2; i++)
   {
    for (j = 1; j <= n; j++)
    {
     fscanf (fp,"%f",&edgewts1[i][j]);
    }
   }
  */
  // Create nodes and edges (from bottom to top)
  for (i = 0; i < m2; i++)
  {
    if (marked[ver[i].edge1 - 1] == 0)
    {
      marked[ver[i].edge1 - 1] = 1;
      nodes[k] = g -> add_node();
      g -> set_tweights(nodes[k], ver[i].bsw, ver[i].fsw);
      k = k + 1;
    }

    if (marked[ver[i].edge2 - 1] == 0)
    {
      marked[ver[i].edge2 - 1] = 1;
      nodes[k] = g -> add_node();
      g -> set_tweights(nodes[k], ver[i].bsw, ver[i].fsw);
      k = k + 1;
    }

    g -> add_edge(nodes[ver[i].edge2 - 1], 
                  nodes[ver[i].edge1 - 1], 
                  ver[i].weight, 
                  ver[i].weight);
//  if( i % 10000 == 0 )
//    printf("doing ver:%d; k=%d\n", i, k);
  }

  /* for (i = 1; i <= m3; i++)
   {
    for (j = 1; j <= n; j++)
    {
     fscanf (fp,"%f",&edgewts1[i][j]);
    }
   }
  */
  // Create nodes and edges (from back to front)
  for (i = 0; i < m3; i++)
  {
    if (marked[tra[i].edge1 - 1] == 0)
    {
      marked[tra[i].edge1 - 1] = 1;
      nodes[k] = g -> add_node();
      g -> set_tweights(nodes[k], tra[i].bsw, tra[i].fsw);
      k = k + 1;
    }

    if (marked[tra[i].edge2 - 1] == 0)
    {
      marked[tra[i].edge2 - 1] = 1;
      nodes[k] = g -> add_node();
      g -> set_tweights(nodes[k], tra[i].bsw, tra[i].fsw);
      k = k + 1;
    }

    g -> add_edge(nodes[tra[i].edge2 - 1], 
                  nodes[tra[i].edge1 - 1], 
                  tra[i].weight, 
                  tra[i].weight);
//  if( i % 10000 == 0 )
//    printf("doing tra:%d; k=%d\n", i, k);
  }


  //Graph::flowtype flow = g -> maxflow();
  printf("now doing maxflow, be patient...\n");
  g -> maxflow();

  // Write the output of mincut to a file
  //fp = fopen("output.txt","w");
  //for (i = 0; i < k; i++)
  // fprintf(fp,"%d \n", g->what_segment(nodes[i]));

  i = 0;
  for (int z = z_start; z <= z_end; z++)
  {
    for (int y = y_start; y <= y_end; y++)
    {
      for (int x = x_start; x <= x_end; x++)
      {
        im_gcut[z][y][x] = g->what_segment(nodes[i]);
        i++;
      }
    }
  }

  delete[] marked;
  delete g;
  delete[] nodes;
}

void decide_bound(unsigned char ***image, double threshold, 
                  int xVol, int yVol, int zVol, 
                  int & x_start, int & x_end, 
                  int & y_start, int & y_end, 
                  int & z_start, int & z_end)
{
  x_start = y_start = z_start = 0;
  x_end = xVol - 1;
  y_end = yVol - 1;
  z_end = zVol - 1;

  double factor = 1.2;
  //z_start
  int bExit = 0;
  for (int z = 0; z < zVol && bExit == 0; z++)
  {
    for (int y = 0; y < yVol && bExit == 0; y++)
    {
      for (int x = 0; x < xVol && bExit == 0; x++)
      {
        if ( image[z][y][x] > threshold*factor )
        {
          z_start = z;
          bExit = 1;
          break;
        }
      }
    }
  }
  //z_end
  bExit = 0;
  for (int z = zVol - 1; z > z_start && bExit == 0; z--)
  {
    for (int y = 0; y < yVol && bExit == 0; y++)
    {
      for (int x = 0; x < xVol && bExit == 0; x++)
      {
        if ( image[z][y][x] > threshold*factor )
        {
          z_end = z;
          bExit = 1;
          break;
        }
      }
    }
  }
  //y_start
  bExit = 0;
  for (int y = 0; y < yVol && bExit == 0; y++)
  {
    for (int z = 0; z < zVol && bExit == 0; z++)
    {
      for (int x = 0; x < xVol && bExit == 0; x++)
      {
        if ( image[z][y][x] > threshold*factor )
        {
          y_start = y;
          bExit = 1;
          break;
        }
      }
    }
  }
  //y_end
  bExit = 0;
  for (int y = yVol - 1; y > y_start && bExit == 0; y--)
  {
    for (int z = 0; z < zVol && bExit == 0; z++)
    {
      for (int x = 0; x < xVol && bExit == 0; x++)
      {
        if ( image[z][y][x] > threshold*factor )
        {
          y_end = y;
          bExit = 1;
          break;
        }
      }
    }
  }
  //x_start
  bExit = 0;
  for (int x = 0; x < xVol && bExit == 0; x++)
  {
    for (int z = 0; z < zVol && bExit == 0; z++)
    {
      for (int y = 0; y < yVol && bExit == 0; y++)
      {
        if ( image[z][y][x] > threshold*factor )
        {
          x_start = x;
          bExit = 1;
          break;
        }
      }
    }
  }
  //x_end
  bExit = 0;
  for (int x = xVol - 1; x > x_start && bExit == 0; x--)
  {
    for (int z = 0; z < zVol && bExit == 0; z++)
    {
      for (int y = 0; y < yVol && bExit == 0; y++)
      {
        if ( image[z][y][x] >= threshold*factor )
        {
          x_end = x;
          bExit = 1;
          break;
        }
      }
    }
  }
}

double graphcut(unsigned char ***image, unsigned char ***label, 
                int ***im_gcut, int ***foreSW, int ***backSW, 
                int xVol, int yVol, int zVol, 
                double kval, double threshold, double whitemean)
{
  //city block matrix
  int ***cityblock;
  matrix_alloc(&cityblock, zVol, yVol, xVol);

  //fore groud seed & back ground seed
  for (int z = 0; z < zVol; z++)
  {
    for (int y = 0; y < yVol; y++)
    {
      for (int x = 0; x < xVol; x++)
      {
        if ( label[z][y][x] == 1 )
        {
          foreSW[z][y][x] = 4000;
        }
        if ( image[z][y][x] <= 0 )
        {
          backSW[z][y][x] = 4000;
        }
        if ( image[z][y][x] > 0 )
        {
          cityblock[z][y][x] = -1;
        }
      }
    }
  }

  do_cityblock(cityblock, zVol, yVol, xVol);

  int x_start, x_end, y_start, y_end, z_start, z_end;
  decide_bound(image, threshold, 
               xVol, yVol, zVol, 
               x_start, x_end, y_start, y_end, z_start, z_end);
  //new range
  int xVol_new = x_end - x_start + 1;
  int yVol_new = y_end - y_start + 1;
  int zVol_new = z_end - z_start + 1;
  //printf("x_new=%d, y_new=%d, z_new=%d", xVol_new, yVol_new, zVol_new);

  double k = kval / (whitemean - threshold);
  //assign memory
  int length_h = (xVol_new-1)*yVol_new*zVol_new;
  int length_v = xVol_new*(yVol_new-1)*zVol_new;
  int length_t = xVol_new*yVol_new*(zVol_new-1);
  edgeW *hor = new edgeW[length_h];
  edgeW *ver = new edgeW[length_v];
  edgeW *tra = new edgeW[length_t];
  memset(hor,0,sizeof(edgeW[length_h]));
  memset(ver,0,sizeof(edgeW[length_v]));
  memset(tra,0,sizeof(edgeW[length_t]));

  int i, j;
  int x, y, z;
  printf("calculating weights...\n");
  //hor
  for (i = 0, j = 0; i < length_h; j++)
  {
    map(j, &x, &y, &z, xVol_new, yVol_new, zVol_new);
    if ( x+1 == xVol_new )
      continue;
    hor[i].edge1 = j + 1;
    hor[i].edge2 = j + 1 + 1;
    //weight
    //map(hor[i].edge1-1, &x, &y, &z, xVol_new, yVol_new, zVol_new);
    x += x_start;
    y += y_start;
    z += z_start;
    hor[i].weight = cityblock[z][y][x] > cityblock[z][y][x+1] ? 
      cityblock[z][y][x] : cityblock[z][y][x+1];
    hor[i].weight = hor[i].weight * hor[i].weight;
    if (hor[i].weight > 1 && hor[i].weight < 6)
      hor[i].weight = 6;
    if (hor[i].weight != 1 && hor[i].weight != 6 && hor[i].weight != 0)
    {
      unsigned char value = image[z][y][x] > image[z][y][x+1] ? 
        image[z][y][x+1] : image[z][y][x];
      hor[i].weight = (int)fabs(hor[i].weight * 
                                (exp(k * (value - threshold)) - 1));
    }
    if (hor[i].weight > 1 && hor[i].weight < 6)
      hor[i].weight = 6;
    if (hor[i].weight > 0 && hor[i].weight < 1)
      hor[i].weight = 1;
    if (hor[i].weight == 0)
      hor[i].weight = 1000;
    //foreground seed and background seed
    hor[i].fsw = foreSW[z][y][x+1];
    hor[i].bsw = backSW[z][y][x+1];
    i++;
  }

  //ver
  for (i = 0, j = 0; i < length_v; j++)
  {
    map(j, &x, &y, &z, xVol_new, yVol_new, zVol_new);
    if ( y+1 == yVol_new )
      continue;
    ver[i].edge1 = j + 1;
    ver[i].edge2 = j + 1 + xVol_new;
    //weight
    //map(ver[i].edge1-1, &x, &y, &z, xVol_new, yVol_new, zVol_new);
    x += x_start;
    y += y_start;
    z += z_start;
    ver[i].weight = cityblock[z][y][x] > cityblock[z][y+1][x] ? 
      cityblock[z][y][x] : cityblock[z][y+1][x];
    ver[i].weight = ver[i].weight * ver[i].weight;
    if (ver[i].weight > 1 && ver[i].weight < 6)
      ver[i].weight = 6;
    if (ver[i].weight != 1 && ver[i].weight != 6 && ver[i].weight != 0)
    {
      unsigned char value = image[z][y][x] > image[z][y+1][x] ? 
        image[z][y+1][x] : image[z][y][x];
      ver[i].weight = (int)fabs(ver[i].weight * 
                                (exp(k * (value - threshold)) - 1));
    }
    if (ver[i].weight > 1 && ver[i].weight < 6)
      ver[i].weight = 6;
    if (ver[i].weight > 0 && ver[i].weight < 1)
      ver[i].weight = 1;
    if (ver[i].weight == 0)
      ver[i].weight = 1000;
    //foreground seed and background seed
    ver[i].fsw = foreSW[z][y+1][x];
    ver[i].bsw = backSW[z][y+1][x];
    i++;
  }

  //tra
  for (i = 0, j = 0; i < length_t; j++)
  {
    map(j, &x, &y, &z, xVol_new, yVol_new, zVol_new);
    if ( z+1 == zVol_new )
      continue;
    tra[i].edge1 = j + 1;
    tra[i].edge2 = j + 1 + xVol_new * yVol_new;
    //weight
    //map(tra[i].edge1-1, &x, &y, &z, xVol_new, yVol_new, zVol_new);
    x += x_start;
    y += y_start;
    z += z_start;
    tra[i].weight = cityblock[z][y][x] > cityblock[z+1][y][x] ? 
      cityblock[z][y][x] : cityblock[z+1][y][x];
    tra[i].weight = tra[i].weight * tra[i].weight;
    if (tra[i].weight > 1 && tra[i].weight < 6)
      tra[i].weight = 6;
    if (tra[i].weight != 1 && tra[i].weight != 6 && tra[i].weight != 0)
    {
      unsigned char value = image[z][y][x] > image[z+1][y][x] ? 
        image[z+1][y][x] : image[z][y][x];
      tra[i].weight = (int)fabs(tra[i].weight * 
                                (exp(k * (value - threshold)) - 1));
    }
    if (tra[i].weight > 1 && tra[i].weight < 6)
      tra[i].weight = 6;
    if (tra[i].weight > 0 && tra[i].weight < 1)
      tra[i].weight = 1;
    if (tra[i].weight == 0)
      tra[i].weight = 1000;
    //foreground seed and background seed
    tra[i].fsw = foreSW[z+1][y][x];
    tra[i].bsw = backSW[z+1][y][x];
    i++;
  }

  // -- free memory
  matrix_free(cityblock, zVol, yVol, xVol);
  matrix_free(backSW, zVol, yVol, xVol);
  matrix_free(foreSW, zVol, yVol, xVol);
  //mincut
  printf("doing mincut...\n");
  mincut(im_gcut, hor, ver, tra, length_h, length_v, length_t, 
         x_start, y_start, z_start, x_end, y_end, z_end);

  delete[] hor;
  delete[] ver;
  delete[] tra;
  return 0;
}

//conv: size 3
int calculate_conv(int ***image, int xCor, int yCor, int zCor, 
                   int conv_m[3][3][3])
{
  int sum = 0;
  for (int z = -1; z <= 1; z++)
  {
    for (int y = -1; y <= 1; y++)
    {
      for (int x = -1; x <= 1; x++)
      {
        sum = sum + image[zCor + z][yCor + y][xCor + x] * 
          conv_m[1 + z][1 + y][1 + x];
      }
    }
  }
  return sum;
}

void post_processing(unsigned char ***image, unsigned char ***image_thre, 
                     double threshold, int ***im_gcut, int ***im_dilute, 
                     int xVol, int yVol, int zVol)
{
  //for loop
  int ***temp, ***compare;
  temp = new int**[zVol];
  for (int i = 0; i < zVol; i++)
  {
    temp[i] = new int*[yVol];
    for (int j = 0; j < yVol; j++)
    {
      temp[i][j] = new int[xVol];
      for (int k = 0; k < xVol; k++)
      {
        temp[i][j][k] = im_gcut[i][j][k];
      }
    }
  }
  compare = new int**[zVol];
  for (int i = 0; i < zVol; i++)
  {
    compare[i] = new int*[yVol];
    for (int j = 0; j < yVol; j++)
    {
      compare[i][j] = new int[xVol];
      for (int k = 0; k < xVol; k++)
      {
        compare[i][j][k] = im_gcut[i][j][k];
      }
    }
  }

  //convoution matrix
  int conv_m[3][3][3];
  for (int z = 0; z < 3; z++)
  {
    for (int y = 0; y < 3; y++)
    {
      for (int x = 0; x < 3; x++)
      {
        if ( (x == 1 && y == 1) || (y == 1 && z == 1) || (z == 1 && x == 1) )
          conv_m[z][y][x] = 1;
        else
          conv_m[z][y][x] = 0;
      }
    }
  }

  int repetition = 10;
  //dilute
  for (int rep = 0; rep < repetition; rep++)
  {
    //dilution: no need to reset im_dilute to 0
    for (int z = 1; z < zVol - 1; z++)
    {
      for (int y = 1; y < yVol - 1; y++)
      {
        for (int x = 1; x < xVol - 1; x++)
        {
          int sum = calculate_conv(temp, x, y, z, conv_m);
          if (sum >= 1)
            im_dilute[z][y][x] = 1;
        }
      }
    }
    //recast to temp
    for (int z = 1; z < zVol - 1; z++)
    {
      for (int y = 1; y < yVol - 1; y++)
      {
        for (int x = 1; x < xVol - 1; x++)
        {
          temp[z][y][x] = im_dilute[z][y][x];
        }
      }
    }
  }
  //erode
  for (int rep = 0; rep < repetition; rep++)
  {
    //dilution: no need to reset im_dilute to 0
    for (int z = 1; z < zVol - 1; z++)
    {
      for (int y = 1; y < yVol - 1; y++)
      {
        for (int x = 1; x < xVol - 1; x++)
        {
          int sum = calculate_conv(temp, x, y, z, conv_m);
          if (sum == 7)
            im_dilute[z][y][x] = 1;
          else
            im_dilute[z][y][x] = 0;
        }
      }
    }
    //recast to temp
    for (int z = 1; z < zVol - 1; z++)
    {
      for (int y = 1; y < yVol - 1; y++)
      {
        for (int x = 1; x < xVol - 1; x++)
        {
          temp[z][y][x] = im_dilute[z][y][x];
        }
      }
    }
  }
  //add only dark layer
  for (int z = 1; z < zVol - 1; z++)
  {
    for (int y = 1; y < yVol - 1; y++)
    {
      for (int x = 1; x < xVol - 1; x++)
      {
        if ( compare[z][y][x] == 0 && 
             image[z][y][x] > threshold )//only dark layer
          temp[z][y][x] = 0;
      }
    }
  }
  //record compare
  for (int z = 1; z < zVol - 1; z++)
  {
    for (int y = 1; y < yVol - 1; y++)
    {
      for (int x = 1; x < xVol - 1; x++)
      {
        compare[z][y][x] = temp[z][y][x];
      }
    }
  }

  //add 1 at cut: bright
  for (int rep = 0; rep < 1; rep++)
  {
    //dilution: no need to reset im_dilute to 0
    for (int z = 1; z < zVol - 1; z++)
    {
      for (int y = 1; y < yVol - 1; y++)
      {
        for (int x = 1; x < xVol - 1; x++)
        {
          int sum = calculate_conv(temp, x, y, z, conv_m);
          if (sum >= 1)
            im_dilute[z][y][x] = 1;
        }
      }
    }
    //recast to temp
    for (int z = 1; z < zVol - 1; z++)
    {
      for (int y = 1; y < yVol - 1; y++)
      {
        for (int x = 1; x < xVol - 1; x++)
        {
          temp[z][y][x] = im_dilute[z][y][x];
          //compare[z][y][x] = im_dilute[z][y][x];
        }
      }
    }
  }
  //only bright
  for (int z = 1; z < zVol - 1; z++)
  {
    for (int y = 1; y < yVol - 1; y++)
    {
      for (int x = 1; x < xVol - 1; x++)
      {
        if ( compare[z][y][x] == 0 && 
             image[z][y][x] < threshold )//only dark layer
          temp[z][y][x] = 0;
      }
    }
  }
  //record compare
  for (int z = 1; z < zVol - 1; z++)
  {
    for (int y = 1; y < yVol - 1; y++)
    {
      for (int x = 1; x < xVol - 1; x++)
      {
        compare[z][y][x] = temp[z][y][x];
      }
    }
  }

  //add 2 layer: dark
  for (int rep = 0; rep < 2; rep++)
  {
    //dilution: no need to reset im_dilute to 0
    for (int z = 1; z < zVol - 1; z++)
    {
      for (int y = 1; y < yVol - 1; y++)
      {
        for (int x = 1; x < xVol - 1; x++)
        {
          int sum = calculate_conv(temp, x, y, z, conv_m);
          if (sum >= 1)
            im_dilute[z][y][x] = 1;
        }
      }
    }
    //recast to temp
    for (int z = 1; z < zVol - 1; z++)
    {
      for (int y = 1; y < yVol - 1; y++)
      {
        for (int x = 1; x < xVol - 1; x++)
        {
          temp[z][y][x] = im_dilute[z][y][x];
        }
      }
    }
  }
  //add only dark layer
  for (int z = 1; z < zVol - 1; z++)
  {
    for (int y = 1; y < yVol - 1; y++)
    {
      for (int x = 1; x < xVol - 1; x++)
      {
        if ( compare[z][y][x] == 0 && 
             image[z][y][x] > threshold )//only dark layer
          im_dilute[z][y][x] = 0;
      }
    }
  }

  //free
  for (int i = 0; i < zVol; i++)
  {
    for (int j = 0; j < yVol; j++)
    {
      delete[] temp[i][j];
    }
    delete[] temp[i];
  }
  delete[] temp;
  for (int i = 0; i < zVol; i++)
  {
    for (int j = 0; j < yVol; j++)
    {
      delete[] compare[i][j];
    }
    delete[] compare[i];
  }
  delete[] compare;

  return;
}
