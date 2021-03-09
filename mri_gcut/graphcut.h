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

#define NODE_BLOCK_SIZE 512
#define ARC_BLOCK_SIZE 1024
#define NODEPTR_BLOCK_SIZE 128

template <class Type> class Block
{
public:
  /* Constructor. Arguments are the block size and
     (optionally) the pointer to the function which
     will be called if allocation failed; the message
     passed to this function is "Not enough memory!" */
  Block(int size, void (*err_function)(char *) = NULL)
  {
    first = last = NULL;
    block_size = size;
    error_function = err_function;
  }

  /* Destructor. Deallocates all items added so far */
  ~Block()
  {
    while (first)
    {
      block *next = first -> next;
      delete[] first;
      first = next;
    }
  }

  /* Allocates 'num' consecutive items; returns pointer
     to the first item. 'num' cannot be greater than the
     block size since items must fit in one block */
  Type *New(int num = 1)
  {
    Type *t;

    if (!last || last->current + num > last->last)
    {
      if (last && last->next) last = last -> next;
      else
      {
        block *next = 
          (block *) new char [sizeof(block) + (block_size-1)*sizeof(Type)];
        if (!next)
        {
          if (error_function) 
            (*error_function)((char*)"Not enough memory!");
          exit(1);
        }
        if (last) last -> next = next;
        else first = next;
        last = next;
        last -> current = & ( last -> data[0] );
        last -> last = last -> current + block_size;
        last -> next = NULL;
      }
    }

    t = last -> current;
    last -> current += num;
    return t;
  }

  /* Returns the first item (or NULL, if no items were added) */
  Type *ScanFirst()
  {
    scan_current_block = first;
    if (!scan_current_block) return NULL;
    scan_current_data = & ( scan_current_block -> data[0] );
    return scan_current_data ++;
  }

  /* Returns the next item (or NULL, if all items have been read)
     Can be called only if previous ScanFirst() or ScanNext()
     call returned not NULL. */
  Type *ScanNext()
  {
    if (scan_current_data >= scan_current_block -> current)
    {
      scan_current_block = scan_current_block -> next;
      if (!scan_current_block) return NULL;
      scan_current_data = & ( scan_current_block -> data[0] );
    }
    return scan_current_data ++;
  }

  /* Marks all elements as empty */
  void Reset()
  {
    block *b;
    if (!first) return;
    for (b=first; ; b=b->next)
    {
      b -> current = & ( b -> data[0] );
      if (b == last) break;
    }
    last = first;
  }

  /***********************************************************************/

private:

  typedef struct block_st
  {
    Type     *current, *last;
    struct block_st   *next;
    Type     data[1];
  }
  block;

  int  block_size;
  block *first;
  block *last;

  block *scan_current_block;
  Type *scan_current_data;

  void (*error_function)(char *);
};

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

template <class Type> class DBlock
{
public:
  /* Constructor. Arguments are the block size and
     (optionally) the pointer to the function which
     will be called if allocation failed; the message
     passed to this function is "Not enough memory!" */
  DBlock(int size, void (*err_function)(char *) = NULL)
  {
    first = NULL;
    first_free = NULL;
    block_size = size;
    error_function = err_function;
  }

  /* Destructor. Deallocates all items added so far */
  ~DBlock()
  {
    while (first)
    {
      block *next = first -> next;
      delete[] first;
      first = next;
    }
  }

  /* Allocates one item */
  Type *New()
  {
    block_item *item;

    if (!first_free)
    {
      block *next = first;
      first = 
        (block *) new char [sizeof(block) + (block_size-1)*sizeof(block_item)];
      if (!first)
      {
        if (error_function) 
          (*error_function)((char*)"Not enough memory!");
        exit(1);
      }
      first_free = & (first -> data[0] );
      for (item=first_free; item<first_free+block_size-1; item++)
        item -> next_free = item + 1;
      item -> next_free = NULL;
      first -> next = next;
    }

    item = first_free;
    first_free = item -> next_free;
    return (Type *) item;
  }

  /* Deletes an item allocated previously */
  void Delete(Type *t)
  {
    ((block_item *) t) -> next_free = first_free;
    first_free = (block_item *) t;
  }

  /***********************************************************************/

private:

  typedef union block_item_st
  {
    Type   t;
    block_item_st *next_free;
  } block_item;

  typedef struct block_st
  {
    struct block_st   *next;
    block_item    data[1];
  }
  block;

  int   block_size;
  block  *first;
  block_item *first_free;

  void (*error_function)(char *);
};




/*
 Nodes, arcs and pointers to nodes are
 added in blocks for memory and time efficiency.
 Below are numbers of items in blocks
*/
class Graph
{
public:
  typedef enum
  {
    SOURCE = 0,
    SINK = 1
  } termtype; /* terminals */

  /* Type of edge weights.
     Can be changed to char, int, float, double, ... */
  typedef short captype;
  /* Type of total flow */
  typedef int flowtype;

  typedef void * node_id;

  /* interface functions */

  /* Constructor. Optional argument is the pointer to the
     function which will be called if an error occurs;
     an error message is passed to this function. If this
     argument is omitted, exit(1) will be called. */
  Graph(void (*err_function)(char *) = NULL);

  /* Destructor */
  ~Graph();

  /* Adds a node to the graph */
  node_id add_node();

  /* Adds a bidirectional edge between 'from' and 'to'
     with the weights 'cap' and 'rev_cap' */
  void add_edge(node_id from, node_id to, captype cap, captype rev_cap);

  /* Sets the weights of the edges 'SOURCE->i' and 'i->SINK'
     Can be called at most once for each node.
     If this function is not called for 'i', then these edges are not present */
  void set_tweights(node_id i, captype cap_source, captype cap_sink);

  /* After the maxflow is computed, this function returns to which
     segment the node 'i' belongs (Graph::SOURCE or Graph::SINK) */
  termtype what_segment(node_id i);

  /* Computes the maxflow. Can be called only once. */
  flowtype maxflow();

  /***********************************************************************/
  /***********************************************************************/
  /***********************************************************************/

private:
  /* internal variables and functions */

  struct arc_st;

  /* node structure */
  typedef struct node_st
  {
    arc_st   *first;  /* first outcoming arc */

    arc_st   *parent; /* node's parent */
    node_st   *next;  /* pointer to the next active node
                    (or to itself if it is the last node in the list) */
    int    mark_count; /* d is valid if mark_count == ::mark_count */
    int    mark_d;  /* distance to the terminal */
    short   is_sink; /* flag showing whether the node is in the source or in the sink tree */

    captype   tr_cap;  /* if tr_cap > 0 then tr_cap is residual capacity of the arc SOURCE->node
                    otherwise         -tr_cap is residual capacity of the arc node->SINK */
  }
  node;

  /* arc structure */
  typedef struct arc_st
  {
    node_st   *head;  /* node the arc points to */
    arc_st   *next;  /* next arc with the same originating node */
    arc_st   *sister; /* reverse arc */

    captype   r_cap;  /* residual capacity */
  }
  arc;

  /* 'pointer to node' structure */
  typedef struct nodeptr_st
  {
    node_st   *ptr;
    nodeptr_st  *next;
  }
  nodeptr;

  Block<node>   *node_block;
  Block<arc>   *arc_block;
  DBlock<nodeptr>  *nodeptr_block;

  void (*error_function)(char *); /* this function is called if a error occurs,
                 with a corresponding error message
                 (or exit(1) is called if it's NULL) */

  flowtype   flow;  /* total flow */

  /***********************************************************************/

  node    *queue_first[2], *queue_last[2]; /* list of active nodes */
  nodeptr    *orphan_first, *orphan_last;  /* list of pointers to orphans */
  int     mark_count;       /* monotonically increasing global counter */

  /***********************************************************************/

  /* functions for processing active list */
  void set_active(node *i);
  node *next_active();

  void maxflow_init();
  void augment(arc *middle_arc);
  void process_source_orphan(node *i);
  void process_sink_orphan(node *i);
};
