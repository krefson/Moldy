/*
 * Data structures and declarations for linked lists from 
 * Paul Davies, "The Indispensable Guide to C", Addison-Wesley, UK, 1995
 */


#define   NEW(x)      (x *)malloc(sizeof(x))
#define   VALID(x)   ((*x) != NULL)
#define   NUM(x)     ((*x)->num)

/*
 * The ROOT and NODE data structures for a linked list
 */

typedef struct node {
   void         *data;           /* pointer to data at this node */
   struct node *next;            /* pointer to next node in list */
} NODE;                          /* a generic pointer to any type of data */

typedef struct {
  long          num;             /* number of nodes in list */
  NODE          *head;           /* pointer to first node in list */
  NODE          *tail;           /* pointer to last node in list */
} ROOT;

extern ROOT *make_root(void);
extern NODE *make_node(void *data);
extern int insert_data(ROOT **root, void *data, int position);
extern NODE *find_previous_node(ROOT **root, NODE *this_node);
extern int insert_at_position(ROOT **root, NODE *node, void *data, int position);
extern NODE *select_node(ROOT **root, int n);
extern int delete_list(ROOT **root);
