/**************************************************************************************
 * list         Routines for creating and handling linked lists, from                 *
 *              Paul Davies, "The Indispensable Guide to C", Addison-Wesley, UK, 1995 *
 **************************************************************************************/
#include <stdio.h>
#include "list.h"

/*
 * Make a new ROOT and set the head and tail pointers to invalid (NULL) and set num
 * to zero, since there cannot be any nodes in the list. Then return a pointer to the newly
 * created ROOT structure, or return a NULL pointer if it cannot be created.
 */

ROOT *make_root(void)
{
   ROOT *root;

   if((root = NEW(ROOT)) != NULL)
   { 
     root->head = root->tail = NULL;
     root->num = 0;
   } 
   return root;
}

/*
 * Make a new node and point it to the data whose address is supplied in the call.
 * A pointer to the new node is returned or NULL on error.
 */

NODE *make_node(void *data)
{
   NODE		*node;

   if((node = NEW(NODE)) != NULL)
   {
      node->data = data;
      node->next = NULL;
   }
   return node;
} 

/*
 * The following function creates and inserts a node in the list.
 */
int insert_data(ROOT **root, void *data, int position)
{
   NODE		*new;

   if( !VALID(root))
      if((*root = make_root()) == NULL)
          return(-1);

   if((new = make_node(data)) == NULL)
   {
       if(NUM(root) == 0)
       {
          free(*root);       /* release the root structure if list is empty */
          *root = NULL;      /* and set program's root pointer to NULL */
       }
       return(-1);           /* return fail */
   }

/*
 * If this node is the first node to be inserted into the list, then set both
 * head and tail pointers in the ROOT structure to point to this new node
 */

   if(NUM(root) == 0)
      (*root)->head = (*root)->tail = new;

/*
 * Otherwise (because a node must already exist in the list), if position is equal to 0,
 * the newly created node will be inserted at the head of the list.
 */
   else
   {
      if(position == 0)
      {
         new->next = (*root)->head;       /* point new node to existing head node */
         (*root)->head = new;             /* make new node the head node */
      }

/* Else insert the new node at the tail and point the existing tail node to it */

      else
      {
         (*root)->tail->next = new;       /* point existing tail node to new node */
         (*root)->tail = new;             /* make new node the new tail node */ 
      }
   }

   (NUM(root))++;                           /* increment node count for list */
   return 0;
}
/*
 * The following function searches a linked list to obtain and return the address
 * of the node prior to that pointed to by this_node. If this_node cannot be
 * located in the list, then NULL is returned to indicate an error.
 * 
 * If this_node exists in the list, then the address of its previous node is returned,
 * or, in the case where the node is the first node in the list, then the address of
 * that first node is returned.
 */
NODE *find_previous_node(ROOT **root, NODE *this_node)
{
   NODE *previous_node, *current_node;

   current_node = previous_node = (*root)->head;
/*
 * While current_node is not pointing to the node we have asked it to search for
 * and it is not pointing to the last node in the list then move current_node
 * on to the next node, saving the address of the previous node
 */
   while((current_node != this_node) && (current_node->next != NULL))
   {
      previous_node = current_node;
      current_node = current_node->next;
   }

   return (previous_node);
}

/*
 * This function inserts a node prior to or subsequent to the node pointed to by
 * 'node'. If position is equal to 0, then a new node is created and inserted
 * prior to the node pointed to by 'node', otherwise it will be inserted after
 * that node. 'data' must point to the object to be inserted into the list.
 */
int insert_at_position(ROOT **root, NODE *node, void *data, int position)
{
   int   flag = 0;
   NODE  *new, *prev;

   if( !VALID(root) )
      return(-1);

   if(((*root)->head == node) && (position == 0))
      flag = insert_data(root, data, 0);
 
   else
      if(((*root)->tail == node) && (position == 1))
         flag = insert_data(root, data, 1);
      else
      {
         if(( prev = find_previous_node(root, node)) == NULL )
             return(-1);

         if(( new = make_node(data)) == NULL )
             return(-1);

         if( position == 1)
         {
            new->next = node->next;
            node->next = new;
         }
         else
         {
            new->next = prev->next;
            prev->next = new;
         }
         NUM(root)++;
      }
   return flag;
}

/*
 * This function returns a pointer to the nth node of the linked list
 */

NODE *select_node(ROOT **root, int n)
{
   NODE		*node;
   int		i;

   if( (VALID(root)) && (n <= NUM(root)) )
   {
      node = (*root)->head;
      for( i = 0; i < n; i++)
         node = node->next;
      return node;         /* return pointer to node */
   }
   else
      return NULL;          /* return fail */
}

/*
 * The following function releases the storage occupied by an entire list. Upon
 * success, the program's 'root' pointer is set to NULL and 0 is returned.
 */

int delete_list(ROOT **root)
{
   NODE		*this;          /* pointer to node to be deleted */
   NODE		*next;          /* pointer to next node to be deleted */

   if( !VALID(root))            /* check validity of list */
       return(-1);

   this = (*root)->head;	/* point 'this' to first node in list */

/* Now traverse the list from start to end releasing each node and its data */

   do {
     next = this->next;
     free(this->data);
     free(this);
     this = next;
   } while( this != NULL);

/* Now free the root node and set the program's root pointer to NULL */

   free(*root);
   *root = NULL;
   return 0;
}
