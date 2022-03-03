#include <stdlib.h>
#include <stdio.h>
#include "node.h"
int main()
{
   // Declare the head node in my list
   node* head = NULL;

   // Set number of nodes and generate a new list
   const int num_nodes = GetNumberOfNodes();
   GenerateList(&head,num_nodes);

   // Print list to screen
   Print(1,head); // foward print
   Print(0,head); // reverse print

// testing
/*
   int GetValue(node* head,int pos);
   printf("Value of pos3: %d\n", GetValue(head,3-1));
   printf("Value of pos5: %d\n", GetValue(head,5-1));

   void Swap2Nodes(node** head, int i, int j);
   Swap2Nodes(&head,5,3);
   printf("after swapping 3 and 5:");
   Print(1,head); // foward print
   Print(0,head); // reverse print
*/

   void SortList(node** head);
   SortList(&head);
   printf("after sorting:");
   Print(1,head); // foward print
   Print(0,head); // reverse print

   // Ask for a key, then search list
   if(num_nodes>0)
   {
      const int key = GetKey();
      SearchList(head,key);
   }

// testing !
/*
   void GenerateList2(node** head);
   const int num_nodes = 7;
   GenerateList2(&head);
*/

   // Delete list (free up memory)
   DeleteList(&head);
}
