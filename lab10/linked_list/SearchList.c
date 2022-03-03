#include <stdio.h>
#include <stdlib.h>
#include "node.h"

void SearchList(const node* head,
		const int key)

{
   if (head->value==key)      
   {
      printf(" Key found at Position: %i\n",
	     head->position);
   }

   if (head->next==NULL)      
   {
      printf("\n"); return;
   }

   SearchList(head->next,key);
}
