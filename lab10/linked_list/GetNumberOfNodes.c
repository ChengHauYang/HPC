#include <stdio.h>
#include <stdlib.h>

int GetNumberOfNodes()
{
   int num_nodes;
   printf("\n Enter the number of nodes: ");
   scanf("%d", &num_nodes);

   while (num_nodes<0)
   {
      printf("\n Invalid input: input must a non-negative integer. ");
      printf("\n Enter the number of nodes: ");
      scanf("%d", &num_nodes);
   }
   
   return num_nodes;
}
