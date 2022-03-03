#include <stdio.h>
#include <stdlib.h>
#include "node.h"
#include <time.h>
#include <math.h>

void GenerateList(node** head, const int num_nodes)
{
   node* temp; srand( time(NULL) );

   for (int i=0; i<num_nodes; i++)
   {
      temp = (node*)malloc(sizeof(node));
      temp->value = rand()%num_nodes; temp->position = 0;
      printf("%4i",temp->value);

      if (*head == NULL)
      {
	 *head = temp;
	 (*head)->next = NULL;
      }
      else
      {
	 temp->next = *head;
	 *head = temp;
      }
   }
   printf("\n");

   node* ptr = *head; int pos = 1;
   while(ptr!=NULL)
   {
      ptr->position = pos;
      pos = pos+1;
      ptr = ptr->next;
   }
}


void GenerateList2(node** head)
{
   node* temp; srand( time(NULL) );

   for (int i=0; i<7; i++)
   {
      temp = (node*)malloc(sizeof(node));
      temp->value = (7-i)*pow(-1,i); temp->position = 0;
      printf("%4i",temp->value);

      if (*head == NULL)
      {
	 *head = temp;
	 (*head)->next = NULL;
      }
      else
      {
	 temp->next = *head;
	 *head = temp;
      }
   }
   printf("\n");

   node* ptr = *head; int pos = 1;
   while(ptr!=NULL)
   {
      ptr->position = pos;
      pos = pos+1;
      ptr = ptr->next;
   }
}
