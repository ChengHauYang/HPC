#include <stdio.h>
#include <stdlib.h>
#include "node.h"

void Print(const int forward, const node* head)
{
   if (head==NULL)
   { printf(" List is emtpy.\n\n"); return; }
   printf("\n");
   printf(" ---------------------------------------------\n");
   printf(" |Pos:|Val:|    Address:    |      Next:     |\n");
   printf(" ---------------------------------------------\n");   
   switch(forward)
   {
      case 0:
	 ReversePrintList(head);
	 break;
      case 1:
	 PrintList(head);
	 break;
      default:
	 printf("\n Error: forward must be 0 or 1.\n");
	 printf("  forward = %i\n",forward); exit(1);
   }
   printf(" ---------------------------------------------\n");
}

void PrintLine(const int pos,const int val,
	       const node* head,const node* next)
{ printf(" |%3i |%3i |%15p |%15p |\n",pos,val,head,next); }

void PrintList(const node* head)
{
   PrintLine(head->position,head->value,head,head->next);
   if (head->next == NULL)      
   { return; }
   PrintList(head->next);
}

void ReversePrintList(const node* head)
{
   if (head->next == NULL)      
   {
      PrintLine(head->position,head->value,head,head->next);
      return;
   }   
   ReversePrintList(head->next);
   PrintLine(head->position,head->value,head,head->next);
}
