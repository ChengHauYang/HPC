#include <stdio.h>
#include <stdlib.h>
#include "node.h"

/*
switches the node in position i and the node in position j.
*/

void Swap2Nodes(node** head, int i, int j)
{
    if (i == j){
      return;
    }

    /// ---Using next to continue the LOOP---
    node *prevI = NULL, *currI = *head;
    //int k=0;
    while (currI && currI->position != i) { // LOOP to "position i-1"
        //k++;
        //printf("%d\n", k);
        prevI = currI;  // for position i it is previous
        currI = currI->next;  // for position i it is current
    }

    /// ---Using next to continue the LOOP---
    node *prevJ = NULL, *currJ = *head;
    while (currJ && currJ->position != j) { // LOOP to "position j-1"
        prevJ = currJ; // for position j it is previous
        currJ = currJ->next; // for position j it is current
    }

    // NullPtr
    if (currI == NULL || currJ == NULL){
        return;
    }

    //KEY: swap "previous+next=current" = swap current
    //prevI->next = currJ;
    //prevJ->next = currI;
    //currJ = prevI->next;
    //currI = prevJ->next;

    //node* tempcur = currJ;  // CAN NOT DO LIKE THIS because curr=head
    //currJ = currI;
    //currI = tempcur;

    // If i is not head of linked list
       if (prevI){
         //printf("prevI != null");
         prevI->next = currJ;
       } else { // Else make j as new head
         //printf("prevI == null");
         *head = currJ;
       }

       // If j is not head of linked list
       if (prevJ){
         //printf("prevJ != null");
         prevJ->next = currI;
       } else { // Else make i as new head
         //printf("prevJ == null");
         *head = currI;
       }


    // Swap "next" pointers "back" => make the sequence normal
    // If NOT doing this, some lines disappear
    node* temp = currJ->next;
    currJ->next = currI->next;
    currI->next = temp;

    int pos = currJ->position;
    currJ->position = currI->position;
    currI->position = pos;
}
