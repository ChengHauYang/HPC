#include <stdio.h>
#include <stdlib.h>
#include "node.h"

/*
use the quicksort algorithm (see Lab #5) to sort the list from smallest value to largest value.
*/

int GetValue(node* head,int pos)
{
    /*
This function return A[pos-1] : start from "0" to "num_nodes-1"
    */

    int pos_temp=pos+1;
    node *prev = NULL, *curr = head;
    while (curr->position != pos_temp) {
        prev = curr;
        curr = curr->next;
    }

    int value =curr->value;

    return value;
}

void quicksort(node** head,int lo, int hi,int num_nodes)
{
    /*
This function uses the quicksort algorithm
to sort.
    */

    lo=lo-1;
    hi=hi-1;
    //printf("hi:%d\n", hi);

    if (hi==-100){
      hi = num_nodes-1;
    }


    int partition(node** head, int lo, int hi);

    if (lo < hi){
        //printf("lo < hi \n");
        int p = partition(head, lo, hi);
        //printf("p:%d \n",p);
        quicksort(head, lo + 1, p,num_nodes);
        quicksort(head, p + 2, hi + 1,num_nodes);
    } else{
      return;
    }

}

int partition(node** head, int lo, int hi)
{
    //printf("new partition===========\n");
    void Swap2Nodes(node** head, int i, int j);

    int pivot = GetValue(*head,hi);
/*
    printf("pivot:%d\n", pivot);
    printf("lo:%d\n", lo);
*/
    int i = lo; // place for swapping
    for(int j=lo;j<hi;j++){ // lo to hi - 1
      //printf("Value of pos %d: %d\n",j, GetValue(*head,j));
      if(GetValue(*head,j) <=pivot){ // swap A[i] with A[j]
        Swap2Nodes(head,i+1,j+1);
        /* swap writing from 1~numnodes
        python writing from 0~numnodes-1 */
        i++;
      }
    }

/*
    printf("i:%d\n", i);
    printf("hi:%d\n", hi);
    Print(1,*head); // foward print
*/
    // swap A[i] with A[hi]
    Swap2Nodes(head,i+1,hi+1);
    //printf("after swap");
    /* swap writing from 1~numnodes
    python writing from 0~numnodes-1 */

    return i;
}

void SortList(node** head)
{

  /// ---Using next to continue the LOOP---
  int num_nodes = 0;
  node *prev = NULL;
  node *curr = *head;
  //printf("dereference head: %15p\n",*head);
  //printf("current: %15p\n",curr);
  curr=curr->next;
  //printf("current: %15p\n",curr);

  while (curr != NULL) {
      prev = curr;
      curr = curr->next;
      //printf("%15p\n",curr);
      num_nodes++;
  }
  num_nodes++;
  //printf("num_nodes: %d\n", num_nodes);

  void quicksort(node** head,int lo, int hi,int num_nodes);
  quicksort(head, 1, num_nodes,num_nodes);


}
