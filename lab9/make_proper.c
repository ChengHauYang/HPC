#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fraction.h"
#include "properfraction.h"

void make_proper(const fraction* in,
                      properfraction* out)
{
  if (in->numer*in->denom>0){
    if (in->numer>0){
      out->whole=(int)in->numer/in->denom;
      out->numer=in->numer%in->denom;
      out->denom=in->denom;
    }else{
      int Pn = - in->numer;
      int Pd = - in->denom;
      out->whole=(int)Pn/Pd;
      out->numer=Pn%Pd;
      out->denom=Pd;
    }
  } else{
    if (-(int)in->numer/in->denom>0){
      if (in->numer>0){
        int Pd = - in->denom;
        out->whole=-(int)in->numer/Pd;
        out->numer=in->numer%Pd;
        out->denom=Pd;
      }else{
        int Pn = - in->numer;
        out->whole = - (int)Pn/in->denom;
        out->numer = Pn%in->denom;
        out->denom = in->denom;
      }
    }else{
      out->whole=0;
      if (in->denom<0){
        out->numer = - in->numer;
        out->denom = - in->denom;
      }else{
        out->numer = in->numer;
        out->denom = in->denom;
      }
    }
  }
}
