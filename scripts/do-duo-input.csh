exec awk '
BEGIN{inp=0;out=0;n=0;i=0;imax=0;itermax=500}
{
 if ($1=="(Transcript"){inp=1; next}
 if ($1=="(<---"){inp=0; next}

 if ($1=="Parameters:"){out=1; next}
 if ($1=="Fitted"){out=0;imax=i;i=0; next}
 if ($1=="Iteration"){iter=$3; next}
 if (iter>itermax){next}
 
 if ( inp == 0 && out == 1){
   if ( $3 == "fit" ){
     i++;
     par[i] = $2
     #print i,par[i];
   };
 };

 if ( inp == 1 ){
   n++;
   if ( $3 == "fit" ){
     char1[n]=$1
     char2[n]=$2
     char3[n]=$3
     F[n] = $1
   }
   else {
     char1[n]="none"
     F[n] = $0
   };
   #F[n] = $0
   #print F[n]
 };
};
END{
  k0=0;
  for ( i0=1;i0<n;i0++ ){
    if ( char1[i0]=="none" ){
      print F[i0];
    }
    else{
      k0++;
      printf("%-12s%21.14e   %3s    (%22.14e)  \n",F[i0],par[k0],char3[i0],char2[i0])
    };
  }
  #print F[i0];

  for ( i0=1;i0<k;i0++ ){
    print i0,par[i0];
  }

};
' "$@"

#  printf("%12.5f %19.8f\n",$2,$4)
