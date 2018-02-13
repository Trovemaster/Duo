exec awk '
BEGIN{thresh=8;n=1;rms=0}
{
 if ($1=="Iteration"){iter=$3; print $0;
   print "rms = " sqrt(rms/n) " ;n=",n
   rms=0;n=0;
 }
 if ( $9=="(" && $8>0){rms=rms+$7*$7;n=n+1}
 if ( sqrt($7*$7)>thresh && $8>0  && $9=="("){print $0; next}
};
END{
};
' "$@"

