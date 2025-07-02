rates<-function(Beta_Matrix,q){
    tmp=Beta_Matrix[,1,]
	  p=dim(tmp)[1]
	  replicate=dim(tmp)[2]
    part1=tmp[1:q,]
    part2=tmp[(q+1):p,]
    r1=1-sum(sum(part1==0))/length(part1);
    r2=sum(sum(part2==0))/length(part2);
    right=0;
    for(i in 1:replicate){
        a=part1[,i];
        b=part2[,i];
        if(sum(a==0)==0 && sum(b==0)==p-q){
		   right=right+1
		   }
		   }    
    r3=right/replicate
	C=matrix(c(r1,r2,r3),1,3)
	C
}
