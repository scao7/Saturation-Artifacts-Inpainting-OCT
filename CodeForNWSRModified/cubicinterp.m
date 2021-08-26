function missingcol=cubicinterp(colidx,value,ii,R)
    missingcol=zeros(R,1);    
    for jj=1:R
        Y =value(jj,:);
        X = colidx;
        cs = csapi(X,Y);   %三次样条函数
        index=ii;
        ivalue=ppval(cs,index);
        missingcol(jj,1)=ivalue;
    end
end
