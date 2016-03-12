clc 
clear all 
nc=input('Enter the value of nc:');
% k=2;
l=input('Enter the length of the input code:');
% l=4;
kc=input('Enter the constrain length:');
% kc=3;
ne=input('Enter the number of errors you need to add to encoded codeword:');
in=double(rand([1,l])>0.5); % Generates a random Input Sequence
% in=[1 0 1 1];
display(in);
% g=[ 1 1 1; 1 0 1];
for j=1:nc
    g(j,:)=double(rand([1,kc])>0.5); % Generates a random generator vectors
    cwa(j,:)=rem(conv(g(j,:),in),2);
end 
display(cwa);
cw=cwa;
cw(1,:)=rem(cwa(1,:)+de2bi(2^(ne)-1,l+kc-1),2);
display(cw);
sd=cell(2^(kc-1),2);
for ns=0:(2^(kc-1)-1);
    a=de2bi(ns,kc-1);
    sd{ns+1,1}=[0,a(1,kc-1:-1:2)];
    sd{ns+1,2}=[1,a(1,kc-1:-1:2)];
end 
so=cell(2^(kc-1),2);
for ns=0:(2^(kc-1)-1);
    a=de2bi(ns,kc-1);
    for j=1:nc
    b=rem(sum([0,a(1,kc-1:-1:1)].*g(j,:)),2);
    so{ns+1,1}=[so{ns+1,1},b];
    b=rem(sum([1,a(1,kc-1:-1:1)].*g(j,:)),2);
    so{ns+1,2}=[so{ns+1,2},b];
    end 
end 
bm=cell(2^(kc-1),l);
for n=1:l  
    bf=cw(1:nc,n)';
    for p=1:2^(kc-1)
      bm{p,n}(1)=sum(rem(bf+so{p,1} ,2));
      bm{p,n}(2)=sum(rem(bf+so{p,2} ,2));
    end 
end
pm=cell(2^(kc),l);
td=cell(2^(kc),l);
pst=[1 zeros(1,2^(kc-1)-1)];
nst=[(bi2de(flip(sd{pst(1),1}))+1),(bi2de(flip(sd{pst(1),2}))+1),zeros(1,2^(kc-1)-2)];
nns=2;
pm{(bi2de(flip(sd{pst(1),1}))+1),1}=bm{1,1}(1);
td{(bi2de(flip(sd{pst(1),1}))+1),1}=[1,0];
td{(bi2de(flip(sd{pst(1),2}))+1),1}=[1,1];
pm{(bi2de(flip(sd{pst(1),2}))+1),1}=bm{1,1}(2);
nstb=zeros(1,2^(kc-1));
flags=zeros(2^(kc-1),l);
for n=2:l
    g=1;
    for h=1:nns
        if flags((bi2de(flip(sd{nst(h),1}))+1),n)==0
        td{(bi2de(flip(sd{nst(h),1}))+1),n}...
            =[nst(h),0];
        pm{(bi2de(flip(sd{nst(h),1}))+1),n}...
            =[pm{(bi2de(flip(sd{nst(h),1}))+1),n},bm{nst(h),n}(1)+pm{nst(h),n-1}];
        flags((bi2de(flip(sd{nst(h),1}))+1),n)=1;
        else 
            td{2^(kc-1)+(bi2de(flip(sd{nst(h),1}))+1),n}...
            =[nst(h),0];
            pm{2^(kc-1)+(bi2de(flip(sd{nst(h),1}))+1),n}...
            =[pm{2^(kc-1)+(bi2de(flip(sd{nst(h),1}))+1),n},bm{nst(h),n}(1)+pm{nst(h),n-1}];
        end 
        nstb(g)=(bi2de(flip(sd{nst(h),1}))+1);g=g+1;
        if flags((bi2de(flip(sd{nst(h),2}))+1),n)==0
        td{(bi2de(flip(sd{nst(h),2}))+1),n}...
            =[nst(h),1];
        pm{(bi2de(flip(sd{nst(h),2}))+1),n}=...
            [pm{(bi2de(flip(sd{nst(h),2}))+1),n},bm{nst(h),n}(2)+pm{nst(h),n-1}];
         flags((bi2de(flip(sd{nst(h),2}))+1),n)=1;
        else 
            td{2^(kc-1)+(bi2de(flip(sd{nst(h),2}))+1),n}...
            =[nst(h),1];
        pm{2^(kc-1)+(bi2de(flip(sd{nst(h),2}))+1),n}=...
            [pm{2^(kc-1)+(bi2de(flip(sd{nst(h),2}))+1),n},bm{nst(h),n}(2)+pm{nst(h),n-1}];
        end
        nstb(g)=(bi2de(flip(sd{nst(h),2}))+1);g=g+1;
    end
    for h=1:2^(kc-1)
        if ~isempty(pm{2^(kc-1)+h,n})
       [val,ordr]=min([pm{h,n}(1) pm{2^(kc-1)+h,n}(1)]);
       if ordr==1
        pm{2^(kc-1)+h,n}=[];
        td{2^(kc-1)+h,n}=[];
       else 
         pm{h,n}(1)=pm{2^(kc-1)+h,n}(1);
         td{h,n}=td{2^(kc-1)+h,n};
         pm{2^(kc-1)+h,n}=[];
         td{2^(kc-1)+h,n}=[];
       end 
        end
    end 
    for t=n:-1:2
        rb1=zeros(2^(kc-1),1);
        rb2=zeros(2^(kc-1),1);
           for rr=1:2^(kc-1)
               if ~isempty(td{rr,t})
                rb1(rr)=td{rr,t}(1);
               end 
           end 
           for rr=1:2^(kc-1)
               if ~isempty(td{rr,t-1})
                rb2(rr)=rr;
               end
           end
           for pr=1:2^(kc-1)
               mtch=0;
               for qr=1:2^(kc-1)
                  if rb2(pr)==rb1(qr)
                    mtch=1;
                  end 
               end 
             if mtch==0
             td{pr,t-1}=[];
             end 
           end 
        end      
    pst=nst;
    nst=nstb;
    nns=min(2*nns,2^(kc-1));
end 
gt=zeros(1,2^(kc-1));
for t=1:2^(kc-1)
    gt(t)=pm{t,l};
end 
[mhd,pn]=min(gt);
for  t=1:2^(kc-1)
    if t~=pn
        td{t,l}=[];
    end 
end 
 for t=l:-1:2
        rb1=zeros(2^(kc-1),1);
        rb2=zeros(2^(kc-1),1);
           for rr=1:2^(kc-1)
               if ~isempty(td{rr,t})
                rb1(rr)=td{rr,t}(1);
               end 
           end 
           for rr=1:2^(kc-1)
               if ~isempty(td{rr,t-1})
                rb2(rr)=rr;
               end
           end
           for pr=1:2^(kc-1)
               mtch=0;
               for qr=1:2^(kc-1)
                  if rb2(pr)==rb1(qr)
                    mtch=1;
                  end 
               end 
             if mtch==0
             td{pr,t-1}=[];
             end 
           end 
        end      
dcdout=[];
for n=1:l
    for h=1:2^(kc-1)
        if ~isempty(td{h,n})
            dcdout=[dcdout td{h,n}(2)];
        end 
    end 
end 
display(dcdout);