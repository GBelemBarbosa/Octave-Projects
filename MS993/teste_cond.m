num_iter=1000000;
pos_pure=0;
pos_neg=0;
pos_zer=0;
neg_zer=0;
for i=1:num_iter
  A=rand(3, 3);
  A_zer=A;
  A_neg_zer=A;
  pos_pure+=cond(A);
  if cond(A)>100000
    A
    cond(A)
  endif
  for j=1:3
    for p=1:3
      if rand(1)>0.5
        A(j, p)*=-1;
        A_neg_zer(j, p)*=-1;
      end
      if rand(1)>0.5
        A_zer_copy=A_zer;
        A_zer_copy(j, p)=0;
        if cond(A_zer_copy)<10^10
          A_zer(j, p)=0;
          A_neg_zer(j, p)=0;
        elseif cond(A_zer_copy)<10^14
          A_zer_copy
          cond(A_zer_copy)
        end
      end
    endfor
  endfor
  pos_neg+=cond(A);
  pos_zer+=cond(A_zer);
  neg_zer+=cond(A_neg_zer);
endfor
pos_pure/=num_iter
pos_neg/=num_iter
pos_zer/=num_iter
neg_zer/=num_iter