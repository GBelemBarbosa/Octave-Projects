graphics_toolkit("qt")
function MakeGif(figHandle, filename)
  persistent persistentFilename = [];
  if isempty(filename)
    error('Can''t have an empty filename!');
  endif
  if ~ishandle(figHandle)
    error('Call MakeGif(figHandle, filename); no valid figHandle was passed!');
  endif
  writeMode = 'Append';
  if isempty(persistentFilename)|(filename!=persistentFilename)
    persistentFilename = filename;
    writeMode = 'Overwrite';
  endif
  imstruct = getframe(figHandle);
  imwrite(imstruct.cdata, filename, 'gif', 'WriteMode',writeMode,'DelayTime',0);
endfunction

A=rand(3, 3);
b=[1;2;3];
for j=1:3
  for p=1:3
    if rand(1)>0.5
      A(j, p)*=-1;
    end
    if rand(1)>0.5
      A(j, p)*=10;
    end
    if rand(1)>0.5
      A_copy=A;
      A_copy(j, p)=0;
      if cond(A_copy)<10^10
        A(j, p)=0;
      end
    end
  endfor
endfor
A
##A =[0.720344124743883,   0.825705713616981,   0.986624961269251;
##  -3.309221225441820,                   0,                  0;
##                   0,                   0,   2.881862123322879]

function weigth=minimize(A, e)
  E=[e(1:3);e(4:6);e(7:9)];
  weigth=100/cond(A+E)+norm(E);
end

if cond(A)<5
  aux=0.1*rand(1,9);

  e=fminunc(@(e) minimize(A, e), aux);
  E=[e(1:3);e(4:6);e(7:9)]
  cond_AE=cond(A+E)
  A_E=A+E
  norm_E=norm(E)

  vects=eye(3);

  f=figure
  quiver3(0, 0, 0, 0, 0, 0)
  hold on;
  for i=1:3
    vect=vects(:,i);
    h=quiver3(0, 0, 0, vect(1), vect(2), vect(3), '--');
    set (h, "color", vects(i,:));
    vect=A*vect;
    h=quiver3(0, 0, 0, vect(1), vect(2), vect(3));
    set (h, "maxheadsize", 0.1/norm(vect));
    set (h, "color", vects(i,:));
  end
  h=patch('Faces',1:3,'Vertices',[A(:,2).';A(:,1).';zeros(1,3)]);
  set(h,'FaceColor','b','FaceAlpha',0.1)
  h=patch('Faces',1:3,'Vertices',[A(:,2).';A(:,3).';zeros(1,3)]);
  set(h,'FaceColor','b','FaceAlpha',0.1)
  h=patch('Faces',1:3,'Vertices',[A(:,1).';A(:,3).';zeros(1,3)]);
  set(h,'FaceColor','b','FaceAlpha',0.1)
  axis("equal")
  hold off;

  inv_norm_inv=1/norm(inv(A))
  inv_norm=1/norm(A)
  cond_A=cond(A)

  c=1/20;
  for j=0:c:1
    f=figure;
    quiver3(0, 0, 0, 0, 0, 0)
    hold on;
    Aux=A+j*E;
    for i=1:3
      vect=vects(:,i);
      h=quiver3(0, 0, 0, vect(1), vect(2), vect(3), '--');
      set (h, "color", vects(i,:));
      vect=Aux*vect;
      h=quiver3(0, 0, 0, vect(1), vect(2), vect(3));
      set (h, "maxheadsize", 0.1/norm(vect));
      set (h, "color", vects(i,:));
    end
    h=patch('Faces',1:3,'Vertices',[Aux(:,2).';Aux(:,1).';zeros(1,3)]);
    set(h,'FaceColor','b','FaceAlpha',0.1)
    h=patch('Faces',1:3,'Vertices',[Aux(:,2).';Aux(:,3).';zeros(1,3)]);
    set(h,'FaceColor','b','FaceAlpha',0.1)
    h=patch('Faces',1:3,'Vertices',[Aux(:,1).';Aux(:,3).';zeros(1,3)]);
    set(h,'FaceColor','b','FaceAlpha',0.1)
    axis("equal")
    MakeGif(f, strcat(num2str(cond_A), '_3d.gif'));
    close(f)
  end
  hold off;

  sig=@(x) 1/(1+exp(-x))

  f=figure
  line([0, 1], [inv_norm_inv, inv_norm_inv], "color", 'r')
  hold on;
  line([0, 0], [0, 0], "color", 'g')
  line([0, 0], [0, 0], "color", 'b')
  line([0, 0], [0, 0], "color", 'y')
  line([0, 1], [0, 0], "color", 'k')
  legend("1/||A^-1||", "||E||", "1/k(A+E)", "||x||/||x-z||", 'AutoUpdate', 'off')
  title(strcat("||A||=", num2str(norm(A)), ", 1/||A^-1||=", num2str(inv_norm_inv), ", k(A)=", num2str(cond_A)))
  x=A\b;
  last=0;
  inv_cond_A=1/cond_A;
  last2=inv_cond_A;
  last3=1;
  for j=c:c:1
    Aux=j*E;
    atual=norm(Aux);
    atual2=1/cond(A+Aux);
    z=(A+Aux)\b;
    atual3=norm(x)/norm(x-z);
    line([j-c, j], [last, atual], "color", 'g')
    line([j-c, j], [last2, atual2], "color", 'b')
    line([j-c, j], [last3, atual3], "color", 'y')
    last=atual;
    last2=atual2;
    last3=atual3;
    ylim([-0.1 max(inv_norm_inv, inv_cond_A)+0.1])
    MakeGif(f, strcat(num2str(cond_A), '_cu.gif'));
  end
  hold off;
  close(f)
end