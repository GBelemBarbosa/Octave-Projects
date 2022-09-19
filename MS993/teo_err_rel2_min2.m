a=2
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
for j=1:3
  for p=1:3
    if rand(1)>0.5
      A(j, p)*=-1;
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
inv_norm_inv=1/norm(inv(A))
cond_A=cond(A)

function weigth=minimize(A, e)
  E=[e(1:3);e(4:6);e(7:9)];
  weigth=10/cond(A+E)+norm(E);
end

aux=0.1*rand(1,9);

e=fminunc(@(e) minimize(A, e), aux);
E=[e(1:3);e(4:6);e(7:9)]
cond_AE=cond(A+E)
A_E=A+E
norm_E=norm(E)

vects=eye(3);

figure
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

for j=0:0.1:1
  f=figure
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
  MakeGif(f, 'test.gif');
  close(f)
end
hold off;

f=figure
line([0, 1], [inv_norm_inv, inv_norm_inv], "color", 'r')
line([0, 1], [0, 0], "color", 'k')
hold on;
last=0;
last2=det(A);
c=0.1;
for j=c:c:1
  Aux=j*E;
  atual=norm(Aux);
  atual2=det(A+Aux);
  line([j-c, j], [last, atual])
  line([j-c, j], [last2, atual2])
  last=atual;
  last2=atual2;
  MakeGif(f, 'norm.gif');
end
hold off;
