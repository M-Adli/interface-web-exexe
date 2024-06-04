eq(X,X).


p(a).
p(b).
q(b).


app([],Ys,Ys).
app([X|Xs],Ys,[X|Zs]) :- app(Xs,Ys,Zs).

    type_rec(app,1,lst,true).
    type_rec(app,2,term,false).
    type_rec(app,3,term,true).

leng([],0).
leng([_X|Xs],s(N)) :- leng(Xs,N).

    type_rec(leng,1,lst,true).
    type_rec(leng,2,nat,true).

lst([]).
lst([_|Xs]) :- lst(Xs).
    
    type_rec(lst,1,lst,true).

nat(0).
nat(s(X)) :- nat(X).

    type_rec(nat,1,nat,true).
    
rev([],[]).
rev([X|Xs],Ys) :- rev(Xs,Zs), app(Zs,[X],Ys).

    type_rec(rev,1,lst,true).
    type_rec(rev,2,lst,true).


add(0,Y,Y).
add(s(X),Y,s(Z)) :- add(X,Y,Z).

    type_rec(add,1,nat,true).
    type_rec(add,2,term,false).
    type_rec(add,3,term,false).

infe(0,_X).
infe(s(X),s(Y)):- infe(X,Y).

     type_rec(infe,1,nat,true).
     type_rec(infe,2,nat,true).

inf(0,s(_X)).
inf(s(X),s(Y)):- inf(X,Y).

     type_rec(inf,1,nat,true).
     type_rec(inf,2,nat,true).

double(0,0).
double(s(X),s(s(Y))) :- double(X,Y).

     type_rec(double,1,nat,true).
     type_rec(double,2,nat,true).

half(0,0).
half(s(0),0).
half(s(s(X)),s(Y)) :- half(X,Y).

    type_rec(half,1,nat,true).
    type_rec(half,2,nat,true).

mul(0,_Y,0).
mul(s(X),Y,Z) :- mul(X,Y,T), add(Y,T,Z).

    type_rec(mul,1,nat,true).
    type_rec(mul,3,nat,false).


fib(0,s(0)).
fib(s(0),s(0)).
fib(s(s(X)),Y) :- fib(s(X),Z),fib(X,T),add(Z,T,Y).

    type_rec(fib,1,nat,true).
    type_rec(fib,2,nat,false).
    
    
perm([],[]).
perm([A|Xs],Ys):- perm(Xs,Zs), insert(A,Zs,Ys).

    type_rec(perm,1,lst,true).
    type_rec(perm,2,lst,false).

ord([]).
ord([_]).
ord([A,B|Xs]):- 
    A =< B,
    ord([B|Xs]).

    type_rec(ord,1,lst,true).

insert(A,As,[A|As]).
insert(A,[B|Xs],[B|Ys]):- insert(A,Xs,Ys).

    type_rec(insert,1,term,false).
    type_rec(insert,2,lst,true).
    type_rec(insert,3,lst,true).

leq_or_gt(X,Y) :- X =< Y.
leq_or_gt(X,Y) :- X > Y.

memb(X,[X|_]).
memb(X,[_|L]) :- memb(X,L).

    type_rec(memb,2,lst,true).