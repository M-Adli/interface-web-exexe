
:- set_prolog_flag(occurs_check, error).
:- use_module(library(lists)).
:- use_module(library(os)).
%:- discontiguous(type_rec/4).

:- dynamic(nat/1).
:- dynamic(lst/1).
:- dynamic(app/3).
:- dynamic(add/3).
:- dynamic(leng/2).
:- dynamic(infe/2).
:- dynamic(double/2).
:- dynamic(half/2).

write_exexe(X) :- % modification
   open('output.txt', append, Stream), % Ouvre le fichier en mode append
    write(Stream, X),nl,
    nl(Stream), 
    close(Stream).
write_exexe1(X) :- % modification
   open('output.txt', append, Stream), % Ouvre le fichier en mode append
    write(Stream, X),
    close(Stream).

% Déclaration d'une variable globale pour stocker le flux de fichier
:- dynamic(file_stream/1). % modification

% Prédicat pour lire le contenu du fichier ligne par ligne et retourner le résultat dans X
read_exexe(X) :-  % modification
    % Si le flux de fichier n'est pas déjà ouvert, ouvrez-le
    (   \+ file_stream(_),
        open('input.txt', read, Stream),
        asserta(file_stream(Stream))
    ;   true
    ),
    % Lire une ligne à partir du flux de fichier
    read_line(Stream, Line),nl,
    % Fermer le flux de fichier s'il est à la fin
    (   Line == end_of_file ->
        close_file_stream,
        delete_file('input.txt') % Supprimer le fichier
    ;   X = Line , \+ \+ (numbervars(X,23,N), % Attribuer la ligne lue à X
       
       write_exexe(X))
    ).

% Prédicat pour lire une ligne à partir du flux de fichier
read_line(Stream, Line) :-
    file_stream(Stream), % Vérifier si le flux de fichier est ouvert
    read(Stream, Line).

% Prédicat pour fermer le flux de fichier
close_file_stream :-
    retract(file_stream(Stream)),
    close(Stream).








app([],Ys,Ys).
app([X|Xs],Ys,[X|Zs]) :- app(Xs,Ys,Zs).

    type_rec(app,1,lst,true).
    type_rec(app,2,term,false).
    type_rec(app,3,term,false).

leng([],0).
leng([_X|Xs],s(N)) :- leng(Xs,N).

    type_rec(leng,1,lst,true).
    type_rec(leng,2,nat,true).

length([], 0).
length([_|T], N) :- length(T, N1), N is N1 + 1.    

lst([]).
lst([_|Xs]) :- lst(Xs).

    type_rec(lst,1,lst,true).

    nat(0).
    nat(s(X)) :- nat(X).

% ok : leng(X,?(L)).
% ok : leng(?(X),L).
% ok : app(X,Y,?(Z))&leng(?(Z),?(T)).
% ok : leng(X,?(L))&leng(?(Y),?(L)).   % preuve puis tronquer Y -> genere same_length
% ok : app(X,Y,Z) => leng(Z,?(T)).
% ok : app(X,Y,Z) => lst(X).
% app(X,Y,Z) => (lst(Y)=> leng(X,?(Lx)) & leng(Y,?(Ly)) & leng(Z,?(Lz)) & add(?(Lx),?(Ly),?(Lz))).
% ok : app(X,Y,Z) =>(lst(Y) => lst(Z)).
% ok : app(X,Y,Z) => lst(X).
% app(X,Y,Z) =>(lst(Y) => lst(Z)).
% app(X,Y,Z) =>(lst(Z) => lst(Y)).

add(0,Y,Y).
add(s(X),Y,s(Z)) :- add(X,Y,Z).

    type_rec(add,1,nat,true).
    type_rec(add,2,term,false).
    type_rec(add,3,term,false).

infe(0,_X).
infe(s(X),s(Y)):- infe(X,Y).

     type_rec(infe,1,nat,true).
     type_rec(infe,2,nat,true).


double(0,0).
double(s(X),s(s(Y))) :- double(X,Y).

     type_rec(double,1,nat,true).
     type_rec(double,2,nat,true).

half(0,0).
half(s(0),0).
half(s(s(X)),s(Y)) :- half(X,Y).

    type_rec(half,1,nat,true).
    type_rec(half,2,nat,true).

fib(0,s(0)).
fib(s(0),s(0)).
fib(s(s(X)),Y) :- fib(s(X),Z),fib(X,T),add(Z,T,Y).

    type_rec(fib,1,nat,true).
    type_rec(fib,2,nat,false).
    



% ok : add(X,Y,?(Z)).
% ok : double(X,?(Y)).
% ok : add(X,0,X).
% double(X,Y) => add(X,X,Y).    % comput, puis necessite commutativite de add
% double(X,Y) => add(X,?(Z),Y). % comput ...
% add(X,X,Y) => double(X,Y).    % echec, il faut var toutes diff sur H
% ok : add(X,0,X).
% add(X,Y,Z) => add(X,s(Y),s(Z)).
% add(s(X),Y,Z)=>add(X,s(Y),Z).
% add(X,s(Y),Z)=>add(s(X),Y,Z).
% add(s(X),Y,Z)=> add(X,s(Y),?(T)).  % ?
% add(X,Y,Z) => add(Y,X,Z).
% ok : infe(X,Y)=>add(X,?(Z),Y).  synthetise la soustraction : Z = Y-X
% half(X,Y) => infe(Y,X). % preuve ok mais bug depl
% ok : double(X,Y) => half(Y,X).

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



% perm(X,?(Y)) & ord(?(Y)).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

:- discontiguous(extex/2).
:- discontiguous(handle_ctnue/9).
:- discontiguous(menu_pair_occ/7).
:- discontiguous(extex/8).
:- discontiguous(handle_end/5).
%:- discontiguous(extex_fail/8).
%:- discontiguous(next/6).


%%%%

:- op(100,xfy,&).
:- op(100,xfy,=>).
:- op(50,xfx,:).

my_name(Atom,Codes) :- %name(Atom,Codes).
    number(Atom), !, number_codes(Atom,Codes).
my_name(Atom,Codes) :- %name(Atom,Codes).
    atom(Atom), !, atom_codes(Atom,Codes).
my_name(Atom,Codes) :- %name(Atom,Codes).
    atom_codes(Atom,Codes).

exexe_system(Cmd) :- write(Cmd),nl. %,shell(Cmd).


%%%%%


mb(X, [X|_]).
mb(X, [_|Z]) :-
    mb(X,Z).

app_pgm([],Y,Y).
app_pgm([A|X],Y,[A|Z]) :-
    app_pgm(X, Y, Z).

rev_pgm(X,Y) :- rev1(X,[],Y).

    rev1([],Z,Z).
    rev1([A|X],Z,Y) :- rev1(X, [A|Z],Y).

% assoc est utilise pour associer une cle avec sa valeur dans une liste
% d'associations representee sous la forme d'une liste de paires (cle,valeur).

assoc(X,[(Y,Z)|L],Z,L) :- X == Y, !.
assoc(X,[(Y,Z)|L],T, [(Y, Z)|Rmd]) :-
    assoc(X,L,T,Rmd).

/*?- assoc(color, [(color, red), (shape, circle), (size, large)], Value,Rmd).Value = red,
Rmd = [(shape, circle), (size, large)].

?- assoc(size, [(color, red), (shape, circle), (size, large)], Value,Rmd).
Value = large,
Rmd = [(color, red), (shape, circle)]
*/

merge([],Y,Y).
merge([A|X],Y,Z) :-
    ins(A,Y,Y1),
    merge(X, Y1 ,Z).

ins(A,[B|X],[B|X]) :-
    A == B , !.
ins(A,[B|X], [B|Y]) :-
    ins(A, X,Y) .
ins(A,[],[A]).

/*
merge([1,6],[1,2,3,4],Result).
Result = [1, 2, 3, 4, 6]

ins(4, [1, 2, 4 , 5, 6], Result).
Result = [1, 2, 4, 5, 6] ;
false.
*/


copy(A,B) :- copy_term(A,B).

/* copy(2,X).
X = 2 */

same_functor(X,Y,F,N) :-
    functor(X,F,N),
    functor(Y,F,N).

/*exemple same_functor*/
/*
?- same_functor(f(a, b), f(c, d), F, N).
F = f,
N = 2.

?-  same_functor(f(a, b), g(c, d), F, N).
false.
*/






% variables

% varexist(X, Y): X est la variable existentielle ?(Y)

varexist(X, ArgX) :-
    nonvar(X),
    functor(X, ?, 1),
    arg(1, X, ArgX),
    var(ArgX).

/*?- varexist(?(y),ArgX).
false.

?- varexist(?(Y),ArgX).
Y = ArgX.
%complexe et existentielle de premier terme Y a cause de ?
?- varexist((Y),ArgX).
false.

?- varexist(f(a,b),ArgX).
false.
compplexe mais pa existentielle
*/





%toutes les variables d'une formule
variables(X,[X]) :-
    var(X), !.
variables(B, Vars) :-
    functor(B, _, N),
    variables_arg(N, B, Vars).

variables_arg(0, _, []).
variables_arg(N, B, Vars) :-
    N > 0,
    arg(N, B, ArgB),
    variables(ArgB, Vars1),
    N1 is N-1,
    variables_arg(N1, B, Vars2),
    merge(Vars1, Vars2, Vars).


/*?- variables(X,Res).
Res = [X].

?- variables(f(a,b),Res).
Res = [].

?- variables(struct(f(a,X,Y,b),g(Z,4,T,9)),Res).

Res = [X, Y, Z, T].
*/



member_var(X, [Y|L], L) :-
    X==Y,!.
member_var(X, [Y|L], [Y|R]) :-
    member_var(X,L,R).

/*?- member_var(X, [a, X, b, Y, c], R).
R = [a, b, Y, c].

?- member_var(X,[X,5,Y,Z],Resultat).
Resultat = [5, Y, Z].

?- member_var(X,[A,5,Y,Z],Resultat).
false.
*/

% variables existentielles et universelles d'une formule

vars(X, [X],[]) :-
    varexist(X,_), !.
vars(X, [], [X]) :-
    var(X), !.
vars(A,Exist, Univ) :-
    functor(A,_,N),
    vars_arg(N, A, Exist, Univ).

vars_arg(0,_, [],[]).
vars_arg(N,B,Exist, Univ) :-
    N > 0,
    arg(N,B, ArgB),
    vars(ArgB, Exist1, Univ1),
    N1 is N-1,
    vars_arg(N1, B, Exist2, Univ2),
    merge(Exist1, Exist2, Exist),
    merge(Univ1, Univ2, Univ).




% variables internes d'une clause

internal_vars(H, B, Vars) :-
    variables(H, VarsH),
    variables(B, VarsB),
    minus1(VarsB, VarsH, Vars).

% minus: L = L1\L2
minus(L1, [X|L2],L) :-
    member_var(X, L1,L3),
    !,
    minus(L3, L2,L).
minus(L1, [X|L2], [X|L]) :-
    minus(L1,L2,L).
minus(L,[],L).

minus1([X|L],M,N) :-
    member_var(X,M,Rmd),
    !,
    minus1(L, Rmd, N).
minus1([X|L],M, [X|N]) :-
    minus1(L,M,N).
minus1([],_,[]).

/* ?- internal_vars([1,2,3,4,5,X,Y],[1,2,3,X],Vars).
Vars = [] ;
false.
?- internal_vars([Y,2,3,4,5,X],[Z],Vars).
Vars = [Z] ;
false.
?- internal_vars([Y,2,3,4,5,X],[Z,Y,T],Vars).
Vars = [Z, T]*/



% liste->clause

list_to_clause([H], H) :- !.
list_to_clause([H|T], (H:-T1)) :-
    %list_to_set(T,T0),
    T0 = T,
    list_to_body(T0,T1).

    list_to_body([A],A) :- !.
    list_to_body([A,B|Cs],(A,Ds)) :-
        list_to_body([B|Cs],Ds).

/*?- list_to_clause([b],H).
H = b.

?- list_to_clause([a,b],H).
H = (a:-b).

?- list_to_clause([a,b,c],H).
H = (a:-b, c).*/

% clause->liste

/*
clause_to_list((H :- Bs), [H | Bl]) :- !,  body_to_list(Bs, Bl).
clause_to_list(H, [H]).

    body_to_list((A,As),[A|L]) :- !, body_to_list(As,L).
    body_to_list(A,[A]).
*/

/*?- ?- clause_to_list(h,X).
X = [h].

?- clause_to_list((h:-s),X).
X = [h, s].

?- clause_to_list((a:-b, c),X).
X = [a, b, c].
*/



%ensemble ->list

set_to_list(true,[]) :- !.
set_to_list((B,Bs),[B|B1s]) :- !, set_to_list(Bs, B1s).
set_to_list(B, [B]) :-  !.

/*?- set_to_list(true,X).
X = [].

?- set_to_list((a,b,c,d),X).
X = [a, b, c, d].

?- set_to_list((A),X).
A = true,
X = [].

?- set_to_list((a),X).
X = [a].
*/


% liste -> ensemble
exexe_list_to_set([], []) :- !.
%exexe_list_to_set([T], T) :- !.
exexe_list_to_set([T|T1], Ts) :- exexe_member_eq(T, T1), !, exexe_list_to_set(T1, Ts).
exexe_list_to_set([T|Tl], [T|Ts]) :- exexe_list_to_set(Tl, Ts).

exexe_member_eq(T, [A|_As]) :- T == A, !.
exexe_member_eq(T, [_A|As]) :- exexe_member_eq(T,As).


/*?- my_list_to_set([],T).
T = [].

?- my_list_to_set([1,2,x,y,5],T).
T = (1, 2, x, y, 5).

?- my_list_to_set([1, 2, 2, 3, 3, 3], Set).
Set = (1, 2, 3)
*/


% corps d'une clause -> formule conjonctive
body_to_formula((A,B), (A&B1)) :-
    !,
    body_to_formula(B, B1).
body_to_formula(A, A).

/*?- body_to_formula((p(X), q(Y)), Formula).
Formula = p(X)&q(Y)
*/

% (form, liste)-> (liste' => form)

mk_formula(Q,[],Q) :- !.
mk_formula(Q,B, (NewB=>Q)) :-
    mk_formula(B, NewB).

mk_formula([B],B) :- !.
mk_formula([B|Bs], (B&NewBs)) :-
    mk_formula(Bs, NewBs).

/*
 ?- mk_formula(p(X), [q(Y), r(Z)], Formula).
Formula = (q(Y)&r(Z)=>p(X)).

?- mk_formula(p(X), [q(Y)], Formula).
Formula = (q(Y)=>p(X)).

?- mk_formula(p(X), [], Formula).
Formula = p(X)
*/

%Unification

%% unification(X, Y, Sun, Sex): donne la substitution unifiant X et Y, decomposee sur les
%% variables existentielles (Sex) et les variables universelles (Sun)

unification(X, Y, Sun, Sex) :-
    substlist(X, Y, Sun1, Sex1),
    unif(Sun1, Sex1, Sun, Sex).

% subslist(X,Y, Suniv, Sexist): commence l'unification de X et Y,donnant des paires
% variable (existentielle ou universelle) - terme, sans instancier; pas d'occur-check

substlist(X,Y,[],[]) :- X==Y, !.
substlist(X,Y,[], [(X,Y)]) :- varexist(X,_), !.
substlist(X,Y,[],[(Y,X)]) :- varexist(Y,_), !.
substlist(X,Y, [(X,Y)],[]) :- var(X), !.
substlist(X, Y, [(Y,X)],[]) :- var(Y), !.
substlist(X,Y,Suniv, Sexist) :- same_functor(X, Y,_,N), !,

substl_args(N,X, Y, Suniv, Sexist).
substl_args(0,_,_,[],[]).
substl_args(N,X,Y, Suniv, Sexist) :-
    N > 0,
    arg(N,X, ArgX), arg(N,Y, ArgY),
    substlist(ArgX, ArgY, Sun1, Sex1),
    N1 is N-1,
    substl_args(N1, X, Y, Sun2, Sex2),
    merge(Sun1, Sun2, Suniv),
    merge(Sex1, Sex2, Sexist).

/* exemple unification*/

/*
?- unification([f(X,Y), ?(Z)], [ f(c,b),a], Sun, Sex).
Sun = [(Y, b), (X, c)],
Sex = [(?(Z), a)] ;
false.

?- unification(g(Z), a, Sun, Sex).
Sun = [],
Sex = [(g(Z), a)] ;
false.

?- unification(f(X, Y), f(c,b), Sun, Sex).
Sun = [(Y, b), (X, c)],
Sex = [] ;
false.*/


%% exemple subslist
/*
?- substlist(f(a, X), f(a, Y), Sun, Sex).
Sun = [(X, Y)],
Sex = [] ;
false.

?- substlist(?(X), a, Sun, Sex).
Sun = [],
Sex = [(?(X), a)].

?- substlist(f(a, ?(X)), f(a,Y), Sun, Sex).
Sun = [],
Sex = [(?(X), Y)] ;
false.
*/


% substlist_l: substlist applique a une liste de paires

substlist_l([],[],[]).
substlist_l([(X,Y)|S], Sun, Sex) :-  substlist(X, Y, Sun1, Sex1),
substlist_l(S, Sun2, Sex2), merge(Sun1, Sun2, Sun), merge(Sex1, Sex2,
Sex).

/*exemple substlist_l*/
/*
?- substlist_l([(a, X), (Z, Y)], Sun, Sex).
Sun = [(Z, Y), (X, a)],
Sex = [] ;
false.

?- substlist_l([(a, X), (Z, ?(Y))], Sun, Sex).
Sun = [(X, a)],
Sex = [(?(Y), Z)].
*/




%unif(Sun,Sex,Suniv,Sexist):l'algorithme d'unification lui-même
unif(Sun, Sex, Su, Se) :-  unif(Sun,Sex,[],[], Su, Se).
unif([(X,T)|Sun], Sex, Suntp, Sextp, Su, Se) :-
    apply([(X,T)], [Sun, Sex, Suntp, Sextp], [Sun1, Sex1, Suntp1, Sextp1]),
    app_pgm(Sun1, Sex1,S), substlist_l(S,Sun2, Sex2),
    unif(Sun2, Sex2, [(X,T)|Suntp1], Sextp1, Su, Se).

unif([], [(X,T)|Sex], Suntp, Sextp, Su, Se) :-
apply_exist([(X,T)], [Sex, Suntp, Sextp], [Sex1, Suntp1, Sextp1]),
    substlist_l(Sex1,Sun2,Sex2),
    unif(Sun2, Sex2, Suntp1, [(X,T)|Sextp1], Su, Se).
unif([],[],Sun, Sex, Sun, Sex).

/*exemple unif*/
/*?
- unif([(a, b)], [(?(X), Y)], Su, Se).
Su = [(a, b)],
Se = [(?(X), Y)] ;
false.

?- unif([(a, b)], [(f(X), Y)], Su, Se).
Su = [(a, b)],
Se = [(f(X), Y)] ;
false.

?- unif([(f(a), b)], [(f(X), Y)], Su, Se).
Su = [(f(a), b)],
Se = [(f(X), Y)] ;
false.

?- unif([(a, a)], [(g(X), Y)], Su, Se).
Su = [(a, a)],
Se = [(g(X), Y)] ;
false.
*/


% unify_decid(A1, A2, VarsG, Suniv, Sexist): donne la substitution
% existentielle unifiant A1 et A2, i.e la substitution universelle doit
% être un simple renommage

unify_decid(A1, A2, VarsG, Suniv, Sexist) :-
    unification(A1, A2, Suniv, Sexist),
    renaming(Suniv, VarsG).
/*
?- unify_decid([f(X,b),g(Y,c), ?(T)], [f(a,b),g(a,c),d],[W,Z], Sun, Sex).
Sun = [(Y, a), (X, a)],
Sex = [(?(T), d)] ;
false.
*/


% renaming (Suniv, VarsG): verifie que la substitution universelle est un renommage
renaming([],_).
renaming([(X,Y)|_], VarsG) :-
    member_var(X, VarsG,_),
    !,
    var(Y).
renaming([_|Suniv], VarsG) :-
    renaming(Suniv, VarsG).

/*
?- renaming([(X, a), (Y, b)], [Z, W]).
true.
*/



% filtre(X, Y, Subst): Y est une instance de X par la substitution Subst;
% map_filtre(X, Yl, Substl); les Yl sont des instances de X par les
% substitutions Substl

filtre(X, Y, Subst) :-  filtre1(X, Y, [], Subst).
filtre1(X, Y, Subst, Subst) :-  var(X), X == Y, !.
filtre1(X, Y, Subst, Subst) :-  var(X), assoc(X, Subst, Z,_), Z == Y, !.
filtre1(X,_,Subst,_) :- var(X),assoc(X, Subst, _, _), !, fail.
filtre1(X, Y, Subst, [(X,Y)|Subst]) :-  var(X), !.
filtre1(_, Y, _, _) :-  var(Y), !, fail.
filtre1(X, Y, Subst, Substfin) :-  same_functor(X, Y, _, N),
        filtre1_arg(N, X, Y, Subst, Substfin).

filtre1_arg(0, _, _, Subst, Subst).
filtre1_arg(N, X, Y, Subst, Substfin) :-
    N > 0,
    arg(N, X, ArgX),
    arg(N, Y, ArgY),
    filtre1(ArgX,ArgY, Subst, Subst1),
    N1 is N - 1,
    filtre1_arg(N1, X, Y, Subst1, Substfin).

/*exemple filtre et filtre1*/
/*
?- filtre(X, X, [(X, a)]).
false.

?- filtre(X, Y, [(X, a)]).
Y = a.

?- filtre(X, a, [(X, a)]).
true.

?- filtre1(X, a, [(X, a), (Y, b)], Subst).
Subst = [(X, a), (Y, b)].

?- filtre1(X, a, [(Z, c)], Subst).
Subst = [(X, a), (Z, c)].
*/



map_filtre(_, [], []).
map_filtre(X, [Y|Yl], [Subst|Substl]) :-
    filtre(X, Y, Subst),
    map_filtre(X, Yl, Substl).

/*exemple MapFiltre(remarque Subst est un liste de liste)*/
/*
?- map_filtre(X, [a,b], [[(X, a)], [(X, b)]]).
true.

?- map_filtre(X, [a,b,c], Subst).
Subst = [[(X, a)], [(X, b)], [(X, c)]]
.*/




% get_subst_prox(A, H, SubstA, SubstH): unifie A et H: A SubstA = H SubstH,
%avec des substitutions separees, restant proches de A

get_subst_prox(A, H, SubstA, SubstH) :-
    unification(H, A, Subst,[]),
    variables(H, VarsH),
    sep_vars(Subst, VarsH, SubstA1, SubstH1),
    nice_rec(SubstA1, VarsH, SubstA2, SubstH2),
    apply(SubstH2, SubstA2, SubstA),
    app_pgm(SubstH1, SubstH2, SubstH).

/*exemple get_subst_prox*/
/*
?- get_subst_prox(f(X, b), f(a, Y), [(X, a)], [(Y, b)]).
true ;
false.

?- get_subst_prox(f(a,Y), f(X, b),SubstA, SubstB).
SubstA = [(Y, b)],
SubstB = [(X, a)] ;
false.

?- get_subst_prox(f(X,Y), f(a, b),SubstA, SubstB).
SubstA = [(Y, b), (X, a)],
SubstB = [] ;
false.
*/



sep_vars([],_,[],[]).
sep_vars([(X,T)|Subst], VarsH, SubstA, [(X,T)|SubstH]) :-
    member_var(X, VarsH,_), !,
    sep_vars(Subst, VarsH, SubstA, SubstH).
    sep_vars([(X,T)|Subst], VarsH, [(X,T)|SubstA], SubstH) :-
    sep_vars(Subst, VarsH, SubstA, SubstH).
/*exemple sep_vars(separer Subst en deux listes distinctes : SubstA
    qui contient les substitutions sans variables de VarsH, et
    SubstH qui contient les substitutions avec des variables de VarsH.)*/

/* ?-  sep_vars([(x, 1), (y, 2), (z, 3), (u, 4)], [x, z], SubstA, SubstH).
SubstA = [(y, 2), (u, 4)],
SubstH = [(x, 1), (z, 3)].*/



% nice_rec(SubstA, Vars, SubstA1, SubstH)
nice_rec([],_,[],[]).
nice_rec([(A,s(H))|SubstA], Vars, [(A,s(A))|SubstA1], [(H,A)|SubstH]) :-
    var(H),
    member_var(H, Vars, Varsrmd), !,
    nice_rec(SubstA, Varsrmd, SubstA1, SubstH).

nice_rec([(A, [H|C])|SubstA], Vars, [(A, [H|A])|SubstA1], [(C,A)|SubstH]) :-
    var(C),
    member_var(C, Vars, Varsrmd), !,
    nice_rec(SubstA, Varsrmd, SubstA1, SubstH).
nice_rec([(A,H)|SubstA], Vars, [(A, H) | SubstA1], SubstH) :-
    nice_rec(SubstA, Vars, SubstA1, SubstH).

/*exemple nice-recexemple, si la tête de SubstA est sous la forme (A, s(H)),
    cela signifie que H est une variable et s(H) est une substitution pour H.
    Dans ce cas, nous verifions si H est une variable de la liste Vars.
    Si c'est le cas, nous retirons H de Vars, mettons a jour SubstA1
    avec la substitution s(A) pour H, et mettons a jour SubstH avec (H, A)()*/

/*?- nice_rec([(X, s(Y)), (Y, [Z|W])] , [Y, Z], SubstA1, SubstH).
SubstA1 = [(X, s(X)), (Y, [Z|W])],
SubstH = [(Y, X)].

?- nice_rec([(X, T), (Y, [Z|W])] , [X, Y, Z], SubstA1, SubstH).
SubstA1 = [(X, T), (Y, [Z|W])],
SubstH = [].

?- nice_rec([(Y, [Z|W])], [X, W], SubstA1, SubstH).
SubstA1 = [(Y, [Z|Y])],
SubstH = [(W, Y)].*/

% prox_subst(Ser, Sex1, Sun): pour la simplification, substitution
% "proche" de la conclusion

prox_subst([],[],[]).
prox_subst([(?(X),Y)|Sex], [(?(X), X) |Sex1], [(Y,X)|Sun]) :-
    var(Y),
    !,
    prox_subst(Sex, Sex1, Sun).

prox_subst([(?(X), Y)|Sex], [(?(X), Y)|Sex1], Sun) :-
    prox_subst(Sex, Sex1, Sun).

/*exemple prox_subst*/
/*
 * ?- prox_subst([(?(a), b)], Subst1, Subst2).
Subst1 = [(?(a), b)],
Subst2 = [].

?- prox_subst([(?(a), X)], Subst1, Subst2).
Subst1 = [(?(a), a)],
Subst2 = [(X, a)].
*/


%Application d'une substitution

%% application de substitutions universelles et existentielles

% apply(Subst, X,Xsubst)

apply(Subst, X,Y) :-
    var(X),
    assoc(X, Subst, Y,_),
    !.
apply(_,X,X) :-
    var(X),
    !.
apply(Subst, X,Y) :-
    same_functor(X,Y,_,N),
    apply_arg(N, Subst, X,Y).

%Exemple apply
/*?- apply([(X, a), (Y, b)], X, Rmd).
Rmd = a.

?- apply([(Y,a)],X,Z).
X = Z

?- apply([(Y,a)],f(X,Y),Z).
Z = f(X, a) ;
false.

?- apply([(X, a), (Y, b)], X, Y).
Y = a.

?- apply([(X, a), (Y, b)], f(X, Y), Rmd).
Rmd = f(a, b) ;
false.

*/

apply_arg(0,_,_,_).
apply_arg(N, Subst, X, Y) :-
       N > 0,
       arg(N, X, ArgX),
       apply(Subst, ArgX, ArgXsubst),
       arg(N, Y, ArgXsubst),
       N1 is N-1,
       apply_arg(N1, Subst, X, Y).

%exemple apply_arg
/*
?- apply_arg(2,[(Z,a),(Y,b)],f(a,b),f(X,Y)).
Y = b,
X = a ;
false.

?- apply_arg(1,[(Z,a),(Y,b)],f(a,b),f(X,Y)).
X = a ;
false.
?- apply_arg(1,[(Z,a),(Y,b)],f(c),f(X)).
X = c ;
false.
X*/

% apply_l(Subst, Xl, Xlsubstl)
apply_l(_,[], []).
apply_l(Subst, [A|L], [Asubst|Lsubst]) :-
    apply(Subst, A, Asubst),
    apply_l(Subst, L, Lsubst).
%apply_l:permet d'appliquer une substitution a une liste d'elements. Voici son fonctionnement :
%Exemple apply_l
/*?- apply_l([(X, a), (Y, b)],[f(X, Y), g(Y, X)],L).
L = [f(a, b), g(b, a)] ;
false.
?- apply_l([(X, a), (Y, b)],[f(X, Y), g(Z, X)],L).
L = [f(a, b), g(Z, a)] ;
false.
*/



% map_apply(Substl, X, Xsubstl)
%map_apply:permet d'appliquer une liste de substitutions a un terme donne.
map_apply([],_,[]).
map_apply([Subst|SubstL], Q, [Qsubst|QsubstL]) :-
    apply(Subst,Q,Qsubst),
    map_apply(SubstL,Q,QsubstL).
%Exempl
/*?- map_apply([[(X,a)],[(Y,b),(Z,c)]],f(X,Y,Z),Result).
Result = [f(a, Y, Z), f(X, b, c)] .

*/

% apply_exist(Sexist, X, Xsubst)
apply_exist(L,X,Y) :-  varexist(X,_), assoc(X,L,Y,_), !.
apply_exist(_,X,X) :-  varexist(X,_), !.
apply_exist(_,X,X) :-  var(X), !.
apply_exist(L,X,Y) :-  same_functor(X,Y,_,N), apply_ex_arg(N,L,X,Y).

apply_ex_arg(0,_,_,_).
apply_ex_arg(N,L,X,Y) :-
    N > 0,
    arg(N,X,ArgX),
    apply_exist(L, ArgX, ArgXsubst), arg(N,Y, ArgXsubst),
    N1 is N - 1,
    apply_ex_arg(N1,L,X,Y).
%Exemple apply_exist
/*?- apply_exist([(?(Z), a), (Y, b)], f(?(Y), g(?(Z), ?(W))), Y_result).
Y_result = f(?(Y), g(a, ?(W))) .

?- apply_exist([(?(Z), a), (Y, b)], f(Y, g(?(Z), ?(W))), Y_result).
Y_result = f(Y, g(a, ?(W))) .

?- apply_exist([(?(Z), a), (?(W), b)], f(Y, g(?(Z), ?(W))), Y_result).
Y_result = f(Y, g(a, b)) .
*/


/*
% apply_exist_l(Sexist, Xl, Xsubstl)
apply_exist_l(_,[],[]).
apply_exist_l(Subst,[A|L],[Asubst|Lsubst]) :-  apply_exist(Subst, A, Asubst),
    apply_exist_l(Subst, L, Lsubst).
*/
%Exemple apply_exist_l
/*?- apply_exist_l([(?(Z), a), (?(W), b)], [f(Y, g(?(Z), ?(W)))], Y_result).Y_result = [f(Y, g(a, b))] .

?- apply_exist_l([(?(Z), a), (?(W), b)], [f(Y, g(?(Z), ?(W))),h(?(W))], Y_result).
Y_result = [f(Y, g(a, b)), h(b)] .
*/

%map_apply_exist(Sexistl, Xl, Xlsubstl)

map_apply_exist([],[],[]).
map_apply_exist([Subst|SubstL], [Q|Ql], [Qsubst|QsubstL]) :-
    apply_exist(Subst, Q, Qsubst), map_apply_exist(SubstL, Ql,
    QsubstL).
%exemple map_apply_exist
/*?- map_apply_exist([[(?(X),a)],[(?(Y),b)]],[[f(?(X),?(Y))],[g(?(Y))]],Xsubst).Xsubst = [[f(a, ?(Y))], [g(b)]] .

?- map_apply_exist([[(?(X),a)],[(?(Y),b)]],[[f(?(X))],[g(?(Y))]],Xsubst).
Xsubst = [[f(a)], [g(b)]] .
*/

%A.2 Regles d'execution etendue et extraction associee

%Remplacement et reduction

replace(Occ,F,G,Grepl) :-
    replace1(Occ,F,G, Gnew),
    reduce(Gnew, Grepl).

%% replace1(Occ, B, G, Gnew): remplacement dans une formule G a une occurrence Occ

%% par une formule B, pour obtenir Gnew

replace1([],B,_,B1) :-  B == B1, !.
replace1([],B,_,B1) :-  var(B1), B = B1.
replace1([],B,_,B1) :-  var(B), B = B1.
replace1([1|Occl],B,(G1 & G2), (Gnew & G2)) :- replace1(Occl,B, G1, Gnew).
replace1([2|Occl],B,(G1 & G2), (G1 & Gnew)) :- replace1(Occl, B, G2, Gnew).
replace1([1|Occl],B,(G1 => G2), (Gnew => G2)) :- replace1(Occl, B, G1, Gnew).
replace1([2|Occl],B,(G1 => G2), (G1 => Gnew)) :- replace1(Occl, B, G2, Gnew).

%exemple replace et replace1
/*
?- replace1([1], true,(G1 & G2), Result).
Result = true&G2 .

?- replace1([2], false,(G1 => G2), Result).
Result = G1=>false .
?- replace([2], false, ((a & true) => b), Grepl).
Grepl = -a .

?- replace([1], false, ((a & true) => b), Grepl).
Grepl = true .
*/



% get subformula(G,SubG): donne une sous-formule SubG de G

%get_subformula(G, SubG) :- replace1(_,SubG,G,G).
%exemple get_subformula
/*?- get_subformula(a&b,SubG).
SubG = a&b ;
SubG = a ;
SubG = b ;
false.

?- get_subformula((a & (b => c)), SubG).
SubG = a&b=>c ;
SubG = a ;
SubG = b=>c ;
SubG = b ;
SubG = c ;
false.
*/



% get_occ(G,SubG,Occ): donne l'occurence de SubG dans G

get_occ(G,SubG,Occ) :- replace1(Occ, SubG,G,G).
%exemple
/*?- get_occ((a & (b => c)), b, Occ).
Occ = [2, 1] ;
*/

is_atom((_ => _)) :-  !, fail.
is_atom((_ & _ )) :- !, fail.
is_atom(-(_)) :- !, fail.
is_atom(F) :-  functor(F,_,N), N > 0.
%exemple is_atom
/*
?- is_atom(_).
false.

?- is_atom(a).
false.

?- is_atom(A).
false.

?- is_atom(f(a,b)).
true.

?- is_atom(f(a)).
true.
*/

% reduce1(X, Y): une passe de reduction par rapport a true et false de X
% a Y par les regles;

% reduce(X, Y): reductions successives jusqu'a stabilite

reduce(X,Y) :-  reduce1(X,Y), X == Y , !.
reduce(X, Y) :-  reduce1(X, X1), reduce(X1,Y).


reduce1(X,X) :- var(X), !.
reduce1(?(X),?(X)) :-  !.
reduce1((X & Y),G) :- X == true, !, reduce(Y,G).
reduce1((X & Y),G) :-  Y == true, !, reduce(X,G).
reduce1((X & _),false) :-  X == false, !.
reduce1((_ & Y), false) :-  Y == false, !.
reduce1((X & Y), (X1 & Y1)) :- !, reduce1(X,X1), reduce1(Y,Y1).
reduce1((X => Y),G) :-  X == true, !, reduce(Y,G).
reduce1((_ => Y), true) :-  Y == true, !.
reduce1((X => _),true) :-  X == false, !.
reduce1((X => Y), -(G)) :- Y == false, !, reduce(X,G).
reduce1((X => Y), (X1 => Y1)) :-!, reduce1(X,X1), reduce1(Y,Y1).
reduce1(X,X).
%exemple reduce1
/*?- reduce1(true & false, Result).
Result = false.

?- reduce1(true => false, Result).
Result = false.

?- reduce1(true => true, Result).
Result = true.

?- reduce1(?(X) & true, Result).
Result = ?(X).

?- reduce1(X, Result).
X = Result.

?- reduce1(?(X), Result).
Result = ?(X).

?- reduce1(X => true, Result).
Result = true.

?- reduce1(X => false, Result).
Result = -X.
?- reduce1((true => false) & (?(Y) => true) & (?(Z) & true), Result).
Result = false&true& ?(Z).
*/









%Inference de clause definie


%% dci (G,SubG, Occ, H, Subst, Newvars): dans G, a SubG et Occ, donnant
%% H et Subst

% dci etendue: avec toutes les clauses

dci_ext(F,SubF, Occ, Hl, Substl, Intvarsl) :- setof((H,Subst, Intvars),
     dci_cl(F,SubF, Occ, H, Subst, Intvars),S),
     separate(S, Hl, Substl, Intvarsl).

% dci_cl: avec une clause

dci_cl(G, SubG, Occ, H, Subst, Internal_vars) :-  same_functor(SubG, Head,_,_),
    copy(Head, Head1), clause(Head1, Body),
    dci(G,SubG,Occ, Head1, Body, H, Sexist, Internal_vars),
    normalize(Sexist, Subst).

% dci_clause: laisse le choix de la clause

dci_clause(G,SubG, Occ, H, Subst, Internal_vars) :-  rm_exist(SubG, Head),
    choose_clause(Inhibit, Head, Head_inst, Body), Inhibit == ok, !,
    dci(G,SubG, Occ, Head_inst, Body, H, Sexist, Internal_vars),
    normalize(Sexist, Subst).

%dci: une fois l'atome et la clause choisis

dci(G,SubG,Occ, Head, Body, H, Sexist2, Internal_vars) :-  %gtrace,
    internal_vars(Head, Body, Internal_vars),
    body_to_formula(Body, Body1),
    variables(G, VarsG),
    unify_decid(Head, SubG, VarsG, Suniv, Sexist),
    apply(Suniv, Sexist, Sexist1),
    subst_new(SubG, Sexist1, Sexist2, Newvars, Sbody),
    apply(Sbody, Body1, Body2),
    app_pgm(Internal_vars, Newvars, Varsl),
    mk_subst(Varsl, Subst),
    apply_exist(Sexist2,G,G1),
    replace(Occ, Body2, G1, G2),
    apply(Subst, G2, G3),
    apply(Suniv, G3,H).


mk_subst([X|L], [(X,?(X))|S]) :-  mk_subst(L,S).
mk_subst([],[]).
%exemple mk_subst
/*
?- mk_subst([X, Y, Z], Subst).
Subst = [(X, ?(X)), (Y, ?(Y)), (Z, ?(Z))].

*/



subst_new(F, Sexist, NewSexist, Newvars, Sbody) :-  variables(F, VarsF),
    vars2(Sexist,VarsS1), minus1(VarsS1, VarsF, Newvars),
    change(Sexist, Newvars, Sex1, Sbody), apply(Sbody, Sex1, NewSexist).

vars2([],[]).
vars2([(X,Y)|S], Vars) :-  varexist(X,_), vars2(S,Vars1), variables(Y, VarsY),
    merge(Vars1, VarsY, Vars).
%exemple Vars2
/*
?- vars2([(?(X),Y),(?(Z),T)], Vars).
Vars = [Y, T].

?- vars2([(?(X),Y),(?(Z),T),(?(D),i)], Vars).
Vars = [Y, T] ;
?- vars2([(?(X),p(a,Y))], Vars).
Vars = [Y] ;
false.

?- vars2([(?(X),p(X,f(W,Z)))], Vars).
Vars = [X, W, Z] ;

*/


change(S,[],S,[]).
change(S, [X|V], Snew, Sbody) :-  change1(S,X, S1, Sbody1),
    change(S1, V, Snew, Sbody2), merge(Sbody1, Sbody2, Sbody).

change1([(X,Y)|S],V,S,[(V,X)]) :-  V == Y, !.
change1([(X,Y)|S],V,[(X, Ynew)|S], [(V,X)]) :-  variables(Y, [VarY]),
    V == VarY, !, apply([(V,X)], Y, Ynew).
change1([(X, [A|Y])|S],V,[(X,[A|X])|S], [(V,X)]) :-  V == Y, !.
change1([(X,Y)|S], V, [(X,Y)|Snew], Sbody) :- change1(S,V, Snew, Sbody).
change1([],_,[],[]).

%%exemple change et change1
/*
?- change1([(X, [a|b]), (T, [d,e,f,g]), (Z, [h,i,j,k])],b, Snew, Sbody).
Snew = [(X, [a|X]), (T, [d, e, f, g]), (Z, [h, i, j, k])],
Sbody = [(b, X)].

?- change([(Z,c),(T, ?(Z)), (Y, a),(?(Y),X)], [Y,c], Snew, Sbody).
Snew = [(T, ?(Z)), (Y, a), (?(Y), X)],
Sbody = [(c, Z)] .

?- change([(Z,c),(T, ?(Z)), (a, Y),(?(Y),W)], [Y,c,W], Snew, Sbody).
Snew = [(T, ?(Z))],
Sbody = [(W, ?(Y)), (c, Z), (Y, a)] .

?- change1([(X, [a|b]), (T, [d,e,f,g]), (Z, [h,i,j,k])],b, Snew, Sbody).
Snew = [(X, [a|X]), (T, [d, e, f, g]), (Z, [h, i, j, k])],
Sbody = [(b, X)].

?- change1([(Y,Z),(T,a)],Z, Snew, Sbody).
Snew = [(T, a)],
Sbody = [(Z, Y)].

?- change1([(Z,Y),(T,a)], Z, Snew, Sbody).
Snew = [(Z, Y), (T, a)],
Sbody = []

*/




%exist -> normal
/*
normalize(X,X) :-  var(X), !.
normalize(X,Y) :-  varexist(X,Y), !.
normalize(X,Y) :-  same_functor(X,Y,_,N), norm_args(N,X,Y).
    norm_args(0,_,_).
    norm_args(N,X,Y) :-
    N > 0,
    arg(N, X, ArgX),
    normalize(ArgX, ArgY),
    arg(N, Y, ArgY),
    N1 is N - 1,
    norm_args(N1, X,Y).*/

normalize(X, X) :- var(X), !.
normalize(X, Y) :- varexist(X, Y), !.
normalize(X, Y) :- same_functor(X, Y, _, N), norm_args(N, X, Y).

norm_args(0, _, _).
norm_args(N, X, Y) :-
    N > 0,
    arg(N, X, ArgX),
    normalize(ArgX, ArgY),
    arg(N, Y, ArgY),
    N1 is N - 1,
    norm_args(N1, X, Y).
%exemple normalize et norm_args
/*?- normalize(?(X),R).
X = R.

?- normalize(X,R).
X = R.

?- normalize(f(X,g(Y,a)),R).
R = f(X, g(Y, a)) .

?- norm_args(2,f(X,Y),f(c,d)).
X = c,
Y = d .

?- norm_args(2,f(a,b),f(_,_)).
true .*/


%Simplification
%% simplification(G,SubG1, Occ1, SubG2, Occ2, H, Subst): un atome
%% positif SubG1 a Occ1
%% et un atome negatif SubG2 a Occ2 sont unifies;
%% les deux atomes sont supprimes;
%% simplification "un pour un"

simplification(G, SubG1, Occ1, SubG2, Occ2, H, Suniv) :-
    unification(SubG1, SubG2,[],Sex), prox_subst(Sex, Sex1, Suniv),
    apply(Suniv, G, G_1), apply_exist(Sex1, G_1, G_2),
    simplif(Occ1, Occ2,G_2, true, H).

% simplif(Occ1, Occ2,G_2, false, H2): H2 = true pour des buts implicatifs
simplif(Occ1, Occ2, G, Bool, H) :-  replace1(Occ1, Bool, G, G1),
    replace(Occ2, Bool, G1, H).

%% simplification gardant l''atome negatif

simplification_kp(G,SubG1, Occ1, SubG2, _, H, Suniv) :-
    unification(SubG1, SubG2,[],Sex), prox_subst(Sex, Sex1, Suniv),
    apply(Suniv, G, G_1), apply_exist(Sex1, G_1, G_2),
    replace(Occ1, true, G_2,H).


% Induction par point fixe

%% induction (A,Q, H_l, Tau_l, Thetal_l): donne les lemmes H_l de la
%% preuve de A => Q
%% par induction par point fire, en backtrackant sur one induction

induction(A, Q, H_l, Tau_l, Thetal_l) :-
    bagof((H, T, Thl), one_induction(A, Q, H, T, Thl), L),
    sep_inv(L, H_l, Tau_l, Thetal_l).

%% one_induction(A,Q, H, Tau, Thetal): l'un des lemmes

one_induction(A,Q, H, Tau, Thl) :-  variables(A, VarsA),
    get_cst_induc(Q, VarsA, List),
    one_induc_gen(A, Q, Tau, Thl1, Qtau, Qthetal, Nonrec),
    chge_cst(Thl1, List, Thl, Templ),
    map_apply_exist(Templ, Qthetal, Qthetal1),
    app_pgm(Qthetal1, Nonrec, Inducbody),
    mk_formula(Qtau, Inducbody, H1), rm_cst(H1, H).

varcst(X, Y) :-  nonvar(X), functor(X, *, 1), arg(1, X, Y).
get_cst_induc(X, V, []) :-  var(X), member_var(X,V,_),!.
get_cst_induc(X, _, [X]) :-  var(X), !.
get_cst_induc(X, _, []) :-  varexist(X, _), !.
get_cst_induc(X, _, []) :-  functor(X, _, 0), !.
get_cst_induc(X, V, L) :- functor(X, _, N), gci(N, X, V, L).

gci(0, _, _, []).
gci(N, X, V, L) :-
    N> 0,
    arg(N, X, ArgX),
    get_cst_induc(ArgX, V, L1),
    N1 is N - 1,
    gci(N1, X, V, L2),
    merge(L1, L2, L).

chge_cst([Th|Thl],L, [Th1|Thl1], [Temp|Templ]) :- chge_c(Th, L, Th1, Temp),
       chge_cst(Thl, L, Thl1, Templ).

chge_cst([], _, [], []).
chge_c([(X, Y)|Thl],L,Thl1, [(Y, X)|Temp]) :-  member_var(X, L, _), !,
    chge_c(Thl, L, Thl1, Temp).
chge_c([C|Thl], L, [C|Thl1], Temp) :- chge_c(Thl, L, Thl1, Temp).
chge_c([], _, [], []).

rm_cst(X, Y) :-  varcst(X, Y),!.
rm_cst(X, X) :-  var(X), !.
rm_cst(X, X) :-  varexist(X, _), !.
rm_cst(X, Y) :-  same_functor(X, Y, _, N), rm_cst(N, X, Y).
rm_cst(N, X, Y) :-
    N > 0,
    arg(N,X, ArgX), rm_cst(ArgX, ArgY),
    arg(N, Y, ArgY),
    N1 is N - 1, rm_cst(N1,X,Y).
    rm_cst(0,_,_).

one_induc_gen(A, Q, Sa, Thetal, Qinst, Qthetal, Nonrec) :-
    same_functor(A, B, F, N), clause(B, Body_set),
    set_to_list(Body_set, Body), get_subst_prox(A,B,Sa,Sb),
    apply_l(Sb,Body, Bodyinst), decomp(Bodyinst, F, N, Rec, Nonrec),
    map_filtre(A, Rec, Thetalist),
    map_new_unifier(Thetalist, Q, A, Thetal1, Thetaexistl),
    map_normalize(Thetaexistl, Thetal2), map_app(Thetal1, Thetal2, Thetal),
    apply(Sa, Q, Qinst), map_apply(Thetal1, Q, Qthetal1),
    map_apply_exist(Thetaexistl, Qthetal1, Qthetal).

map_normalize([], []).
map_normalize([A|X], [B|Y]) :-  normalize(A, B), map_normalize(X, Y).

map_app([], [], []).
map_app([A|X], [B|Y], [C|Z]) :-  app_pgm(A, B, C), map_app(X, Y, Z).

sep_inv(L, Hl, Tl, Thl) :-  rev_pgm(L,L1), separate(L1, Hl, Tl, Thl).

    separate([], [], [],[]).
    separate([(H, T, Th)|L], [H|Hl], [T|Tl], [Th|Thl]) :-  separate(L, Hl, Tl, Thl).

% decomp(Body, F, N, Rec, Nonrec): donne la partie recursive Rec (ie de foncteur F,N),
% et le reste Nonrec du corps Body

decomp([],_,_,[],[]).
decomp([B|Bl], F, N, [B|Rec], Nonrec) :-  functor(B, F, N), !,
    decomp(Bl, F, N, Rec, Nonrec).

decomp([B|Bl], F, N, Rec, [B|Nonrec]) :-  decomp(Bl, F, N, Rec, Nonrec).

%% new_unifier(Subst, Q, A, Newsubst, Substexist): les variables
%% existentielles de Q-A
%% sont remplacees par de nouvelles variables (Newsubst),
%% et les variables universelles de Q-A deviennent ezistentielles
%% (Substexist)

new_unifier(Subst, Q, A, Newsubst, Substexist) :-  varsout((Q => A), A, Exist, Univ),
    renew(Subst, Univ, Newsubst), create(Exist, Substexist).

varsout(Q, A, Exist, Univ) :-  vars(Q, ExistQ, UnivQ), vars(A,ExistA, UnivA),
    minus(ExistA, ExistQ, Exist), minus(UnivA, UnivQ, Univ).

renew(Subst,[], Subst).
renew(Subst, [X|Univ], [(X,?(_))|Newsubst]) :-  renew(Subst, Univ, Newsubst).

create([],[]).
create([X|Exist], [(X,_)|Subst]) :-  create(Exist, Subst).

map_new_unifier([], _, _, [], []).
map_new_unifier([Subst|SubstL], Q, A, [Newsubst|NewsubstL], [Exist|ExistL]) :-
    new_unifier(Subst, Q, A, Newsubst, Exist),
    map_new_unifier(SubstL, Q, A, NewsubstL, ExistL).

%Induction structurelle

%% structural (F, X, Type, G_l, ConcL_l, HypL_l)
%% donne les lemmes G_L de la preuve par induction structurelle sur X:Type

structural(F, X, tree(X),
            [Gb, ((Ghyp1&Ghyp2)=>Gconcl)],
            [[(X, leaf(X))], [(X, branch(Y1,Y2))]],[[],
            [[(X,Y1)|Subst1], [(X,Y2)|Subst2]]]) :-
    !,
    get_univ(F, Csts1),
    get_univ(F, Csts2),
    apply_exist(Csts1, F, F1),
    apply_exist(Csts2, F, F2),
    rm_exist(Csts1, Subst1),
    rm_exist(Csts2, Subst2),
    apply([(X,Y1)], F1, Ghyp1),
    apply([(X,Y2)], F2, Ghyp2),
    apply([(X, branch(Y1, Y2))],F,Gconcl),
    apply([(X, leaf(X))],F,Gb).

structural(F, X, Type, G_l, Concl_l, Hypl_l) :-
    clauses(Type, Rec, Base_l),
    struct_base(F, X, Base_l, Gb_l, Conclb_l, Hyplb_l),
    struct_rec(F, X, Rec, Grec, Conclrec, Hyplrec),
    app_pgm(Gb_l, [Grec], G_l),
    app_pgm(Conclb_l, [Conclrec], Concl_l),
    app_pgm(Hyplb_l, [Hyplrec], Hypl_l).

clauses(Head, Rec, Base_l) :-  setof((Head, Body), clause(Head, Body),S),
    mk_clauses(S, Rec, Base_l).

  mk_clauses([(Head, true)|S], Rec, [Head|Base_l]) :-  !,
    mk_clauses(S, Rec, Base_l).

  mk_clauses([(Head, Body)|S], [(Head :- Body)|Rec], Base_l) :-
    mk_clauses(S, Rec, Base_l).
    mk_clauses([],[], []).

struct_base(_, _, [], [], [], []).
struct_base(F, X, [Base|Base_l], [Gb|Gb_l], [[(X,Xb)]|Conclb_l],[[]|Hyplb_l]) :-
    arg(1, Base, Xb), apply([(X,Xb)], F, Gb),
    struct_base(F, X, Base_l, Gb_l, Conclb_l, Hyplb_l).

struct_rec(F, X, [(Concl :- Hyp)], (Ghyp => Gconcl),[(X,Xc)], [[(X,Xh)|Subst]]) :-
    arg_nice(1, Concl,X,Xc), arg_nice(1, Hyp, X, Xh), get_univ(F,Csts),
    apply_exist(Csts, F, Fc), rm_exist(Csts, Subst),
    apply([(X, Xh)], Fc, Ghyp),
    apply([(X, Xc)], F, Gconcl).

   arg_nice(N, F, X,YX) :-  arg(N, F, Y), nice(X,Y,YX).

   nice(X,Y,X) :-  var(Y), !.
   nice(X,s(Y),s(X)) :- var(Y), !.
   nice(X,s(s(Y)), s(s(X))) :- var(Y), !.
   nice(X, [A|Y], [A|X]) :-  var(Y), !.
   nice(X, [A,B|Y], [A,B|X]) :-  var(Y), !.

 get_univ(X, [(X, _)]) :-  varexist(X, _), !.
 get_univ(X, []) :-  var(X), !.
 get_univ(X,[]) :-  functor(X, _, 0), !.
 get_univ(F,Csts) :- functor(F, _, N), getunivargs(N,F,Csts).

   getunivargs(0, _, [] ).
   getunivargs(N,F,Csts) :-
      N > 0,
      arg(N, F, ArgF), get_univ(ArgF, Csts1),
      N1 is N - 1,
      getunivargs(N1,F,Csts2), merge_s1(Csts1, Csts2, Csts).

   merge_s1([], S, S).
   merge_s1([(X, Y)|S1], S2, [(X,Y)|S]) :-
       assoc(X,S2,_,  S2rmd), !,
       merge_s1(S1, S2rmd,S).
   merge_s1([(X, Y)|S1], S2, [(X,Y)|S]) :-  merge_s1(S1,S2,S).

/*
?- trace,struct_rec(double(_32562, _32564)=>add(_32562, ?(_32568), _32564), _32562, [(nat(s(_38744)):-nat(_38744))], _41554, _41556, _41558).
   Call: (13) struct_rec(double(_14268, _14270)=>add(_14268, ?(_14274), _14270), _14268, [(nat(s(_14294)):-nat(_14294))], _14324, _14326, _14328) ? creep
   Call: (14) arg_nice(1, nat(s(_14294)), _14268, _16308) ? skip
   Exit: (14) arg_nice(1, nat(s(_14294)), _14268, s(_14268)) ? creep
   Call: (14) arg_nice(1, nat(_14294), _14268, _16326) ? skip
   Exit: (14) arg_nice(1, nat(_14294), _14268, _14268) ? creep
   Call: (14) get_univ(double(_14268, _14270)=>add(_14268, ?(_14274), _14270), _19744) ? skip
   Exit: (14) get_univ(double(_14268, _14270)=>add(_14268, ?(_14274), _14270), [(?(_14274), _20626)]) ? creep
   Call: (14) apply_exist([(?(_14274), _20626)], double(_14268, _14270)=>add(_14268, ?(_14274), _14270), _21576) ? skip
   Exit: (14) apply_exist([(?(_14274), _20626)], double(_14268, _14270)=>add(_14268, ?(_14274), _14270), double(_14268, _14270)=>add(_14268, _20626, _14270)) ? creep
   Call: (14) rm_exist([(?(_14274), _20626)], _16320) ? skip
   Exit: (14) rm_exist([(?(_14274), _20626)], [(_14274, _20626)]) ? creep
   Call: (14) apply([(_14268, _14268)], double(_14268, _14270)=>add(_14268, _20626, _14270), _16294) ? skip
   Exit: (14) apply([(_14268, _14268)], double(_14268, _14270)=>add(_14268, _20626, _14270), double(_14268, _14270)=>add(_14268, _20626, _14270)) ? creep
   Call: (14) apply([(_14268, s(_14268))], double(_14268, _14270)=>add(_14268, ?(_14274), _14270), _16296) ? skip
   Exit: (14) apply([(_14268, s(_14268))], double(_14268, _14270)=>add(_14268, ?(_14274), _14270), double(s(_14268), _14270)=>add(s(_14268), ?(_14274), _14270)) ? creep
   Exit: (13) struct_rec(double(_14268, _14270)=>add(_14268, ?(_14274), _14270), _14268, [(nat(s(_14294)):-nat(_14294))], (double(_14268, _14270)=>add(_14268, _20626, _14270))=>double(s(_14268), _14270)=>add(s(_14268), ?(_14274), _14270), [(_14268, s(_14268))], [[(_14268, _14268), (_14274, _20626)]]) ? creep
_41554 = (double(_32562, _32564)=>add(_32562, _A, _32564))=>double(s(_32562), _32564)=>add(s(_32562), ?(_32568), _32564),
_41556 = [(_32562, s(_32562))],
_41558 = [[(_32562, _32562), (_32568, _A)]] .
*/


%A.3 Boucle interactive

%Entree dans la boucle

exexe :-
   %write_exexe('?- exexe.'),
   give_pred_types(File), nl,
    %write('>>>> Taper: g. ou: guide. pour obtenir de l''aide.'), nl, nl,
    handle_begin([File]).

%% handle_begin

handle_begin(File_l) :-  give_rule(Rule),extex(Rule, File_l).

give_pred_types(File) :-
    %write('>>>> Le fichier "pred" contient des definitions de predicats et des'), nl,
    %write('>>>> informations sur leurs arguments (types, definition recursive).'), nl,
    %write('>>>> Preferez-vous charger un autre fichier (o/n)? '),
    %read(Rep), other_file(Rep, File).
    File=File, true.


%other_file(o, File) :- !, write('>>>> Nom de fichier: '), read(File), nl,
%    name(File, Filen), app_pgm(Filen," > tempexex",S1),
%    app_pgm("cat",S1,S), system(S), system("cp tempexex tempexex1"),
%    %reconsult(File).
%    true.
%other_file(_,pred) :-  nl,
%    %system("copy pred tempexex"),
%    %system("copy tempexex tempexex1"),
%    %reconsult(pred).
%    true.


% guide

extex(guide, File_l) :-
    write_exexe('>>>> Regles disponibles: '), nl,
    write_exexe('charge (F)   : charge le fichier F'), nl,
    write_exexe('charge       : donne les noms des fichiers charges'), nl,
    write_exexe('def(P)       : donne les clauses definissant le predicat P'), nl,
    write_exexe('info (P)     : donne les clauses concernant le predicat P'), nl,
    write_exexe('info(P,F)    : donne les clauses du fichier F concernant le predicat P'), nl,
    write_exexe('exemple      : montre un exemple de session'), nl,
    write_exexe('init         : demande le but initial a prouver'), nl,
    write_exexe('f, fin       : fin'), nl, nl, !, handle_begin(File_l).
extex(g,File_l) :-  extex(guide, File_l).

% info

extex(info(P), File_l) :-  system_info(P,tempexex), !, handle_begin(File_l).
extex(info(P,File), File_l) :-  system_info(P,File), !, handle_begin(File_l).
extex(info(P,File), File_l) :-
    write_exexe(' Il n y a aucune clause dans '), write_exexe(File),
    write_exexe(' concernant le predicat'), write_exexe(P), nl, nl, !,
    handle_begin(File_l).

system_info(P,File) :-  P=P, File=File,%name(P, Pname), name(File, Filename),
    true.
    %app_pgm(Pname,"' | egrep 'type rec'", F1), app_pgm(" | egrep '", F1,F2),
    %app_pgm(Filename, F2,F3), app_pgm("'; more ",F3,F4),
    %app_pgm(Pname, F4, F5), app_pgm(" | egrep '^", F5, F6),
    %app_pgm(Filename, F6,F7), app_pgm("more ",F7, F), system(F), nl.
% exemple

extex(exemple, File_l) :-  exexe_system("more session"), !, handle_begin(File_l).

% charge(F)

extex(charge(F), File_l) :-
    exexe_system("mv tempexex tempexex1"),
    my_name(F, Fn),
    app_pgm(Fn," tempexex1 > tempexex", S1), app_pgm("cat",S1,S),
    exexe_system(S), reconsult(tempexex), !, nl, app_pgm(File_l, [F], File_l1),
    handle_begin(File_l1).

% def

extex(def(P), File_l) :-  handle_def(P), !, handle_begin(File_l).

handle_def(P) :-  bagof(Cl, is_clause(P,Cl), Cl_1), rm_true(Cl_1,Cl_l1),
    write_cl(Cl_l1), !.

handle_def(P) :-
    write_exexe('Il n''y a pas de clause chargee definissant le predicat '),
    write_exexe(P), nl, nl.

is_clause(P, (Head:-Body)) :-
    mb(N, [1,2,3,4,5]),
    functor(Head,P,N),
    clause(Head, Body).

rm_true([],[]).
rm_true([(H:-true)|C1], [H|Cl1]) :-  !, rm_true(C1,Cl1).
rm_true([C|C1], [C|Cl1]) :-  rm_true(C1, Cl1).

write_cl([]) :-  nl.
write_cl([C|Cl]) :-
    variables(C,Vars), write_varnames(C,Vars), nl,
    write_cl(Cl).

% charge

extex(charge, File_l) :-  write_exexe('fichiers charges: '), write_l(File_l), !,
    handle_begin(File_l).
write_l([F]) :-  !, write_exexe(F), write_exexe('.'), nl, nl.
write_l([F|Fl]) :-  write_exexe(F), write_exexe(', '), write_l(Fl).

% fin

extex(fin,_) :-  exexe_system("rm tempexex tempexex1"), !, write_exexe('ciao.'), !.

extex(fin,_) :-  write_exexe('sortie de secours'), !.

extex(f,_) :-  extex(fin,_).

% init

extex(init, File_l) :-  give_formula(F), variables(F, Vars),
    test_echo(F, Vars), !, extex(F, ['F'],0,Vars,[],[],F,File_l).

test_echo((H=>C), Vars) :-  !, vars(H,[],_), write_exexe(' but initial:'),
    write_ls_f([((H=>C), ['F'])], Vars).
test_echo(F, Vars) :-  write_exexe(' but initial:'),
    write_ls_f([(F,['F'])], Vars).

% give formula

give_formula(Formula) :-
    write_exexe1('>>>>>>>>>> BUT:'), nl, nl, read_exexe(Formula), nl.
% erreur

%extex(_, File_l) :-
    write_exexe('>>>> Utiliser: init, f/fin, g/guide, exemple, charge (F), charge,'), nl,
    write_exexe('>>>> def(P), info(P), info(P,F).'), nl, nl,
    handle_begin(File_l).

%%Corps de la boucle

extex(F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    F\==nil,
    give_rule(Rule),
    handle_ctnue(Rule, F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l).
   %give_rule1(Rule) :-  write_exexe1('>>>>>>>>>> REGLE: '), read_exexe(Rule), nl.
% get a rule

give_rule(Rule) :-  write_exexe1('>>>>>>>>>> REGLE: '), read_exexe(Rule), nl.

%%Gestion des regles: divers

%% handle_ctnue(Rule, F, Name, Nmax, Vars, Flist, Plist, Fdep)

% guide

handle_ctnue(guide, F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    write_exexe('>>>> Regles disponibles : '), nl,
    write_exexe('charge     : donne les noms des fichiers charges'), nl,
    write_exexe('def(P)     : donne les clauses definissant le predicat P'), nl,
    write_exexe('info(P)    : donne les clauses concernant le predicat P'), nl,
    write_exexe('info (P,F) : donne les clauses du fichier F concernant le predicat P'), nl,
    write_exexe('exemple    : montre un exemple de session'), nl,
    write_exexe('exist      : existentialise une variable de la conclusion'), nl,
    write_exexe('univ       : universalise une variable de la conclusion'), nl,
    write_exexe('transfert  : transfere un atome dans la conclusion'), nl,
    write_exexe('struct     : regle d''induction structurelle'), nl,
    write_exexe('comput     : regle d''induction par point fixe'), nl,
    write_exexe('dci        : regle de \'definite clause inference\' non-deterministe'), nl,
    write_exexe('dcir       : regle de \'definite claune inference\' deterministe'), nl,
    write_exexe('simpl      : regle de simplification'), nl,
    write_exexe('simpl_g    : regle de simplification gardant 1''hypothese '), nl,
    write_exexe('postul     : regle de postulat '), nl,
    write_exexe('b_cour     : affiche le but courant'), nl,
    write_exexe('b_init     : affiche le but initial'), nl,
    write_exexe('tout_es    : affiche la liste de toutes les clauses generees'), nl,
    write_exexe('depl_es    : affiche la liste des clauses generees, apres depliages'),nl,
    write_exexe('efface     : fin de la preuve, on peut alors en commencer une autre'), nl,
    write_exexe('recomm     : recommence la preuve, avec le meme but initial'), nl,
    write_exexe('f, fin     : fin'), nl, nl, !,
    extex(F,Name, Nmax, Vars, Flist, Plist, Fdep, File_l).




handle_ctnue(g, F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    handle_ctnue(guide,F,Name,Nmax,Vars,Flist,Plist,Fdep,File_l).

% def(P)

handle_ctnue(def(P),F, Name,Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    handle_def(P), ! , extex(F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l).

% info(P)

handle_ctnue(info(P), F, Name,Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    system_info(P,temperex), !,
    extex(F,Name,Nmax,Vars,Flist,Plist,Fdep,File_l).

% charge

handle_ctnue(charge,F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    write_exexe('fichier(s) charge(s): '), write_l(File_l), !,
    extex(F,Name,Nmax,Vars,Flist,Plist,Fdep,File_l).

% exemple d''une session

handle_ctnue(exemple,F,Name,Nmax,Vars,Flist,Plist,Fdep,File_l) :-
    exexe_system("more session"), !,
    extex(F,Name, Nmax, Vars, Flist, Plist, Fdep, File_l).

% fin
handle_ctnue(fin,_,_,_,_,_,_,_,_) :-  extex(fin,_).
handle_ctnue(f,_,_,_,_,_,_,_,_) :-  extex(fin,_).


%but courant
handle_ctnue(b_cour,F,Name,Nmax,Vars,Flist,Plist,Fdep,File_l) :-
      write_exexe('    but courant: '), write_ls_f([(F,Name)],Vars), !,
      extex(F,Name,Nmax,Vars,Flist,Plist,Fdep,File_l).


% but initial
handle_ctnue(b_init, F, Name, Nmax, Vars,Flist, Plist,Fdep,File_l) :-
    write_exexe('    but initial: '), write_ls_f([(Fdep, ['F'])],Vars), !,
    extex(F,Name,Nmax,Vars,Flist,Plist,Fdep,File_l).

% tout_es

handle_ctnue(tout_es, F,Name,Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    write_ls_es(Plist,Vars), !,
    extex(F, Name, Nmax,Vars,Flist, Plist, Fdep, File_l).

% depl_es

handle_ctnue(depl_es, F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    write_exexe(' but prouve:'), write_ls_f([(Fdep, ['F'])], Vars),
    write_exexe(' predicat associe: '), es_pred(Fdep, ['F'], Pred), nl,
    %gtrace,
    write_exexe('     '), write_varnames(Pred, Vars), nl, nl,
    write_exexe(' clauses d''entree-sortie apres depliages: '), nl,
    simpl_by_unfolding(Plist, Plist1),
    update_vars(Plist1, Vars, Vars1),
    write_cls(Plist1, Vars1), !,
    extex(F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l).

% efface

handle_ctnue(efface,_,_,_,_,_,_,_,File_l) :-
    write_exexe('>>>> Interruption de la preuve. '), nl,
    write_exexe('>>>> g/guide pour obtenir de l''aide.'), nl, nl,
    handle_begin(File_l).

% recommence

handle_ctnue(recomm,_,_,_, Vars,_,_,Fdep, File_l) :-
    write_exexe('>>>> Nouvelle preuve, avec comme but initial: '), nl,
    write_ls_f([(Fdep, ['F'])], Vars), variables(Fdep, Vars1), !,
    extex(Fdep, ['F'], 0, Vars1,[],[], Fdep, File_l).

% remove

handle_ctnue(remove,_,_, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    next([], Vars, Flist, Fnext, Name, Flistrmd), !,
    extex(Fnext, Name, Nmax, Vars, Flistrmd, Plist, Fdep, File_l).

%Gestion des regles d''execution etendue

% postulat

handle_ctnue(postul, A, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    assume_es(A, Name, Nmax, Vars, Plist, Plist1),
    next([(true,_)], Vars,Flist, Fnext, Name1, Flistrmd),!,
    extex(Fnext, Name1, Nmax, Vars, Flistrmd, Plist1, Fdep,File_l).

% induction par point fize

handle_ctnue(comput, (A => Q), Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    %trace,
    induction(A,Q,Hl, Tal, Thll),
    update_vars([Tal|Thll], Vars, Vars1),
    new_names(Nmax, Hl, Nmax1, Hlind),
    write_ls_form(Hlind, Vars1), !,
    ind_es((A => Q), Name, Vars1, Hlind, Tal, Thll, Plist, Plist1),
    next(Hlind, Vars1, Flist, Fnext, Name1, Flistrmd), !,
    extex(Fnext, Name1, Nmax1, Vars1, Flistrmd, Plist1, Fdep, File_l).

handle_ctnue(comput, (A => Q), Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    extex_fail((A => Q), Name, Nmax, Vars, Flist, Plist, Fdep, File_l), !.

% induction structurelle

handle_ctnue(struct, F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    choose_var(F, Var, Type, Vars),
    %gtrace, % FM
    structural(F, Var, Type, Gl, Concl, Hypll),
    update_vars([Concl|Hypll], Vars, Vars1),
    new_names(Nmax, Gl, Nmax1, Glind),
    write_ls_form(Glind, Vars1),
    ind_es(F, Name, Vars1, Glind, Concl, Hypll, Plist, Plist1),
    next(Glind, Vars1, Flist, Fnext, Name1, Flistrmd), !,
    extex(Fnext, Name1, Nmax1, Vars1, Flistrmd, Plist1, Fdep, File_l).

handle_ctnue(struct, F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    extex_fail(F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l), !.

% existentialisation: transforme des variables universelles en
% existentielles avant induction

handle_ctnue(exist, (A => B), Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    %trace, % FM
    mkexist(A,B, Vars, B, B1), !,
    write_ls_form([((A=>B1), Name)], Vars), !,
    extex((A=>B1), Name, Nmax, Vars, Flist, Plist, Fdep, File_l).

handle_ctnue(exist, C, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    mkexist(true, C, Vars, C, C1),
    write_ls_form([(C1, Name)], Vars), !,
    extex(C1, Name, Nmax, Vars, Flist, Plist, Fdep, File_l).
handle_ctnue(exist,F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    extex_fail(F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l), !.

mkexist(H,C,Vars, A, Cnew) :-
    vars_to_change(A, List, H, Vars),
    chg_exist(List, C, Cnew).

% universalisation: transforme des variables existentielles
% en universelles avant induction

handle_ctnue(univ, (A => B), Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    mkuniv(B, B1, Vars), !,
    write_ls_form([((A=>B1), Name)], Vars), !,
    extex((A=>B1), Name, Nmax, Vars, Flist, Plist, Fdep, File_l).

handle_ctnue(univ, C, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    mkuniv(C, C1, Vars),
    write_ls_form([(C1, Name)], Vars), !,
    extex(C1, Name, Nmax, Vars, Flist, Plist, Fdep, File_l).

handle_ctnue(univ, F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    extex_fail(F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l), !.

mkuniv(B,Bnew, Vars) :-  vars(B, ExistB,_),
    choose_vars(ExistB, Vars, EVars, 'universelles:'),
    chg_univ(EVars, B, Bnew).

% transfert d'un atome de l'hypothese vers la conclusion

handle_ctnue(transfert, (H=> C), Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    choose_atom(H,A, Hnew, Vars), mkexist(Hnew, (A&C), Vars, A, Cnew),
    mk_form(Hnew, Cnew, NewF), write_ls_form([(NewF, Name)], Vars), !,
    extex(NewF, Name, Nmax, Vars, Flist, Plist, Fdep, File_l).

handle_ctnue(transfert, F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    extex_fail(F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l), !.

vars_to_change(A, List, Hnew, Vars) :-
    vars_univ_corr(A, Hnew, VarsA),
    choose_vars(VarsA, Vars, List, 'existentielles :').

vars_univ_corr(A,H,VarsA) :-
    vars(A,_, UnivA),
    vars(H,_, UnivH),
    minus1(UnivA, UnivH, VarsA).

choose_vars([],_,[],_) :-  !.
choose_vars([X|V], Vars, List, Specif) :-
    write_exexe('>> indiquez les variables qui deviennent '), write_exexe(Specif), nl,
    propose_vars([X|V], Vars, List).

propose_vars([X|V], Vars, List) :-  write_varnames(X,Vars),
    write_exexe(' (o|n): '), read_exexe(Rep), propaux(Rep, X, List, List1),
    propose_vars(V, Vars, List1).

propose_vars([],_,[]) :-  nl.
propaux(o,X, [X|List], List).
propaux(n,_, List, List).

chg_exist(_,X,X) :-  varexist(X,_), !.
chg_exist(Vars, X,?(X)) :-  var(X), member_var(X, Vars,_), !.
chg_exist(_,X,X) :-  var(X), !.
chg_exist(Vars, F, Fex) :-  same_functor(F, Fex,_,N),
    chg_exist_arg(N, Vars, F, Fex).

chg_exist_arg(0,_,_,_).
chg_exist_arg(N, Vars, F, Fex) :-
    N >0, arg(N,F,ArgF),
    chg_exist(Vars, ArgF, ArgFex), arg(N, Fex, ArgFex),
    N1 is N-1,
    chg_exist_arg(N1, Vars, F, Fex).

chg_univ(_,X,X) :- var(X), !.
chg_univ(EVars, X,Y) :-  varexist(X,Y), member_var(X, EVars,_), !.
chg_univ(_,X,X) :-  varexist(X,_), !.
chg_univ(EVars, F, Fex) :-  same_functor(F, Fex,_,N),
    chg_univ_arg(N, EVars, F, Fex).

chg_univ_arg(0,_,_,_).
chg_univ_arg(N, EVars, F, Fex) :-
    N > 0, arg(N,F,ArgF),
    chg_univ(EVars, ArgF, ArgFex), arg(N, Fex, ArgFex),
    N1 is N-1,
    chg_univ_arg(N1, EVars, F, Fex).


mk_form(true, C,C) :-  !.
mk_form(H,C,C) :-  rm_true_f(H, true),!.
mk_form(H,C, (H1 => C)) :-  rm_true_f(H,H1).

rm_true_f((true&B), B) :- !.
rm_true_f((A&true), A) :-  !.
rm_true_f((A&B), (A1&B1)) :-!, rm_true_f(A,A1), rm_true_f(B,B1).
rm_true_f(A,A).

% deci = dci non-deterministe (etendue)

handle_ctnue(dci, F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    %trace,
    choose_atom_occ_dci(F,SubF, Occ, Vars), !,
    dci_fail(F,SubF, Occ, Vars, Name, Nmax, Flist, Plist, Fdep, File_l).

dci_fail(F, SubF, Occ, Vars, Name, Nmax, Flist, Plist, Fdep, File_l) :-
    dci_ext(F,SubF, Occ, Hl, Substl, Intvarsl), !,
    update_vars([Substl, Intvarsl], Vars, Vars1),
    new_names(Nmax, Hl, Nmax1, Hnl),
    write_ls_form(Hnl, Vars1),
    map_one_es(F, Name, Vars1, Hnl, Substl, Plist, Plist1),
    next(Hnl, Vars1, Flist, Fnext, Name2, Flistrmd), !,
    extex(Fnext, Name2, Nmax1, Vars1, Flistrmd, Plist1, Fdep, File_l).
dci_fail(F,_,_, Vars, Name, Nmax, Flist, Plist, Fdep, File_l) :-
    extex_fail(F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l), !.

handle_ctnue(dci,F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    handle_ctnue(dci_atom, F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l).



% dci_atom permet de choisir un atome non-instancie
handle_ctnue(dci_atom, F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    choose_atom_occ(F,pos, SubF, Occ, Vars),
    dci_ext(F, SubF, Occ, Hl, Substl, Intvarsl),
    update_vars([Substl, Intvarsl], Vars, Vars1),
    new_names(Nmax, Hl, Nmax1, Hnl),
    write_ls_form(Hnl, Vars1),
    map_one_es(F, Name, Vars1, Hnl, Substl, Plist, Plist1),
    next(Hnl, Vars1, Flist, Fnext, Name2, Flistrmd), !,
    extex(Fnext, Name2, Nmax1, Vars1, Flistrmd, Plist1, Fdep, File_l).
handle_ctnue(dci_atom, F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    extex_fail(F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l), !.

% dcir choisit la clause et restreint le choix de l''atome

handle_ctnue(dcir, F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    choose_atom_occ_dci(F,SubF,Occ, Vars), !,
    dci_cl_fail(F, SubF, Occ, Vars, Name, Nmax, Flist, Plist, Fdep, File_l).

dci_cl_fail(F,SubF, Occ, Vars, Name, Nmax, Flist, Plist, Fdep, File_l) :-
    dci_clause(F,SubF, Occ,H,Subst, Intvars),
    update_vars([Subst|Intvars], Vars, Vars1),
    new_names(Nmax, [H], Nmax1, [(H, Name1)]),
    write_ls_form([(H,Name1)],Vars1),
    one_es(F, Name, Name1, Vars1, H, Subst, Plist, Plist1),
    next([(H,Name1)], Vars1, Flist, Fnext, Name2, Flistrmd), !,
    extex(Fnext, Name2, Nmax1, Vars1, Flistrmd, Plist1, Fdep, File_l).
dci_cl_fail(F,_,_, Vars, Name, Nmax, Flist, Plist, Fdep, File_l) :-
    extex_fail(F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l), !.

handle_ctnue(dcir, F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    handle_ctnue(dcir_atom, F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l).

% dcir_atom choisit la clause mais permet un atome non-instancić

handle_ctnue(dcir_atom, F,Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    choose_atom_occ(F,pos,SubF,Occ, Vars),
    dci_clause(F, SubF, Occ, H, Subst, Intvars),
    update_vars([Subst|Intvars], Vars, Vars1),
    new_names(Nmax, [H], Nmax1, [(H,Name1)]),
    write_ls_form([(H,Name1)], Vars1),
    one_es(F, Name, Name1, Vars1, H, Subst, Plist, Plist1),
    next([(H,Name1)], Vars1, Flist, Fnext, Name2, Flistrmd), !,
    extex(Fnext, Name2, Nmax1, Vars1, Flistrmd, Plist1, Fdep, File_l).
handle_ctnue(dcir_atom, F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    extex_fail(F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l), !.

% simplification (un atome pour un atome)

handle_ctnue(simpl, F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    %gtrace,
    choose_pair_atom_occ(F, Inhibit, SubF2, Occ2, SubF1, Occ1, Vars),
    simpl_inhib(Inhibit, F, SubF1, Occ1, SubF2, Occ2, Name, Nmax, Vars, Flist, Plist, Fdep, File_l).

simpl_inhib(no, F,_,_,_,_, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-  !,
    extex(F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l).

simpl_inhib(_,F,SubF1, Occ1, SubF2, Occ2, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    %trace,
    simplification(F, SubF1, Occ1, SubF2,Occ2, H, Subst),
    update_vars(Subst, Vars, Vars1),
    new_names(Nmax, [H], Nmax1, [(H,Name1)]),
    write_ls_form([(H, Name1)], Vars1),
    one_es(F, Name, Name1, Vars1, H, Subst, Plist, Plist1),
    next([(H, Name1)], Vars1, Flist, Fnext, Name2, Flistrmd), !,
    extex(Fnext, Name2, Nmax1, Vars1, Flistrmd, Plist1, Fdep, File_l).

handle_ctnue(simpl, F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    extex_fail(F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l), !.

% simplification gardant l''atome negatif

handle_ctnue(simpl_g, F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    choose_pair_atom_occ(F,Inhibit, SubF2, Occ2, SubF1,Occ1,Vars),
    s_kp_inhib(Inhibit, F, SubF1, Occ1, SubF2, Occ2, Name, Nmax, Vars, Flist,
               Plist, Fdep, File_l).

s_kp_inhib(no, F,_,_,_,_, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-  !,
    extex(F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l).

s_kp_inhib(_,F, SubF1, Occ1, SubF2, Occ2, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    simplification_kp(F, SubF1, Occ1, SubF2,Occ2,H,Subst),
    update_vars(Subst, Vars, Vars1),
    new_names(Nmax, [H], Nmax1, [(H,Name1)]),
    write_ls_form([(H, Name1)], Vars1),
    one_es(F, Name, Name1, Vars1, H, Subst, Plist, Plist1),
    next([(H,Name1)], Vars1, Flist, Fnext, Name2, Flistrmd), !,
    extex(Fnext, Name2, Nmax1, Vars1, Flistrmd, Plist1, Fdep, File_l). %flag1

handle_ctnue(simpl_g, F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    extex_fail(F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l), !.

% rattrapages

extex_fail(F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    write_exexe('>>>> Cette regle echoue ... '), nl, nl,
    extex(F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l).

% problemes divers

handle_ctnue(_, F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l) :-
    write_exexe('>>>> Probleme ! (g/guide pour obtenir de l''aide)'), nl, nl,
    extex(F, Name, Nmax, Vars, Flist, Plist, Fdep, File_l).

%Menus, affichage

%% next(Flist, Vars, Flist1, F, Name, Flistrmd)

next(Flist, Vars, Flist1, F, Name, Flistrmd) :-
    app_true(Flist, Flist1, Flist2),
    %gtrace,
    menu_form(Flist2, Vars, 'quel but ?', 'nouveau but courant: ',
              (F,Name), Flistrmd).

% supprime "true" de la pile

app_true(Fl, Flist, NewFlist) :-  success(Fl, 0, Flist, 0, NewFl),
    app_pgm(NewFl, Flist, NewFlist).
success([(true,_)|Fl], N, Flist, Flag, NewFl) :-  !,
    N1 is N + 1,
    success(Fl, N1, Flist, Flag, NewFl).
success([F|Fl],N,Flist,_, [F|NewFl]) :-  success(Fl, N, Flist,1,NewFl).
success([],_,[],0,[]) :-  !, write_exexe('>>>>>>>>>> succes total!'), nl, nl.
success([],0,_,_,[]) :- !.
success(0,1,_,_,[]) :-  !, write_exexe('>>>> 1'),
    write_exexe(' succes partiel'), nl, nl.
success([],N,_,_,[]) :-  write_exexe1('>>>>'), write_exexe1(N),
    write_exexe1(' succes partiels '), nl, nl.

%% menu(List, Vars, Text1, Text2, Choice, Rest): choisit dans List un
%% element Choice,
%%laissant Rest; affiche Text1 (general) ou Text2 (une seule possibilite)

menu_cl(no,[],_,nil,[]) :-  !.
menu_cl(ok,[F],_, F,[]) :-  !.
menu_cl(Inhibit,List, Vars, Choice, Rest) :-
    write_exexe('>> quelle clause (a. pour abandon)?'), nl,
    display_cl(List, Vars, 1), read(Rep), nl,
    inhibit(Rep, Inhibit, List, Rest,Choice).

menu_occ([],_,_,_,([],nil),[]) :-  !.
menu_occ([(Occ,F)],Vars,_, Text2, (Occ,F),[]) :-  !, write_exexe(' '),
    write_exexe(Text2), nl, write_varnames(F,Vars), nl, nl.
menu_occ(Occlist, Vars, Text,_, Choice, Rest) :-
    write_exexe('>>'), write_exexe(Text), nl, display_occ(Occlist, Vars, 1),
    read_exexe(N), nl, remove(Occlist, N, Rest, Choice).

menu_pair_occ(no,[],_,_,_,([],nil,[],nil),[]) :-  !.
menu_pair_occ(Inhibit, [(Oc1, A1, Oc2, A2)], Vars,_, Text2, (Oc1, A1, Oc2, A2),[]) :-  !,
    write_exexe('   '), write_exexe(Text2), nl, write_varnames((A1,A2), Vars), nl, nl,
    write_exexe1('>> confirmation (o/n):'), read_exexe(Rep), nl, inhibit(Rep, Inhibit).
menu_pair_occ(Inhibit, Pairlist, Vars, Text,_, Choice, Rest) :-
    write_exexe('>>'), write_exexe(Text), nl, display_pair(Pairlist, Vars,1),
    read_exexe(N), nl, inhibit(N, Inhibit, Pairlist, Rest, Choice).

inhibit(n,no) :-  !.
inhibit(_,ok).

inhibit(a,no,_,_,_) :-  !.
inhibit(N,ok, Choicelist, Rest, Choice) :-  remove(Choicelist, N, Rest, Choice).

menu_form([],_,_,_,(nil,_),[]) :-  !.
menu_form([(F,Name)], Vars,_, Text2, (F, Name),[] ) :-  !, write_exexe('  '),
    write_exexe(Text2), mk_form_name(Name,_), write_ls_f([(F,Name)], Vars).


menu_form(List, Vars, Text,_, Choice, Rest) :-  write_exexe1('>>'), write_exexe(Text), nl,
    display_form(List, Vars),
    %get_code(_), % Correction % modification
    get_l(Name), nl,
    remove_form(List, Name, Rest, Choice).


get_l(X) :-  % modification
    open('getg.txt', read, Stream),
    get_code(Stream, F),
    get_code(Stream,N),
    
    close(Stream), get_l1(F,N,X).

get_l1(F,46,X) :-  !,
    get_code(10),
    my_name(X, [F]).


% Prédicat pour vérifier la séquence de caractères dans le flux
get_l1(F, N, X) :- %modification
    open('getg1.txt', read, Stream), % Réouvrir le flux pour lire
    % Lit les caractères et vérifie le code ASCII
    get_code(Stream, Code1), % Lit le prochain caractère
    (Code1 = 46 -> % Vérifie si le caractère est un point (46)
        % Lit le prochain caractère
        get_code(Stream, Code2),
        (Code2 = 10 -> % Vérifie si le caractère est une nouvelle ligne (10)
            my_name(X, [F, N]) % Appelle my_name/2
        ;   % Si le deuxième caractère n'est pas une nouvelle ligne, affiche une erreur
            write_exexe('Erreur : caractère attendu \\n non trouvé.'),
            fail % Échec de la lecture
        )
    ;   % Si le premier caractère n'est pas un point, affiche une erreur
        write_exexe('Erreur : caractère attendu . non trouvé.'),
        fail % Échec de la lecture
    ),
    close(Stream),% Ferme le flux
    write_exexe(X). 



    

menu_var([],_,nil,[]) :-  !.
menu_var([F], Vars, F,[]) :-  !,
    write_exexe(' variable d''induction structurelle choisie: '), nl,
    write_varnames(F, Vars), nl, nl.

menu_var(List, Vars, Choice, Rest) :-
    write_exexe('>> choisissez une variable d''induction structurelle: '), nl,
    display_var(List,Vars,1), read_exexe(N), nl, remove(List, N, Rest, Choice).

display_cl([],_,_) :-  nl.
display_cl([C|Clist], Vars,N) :-  write(N), put_code(32),  put_code(32),
    write_varnames(C, Vars), nl,
    N1 is N + 1,
    display_cl(Clist, Vars, N1).

display_occ([],_,_) :-  nl.
display_occ([(_,P)|Plist], Vars, N) :-  write(N),  put_code(32),  put_code(32),
    write_varnames(P,Vars), nl,
    N1 is N + 1,
    display_occ(Plist, Vars, N1).

display_pair([],_,_) :-  nl.
display_pair([(_,A1,_, A2)|Plist], Vars, N) :-  write(N),  put_code(32),  put_code(32),
    write_varnames(A1, Vars), write_exexe(','), write_varnames(A2, Vars), nl,
    N1 is N + 1, display_pair(Plist, Vars, N1).

display_var([],_,_) :-  nl.
display_var([(V: Type)|Vlist], Vars, N) :-  write(N), put_code(32), put_code(32),
    write_varnames(V,Vars), write_exexe1(' : '), write_varnames(Type, Vars),nl,
    N1 is N + 1, display_var(Vlist, Vars, N1).

display_form([],_) :-  nl.
display_form([(F, Name)|Flist], Vars) :-  length(Name, 1), !,
    mk_form_name(Name, Str), write_exexe1(Str), write_exexe1(': '),
    write_varnames(F, Vars), nl, display_form(Flist, Vars).
display_form([(F, Name)|Flist], Vars) :-  length(Name, 2), !,
    mk_form_name(Name, Str), write_exexe1(Str), write_exexe1(': '),
    write_varnames(F, Vars), nl, display_form(Flist, Vars).

remove([F|Flist], 1, Flist, F) :-  !.
remove([F|Flist], N, [F|Flist1], Choice) :-
    N1 is N - 1,
    remove(Flist, N1, Flist1, Choice).

%remove_var([(V : T)|Vlist], X, Vlist, (V : T)) :-  V == X, !.
%remove_var([V|Vlist], X, [V|Vlist1], Choice) :-
%    remove_var(Vlist, X, Vlist1, Choice).

remove_form([(F, Name) | Flist], String, Flist, (F, Name)) :-
    mk_form_name(Name, Str), Str == String, !.
remove_form([F|Flist], String, [F|Flist1], Choice) :-
    remove_form(Flist, String, Flist1, Choice).

%% choix

% choose_clause(Inhibit, Head, Head_inst, Body)

choose_clause(Inhibit, Head, Head_inst, Body) :-  variables(Head, Vars),
    setof((Head:-B), Vars^clause(Head, B), Cl), variables(Cl, Vars1),
    menu_cl(Inhibit, Cl, Vars1, (Head_inst:-Body),_).

/*sous-but positif/negatif: SH est un sous-but positif/negatif de H a
% Occ*/

positive(SH,H,[]) :-  SH == H.
positive(SH, (_ => H2), Occ) :-  !, Occ =[2|Occ1], get_occ(H2,SH, Occ1).
positive(SH,H,Occ) :-  get_occ(H, SH, Occ).
negative(SH, (H1=> _), [1|Occ]) :-  !, get_occ(H1, SH,Occ).

% choose var (variable d''induction)

choose_var(F, Var, Type, Vars) :-
    setof((Var: Type), var_type(F, Var, Type), VTl),
    menu_var(VTl, Vars, (Var: Type),_).

var_type(F, Var, Type) :-
    functor(F, Pred,_),
    type_rec(Pred, N, Type1, _),
    arg(N,F,Arg),
    ctns_univ(Arg, Var,_),
    functor(Type, Type1,1),
    arg(1, Type, Var).

var_type((A =>_), Var, Type) :-  var_type(A,Var, Type).
var_type((_ => B), Var, Type) :-  var_type(B, Var, Type).
var_type((A&_), Var, Type) :-  var_type(A, Var, Type).
var_type((_&B), Var, Type) :-  var_type(B, Var, Type).
var_type(-(A), Var, Type) :-  var_type(A, Var, Type).

ctns_univ(X,X,N) :-  var(X), !, N = 1.
ctns_univ(X,_,_) :-  ground(X), !, fail.  % FM
ctns_univ(X,_,_) :-  varexist(X,_), !, fail.
ctns_univ(X,Y,_) :-  functor(X,_,M), arg(M,X,Y), var(Y), !.   % FM, for [X|Xs] ?

/*
% choose subf_occ(F,Sign, SubF, Occ, Vars)

choose_subf_occ(F,pos, SubF, Occ, Vars) :-
    setof((Oc,G), (get_occ(F,G,Oc), positive(G,F,Oc)),Gl),
    menu_occ(Gl, Vars, 'a quel sous-but positif? ',
             'sous-but positif choisi: ', (Occ, SubF),_).

choose_subf_occ(F,neg, SubF, Occ, Vars) :-
    setof((Oc,G), (get_occ(F,G,Oc), negative(G,F,Oc)),Gl),
    menu_occ(Gl, Vars, 'a quel sous-but negatif?',
             'sous-but negatif choisi: ', (Occ,SubF),_).
*/

% choose atom (transfert)

choose_atom(Hyp, Atom, Hypnew, Vars) :-
    setof((Hn, A), atom_hyp(Hyp, A,Hn), AHl),
    menu_occ(AHl, Vars, 'a quel atome negatif?',
             'atome negatif choisi: ', (Hypnew, Atom),_).

atom_hyp(A,A,true) :-  is_atom(A).
atom_hyp(H,A,Hn) :-  atom_hyp1(H,A,Hn).
atom_hyp1((A&H),A,H) :- is_atom(A).
atom_hyp1((H&A),A,H) :-  is_atom(A).
atom_hyp1((H1&H2), A, (H1n&H2)) :- atom_hyp1(H1, A, H1n).
atom_hyp1((H1&H2),A, (H1&H2n)) :- atom_hyp1(H2,A,H2n).

% choose_atom_occ(F,Sign, Atom, Occ, Vars)

choose_atom_occ(F,pos, Atom, Occ, Vars) :-
    setof((Oc,G), atom_occ(F,G,Oc,pos),Gl),
    menu_occ(Gl,Vars, 'a quel atome positif ? ',
             'atome positif choisi : ', (Occ,Atom),_).

choose_atom_occ(F,neg, Atom,Occ, Vars) :-
    setof((Oc,G), atom_occ(F,G,Oc,neg),Gl),
    menu_occ(Gl, Vars, 'a quel atome negatif ? ',
             'atome negatif choisi : ', (Occ, Atom),_).

atom_occ(F,G,Oc,pos) :-  get_occ(F,G,Oc), is_atom(G), positive(G,F,Oc).
atom_occ(F,G,Oc,neg) :-  get_occ(F,G,Oc), is_atom(G), negative(G,F,Oc).

%choose pair_atom_occ

choose_pair_atom_occ(F, Inhibit, Atomneg, Occneg, Atompos, Occpos, Vars) :-
    setof( (Oc1,A1,Oc2,A2), pair_atom_occ(F,Oc1,A1,Oc2,A2), P1),
    menu_pair_occ(Inhibit, P1, Vars, 'a quelle paire d atomes (a. pour abandon) ?',
        'paire d atomes choisie : ',
        (Occneg, Atomneg, Occpos, Atompos),_).

pair_atom_occ(F,Oc1,A1,Oc2,A2) :-  atom_occ(F,A1,Oc1,neg),
    atom_occ(F,A2,Oc2,pos), unification(A1,A2,[],_). %Suniv = []

%choose_atom_occ_determ

choose_atom_occ_determ(F,Atom,Occ, Vars) :-
    setof((Oc,G),atom_occ(F,G,Oc,pos),Gl),
    menu_occ(Gl, Vars, 'a quel atome positif ? ', 'atome positif choisi : ',(Occ,Atom),_).

% choose_atom_occ_dci

choose_atom_occ_dci(F,Atom, Occ, Vars) :-
    setof((Oc,G),dci_atom_occ(F,Oc,G),Gl),
    menu_occ(Gl, Vars, 'a quel atome positif ? ', 'atome positif choisi : ',(Occ,Atom),_).

dci_atom_occ(F,Oc,G) :-
    get_occ(F,G,Oc),
    is_atom(G),
    positive(G,F,Oc),
    instantiated_type(G).

instantiated_type(G) :-
    functor(G,F,N),
    inst_arg(F,N,G).

inst_arg(F,N,G) :-
    N > 0,
    arg(N,G,ArgG), type_rec(F,N,_, true),
    inst(ArgG), !.

inst_arg(F,N,G) :-
    N > 0 , N1 is N - 1,
    inst_arg(F,N1,G).

inst(X) :-  var(X), !, fail.
inst(X) :-  varexist(X,_), !, fail.
inst(_).

%% affichage des buts

/*write_ls_form([(F,N)],Vars) :-  !, write(' nouveau but:'),
   %gtrace,
   write_ls_f([(F,N)], Vars).*/

write_ls_form([FN|Lind], Vars) :-
    write_exexe(' nouveaux buts:'),
    write_ls_f([FN|Lind], Vars).

write_ls_f([],_) :-  nl.
write_ls_f([(true,_)|Fs], Vars) :-
    !, nl, write_exexe(' vrai'), nl,
    write_ls_f(Fs, Vars).
write_ls_f([(F, Name)|Fs], Vars) :-
    length(Name, 1), !,
    mk_form_name(Name, String), nl,
    write_exexe1('('), write_exexe1(String),
    write_exexe1(')'),
    write_varnames(F, Vars), nl,
    write_ls_f(Fs, Vars).
write_ls_f([(F,Name)|Fs], Vars) :-
    length(Name, 2), !,
    mk_form_name(Name, String),
    nl, write_exexe1('('), write_exexe1(String),
    write_exexe1(')'), nl, write_varnames(F, Vars), nl,
    write_ls_f(Fs,Vars).

write_cls([F|Fs], Vars) :-
    write_varnames(F, Vars), %write_exexe1('.'), 
    nl,
    write_cls(Fs, Vars).
write_cls([],_) :-  nl.

write_varnames(F, Vars) :-  lettervars(F, Vars, F1), write_exexe(F1).

lettervars(X, Vars, Y) :-  var(X), !, names(Names), get_name(X, Vars, Names, Y).
lettervars(F, Vars, G) :-  same_functor(F,G,_,N), lettervars_arg(N,F,Vars,G).

lettervars_arg(0,_,_,_).
lettervars_arg(N,F,Vars,G) :-   N > 0, arg(N,F,ArgF),
    lettervars(ArgF, Vars,ArgG), arg(N,G,ArgG),N1 is N - 1,
    lettervars_arg(N1,F,Vars,G).

get_name(X, [V|_],[Y|_],Y) :-  X == V, !.
get_name(X, [_|Vars], [_|Names], Y) :-  get_name(X,Vars, Names,Y).

names(['X','Y', 'Z', 'T', 'U', 'V', 'W', 'A', 'B','C' ,'D','E','F',
       'X1','Y1', 'Z1', 'T1', 'U1', 'V1', 'W1', 'A1', 'B1','C1' ,'D1','E1','F1',
       'X2','Y2', 'Z2', 'T2', 'U2', 'V2', 'W2', 'A2', 'B2','C2' ,'D2','E2','F2',
       'X3','Y3', 'Z3', 'T3', 'U3', 'V3', 'W3', 'A3', 'B3', 'C3', 'D3', 'E3', 'F'
       ]).
update_vars([], Vars, Vars).
update_vars([NewF|Fl], Vars, Varsnew) :-  variables(NewF, Vars1),
    merge(Vars1, Vars, Vars2), update_vars(Fl, Vars2, Varsnew).

%%Sortie de la boucle

extex(nil,_,_, Vars,[], Plist, Fdep, File_l) :-  !,
    write_exexe('>>>> Voulez-vous voir toutes les clauses d''entree-sortie (tout_es)'), nl,
    write_exexe('>>>> ou seulement leur simplification par depliages (depl_es) ? '), nl, nl,
    give_rule(Rule), handle_end(Rule, Vars, Plist, Fdep, File_l).

%% handle_end

% tout_es

handle_end(tout_es, Vars, Plist, Fdep, File_l) :-  write_exexe(' but prouve:'),
    write_ls_f([(Fdep, ['F'])], Vars), write_exexe('   '), write_exexe(' predicat associe: '),%modification
    es_pred(Fdep, ['F'], Pred), nl,
    write_varnames(Pred, Vars), nl, nl, write_exexe('   '),write_exexe('clauses collectees: '), nl,%modification
    write_cls(Plist, Vars),
    give_rule(Rule),
    handle_end(Rule, Vars, Plist, Fdep,File_l).


% deples

handle_end(depl_es, Vars, Plist, Fdep, File_l) :-  write_exexe('but prouve:'),
    write_ls_f([(Fdep, ['F'])], Vars),  write_exexe('     '), write_exexe('predicat associe: '),%modification
    es_pred(Fdep, ['F'], Pred), nl,
    write_varnames(Pred, Vars), nl, nl,write_exexe('     '),%modification
    write_exexe(' clauses d''entree-sortie apres depliages: '), nl,
    simpl_by_unfolding(Plist, Plist1),
    update_vars(Plist1, Vars, Vars1),
    write_cls(Plist1, Vars1),
    try_truncate(Fdep, Plist1,_, Vars1),nl,nl,
    write_exexe('>>>> Fin de la synthese.'),  nl,
    write_exexe('>>>> Vous pouvez commencer une autre preuve (g. pour etre aide).'),
    nl, nl,
    handle_begin(File_l).

% fin

handle_end(f,_,_,_,_) :-  extex(fin,_).

% efface

handle_end(efface,_,_,_, File_l) :-  handle_ctnue(efface,_,_,_,_,_,_,_, File_l).

% recommencement

handle_end(recomm,_,_, Fdep, File_l) :-
    handle_ctnue(again,_,_,_,_,_,_, Fdep, File_l).

% erreur

handle_end(_,V,Pl,F,File_l) :-
    write_exexe('>>>> Utiliser: tout_es, depl_es, f/fin.'), nl,
    give_rule(Rule), handle_end(Rule, V, Pl, F, File_l).

%% try_truncate

try_truncate(Fdep, P,Ptr, Vars) :-  quest(o), !, give_truncs(Fdep, Argl), nl,
    truncate(P, [('es_F', Argl)], Ptr),
    write_exexe(' clauses d''entree-sortie tronquees: '), nl, write_cls(Ptr, Vars).
try_truncate(_,_,_,_) :-  nl.

quest(X) :-  write_exexe('>>>> Voulez-vous tronquer es_F (o/n) ? '), read_exexe(X), nl.

%% get trunc_list

give_truncs(Fdep, Argl) :-
    multiple_args(Fdep, MArgs, SArgs),
    MArgs \== [],!,
    propose(MArgs, SArgs, Argl).

give_truncs(Fdep, Argl) :-  multiple_args(Fdep,[], SArgs), prop(SArgs, Argl,_).
multiple_args(F, MArgs, SArgs) :-  vars_mult(F, Vars),
    sep_mult(1, Vars, MArgs, SArgs).

vars_mult(X, [(X,1)]) :- var(X), !.
vars_mult(B, Vars) :-  functor(B,_,N), vars_mult_arg(N,B, Vars).
vars_mult_arg(0,_,[]).
vars_mult_arg(N,B,Vars) :-  N > 0, arg(N,B,ArgB), vars_mult(ArgB, Vars1),
    N1 is N - 1, vars_mult_arg(N1,B, Vars2), merge_m(Vars1, Vars2, Vars).

merge_m([],B,B).
merge_m([(X,N)|A],B,C) :-  insert_m(X,N,B,B1), merge_m(A,B1,C).

insert_m(X,N,[(Y,M)|A], [(X,P)|A]) :-  X == Y, !, P is N + M.
insert_m(X,N,[(Y,M)|A], [(Y,M)|B]) :-  insert_m(X,N,A,B).
insert_m(X,N,[],[(X,N)]).

sep_mult(N, [(_,M)|Vars], MArgs, [N|SArgs]) :-  M == 1, !, N1 is N + 1,
    sep_mult(N1, Vars, MArgs, SArgs).
sep_mult(N, [(_,_)|Vars], [N|MArgs], SArgs) :-  N1 is N + 1,
    sep_mult(N1, Vars, MArgs, SArgs).
sep_mult(_,[],[],[]).

propose(MArgs, SArgs, Argl) :-  prop(MArgs, Argl1, Rep),
    propctnue(Rep, SArgs, Argl1, Argl).

propctnue(fin,_, Argl, Argl) :- !.
propctnue(ctnue, SArgs, Argl1, Argl) :- prop(SArgs, Argl2,_), !,
    app_pgm(Argl1,Argl2, Argl).

prop(Args, Argl, Rep) :-  prop1(Args,[],_,Argl,Rep).
prop1([], Argl, Rep, Argl, Rep).
prop1([N|Args],Atemp,_, Argl, Rep) :-  quest(N, R, Ctnue),
    prop1ctnue(N, Args, Atemp, Ctnue, Argl, Rep, R).
prop1ctnue(N, Args, Atemp, ctnue, Argl, Rep,o) :-
    prop1(Args, [N|Atemp], ctnue, Argl, Rep).
prop1ctnue(_,Args, Atemp, ctnue, Argl, Rep,n) :-
    prop1(Args, Atemp, ctnus, Argl, Rep).
prop1ctnue(N,_,Atemp,fin, [N|Atemp],fin,o).
prop1ctnue(_,_, Argl, fin, Argl, fin,n).

quest(N,R,C) :-
    write_exexe1('>>>> par rapport a l''argument numero '),
    write_exexe1(N), write_exexe1(' (o/n/ofin/nfin) ? '), read_exexe(Rep), dispatch(Rep, R,C).
dispatch(ofin,o,fin) :-  !.
dispatch(o,o,ctnue) :- !.
dispatch(nfin,n,fin) :-!.
dispatch(n,n,ctnue) :-!.

%A.4 Predicats d'entree-sortie
% noms des predicats d'entree-sortie
% next_name (Nmax, Name)

next_name(Nmax, Name) :-  form_names(Names), next_n(Nmax, Names, Name).
form_names(['F','G','H','I','J','K','L','M','N','P','Q','R','S','T',
            'U', 'V','W', 'X','Y', 'Z', 'A','B', 'C','D','E','O']).
next_n(0,[Name|_], [Name]) :- !.
next_n(s(N), [_ |Names], Name) :-  next_n(N, Names, Name).

%% new_names (Nmax, Hl, Nmax, Hnl): manipulation d'indices (induction,
%% dci)

new_names(Nmax, Hl, Nmax, Hnl) :-  leng_true(Hl,0), !,
    next_name(Nmax, Name), mk_n_l(Hl, Name, Hnl).
new_names(Nmax, Hl,s(Nmax), Hnl) :-  leng_true(Hl,1), !,
    next_name(s(Nmax), Name), mk_n_l(Hl, Name, Hnl).
new_names(Nmax, Hl,s(Nmax), Hnl) :-  next_name(s(Nmax), Name),
    mk_n_l(Hl, Name, 1, Hnl).

leng_true([],0).
leng_true([true|L],N) :-  !, leng_true(L,N).
leng_true([_|L],N) :-  leng_true(L,N1), N is N1 + 1.

mk_n_l([true|Hl], Name, [(true,[])|Hnl]) :-  !, mk_n_l(Hl, Name, Hnl).
mk_n_l([H|Hl], Name, [(H,Name)|Hnl]) :-  mk_n_l(Hl, Name, Hnl).
mk_n_l([],_,[]).
mk_n_l([true|Hl], Name, N, [(true,[] )|Hnl]) :-  !, mk_n_l(Hl, Name, N, Hnl).
mk_n_l([H|Hl], Name, N, [(H, Name1)|Hnl]) :-  app_pgm(Name, [N], Name1),
    N1 is N + 1, mk_n_l(Hl, Name, N1, Hnl).
mk_n_l([],_,_,[]).

% predicat d'entree-sortie associe a une formule: es_pred(F, Name, Pred)

es_pred(F,Name, Pred) :-
    variables(F,Vars),
    mk_pred_name(Name, Str),
    mk_pred(Str, Vars, Pred).
mk_pred_name(Name, Str) :-
    char_to_num(['e', 's','_'|Name],L),
    my_name(Str,L).
  
mk_form_name(Name, Str) :-
    char_to_num(Name, L),
    my_name(Str,L).
    

char_to_num([],[]).
char_to_num([A|X], [B|Y]) :-
    my_name(A, [B]),
    char_to_num(X,Y).

mk_pred(Str,Vars, Pred) :-
    length(Vars,N),
    functor(Pred, Str,N),
    mk_args(1, Pred, Vars).

mk_args(_,_,[]).
mk_args(N, Pred, [X|Vars]) :-
    arg(N, Pred, X),
    N1 is N + 1,
    mk_args(N1, Pred, Vars).

%% manipulation des predicats d'entree-sortie

% postulat

assume_es(A, Name,_, Vars, Plist, Plist1) :-  es_pred(A, Name, Pred),
    rm_exist(A, A1), imp_and_to_cl(A1, Body),
    write_ls_cl([(Pred:- Body)], Vars),
    app_pgm(Plist, [(Pred:- Body)], Plist1).
imp_and_to_cl((A & B), (A1,B1)) :-  !, imp_and_to_cl(A,A1),
    imp_and_to_cl(B,B1).
imp_and_to_cl((A => B), (A1,B1)) :-  !, imp_and_to_cl(A,A1),
    imp_and_to_cl(B,B1).
imp_and_to_cl(A,A).

rm_exist(X,Y) :-  varexist(X,Y), !.
rm_exist(X,X) :-  var(X), !.
rm_exist(X,Y) :-  same_functor(X,Y,_,N), rm_exist_arg(N,X,Y).

rm_exist_arg(0,_,_).
rm_exist_arg(N,X,Y) :-  N > 0, arg(N,X, ArgX), rm_exist(ArgX, ArgY),
    arg(N,Y,ArgY), N1 is N - 1, rm_exist_arg(N1,X,Y).

% induction, dci

one_es(F, Name,_, Vars, true, Subst, Plist, Plist1) :-  !, es_pred(F, Name, Pred),
    apply(Subst, Pred, Head), write_ls_cl([Head], Vars),
    app_pgm(Plist, [Head], Plist1).

one_es(F, Name, Name1, Vars, H, Subst, Plist, Plist1) :-  es_pred(F, Name, Pred),
    apply(Subst, Pred, Head), es_pred(H, Name1, Body),
    write_ls_cl([(Head:- Body)], Vars),
    app_pgm(Plist, [(Head:- Body)], Plist1).

map_one_es(F, Name, Vars, Hnl, Substl, Plist, Plist1) :-
    map_es_rec(F, Name, Hnl, Substl, Clauses),
    app_pgm(Plist, Clauses, Plist1), write_ls_cl(Clauses, Vars).

map_es_rec(F, Name, [(H, Name1)|Hnl], [Subst|Substl], [Clause|Clauses]) :-
    es1(F, Name, Name1, H, Subst, Clause),
    map_es_rec(F,Name,Hnl,Substl,Clauses).
map_es_rec(_,_,[],[],[]).

es1(F, Name,_, true, Subst, Head) :-
    !,
    es_pred(F, Name, Pred),
   apply(Subst, Pred, Head).
es1(F, Name, Name1, H, Subst, (Head:- Body)) :-
    es_pred(F, Name, Pred),
    apply(Subst, Pred, Head),
    es_pred(H, Name1, Body).

ind_es(F, Name, Vars, H_l, Tau_l, Thetal_l, Plist, Plist1) :-
    ind_es_rec(F, Name, H_l, Tau_l, Thetal_l, Clauses),
    app_pgm(Plist, Clauses, Plist1),
    write_ls_cl(Clauses, Vars).

ind_es_rec(_,_,[],[],[],[]).
ind_es_rec(F, Name, [Hind|Hl], [Ta|Tal], [Thl|Thll], [Clause|Clauses]) :-
    one_ind_es(F, Name, Hind, Ta, Thl, Clause),
    ind_es_rec(F, Name, Hl, Tal, Thll, Clauses).

one_ind_es(F, Name, (H, NameH), Ta, Thl, (Head:- Body)) :-
    es_pred(F, Name, Pred),
    apply(Ta, Pred, Head),
    map_apply(Thl, Pred, Body_l1),
    es_pred(H, NameH, PredH),
    app_pgm(Body_l1, [PredH], Body_l),
    %list_to_set(Body_l,Body).
    exexe_list_to_set(Body_l,Body_s),
    list_roundlist(Body_s,Body). % FM

    list_roundlist([A],A) :- !.
    list_roundlist([A,B|Cs],(A,RL)) :- list_roundlist([B|Cs],RL).


% affichage des clauses d'entree-sortie

write_ls_cl([C], Vars) :-
    !,
    write_exexe(' clause associee:'),
    nl,
    write_cls([C], Vars).
write_ls_cl(L,Vars) :-
    write_exexe(' clauses associees:'),
    nl,
    write_cls(L,Vars).

write_ls_es([C], Vars) :-
    !,
    write_exexe(' clause:'),
    nl,
    write_cls([C], Vars).
write_ls_es(L,Vars) :-
    write_exexe('  clauses:'),
    nl,
    write_cls(L,Vars).

%% simplification par depliages

% unfold_pred(Clause, BridgeClause, UnfClause)

unfold((H:-Bs), (Haux:-Bsaux), Cl) :-
    !,
    set_to_list(Bs, Bl),
    set_to_list(Bsaux, Blaux),
    head_of(Haux, F),
    collect_atom(Bl, F, Befl, Bi, Aftl),
    unf_synt(H, Befl, Bi, Aftl, [Haux|Blaux], Cll),
    list_to_clause(Cll,Cl).
unfold((H:-Bs), Haux, Cl) :-
    !,
    set_to_list(Bs, Bl),
    head_of(Haux, F),
    collect_atom(Bl, F, Befl, Bi, Aftl),
    unf_synt(H,Befl, Bi, Aftl, [Haux], Cll),
    list_to_clause(Cll,Cl).


unf_synt(H, Befl, Bi, Aftl, [Haux|Aux], [Hi|Body]) :-
    unification(Bi, Haux, Suniv,[]),
    apply_l(Suniv, [H|Befl], [Hi|Befli]),
    apply_l(Suniv, Aftl, Aftli),
    app_pgm(Aux, Aftli, Bodyend),
    app_pgm(Befli, Bodyend, Body).

head_of((H :- _),F) :-
    !,
    functor(H,F,_).
head_of(H,F) :- functor(H,F,_).

collect_atom([B|Bl],F,[],B,Bl) :- head_of(B,F), !.
collect_atom([B|Bl],F,[B|Befl], Bi, Aftl) :- collect_atom(Bl,F,Befl, Bi, Aftl).


%simpl_by_unfolding(Plist, Prog)

simpl_by_unfolding(Plist, Prog) :-
    divide(Plist, Main, Others),
    partialev(Main, Others, [], Prog).

divide([P|Plist], [P|Main], Others) :-
    head_of(P, 'es_F'), !,
    divide(Plist, Main, Others).
divide([P|Plist], Main, [P|Others]) :-
    divide(Plist, Main, Others).
divide([],[],[]).


% partialev(Main, Others, [], Prog).
% ex : partialev([(a:-b),(a:-a,c)],[b,c],[],P).
% partialev([(es_F(0,Y,Y):-es_G(Y))],[(es_G(0):-es_H1),(es_G(s(Y)):-es_G(Y),es_H2(Y)),  es_H1, (es_H2(Y):-es_I(Y)),es_I(Y)],[],P).

partialev([],_, Prog,Prog).
partialev([Cl|Main], Others, P, Prog) :-
    one_partialev(Cl, Main, Others, P, Main1, Others1, P1),
    partialev(Main1, Others1, P1, Prog).
    %partialev(Main, Others, P1, Prog). %FM

one_partialev(Cl, Main, Others, P, Main1, Others1, P1) :-
    try_unf(Cl, Others, Clunf, Others_rmd), !,
    one_partialev(Clunf, Main, Others_rmd, P, Main1, Others1, P1).
one_partialev(Cl, Main, Others, P, Main1, Others_rmd, P1) :-
    collect_preds(Cl,Preds),
    collect_clauses(Others, Preds, Clauses, Others_rmd),
    app_pgm(Main, Clauses, Main1),
    app_pgm(P, [Cl],P1).

collect_preds((H:-B), Preds) :-
    !,
    head_of(H,F),
    set_to_list(B,Bl),
    collect_p_l(Bl, F, Preds).
collect_preds(_,[]).

collect_p_l([],_,[]).
collect_p_l([B|Bl], F, Preds) :-
    head_of(B,F), !, collect_p_l(Bl,F,Preds).
collect_p_l([B|Bl], F, [G|Preds]) :-
    head_of(B,G), collect_p_l(Bl,F,Preds).

collect_clauses([Cl|Others], Preds, [Cl|Clauses], Others_rmd) :-
    head_of(Cl,F), mb(F,Preds), !,
    collect_clauses(Others, Preds, Clauses, Others_rmd).
collect_clauses([Cl|Others], Preds, Clauses,[Cl|Others_rmd]) :-
    collect_clauses(Others, Preds, Clauses, Others_rmd).
collect_clauses([],_,[],[]).

% try_unf

try_unf(Cl,Others, Clunf, Rmd) :-
    multiple_preds(Others, Preds),
    try_unf(Cl, Others, Preds, Clunf, Rmd).
try_unf(Cl, [Cl1|Others], Preds, Clunf, [Cl1|Others_rmd]) :-
    head_of(Cl1,F),
    mb(F,Preds), !,
    try_unf(Cl, Others, Preds, Clunf, Others_rmd).% recursive pred
try_unf(Cl, [Cl1|Others],_, Clunf, Others) :-
    unfold(Cl, Cl1, Clunf), !.% unfolding here
try_unf(Cl, [Cl1|Others], Preds, Clunf, [Cl1|Others1]) :-
    try_unf(Cl, Others, Preds, Clunf,Others1). %irrelevant clause

multiple_preds([],[]).
multiple_preds([Cl|Cls], Preds) :-
    head_of(Cl,F),
    collect_clauses(Cls,[F],Coll,_),
    Coll = [], !,
    multiple_preds(Cls, Preds).
multiple_preds([Cl|Cls], [F|Preds]) :-
    head_of(Cl,F),
    collect_clauses(Cls, [F],_, Rmd),
    multiple_preds(Rmd, Preds).

%A.5 Troncature

%% truncate(Pgm, P_Arg_list, Pgmtr)

truncate(Pgm, [(F,A)|FAl], Pgmtr) :-
    up(A, Argl),
    trunc(Pgm,F,Argl,Auxl,Pgm1),
    app_pgm(FAl, Auxl, FAl1),
    truncate(Pgm1, FAl1, Pgmtr).
truncate(Pgm,[],Pgm).

up([],[]).
up([A|L], [M|Ml]) :-
    max(A,L,M,L1),
    up(L1, Ml).

max(A,[B|L], M,[A|L1]) :-
    B >= A,
    max(B,L,M,L1).
max(A, [B|L],M, [B|L1]) :-
    B < A,
    max(A,L,M,L1).
max(A,[],A,[]).

%% trunc(Pgm, F, Argl, Auxl, Pgmtr)

trunc([],_,_,[],[]).
trunc([(H:-Bs)|Cls], F, Argl, Auxl, [Cltr|Clstr]) :-
    functor(H,F,N), !,
    trunc_get(H,Bs, F, N, Argl, Auxl1, Cltr),
    trunc(Cls, F, Argl, Auxl2, Clstr),
    app_pgm(Auxl1, Auxl2, Auxl).
trunc([(H:-Bs)|Cls], F, Argl, Auxl, [Cltr|Clstr]) :-
    !, set_to_list(Bs, Bl),
    trunc_cl([H|Bl], F, Argl, Cl),
    list_to_clause(Cl,Cltr),
    trunc(Cls, F, Argl, Auxl, Clstr).
trunc([H|Cls],F,Argl, Auxl, [Htr|Clstr]) :-
    !, trunc_cl([H], F, Argl, [Htr]),
    trunc(Cls, F, Argl, Auxl, Clstr).

trunc_get(H,Bs, F, N, Argl, Auxl, Cltr) :-
    set_to_list(Bs, Bl),
    decomp(Bl, F, N, Recl, Nonrecl),
    varscorresp([H|Recl], Argl, Vars),
    argscorresp(Nonrecl, Vars, Auxl),
    trunc_cl([H|Bl],F,Argl,Cl),
    list_to_clause(Cl, Cltr).

varscorresp([],_,[]).
varscorresp([P|Pl], Argl, Vars) :-
    varscorr(P,Argl,[],VarsP),
    varscorresp(Pl, Argl, Varsl),
    merge(VarsP, Varsl, Vars).
varscorr(P,[N|Argl], Vars, VarsP) :-
    arg(N,P,V),
    merge([V], Vars, Vars1),
    varscorr(P,Argl, Vars1, VarsP).


varscorr(_,[], Vars, Vars).
argscorresp([],_,[]).
argscorresp([P|Pl], Vars, [(F, Argl)|PAl]) :-
    functor(P,F,N),
    is_es(F),
    argscorr(N,P, Vars, Argl),
    Argl \== [], !,
    argscorresp(Pl, Vars, PAl).
argscorresp([_|Pl], Vars, PAl) :-
    argscorresp(Pl, Vars, PAl).

argscorr(N,P,Vars, [N|Argl]) :-
    N > 0, arg(N,P,ArgP),
    member_var(ArgP, Vars,_),
    !,
    N1 is N - 1,
    argscorr(N1, P, Vars, Argl).
argscorr(N, P, Vars, Argl) :-
    N > 0, N1 is N - 1,
    argscorr(N1,P,Vars, Argl).
argscorr(0,_,_,[]).

is_es(F) :-  mb(F, [
'es_F','es_G','es_H', 'es_I', 'es_J', 'es_K', 'es_L', 'es_M',
'es_N', 'es_P', 'es_Q', 'es_R', 'es_S', 'es_T', 'es_U', 'es_V',
'es_G1','es_H1','es_I1','es_J1','es_K1', 'es_L1','es_M1',
'es_N1','es_P1','es_Q1', 'es_R1', 'es_S1', 'es_T1','es_U1',
'es_G2','es_H2', 'es_I2','es_J2', 'es_K2', 'es_L2', 'es M2',
'es N2', 'es P2','es_Q2','es_R2', 'es_S2', 'es_T2','es_U2'
]).

%% trunc_cl

trunc_cl([],_,_,[]).
trunc_cl([P|Pl],F,L, [P1|Pl1]) :-
    functor(P,F,_),
    !,
    trunc_pred(P,L,P1),
    trunc_cl(Pl,F,L,Pl1).
trunc_cl([P|Pl], F,L, [P|Pl1]) :-
    trunc_cl(Pl,F,L,Pl1).

trunc_pred(P,[],P).
trunc_pred(P, [N|L], Ptr) :-
    trunc1(P,N,P1),
    trunc_pred(P1,L,Ptr).

trunc1(P,N,Q) :-
    functor(P,F,M),
    M1 is M - 1,
    functor(Q,F,M1),
    args(M,P, ArgsP),
    N1 is (M - N) + 1,
    remove(ArgsP, N1,ArgsQ,_),
    args(M1,Q,ArgsQ).

args(N,P,[Arg|Args]) :-
    N > 0,
    arg(N,P,Arg),
    N1 is N - 1, args(N1,P,Args).
args(0,_,[]).
