

% Preuve puis tronquer Z :
% Synthétise la relation qui calcule
% la longueur de la concatenation de X et Y
% sans calculer explicitement la concatenation
% 
% app(X,Y,?Z)&leng(?Z,?T).


% idem mais avec une preuve differente
%
% app(X,Y,Z) => leng(Z,?T).


% Preuve puis tronquer Y :
% Synthétise la relation same_length/2
% sans calculer la longueur
%
% leng(X,?L)&leng(?Y,?L).  


% Synthetise la soustraction : Z = Y-X
%
% infe(X,Y)=>add(X,?Z,Y).  

  
% Synthetise une relation qui calcule fib(N)
% en évitant un double appel récursif 
%  
% fib(X,?Y) & fib(s(X),?Z).
