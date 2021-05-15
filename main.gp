printmessage(m) = print(concat([Str(x)|x<-m]));
[N,k,c] = readvec("input.txt");
\\k = [b_1,...b_n]
\\c chiffré = sum m_i b_i



/* ********* Travail effectué en binôme (Daphné-Marina) ********* */


/* ****************************************************************************** */
/* ******************************* Explications : ******************************* */
/* ****************************************************************************** */


/*
 * RAPPEL Cryptosystème Sac à Dos :
   * N>sum_{i=1}^{n} (a_i)
   * b_i = k*a_i mod N
   * clef publique b1...bn
   * chiffrement    m |-> sum m_i b_i [N]
   * dechiffrement  k^{-1}*c [N] = sum m_i a_i, puis méthode gloutonne.

 * Pour résoudre le problème du sac à dos, on cherche à appliquer l'algorithme LLL.
   Pour ce fair, on utilisera la fonction qflll(x,{flag=0}):
   qflll(x,{flag=0}): LLL reduction of the vectors forming the matrix x
   (gives the unimodular transformation matrix T such that x*T is LLL-reduced).
   flag is optional, and can be 0: default, 1: assumes x is integral,
   2: assumes x is integral, returns a partially reduced basis,
   4: assumes x is integral, returns [K,T] where K is the integer kernel of x and T the LLL reduced image,
   5: same as 4 but x may have polynomial coefficients, 8: same as 0 but x may have polynomial coefficients.

 * On cherche m_1,..., m_n tels que c = sum m_i b_i.

 * On va construire une matrice M telle que décrite dans le cours,
   puis utiliser qflll dessus, afin de retrouever le vecteur [m1,...,m_n].
 * L'algorithme LLL nous permettra de trouver de petits vecteurs qui satisfont l'équation.
 * On gardera le vecteur [m_1-1/2,...,m_n-1/2] de norme sqrt(n)/2,
   qui correspond à la solution m_1,...,m_n.
 
 
 * Le principe est, après avoir créé la matrice, de trouver une colonne de la
   base réduite telle que tous les coefficients soient égaux à +/- 0,5.
   On calcule alors s= g*k, avec k = [b1,......,bn] et si s==c, on renvoie
   la solution, qui est le texte en clair g.
   
 
 * Notre implémentation ne marchant pas (problème de densité ?),
 * nous avons décidé de nous pencher sur une autre méthode de résolution (merci Marc).
 * Elle consiste à effectuer une disjonction de cas sur m_1 à valeurs dans {0,1}.
 * Dans les deux cas, on met à jour la matrice A, puis on calcule une base réduite
 * puis, pour chaque colonne, on regarde si elle vérifie 
   la propriété de ne posséder que deux coefficients, dont un nul.
* Si c'est le cas, et que le message vérifie l'égalité sac à dos, 
  on renvoie le message.
 
*/



/* ****************************************************************************** */
/* ******************************** Fonctions : ********************************* */
/* ****************************************************************************** */

\\http://j.deroff.free.fr/rapportter.pdf

/* Création de la matrice M */
create_matrix(B,n)={
  my(A);
  A=matid(n);
  A[n,n]=-c*B;
  for(i=1, n-1, A[i,n]=-1/2; A[n,i]=B*k[i+1]);
  return (A);
};





/* Tests pour la première méthode (qui ne marche pas) */
/*
\\Recherche d'une colonne g de G telle que tous les coeffs valent +/- 0.5,
\\sauf le dernier qui vaut 0 
get_right_column(G,n)={
  my(i,j);
  \\i donne la ligne
  \\j donne la colonne
  i=0; j=0;
  for(j=1,n, 
	  for(i=1, n, 
  	  if(abs(G[i,j]) != 0.5, break);
      if(i==n && G[i+1,j]==0 , return (j));  
    );
  );
  return (-1);
};


\\ on ajoute 0.5 à la colonne trouvée et on appelle le résultat g
\\ on calcule s= g*k, avec k = [b1,......,bn]

get_s(G,j,k,n)={
  my(g,s);
  g=vector(n);
  for(i=1, n, g[i] = 0.5 + G[i,j]);
  return ([g*(k~)%N,g]);
};

\\Si s=c, on a trouvé une solution et le texte en clair est le vecteur g.

\\Sinon :
\\On re-définit le vecteur g =0.5−g, 
\\On re-calcule s=gB. Si s=c alors le texte en clair est le vecteur g.
get_s_prime(G,j,k,n)={
  my(g);
  g=vector(n);
  for(i=1, n, g[i] = 0.5 - G[i,j]);
  return ([g*(k~)%N,g]);
};
*/




/* Tests pour la seconde méthode */


/* 
 * On cherche une colonne qui ne posséde que deux coefficients, dont un nul.
 * Si une telle colonne existe, on renvoie l'indice de la colonne
*/
right_coefs(v)={
  my(s);
  s = Set(v);
  return ( (#s==2) && setsearch(s,0) );
};


/* Fonction décalant les coefficients d'un vecteur de  v de val, sauf le dernier */
shifts(v, val)={
	for(i=1,#v-1, v[i]=v[i]+val);
	v;
};


/*
 * Permet de trouver le message en appliquant la seconde méthode décrite plus haut.
*/
find_message(B,n,c)={
	my(i,A,G,g);
	i=1;
	A=create_matrix(B,n);
	while(true,
	  for(l=0, 1,
	    if(l==0,
	      A[n,n]= - B*(c-l+i*N),
	      A[n,n]= - B*(c-k[l]+i*N);
	    );
	    G = A*qflll(A);
	    for(j = 1, n,
              g=G[,j]~;
	      g = shifts(g, 1/2);
	      if(right_coefs(g), 
	        g = concat(Vec(l), g[1..-2]);
	        if((g*(k~))%N == c, printmessage(g); return (g);
	        );
	      );
	    );
	  );
  	  i=i+1;
	  \\print(i);
	);
};








/* ****************************************************************************** */
/* ******************************* Application : ******************************** */
/* ****************************************************************************** */



n = #k; \\n=140
B = ceil(sqrt(n * 2^n));
find_message(B,n,c);
