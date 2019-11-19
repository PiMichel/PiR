##################################################################
# Research subject : Financially Optimize a water supply network #
# Author: Pierre MICHEL                                          #
# Tutor: Konstantin BRENNER                                      #
##################################################################

library(slam) #rappeler la matrice simple 
library(Rglpk) # l'optimisation 

#read inp function
Readinp<- function(fname){
L=NULL
connectivity=NULL
q=NULL
n=NULL
ENopen( paste(fname,".inp",sep=""), "test.rpt")
n=c(ENgetcount("EN_LINKCOUNT"))
for (i in 1:n) {
  L=c(append(L,ENgetlinkvalue(i, "EN_LENGTH")))
  connectivity <- rbind( connectivity, ENgetlinknodes(i) )
  q=c(append(q,ENgetnodevalue(i, "EN_BASEDEMAND")))
}

print(connectivity)

inp<-list("length"=L,"connectivity"=connectivity,"demand"=q)
return(inp)
}
#High level functions
#la matrice des coefficients avant la perte de charge hi
ComputeD <- function( connec, upstream_p, coeff, nDiam, numNoeuds, numConduites ){
  #Assemble inequality constraints matrix
  
  nele = Compute_nele(connec, upstream_p, numNoeuds, nDiam)
  row=c()
  col=c()
  val=c()
  row[1]=numNoeuds
  col[1]=numConduites*nDiam+1
  val[1]=1
  row_pt=matrix(0,1,numNoeuds)
  row_pt[numNoeuds]=1
  
  for(i in 1:numConduites){
    result = set_D_branch(i,row, col,val,row_pt,nele,connec, coeff, nDiam, upstream_p) 
    row = result$row
    col = result$col 
    val = result$val
    row_pt = result$row_pt
  }
  
  return (list("row" = row, "col" = col, "val" = val))
}

#la matrice de la direction d'ecoulement 
ComputeA <- function( connectivity ){
  n = max( connectivity )
  m = n-1
  row = c()
  col = c()
  val = c()
  k=1
  for( j in 1:m){
    for (i in 1:2){
      node=connectivity[j,i]
      if (node <n){ 
        row[k] = node
        col[k] = j
        val[k] = (-1)^(i+1)
        k = k+1
      }
    }
  }
  return( simple_triplet_matrix(row, col, val ) )
}


Compute_nele <- function( connec, upstream_p, numNoeuds, nDiam){
  gen=matrix(-1,1,numNoeuds)
  gen[numNoeuds]=0
  
  for(i in 1:(numNoeuds-1))
  {
    gen = set_gen_branch(i, gen, connec, upstream_p) # use connec as global and modifies gen
  }
  return(nDiam*gen+1)
}




#Low level functions

Reynolds <- function (Q,D)
{
  nu=1e-6 
  Re=Q/(nu*D^2*pi/4)
  print (Re)
}

#les coefficients en utilisant la formule de la perte de charge Hazen-Williams
coefficientsHW <- function( Q, D, c)
{
  m = length(Q)  #le nombre de conduites longs
  p =  length(D) #le nombre de diametres a choisir 
  k=1
  ret = matrix(0,1,m*p) #une matrice creuse pour metre tous les coefficients
  for( i in 1:m)
  {
    for (j in 1:p)
    {
      #la formule de la perte de charge Hazen-Williams
      ret[k] = -10.67*abs(Q[i])^0.852*Q[i]/c[j]^1.852/D[j]^4.871
      k = k+1
    }
  }
  return( ret )
}


colebrook<- function(Re, erel) 
{
  #pécision requise
  tolerance = 10^(-6)
  #initialisation de lambda
  lambda0 = 0.02
  #initialisation d'une variable "différence entre deux approximations successives"
  difference = tolerance + 1.0
  #//boucle avec un critère d'arrêt 
  while (difference > tolerance) {
    lambda = ( -2*log10( erel/3.7 + 2.51/Re/sqrt(lambda0) ) )^(-2)
    difference = abs( lambda - lambda0 )
    lambda0 = lambda
  }
  print(lambda)
}

coefficientsDW <- function( Q, D, eps)
{
  m = length(Q)  #le nombre de conduites longs
  p =  length(D) #le nombre de diametres a choisir 
  k=1
  ret = matrix(0,1,m*p) #vecteur-ligne de taille 1xmp
  for( i in 1:m)
  {
    for (j in 1:p)
    {
      #la formule de la perte de charge Darcy
      erel = eps[j]/D[j]
      Re = Reynolds(Q[i], D[j])
      lambda = colebrook( Re, erel )
      
      ret[k]= lambda*Q[i]^2/D[j]^5*8/pi^2
 
      k = k+1
    }
  }
  return( ret )
}





#la fonction pour mettre les chiffres dans l'ordre de l'amont vers l'aval 
set_flow_direction <- function(connectivity,Q)
{
  a=c()
  n = max( connectivity ) 
  m = n-1
  for(i in 1:m)
  {
    # Q est le debit dans la conduite, si Q est negatif, ca veut dire que la 
    #direction n'est pas dans le bon sens, on change les labels au bout de la conduite 
  if(Q[i]<0)
    {
      a = connectivite[i,1] # le noeud a l'amont de la conduite i
      connectivity[i,1] = connectivity[i,2]
      connectivity[i,2]=a # le noeud a l"aval de la conduite i
    }
  }
  return(connectivity) 
}

#la fonction pour chercher la conduite se trouvant al'amont du noeud i,
#c'est-a-dire le noeud i est la sortie de la conduite 
set_upstream_pipes <- function(connectivity){
  upstream_pipes=c() # le vecteur pour mettre le numero de conduite se trouvant a l'amont du noeud i
  n = max( connectivity )
  m = n-1 #le nombre total de conduites
  for(j in 1:m)
  {
    #connectivity[j,2] est le noeuds i se trouvent a l'aval de la conduite j
    upstream_pipes[connectivity[j,2]]=j
  }
  return (upstream_pipes)
}




set_D_branch <- function(i, row, col,val,row_pt,nele,connec, coeff,nDiam, upstream_pipes)
{
  if(row_pt[i]==0)
  {
    upstream = connec[upstream_pipes[i]][1]
    if(row_pt[upstream]==0)
    {
      result = set_D_branch(upstream,row, col,val,row_pt,nele,connec, coeff, nDiam, upstream_pipes) 
      row = result$row
      col = result$col 
      val = result$val
      row_pt = result$row_pt
    }
    
    upstream = connec[upstream_pipes[i]][1]
    upstream_pt = row_pt[upstream]
    
    row=append(row,matrix(i,1,nele[i]))
    
    #copy upstream indices and values 
    col=append(col,col[upstream_pt:(upstream_pt+nele[upstream]-1)])
    val=append(val,val[upstream_pt:(upstream_pt+nele[upstream]-1)])
    
    #add new indices and values 
    col=append(col,( nDiam*(upstream_pipes[i]-1)+1 ):(nDiam*upstream_pipes[i]) )
    val=append(val,coeff[ ( (upstream_pipes[i]-1)*nDiam+1 ):(upstream_pipes[i]*nDiam) ])
    
    row_pt[i] =length(row)-nele[i]+1
    
  }
  return (list("row" = row, "col" = col, "val" = val, "row_pt" = row_pt))
}

set_gen_branch <- function(i, gen, connec, upstream_pipes){
  
  if (gen[i]==-1){   
    upstream = connec[upstream_pipes[i]][1]
    
    if(gen[upstream]==-1) {
      gen = set_gen_branch( upstream, gen, connec, upstream_pipes )
    }
    gen[i]=gen[upstream]+1
  }
  return( gen )
}


ComputeE <- function(connectivity,p){
  n = max( connectivity )
  m = n-1
  row2 = c()
  col2 = c()
  val2 = c()
  k=1
  for( j in 1:m){
    for (i in 1:p){
      row2[k] = j
      col2[k] = (j-1)*p+i
      val2[k] = 1
      k = k+1
    }
  }
  if (is_fixed_level ==1){ 
    row2[k] = m+1
    col2[k] = m*p+1
    val2[k] = 1
    k = k+1
  }
  if (is_fixed_level !=1){ 
    row2[k] = m+1
    col2[k] = m*p+1
    val2[k] = 0
    k = k+1
  }
  return( list( "row" = row2, "col" = col2, "val" = val2 ) )
  
}

#la fonction pour installer les couts pour differentes conduites
set_sigma <- function( nPipe, cout_conduite, cout_electri){
  sigma=c( ) #le vecteur pour mettre le cout 

#pour chaque conduite, il existe plusieurs petits conduites de 
#differents couts  
  for (i in 1:nPipe)
  {
    sigma = c(sigma, cout_conduite)
  }
#ajouter le cout d'electricite   
  if (is_fixed_level == 1)
  {
    sigma = c(sigma, 0)                   
  }else{
      sigma = c(sigma, cout_electri)
  }
}

compute_flow <- function(connectivity, q){
sp_matA <- ComputeA(connectivity)
return( -solve(sp_matA, matrix(q)) )
}

setup_upstream_pipes <- function(connectivity,Q)
{
connectivity <-set_flow_direction(connectivity,Q)
return( set_upstream_pipes(connectivity) )
}
# 
# 
# 
optimize <- function(connectivity,q,D,Q,friction_coef,cout,cout_elec,is_fixed_level,L,zn,z,p_min)
{
  #connectivity: topologie du reseau
  #q: debit de base
  #D: liste de diam disponible
  # Q: debit par lien (arc)
  # friction_coef: soit coef de HW, soit rugosite en fonction de la formule de PDC choisit
  # cout: cout d'un metre de conduite en fonction de diametre
  # cout_elec:
  # a terminer...


  #Script starts here
  m = dim(connectivity)[1] #le nombre des conduites
  n = m + 1 #le nombre des noeuds
  if( n != max(connectivity) )
  { print("Error: Maximal number of nodes != maximal nodes number"); break; }



  coeff <- coefficientsHW( Q, D, friction_coef) #coefficient DH=coeff*L


  connectivity <-set_flow_direction(connectivity,Q)
  upstream_pipes <- set_upstream_pipes(connectivity)


  nDiam = length(D)
  E <- ComputeE( connectivity, nDiam )
  sp_matE <- simple_triplet_matrix( E$row, E$col, E$val )
  mat_eq <-sp_matE


  nNode =length(q)+1
  nPipe = nNode-1 #idem

  
  
  sigma <- set_sigma( nPipe, cout, cout_elec)

  if (is_fixed_level == 1){
    rhs_eq = c(matrix(L), zn)
  } else
  {
    rhs_eq =c(matrix(L),0)
  }

  rhs_ineq <- z+p_min
  sp_D <- ComputeD( connectivity, upstream_pipes, coeff, nDiam, n, m )
  mat_ineq <- simple_triplet_matrix(sp_D$row,sp_D$col, sp_D$val)
  S <- rbind(mat_ineq, mat_eq )##combinaison de D et E

  k1 <- length(rhs_ineq)
  k2 <- length(rhs_eq)
  dir <- matrix(NA,1,k1+k2)
  for ( i in 1:(k1+k2)){
    if (i <= k1){
      dir[i]=">="
    }
    else{
      dir[i]="=="
    }
  }
  
  
  obj <- sigma
  mat <- S
  rhs <- c(rhs_ineq,rhs_eq) ##combinaison de z et L
  max <- FALSE

  
  return( Rglpk_solve_LP(obj, mat, dir, rhs,max=max) )
}


