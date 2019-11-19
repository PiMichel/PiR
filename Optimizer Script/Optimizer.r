##################################################################
# Research subject : Financially Optimize a water supply network #
# Author: Pierre MICHEL                                          #
# Tutor: Konstantin BRENNER                                      #
##################################################################

# Source de mes fonctions
source('functions.R')
library(epanet2toolkit)

p_min = 4 #la pression minimale 

# Diametres disponibles
Dmin = 10
Dmax = 200
nD = 9
dD = (Dmax - Dmin)/(nD-1)
D = (Dmin + dD*(0:(nD-1)))*(10^(-3))

cout = 20 + 0.1*D # tableau prix par diametre [euro/m]


 data <- Readinp("test")
 connectivity <- data$connectivity
 q=data$demand
 L=data$length

#la topographie
z = matrix(25, byrow=T, nrow=6)
zn=60 #charge totale du reservoir 
is_fixed_level =1 #numero de de la bache

# Hazen-Williams coef ou rugosite si Darcy est utilise
friction_coef = matrix(100,1,nD) 

# les diametres disponibles pour les conduites
Q<-compute_flow(connectivity, q)

  resolution <- optimize(connectivity,q,D,Q,friction_coef,cout,cout_elec,is_fixed_level,L,zn,z,p_min)
  if( resolution$status != 0 )
  {
    print( "ERROR: Optimal solution has not been found!" )
    break
  }
 
#-----------------------------------------------------------------------------------
#Print des résultats
cat("Diamètres disponibles :", D,"\n")
for (j in 0:(length(L)-1)){
cat("Optimisation de la conduite numéro :",j+1)
j=j*nD
print(round(resolution$solution[(1+j):(j+nD)],0)) #arrondit au centimère pres pour la mise en forme du résultat
}

#------------------------------------------------------------------------------
#Calcul des diametres "moyens" pour le nouveau fichier .inp
Dopti=NULL
for (j in 0:(length(L)-1)){
i=j+1
j=j*nD
Dopti=c(append(Dopti,sqrt(sum(resolution$solution[(1+j):(j+nD)]*(D[1:9]^2))/L[i])))
}
cat("Diamètre optimisé pour EPANET : ", Dopti, "\n")

#Ecriture dans un second fichier inp avec les nouveaux diamètres
for (i in 1:(length(L))){
ENgetlinkvalue(i, "EN_DIAMETER")
ENsetlinkvalue(i, "EN_DIAMETER", Dopti[i]) 
}
ENsaveinpfile("test2.inp")
ENclose()