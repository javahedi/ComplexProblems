!11juin18:  write Hamiltonian HK for the k vector read from the file KVECTOR

MODULE precision_reel
  integer, parameter :: prec = kind(1.d0)        ! precision des reels 
  real(kind=prec),parameter :: zero = 1.D-3
end MODULE precision_reel 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE constantes_physiques
USE precision_reel
  real(kind=prec),parameter :: pi = 3.141592654d0
  real(kind=prec),parameter :: twopi = 6.28318530717958648d0
  real(kind=prec),parameter :: hh = 6.62620d-34, hbar = hh / twopi,hbar2 = hbar*hbar
  real(kind=prec),parameter :: e_electron = 1.60219d-19
  real(kind=prec),parameter :: ry_ev = 13.6058d0, ry_j = 13.6058d0 * e_electron
  real(kind=prec),parameter :: ua = 0.529177d0
  real(kind=prec),parameter :: clumiere = 2.99792458d8

end MODULE constantes_physiques 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE constantes

  integer, parameter :: dimension_reseau = 3
  integer, parameter :: dimension_pos_atom = 3
  integer, parameter :: nspin = 1  ! 1 pas de spin, 2 avec spin (le nb d'orbitale est doublé)

  ! tailles max pour les  tableaux :

  integer, parameter :: natomes_max = 12000 ! nb max d'atomes

  integer, parameter :: ntypes_atomes_max =  2
                        !                =  nb max de types d'atome

  integer, parameter :: lm_max = 1  ! nb obitale par atome max orbitale s p et/ou d uniquement

  integer, parameter :: ndim_max =  natomes_max*lm_max !  (sans spin)
                        !        = nb max d'orbitale dans la maille 

  integer, parameter :: npas_dos_max = 2001

  integer, parameter :: nb_print_max = 10  ! nb max d'atomes affiches...
                                            ! pour limiter le fichier de sortie
 

  ! nom des fichier d'entree:
  character*30, parameter :: nom_structure="INPUT_structure"
  character*30, parameter :: nom_atomes="INPUT_atomes"
  character*30, parameter :: nom_SlaterKoster="INPUT_SlaterKoster"
  character*30, parameter :: nom_DOS_transp="INPUT_DOS_transp"
  character*30, parameter :: nom_bnds="INPUT_bnds"
  character*30, parameter :: nom_fichier_DOS_all="dos_all"
  character*30, parameter :: nom_fichier_DOS_all_at="dos_all_at"
  character*30, parameter :: nom_fichier_part_ratio="part_ratio"
  character*30, parameter :: nom_fichier_part_ratio_m="part_ratio_moy"
  character*30, parameter :: nom_fichier_DOS_part="DOS_part"
  character*30, parameter :: nom_fichier_bnds="bnds"
  character*30, parameter :: nom_fichier_encours="EN_COURS"
  character*30, parameter :: nom_fichier_spin="INPUT_SO"

end MODULE constantes
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE VarReseau           ! Pour cristal  3D 
  USE precision_reel 
  USE constantes

  real(kind=prec) :: a_lattice       ! parametre de maille 
  real(kind=prec),dimension(1:dimension_reseau,1:dimension_pos_atom) :: a,b 
  !    a et b : vecteurs de maille du reseau reel et du reseau reciproque 
  integer :: natomes ! nb d'atomes dans une maille

  integer, dimension(0:natomes_max) :: type_atome 
            ! type_atome(0) = nb de type differents
            ! type_atome(i) = type de l'atome i 
  
  real(kind=prec),dimension(1:natomes_max,1:dimension_pos_atom) :: vec_Rat !pos. at en unite de a
  integer, parameter :: nbvoisins_max = 100
  integer,dimension(1:natomes_max) :: nbvoisins ! Nb de voisins du site i  
  integer,dimension(1:natomes_max,1:nbvoisins_max,0:dimension_reseau) :: voisins 
      ! voisins(i,j,0) est le jieme voisin du site i
      ! voisins(i,j,1:3) indique dans quelle maille est ce voisin (vec R)

  real(kind=prec) dist_voisin(1:ntypes_atomes_max,1:ntypes_atomes_max)

  
end MODULE VarReseau

MODULE VarHamil          
  USE precision_reel 
  USE constantes
  USE VarReseau

  integer ndim     ! nb total d'orbitale dans la maille (sortie) 
  integer, dimension(1:ntypes_atomes_max,0:lm_max) :: orbitales_lm
       ! orbitale_lm(i,0) nb d'orbitales sur l'atome de type i
       ! orbitale_lm(i,p) valeur lm de l'orbitale p de l'atome de type i 
  integer, dimension (1:natomes_max,1:lm_max) :: indice_orb
       !  indice_orb(i,j) indice de 1 a ndim de la jeme sur l'atome i
  complex(kind=prec),  & 
       dimension(1:natomes_max,0:nbvoisins_max,1:lm_max,1:lm_max,1:nspin,1:nspin) :: g,gvx
       ! g(i1,j,lm1,lm2) saut entre lm1 de l'atome i et lm2 du voisin numero j de i

  ! estimation rapide des bords de bandes
  real(kind=prec) ee_min,ee_max

end MODULE VarHamil


program writeHk


USE precision_reel
USE constantes
USE VarReseau
USE VarHamil
implicit none

  integer nb_total_pk
  real(kind=prec) vec_k(1:dimension_reseau),zero2
  complex(kind=prec),dimension(1:ndim_max*nspin,1:ndim_max*nspin) :: Hk
  integer i,j

zero2 = 1.0d-10

! Entree reseau et ecriture tableau voisins
call Reseau_atome

! Entree Hamiltonien et ecriture du tab. g(tau avec voisins)
call Hamil_Vx_entree
!print*,'Hamil_Vx_entre OK'

open(18,file='KVECTOR') ! Fichier d'entree
 read(18,*) vec_k(1:dimension_reseau)
 nb_total_pk = 1 
close(18)
print*,''
print*,'Calculation of Hk for k = ',vec_k(1:dimension_reseau)

! Ecriture de Hk et diagonalisation
call Ecriture_Hk(vec_k,Hk)  ! sortie: Hk 
!print*,'ecriture tableau Hk OK'
!print*,'ndim = ', ndim

open(18,file='HK')
!write(18,*) ndim*nspin
do i=1,ndim*nspin 
   do j=1,i ! ndim*nspin
   !  print*,'i,j,hij = ',i,j, Hk(i,j)
   if (sqrt(real(Hk(i,j)*conjg(Hk(i,j)))).gt.zero2) then
	
      write(18,*) i, j,  real(Hk(i,j)),  aimag(Hk(i,j))
   endif
   enddo
enddo
close(18)

!call sleep(10)

end program writeHk




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Reseau_atome
!
! Positions des atomes et recherche des voisins
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE precision_reel
USE constantes
USE VarReseau
IMPLICIT  NONE
real(kind=prec) R(1:dimension_pos_atom),rij(1:dimension_pos_atom),dij,dij_min
integer i,j,n_a(1:3),n1,n2,n3,index_clas(1:nbvoisins_max)
real(kind=prec) vecteur_c(1:3), prod_mix,distij(1:nbvoisins_max,0:3)


character nom*30

print*,''
print*,'Sub: Reseau_atome :    Pour cristal 3D '
print*,''

open(18,file=nom_structure,action='read')
!print*,'KL'
! Vecteur du reseau
read(18,*) a_lattice
print*,'Parametre de maille a (Ang.) = ',a_lattice
do i=1,dimension_reseau
   read(18,*) a(i,:)
enddo


print*,' Vecteurs du reseau (unite de a) :'
do i=1,dimension_reseau
   print "('    a(',I1,') =',3F10.5)",i,a(i,:)
enddo

!! Reseau reciproque
!!b(1,:) = (/ 1.0 , 0.0 , 0.0 /)
!!b(2,:) = (/ 0.0 , 1.0 , 0.0 /)
!vecteur_c(1:3) = (/ 0.0 , 0.0 , 1.0 /)

call produit_vec(a(2,:),a(3,:),b(1,:))
prod_mix = DOT_PRODUCT(a(1,:),b(1,:))
b(1,:) = b(1,:)/prod_mix
call produit_vec(a(3,:),a(1,:),b(2,:))
b(2,:) = b(2,:)/prod_mix
call produit_vec(a(1,:),a(2,:),b(3,:))
b(3,:) = b(3,:)/prod_mix

print*,' Vecteurs du reseau reciproque (unite de 2pi/a) :'
do i=1,dimension_reseau
   print "('    b(',I1,') =',3F10.5)",i,b(i,:)
enddo


! Lecture des positions atomiques
! La position des atomes doit etre donne en 3D
read(18,*) natomes ! Nb total d'atomes
if (natomes.gt.natomes_max) stop &
     '!!!!!!!!!! nb atomes > natomes_max (variable.f90)'

read(18,*) type_atome(0)  ! nb d'atome de type different (inequivalent)

if (type_atome(0)>ntypes_atomes_max) stop & 
     '!!!!!!!!!! nb de type atomes > ntypes_atomes_max (variable.f90)'

do i=1,natomes
   read(18,*) type_atome(i),Vec_Rat(i,1:3)
   !print*,'type, ri = ', type_atome(i),Vec_Rat(i,1:3) 
   if (type_atome(i) > type_atome(0))   &
    stop '!!! Sub:Reseau_atome: type_atome > nb type atome (variable.f90)'
enddo

!!read(18,*) type_atome(0)  ! nb d'atome de type different (inequivalent)
!!if (type_atome(0)>ntypes_atomes_max) stop & 
!!     '!!!!!!!!!! nb de type atomes > ntypes_atomes_max'
!!do i=1,natomes
!!   read(18,*) type_atome(i)
!!   if (type_atome(i) > type_atome(0))   &
!!        stop '!!! Sub:Reseau_atome : type_atome > nb de type atome'
!!enddo

print*,''
print "('Nb de atomes dans la maille initiale : natomes =',I7)", natomes
print "('Nb de types atomes : type_atome(0) = ',I3)",type_atome(0)
print*,'Nb atome : Position atome dans la maille (unite de a), type : '
do i=1,natomes
   print "(I7,' : ',F10.5,' ',F10.5,' ',F10.5,' , de type : ',I3)",i,vec_Rat(i,1:3) &
        ,type_atome(i)
enddo
close(18)

! Recherche des voisins:

! Si les distances max dependent de la nature des atomes 
!   ==> elles sont dans le fichier nom_SlaterKoster
  open(19,file=nom_SlaterKoster,action='read')

  read(19,*) nom

  read(19,*) i 
  if (i<type_atome(0)) STOP & 
       '!!!!! INPUT Slater-Koster : nb type atomes incompatible !' 

  do i=1,(type_atome(0)*(type_atome(0)+1)/2)
     read(19,*) nom 
     read(19,*) n1,n2
     read(19,*) dist_voisin(n1,n2)
     dist_voisin(n2,n1) = dist_voisin(n1,n2)
     do j=0,14                                       ! avril17
        read(19,*) nom
     enddo
  enddo
close(19)

! SI la distance max est dans le fichier nom_structure
!read(18,*) dist_voisin       ! distance voisin max (en unite de a)
!close(18)
print*,''
print*,'Voisins jusqu"a dist_voisin :'
print*,'  Type atome, type atome, dist (unite de a), dist (Ang)'
do n1=1,type_atome(0)
   do n2=1,type_atome(0)
        print*,'   ',n1,',  ',n2,',  ',real(dist_voisin(n1,n2)),',  ',real(dist_voisin(n1,n2)*a_lattice)
   enddo
enddo

! Parametres pour la recherche des voisins dans les mailles voisines
! +/- n_a(i) maille dans la direction a(i)
if(natomes.lt.20) then 
 n_a(1) = 20 !100
 n_a(2) = 20 ! 100
 n_a(3) = 0 ! 100
else if(natomes.lt.100) then
 n_a(1) = 10
 n_a(2) = 10
 n_a(3) = 0
else 
 n_a(1) = 2 
 n_a(2) = 2
 n_a(3) = 2
endif


print*,''
print*,'Pour la recherche des voisins de chaque atomes de la maille (0,0,0)'
print*,'   parametres n_a(i) (ie. +/- n_a(i) dans la direction a(i)):',n_a(:)




nbvoisins = 0
voisins = 0

do i=1,natomes
   !4avril18 print*,''
   !4avril18 print "('Voisins de atome nb',I7,' (type ',I3,'),: r_i = ',3(F11.6,' '))" & 
   !4avril18       ,i,type_atome(i),vec_Rat(i,1:3)
   !print "('Voisins de atome nb',I7,' : r_i = ',3(F11.6,' '))" ,i,vec_Rat(i,1:3)
   dij_min = 1000
   do n1=-n_a(1),n_a(1),1
   do n2=-n_a(2),n_a(2),1
   do n3=-n_a(3),n_a(3),1
      do j=1,natomes
         R = a(1,:)*n1 + a(2,:)*n2 + a(3,:)*n3
         rij = R + vec_Rat(j,:) - vec_Rat(i,:)
! if (i.eq.1) then 
!            write(81,*) real(rij(1:3))
! endif
         dij = sqrt(dot_product(rij,rij))
!         print*,dij,i,j,type_atome(i),type_atome(j),n1,n2
!         print*,'rij=',rij


         if (dij<=zero) then
            if (i.ne.j) then 
               print*,'i,j,rij =',i,j,rij
               print*,'vec_Rat(i,:)',vec_Rat(i,:)
               print*,'vec_Rat(j,:)',vec_Rat(j,:)
               stop '!!!!!Sub: Reseau_atome : atomes i et j ont meme position:'
            endif
	!!!! POUR supprimer les paires inter-plans A-A d'atomes juste au dessus l'un de l'autre.
        !!!! else if (dij<=dist_voisin(type_atome(i),type_atome(j)).and.((rij(1)**2+rij(2)**2)>(2*zero**2))) then
	!!!! sinon calcul normal :
         else if (dij<=dist_voisin(type_atome(i),type_atome(j))) then
!            if (abs(n1).eq.n_a(1).or.abs(n2).eq.n_a(2).or.abs(n3).eq.n_a(3)) then
!		print*,'Sub.: Reseau !!!!!! Augmenter n_a(:) = (n max) dans le recherche des voisins:'
!                print*,'  !!!!!!! n1 n2 n3 = ',n1, n2, n3
!	 	print*,'  !!!!!!! n_a(1), n_a(2), n_a(3) = ',n_a(1), n_a(2), n_a(3)
!		stop
!	    endif
            nbvoisins(i) = nbvoisins(i) + 1
            if (nbvoisins(i) > nbvoisins_max) &
               stop'!!!!!Sub: Reseau_atome : nbvoisins(i) > nbvoisins_max'
            voisins(i,nbvoisins(i),0) = j
            voisins(i,nbvoisins(i),1) = n1
            voisins(i,nbvoisins(i),2) = n2
            voisins(i,nbvoisins(i),3) = n3
            distij(nbvoisins(i),0) =  dij
            distij(nbvoisins(i),1:3) =  rij(1:3)
            !4avril18 print "(I7,' (type ',I3,'), dij (a) = ',F8.6,' R = ',3I3,', r_ij = ',3(F10.6,' '))" &
            !4avril18     ,j,type_atome(j),dij,n1,n2,n3,(rij)
            if (dij<dij_min) then 
               dij_min = dij
            endif
         endif
      enddo
   enddo
   enddo
   enddo

   if (i.le.minval((/nb_print_max,natomes /))) then 
       call classement(distij(1:nbvoisins_max,0),index_clas(1:nbvoisins_max),  &
                              nbvoisins_max,nbvoisins(i))
       print*,''
       print "('Voisins de atome nb',I7,' (type ',I3,'),: r_i = ',3(F11.6,' '))" & 
         ,i,type_atome(i),vec_Rat(i,1:3)
       print "('   nb de voisins = ',I4)",nbvoisins(i)
       do n2=1,nbvoisins(i)
            n1 = index_clas(n2)
            print "(I7,' (type ',I3,'), dij (a) = ',F8.6,' R = ',3I3,', r_ij = ',3(F10.6,' '))" &
                 ,voisins(i,n1,0),type_atome(voisins(i,n1,0)),distij(n1,0),voisins(i,n1,1), &
                 voisins(i,n1,2),voisins(i,n1,3),distij(n1,1:3)
       enddo
       print "('nb de voisins de atome i = ',I3,' : ',I3)",i,nbvoisins(i)
       print "('distance minimale ij = ',F10.6,' unite de a = ',F10.6,' Ang.')", & 
        dij_min,dij_min*a_lattice
   end if
enddo

end subroutine Reseau_atome

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE produit_vec(a,b,c)
!
!	c = produit vectoriel a x b             vecteurs a 3 dimension
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE precision_reel
IMPLICIT  NONE
real(kind=prec) a(1:3),b(1:3),c(1:3)

 c(1) = a(2)*b(3) - a(3)*b(2)
 c(2) = a(3)*b(1) - a(1)*b(3)
 c(3) = a(1)*b(2) - a(2)*b(1)


end subroutine produit_vec


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine classement(x,index_clas,dimens_tab,nmax)
!
! entree : x,dimens_tab,nmax
! sortie : index_clas
!          x(index_clas(1:nmax)) sont classes dans l'ordre croissant
!
USE precision_reel 
implicit none
integer dimens_tab
real(kind=prec), dimension(1:dimens_tab) :: x
integer, dimension(1:dimens_tab) :: index_clas
integer nmax,n,i,j

!allocate (x(1:dimens_tab),index_clas(1:dimens_tab))
!print*,'avant classement'
!print*,x(1:nmax)
index_clas = 0
index_clas(1)=1
do n=2,nmax
   i=1
   do while (x(n).ge.x(index_clas(i)).and.i.lt.n)
     i=i+1
  enddo
  !print*,'n,i = ',n,i
  do j=n-1,i,-1
    !if(index_clas(j).ge.i) then 
     index_clas(j+1)= index_clas(j)
    !endif
  enddo
  index_clas(i) = n
  !print*,'index_clas = ',index_clas(1:n)
  !print*,x(index_clas(1:n))
enddo
!print*,'apres classement'
!print*,x(index_clas(1:nmax))
!print*,'index_clas = ',index_clas(1:nmax)
!deallocate(x,index_clas)
end subroutine classement



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Hamil_Vx_entree
!
!    version du 25/03/2008 : 
!    Ecriture du tableau gvx pour calcul de Vx
!
!    Ecriture du tableau g avec terme de saut entre i et ses voisins
!
! 2janv18: avec spin-orbite pour orbitale p et d
! 9janv18: lecture parametre spin-orbite pour orbitale p et d and INPUT_SO
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE precision_reel 
  USE constantes
  USE VarReseau
  USE VarHamil
  USE constantes_physiques
  implicit none
 
  integer ntype,nb_types_at,lm,i,j,orb1,orb2,lm_i,lm_j,i2
  integer i_type,i2_type,ndim_res,lm_1,lm_2
  character nom*30
  real(kind=prec) ee_lm

  real(kind=prec),dimension(1:ntypes_atomes_max,0:lm_max) :: energie_at_type
  real(kind=prec),dimension(1:ntypes_atomes_max,1:ntypes_atomes_max,0:14,1:4) :: &   !17avril17
       SlaterKoster_par

  real(kind=prec) R_i2(1:dimension_pos_atom),dmin,SlaterKoster_AttenD

  character(len=7), dimension(1:9) :: nom_orbitale
!  character*7 nom_orbitale(1:9)

  complex(kind=prec),parameter :: i_imag = (0.0d0,1.0d0)

  real(kind=prec),parameter :: phi0 = hh*clumiere / e_electron
  real(kind=prec) champB,coef_phase_g,dunite

  real(kind=prec) direction_transp(1:dimension_reseau)

!2janv18 
  complex(kind=prec),dimension(1:6,1:6) :: MSOCp  
  complex(kind=prec),dimension(1:10,1:10) :: MSOCd  


  real(kind=prec),dimension(1:2,1:ntypes_atomes_max) :: lambda_SOC  
  ! avec           lambda_SOC(l,ntype) l=1 orbital p  et l=2 orbital d



direction_transp = 0.0d0
direction_transp(1) = 1.0d0
print*,''
if (nspin.eq.1) then
  print*,'Calcul sans spin: nspin = 1.'
else
  print*,'Calcul avec spin et spin-orbit (pour orbitales p et d uniquement'
endif
print*,'!!!! Sub Hamil_Vx_entree: transport pour direction : ',real(direction_transp)

!nom_orbitale(1:9) = (/'s','px','py','pz','dxy','dyz','dzx','dx2-y2','d3z2-r2'/)
nom_orbitale(1:9) = (/'1','2','3','4','5','6','7','8','9'/)

orbitales_lm = 0
energie_at_type = 0.0d0
SlaterKoster_par = 0.0d0


! 19/05/2008   15/05/2008
  dunite = a_lattice ! en Ang. 
  print*,' dunite = ',dunite,' Ang.'

! Lecture du fichier atome (nb d'orbitale par type d'atome)
open(18,file=nom_atomes,action='read')
read(18,*) nb_types_at
if (nb_types_at.ne.type_atome(0)) & 
     stop '!!! nb types atomes dans INPUT_atome differents !'

print*,''
print*,'Sub: Hamil_Vx_entree : lecture des parametres du hamiltonien'
print "(I3,' types atomes differents :')",nb_types_at
do ntype = 1, type_atome(0)
   read(18,*) nom
   read(18,*) orbitales_lm(ntype,0)
   if(orbitales_lm(ntype,0).gt.lm_max) &
    stop '!!!STOP: INPUT_atome : nb orbitals > lm_max (variable.f90)'
   do i=1,orbitales_lm(ntype,0)
      read(18,*) lm,ee_lm
      orbitales_lm(ntype,i) = lm
      energie_at_type(ntype,i) = ee_lm
   enddo
enddo
close(18)
do ntype = 1, nb_types_at
   print "('  atome de type : ',I3,' contient ',I3,' orbitales : ')" & 
        ,ntype,orbitales_lm(ntype,0)
   print "('        lm = ',9(I7,' ;'))", orbitales_lm(ntype,1:orbitales_lm(ntype,0))
   print "('      soit : ',9(A7,' ;'))", nom_orbitale(orbitales_lm(ntype,1:orbitales_lm(ntype,0)))
   print "('        Eii= ',9(F7.3,' ;'))",energie_at_type(ntype,1:orbitales_lm(ntype,0))
enddo
! fin lecture fichier atome (nb d'orbitale par type d'atome)


! Lecture du fichier spin-orbite (uniquement si nspin=2)     ! 9janv18
lambda_SOC = 0.0d0
if(nspin.eq.2) then
 open(18,file=nom_fichier_spin,action='read')
 read(18,*) nb_types_at
 if (nb_types_at.ne.type_atome(0)) & 
     stop '!!! nb types atomes dans INPUT_SO differents !'
 print*,''
 print*,'Sub: Hamil_Vx_entree: nspin=2, lecture des parametres du couplage spin-orbite dans ',nom_fichier_spin
 !print "(I3,' types atomes differents :')",nb_types_at
 do ntype = 1, type_atome(0)
   read(18,*) nom
   read(18,*) lambda_SOC(1,ntype), lambda_SOC(2,ntype)
 enddo
 close(18)
endif
! fin Lecture du fichier spin-orbite (uniquement si nspin=2)


! Lecture des parametres de SlaterKoster et du champ magnetique champB (suivant z) en eV
call Lec_SlaterKoster_par(ntypes_atomes_max,nb_types_at, &
                          nom_SlaterKoster,SlaterKoster_par,champB)

 coef_phase_g = pi * champB / phi0
 print*,''
 print*,'  Champ magnetique B = ',real(champB),' soit coef phase = ',real(coef_phase_g)


!2janv18 Spin-Orbite SO
!
MSOCp= cmplx(0.d0,0.d0) 
MSOCd= cmplx(0.d0,0.d0) 
if (nspin.eq.1) then
 print*,''
 print*,'Sub Hamil_Vx_entree: spin=1 : Calcul sans spin (et donc sans couplage spin-orbite)'
else if (nspin.eq.2) then
 print*, ' '
 print*,'Sub Hamil_Vx_entree: nspin=2 : avec couplage-spin orbite pour orbitale p :'
 ! d'apres Boyer-Richard J. Phys. Chem. Lett. 2016, 7, 3833−3840 page 2825
 ! px+ py+ pz+ px- py- pz-
 ! 1   2   3   4   5   6
 MSOCp(1,2)= -i_imag 
 MSOCp(1,6)= 1.d0
 MSOCp(2,1)= i_imag
 MSOCp(2,6)= -i_imag
 MSOCp(3,4)= -1.0d0
 MSOCp(3,5)= i_imag
 MSOCp(4,3)= -1.d0
 MSOCp(4,5)= i_imag
 MSOCp(5,3)= -i_imag
 MSOCp(5,4)= -i_imag
 MSOCp(6,1)= 1.d0
 MSOCp(6,2)= i_imag
 print*,"d'apres Boyer-Richard J. Phys. Chem. Lett. 2016, 7, 3833−3840 page 2825 :" 
 print*,'Matrice MSOC dans la base px+ py+ pz+ px- py- pz- :'
 do i=1,6
  !print*,MSOCp(i,:) 
  print "(6('(',(F7.2),',',(F7.3),')',1X))", real(MSOCp(i,1)),imag(MSOCp(i,1)), &
     real(MSOCp(i,2)),imag(MSOCp(i,2)),real(MSOCp(i,3)),imag(MSOCp(i,3)),real(MSOCp(i,4)),imag(MSOCp(i,4)), &
     real(MSOCp(i,5)),imag(MSOCp(i,5)),real(MSOCp(i,6)),imag(MSOCp(i,6))
  do j=1,i  !verification hermiticite de MSOCp
    if (MSOCp(i,j).ne.conjg(MSOCp(j,i))) stop 'matrice couplage spin-orbit p non hermicienne'
  enddo
 enddo
 print*,'H_SOC = 0.5 lambda_SOC MSOCp, avec :'
 print*,'  lambda (SOC) pour atome de type 1 orbitales p = ',lambda_SOC(1,1)
 print*,'  lambda (SOC) pour atome de type 2 orbitales p = ',lambda_SOC(1,2)

 print*, ' '
 print*,'Sub Hamil_Vx_entree: nspin=2 : avec couplage-spin orbite pour orbitale d :'
 ! 5 dxy+, 6 dyz+, 7 dzx+, 8 dx2-y2+, 9 d3z2-r2+   dxy-,  dyz-,  dzx-,  dx2-y2-,  d3z2-r2-
 !      1       2       3          4           5      6      7      8         9         10
 MSOCd(2,3)= i_imag !
 MSOCd(6,3)= -i_imag !
 MSOCd(9,3)= -1.0d0 ! 
 MSOCd(10,3)= sqrt(3.0d0) !
 MSOCd(7,8)= -i_imag !
 MSOCd(1,8)= i_imag !
 MSOCd(4,8)= 1.0d0 !
 MSOCd(5,8)= -sqrt(3.0d0) !
 MSOCd(8,5)= -sqrt(3.0d0) !
 MSOCd(7,5)= -i_imag*sqrt(3.0d0)  !
 MSOCd(3,10)= sqrt(3.0d0) !
 MSOCd(2,10)= -i_imag*sqrt(3.0d0) !
 MSOCd(3,2)= -i_imag !
 MSOCd(6,2)= -1.0d0 !
 MSOCd(9,2)= i_imag !
 MSOCd(10,2)= i_imag*sqrt(3.0d0) !
 MSOCd(8,7)= i_imag !
 MSOCd(1,7)= 1.0d0 !
 MSOCd(4,7)= i_imag !
 MSOCd(5,7)= i_imag*sqrt(3.0d0) !
 MSOCd(4,1)= -2.0d0*i_imag !
 MSOCd(8,1)= -i_imag !
 MSOCd(7,1)= 1.0d0 !
 MSOCd(6,9)= -2.0d0*i_imag !
 MSOCd(3,9)= -1.0d0 !
 MSOCd(2,9)= -i_imag !
 MSOCd(9,6)= 2.0d0*i_imag !
 MSOCd(3,6)= i_imag !
 MSOCd(2,6)= -1.0d0 !
 MSOCd(1,4)= 2.0d0*i_imag
 MSOCd(8,4)= 1.0d0 !
 MSOCd(7,4)= -i_imag !
 print*,"d'apres calcul perso :" 
 print*,'Matrice MSOC dans la base dxy+, dyz+, dzx+, dx2-y2+, d3z2-r2+ dxy-, dyz-, dzx-, dx2-y2-, d3z2-r2- :'
 do i=1,10
  !print*,MSOCp(i,:) 
  print "(10('(',(F7.2),',',(F7.3),')',1X))", real(MSOCd(i,1)),imag(MSOCd(i,1)), &
     real(MSOCd(i,2)),imag(MSOCd(i,2)),real(MSOCd(i,3)),imag(MSOCd(i,3)),real(MSOCd(i,4)),imag(MSOCd(i,4)), &
     real(MSOCd(i,5)),imag(MSOCd(i,5)),real(MSOCd(i,6)),imag(MSOCd(i,6)),real(MSOCd(i,7)),imag(MSOCd(i,7)), &
     real(MSOCd(i,8)),imag(MSOCd(i,8)),real(MSOCd(i,9)),imag(MSOCd(i,9)),real(MSOCd(i,10)),imag(MSOCd(i,10))
  do j=1,i  !verification hermiticite de MSOCd
    if (MSOCd(i,j).ne.conjg(MSOCd(j,i))) stop 'matrice couplage spin-orbit d non hermicienne'
  enddo
 enddo
 print*,'H_SOC = 0.5 lambda_SOC MSOCp, avec :'
 print*,'  lambda (SOC) pour atome de type 1 orbitales d = ',lambda_SOC(2,1)
 print*,'  lambda (SOC) pour atome de type 2 orbitales d = ',lambda_SOC(2,2)
else
 print*,''
 print*,'!!!! nspin (variable.f90) doit etre egal a 1 (calcul sans spin)', & 
        'ou 2 (calcule avec spin et spin-orbite).'
 stop
endif


! Ecriture du tableau g avec terme de saut atome i et ses entre voisins
g = 0.0d0
gvx = 0.0d0

! terme non diagonaux
do i=1,natomes
   i_type = type_atome(i)
   do j=1,nbvoisins(i)
      i2 = voisins(i,j,0)  ! numero du jeme voisin de i
      i2_type = type_atome(i2)
      R_i2(:) = vec_Rat(i2,:)  ! position dans la maille 
      do ndim_res = 1,dimension_reseau
        R_i2(:) = R_i2(:) + a(ndim_res,:)*voisins(i,j,ndim_res) ! position reelle
      enddo 
      do orb1=1,orbitales_lm(i_type,0)
         lm_i = orbitales_lm(i_type,orb1)
         do orb2=1,orbitales_lm(i2_type,0)
            lm_j = orbitales_lm(i2_type,orb2)
            !dmin = 0.577350 ! distance 
	    !dmin = 0.0d0
     !!15/05/2008       g(i,j,orb1,orb2) = + SlaterKoster(lm_i,lm_j, &
     !!            vec_Rat(i,1),vec_Rat(i,2),vec_Rat(i,3), &
     !!            R_i2(1),R_i2(2),R_i2(3),  &
     !!            SlaterKoster_par(i_type,i2_type,:)) * & 
!2janv18    g(i,j,orb1,orb2) = + SlaterKoster_AttenD(lm_i,lm_j, &
            g(i,j,orb1,orb2,1,1) = + SlaterKoster_AttenD(lm_i,lm_j, &
                 vec_Rat(i,1),vec_Rat(i,2),vec_Rat(i,3), &
                 R_i2(1),R_i2(2),R_i2(3),  &
                 dunite,SlaterKoster_par(i_type,i2_type,:,:) ) * & 
  !! fin modif 15/5/2008
                 exp(i_imag*coef_phase_g*(R_i2(2)-vec_Rat(i,2))*(R_i2(1)+vec_Rat(i,1))) ! effet de B
            ! 25/03/2008:
!2janv18    gvx(i,j,orb1,orb2) = g(i,j,orb1,orb2) * &
            gvx(i,j,orb1,orb2,1,1) = g(i,j,orb1,orb2,1,1) * &
              DOT_PRODUCT((R_i2(1:dimension_reseau)-vec_Rat(i,1:dimension_reseau)),direction_transp)
	   ! print*,'i,j,i2,R_i,R_i2,=',i,j,i2,vec_Rat(i,:),R_i2(:)
           ! print*,'    g,gvx=',g(i,j,orb1,orb2),gvx(i,j,orb1,orb2)
	   ! print*,'      DXikj=',DOT_PRODUCT((R_i2(1:dimension_reseau)-vec_Rat(i,1:dimension_reseau)),direction_transp)
         enddo
      enddo
   enddo
enddo

! termes diagonaux:
do i=1, natomes
   i_type = type_atome(i)
   do orb1=1,orbitales_lm(i_type,0)
!2janv18      g(i,0,orb1,orb1) = energie_at_type(type_atome(i),orb1)
              g(i,0,orb1,orb1,1,1) = energie_at_type(type_atome(i),orb1)
   enddo
enddo

! 2janv18  spin 1 et spin 2 identiques et non couples (pas encore de couplage spin-orbite)
if (nspin.eq.2) then  
  g(:,:,:,:,nspin,nspin) = g(:,:,:,:,1,1)
  gvx(:,:,:,:,nspin,nspin) = gvx(:,:,:,:,1,1)
endif

! 2janv18: couplage spin-orbite uniquement entre orbitale p et d d'un même atome 
if (nspin.eq.2) then
do i=1, natomes
   i_type = type_atome(i)
   do orb1=1, orbitales_lm(i_type,0)
      lm_1 = orbitales_lm(i_type,orb1)
      do orb2=1, orbitales_lm(i_type,0)
	  lm_2 = orbitales_lm(i_type,orb2)
          !print*,'i,orb1,orb2 = ',i,orb1,orb2
	  if ((lm_1.eq.2.or.lm_1.eq.3.or.lm_1.eq.4).and.(lm_2.eq.2.or.lm_2.eq.3.or.lm_2.eq.4) &  ! orbital p
               .and.(lm_1.ne.lm_2)) then
	    g(i,0,orb1,orb2,1,1) = MSOCp(lm_1-1,lm_2-1) * lambda_SOC(1,i_type) / 2.0d0
            g(i,0,orb1,orb2,1,nspin) = MSOCp(lm_1-1,lm_2+2) * lambda_SOC(1,i_type) / 2.0d0
            g(i,0,orb1,orb2,nspin,1) = MSOCp(lm_1+2,lm_2-1) * lambda_SOC(1,i_type) / 2.0d0
            g(i,0,orb1,orb2,nspin,nspin) = MSOCp(lm_1+2,lm_2+2) * lambda_SOC(1,i_type) / 2.0d0
	  else if (lm_1.ge.5.and.lm_1.le.9.and.lm_2.ge.5.and.lm_2.le.9 &                   ! orbital d
               .and.lm_1.ne.lm_2) then
	    g(i,0,orb1,orb2,1,1) = MSOCd(lm_1-4,lm_2-4) * lambda_SOC(2,i_type) / 2.0d0
            g(i,0,orb1,orb2,1,nspin) = MSOCd(lm_1-4,lm_2+1) * lambda_SOC(2,i_type) / 2.0d0
            g(i,0,orb1,orb2,nspin,1) = MSOCd(lm_1+1,lm_2-4) * lambda_SOC(2,i_type) / 2.0d0
            g(i,0,orb1,orb2,nspin,nspin) = MSOCd(lm_1+1,lm_2+1) * lambda_SOC(2,i_type) / 2.0d0
          endif
      enddo
   enddo
enddo
endif 

gvx = gvx / i_imag

! Ecriture g
!do orb1=1,3
!do orb2=1,3
! print*,'g pour orb1, orb2 =',orb1,orb2
! print*,g(1,0,orb1,orb2,:,:)
!enddo
!enddo
!stop
!print*,'gvx = ',gvx



! Ecriture du tableau indice_orb 
! indice_orb(i,orb1) : indice dans la maille de 1 a ndim 
!                             de l'orbitale numero orb1 de l'atome i
ndim = 0
do i=1,natomes
   i_type = type_atome(i)
!   print*,'i,orbitales_lm(i_type,0) =',i,orbitales_lm(i_type,0)
   do orb1=1,orbitales_lm(i_type,0)
!      print*,'i,orb1 =',i,orb1
      ndim = ndim + 1
      indice_orb(i,orb1) = ndim
   enddo
enddo
print*,''
print "('Ecriture du tableau indice_orb, ndim = ',I5)",ndim
!print*,indice_orb
if (ndim>ndim_max) stop '!!!! Sub: Hamil_entree : ndim > ndim_max'

! estimation rapide des energies de bords de bandes:
ee_min = SUM(energie_at_type) / nb_types_at - &
!ee_min = SUM(energie_at_type(1:nb_types_at)) / nb_types_at - &
         MAXVAL(abs(SlaterKoster_par))*MAXVAL(nbvoisins(1:natomes))
ee_max = SUM(energie_at_type) / nb_types_at + &
!ee_max = SUM(energie_at_type(1:nb_types_at)) / nb_types_at + &
         MAXVAL(abs(SlaterKoster_par))*MAXVAL(nbvoisins(1:natomes))


end SUBROUTINE Hamil_Vx_entree

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Lec_SlaterKoster_par(ntypes_atomes_max, ntypes_atomes, &
                                nom_SlaterKoster,SlaterKoster_par,champB)
!
!    23/05/2008: dernier colonne de SlaterKoster_par(....., 1:3)
!         SlaterKoster_par(....., 1) = Gamma = V
!         SlaterKoster_par(....., 2) = q 
!         SlaterKoster_par(....., 1) = dmin 
!
!    Lecture du fichie INPUT_SlaterKoster 
!    contenant les parametres utile pour Sub: SlaterKoster
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE precision_reel
  
  implicit none

  integer ntypes_atomes_max,ntypes_atomes  ! entree
  character*30 nom_SlaterKoster            ! entree

  real(kind=prec),dimension(1:ntypes_atomes_max,1:ntypes_atomes_max,0:14,1:4) :: &  !avril17
       SlaterKoster_par  ! sortie
  real(kind=prec) champB

  integer i,j,n1,n2
  character nom*30


  SlaterKoster_par = 0.0d0

  open(18,file=nom_SlaterKoster,action='read')

  read(18,*) champB

  read(18,*) i
  if (i<ntypes_atomes) STOP &
       '!!!!! INPUT Slater-Koster : nb type atomes incompatible !'

  do i=1,(ntypes_atomes*(ntypes_atomes+1)/2)
     read(18,*) nom
     read(18,*) n1,n2
     read(18,*) nom ! pour sauter la ligne avec distance max des voisins
     !read(18,*) SlaterKoster_par(n1,n2,0,1)!,SlaterKoster_par(n1,n2,0,2)
     read(18,*) nom
     !SlaterKoster_par(n2,n1,0,1) = SlaterKoster_par(n1,n2,0,1)
     !SlaterKoster_par(n2,n1,0,2) = SlaterKoster_par(n1,n2,0,2)
     do j=1,14
        read(18,*) SlaterKoster_par(n1,n2,j,1:4)
     enddo
     if (n1.ne.n2) then                              ! janv17                             
        SlaterKoster_par(n2,n1,1:14,1:4) = SlaterKoster_par(n1,n2,1:14,1:4)
        SlaterKoster_par(n2,n1,2,1:4) = SlaterKoster_par(n1,n2,3,1:4)  ! hermiticite il faut Vs(n2)-p(n1) = Vp(n1)-s(n2)
        SlaterKoster_par(n2,n1,3,1:4) = SlaterKoster_par(n1,n2,2,1:4)  !                 et  Vp(n2)-s(n1) = Vs(n1)-p(n2)
        SlaterKoster_par(n2,n1,6,1:4) = SlaterKoster_par(n1,n2,7,1:4)  !
        SlaterKoster_par(n2,n1,7,1:4) = SlaterKoster_par(n1,n2,6,1:4)  !                 ...
        SlaterKoster_par(n2,n1,8,1:4) = SlaterKoster_par(n1,n2,9,1:4)  !
        SlaterKoster_par(n2,n1,9,1:4) = SlaterKoster_par(n1,n2,8,1:4)  !
        SlaterKoster_par(n2,n1,10,1:4) = SlaterKoster_par(n1,n2,11,1:4)!
        SlaterKoster_par(n2,n1,11,1:4) = SlaterKoster_par(n1,n2,10,1:4)!
     endif
  enddo

  print*,''
  print*,'Sub: Lec_SlaterKoster_par : Parmetres de Slater-Koster (AttenD version janvier 2018):'
  print*,'     q et dmin dependent de chaque parametre de Slater-Koster'
  do n1=1,ntypes_atomes
     do n2=1,ntypes_atomes
        print "('  entre atome de type ,'I3,' et ',I3)",n1,n2
        print "('    Distance dmin (si <>0 : pour calcul de terme de saut) =',F10.5,' q = ',F10.5)", &
               SlaterKoster_par(n1,n2,0,1)!,SlaterKoster_par(n1,n2,0,2)
        print*,'Pour chaque parametre de S-K: Valeur de V  ,   q     , dmin ,  Rcut:'
        print "('    Vss sigma = ',4F10.5)",SlaterKoster_par(n1,n2,1,:)
        print "('    Vsp sigma = ',4F10.5)",SlaterKoster_par(n1,n2,2,:)
	print "('    Vps sigma = ',4F10.5)",SlaterKoster_par(n1,n2,3,:)
	if (n1.eq.n2.and.(SlaterKoster_par(n1,n2,2,1).ne.SlaterKoster_par(n1,n2,3,1) &
			  .or.SlaterKoster_par(n1,n2,2,2).ne.SlaterKoster_par(n1,n2,3,2) &
			  .or.SlaterKoster_par(n1,n2,2,3).ne.SlaterKoster_par(n1,n2,3,3) &
			  .or.SlaterKoster_par(n1,n2,2,4).ne.SlaterKoster_par(n1,n2,3,4))) stop &
           "STOP: couplage entre atomes identiques Vsp_sigma et Vps_sigma doivent etre egaux !"
        print "('    Vpp sigma = ',4F10.5)",SlaterKoster_par(n1,n2,4,:)
        print "('    Vpp pi =    ',4F10.5)",SlaterKoster_par(n1,n2,5,:)
        print "('    Vsd sigma = ',4F10.5)",SlaterKoster_par(n1,n2,6,:)
        print "('    Vds sigma = ',4F10.5)",SlaterKoster_par(n1,n2,7,:)
	if (n1.eq.n2.and.(SlaterKoster_par(n1,n2,6,1).ne.SlaterKoster_par(n1,n2,7,1) &
			  .or.SlaterKoster_par(n1,n2,6,2).ne.SlaterKoster_par(n1,n2,7,2) &
			  .or.SlaterKoster_par(n1,n2,6,3).ne.SlaterKoster_par(n1,n2,7,3) &
			  .or.SlaterKoster_par(n1,n2,6,4).ne.SlaterKoster_par(n1,n2,7,4))) stop &
           "STOP: couplage entre atomes identiques Vsd_sigma et Vds_sigma doivent etre egaux !"
        print "('    Vpd sigma = ',4F10.5)",SlaterKoster_par(n1,n2,8,:)
        print "('    Vdp sigma = ',4F10.5)",SlaterKoster_par(n1,n2,9,:)
	if (n1.eq.n2.and.(SlaterKoster_par(n1,n2,8,1).ne.SlaterKoster_par(n1,n2,9,1) &
			  .or.SlaterKoster_par(n1,n2,8,2).ne.SlaterKoster_par(n1,n2,9,2) &
			  .or.SlaterKoster_par(n1,n2,8,3).ne.SlaterKoster_par(n1,n2,9,3) &
			  .or.SlaterKoster_par(n1,n2,8,4).ne.SlaterKoster_par(n1,n2,9,4))) stop &
           "STOP: couplage entre atomes identiques Vpd_sigma et Vdp_sigma doivent etre egaux !"
        print "('    Vpd pi =    ',4F10.5)",SlaterKoster_par(n1,n2,10,:)
        print "('    Vdp pi =    ',4F10.5)",SlaterKoster_par(n1,n2,11,:)
	if (n1.eq.n2.and.(SlaterKoster_par(n1,n2,10,1).ne.SlaterKoster_par(n1,n2,11,1) &
			  .or.SlaterKoster_par(n1,n2,10,2).ne.SlaterKoster_par(n1,n2,11,2) &
			  .or.SlaterKoster_par(n1,n2,10,3).ne.SlaterKoster_par(n1,n2,11,3) &
			  .or.SlaterKoster_par(n1,n2,10,4).ne.SlaterKoster_par(n1,n2,11,4))) stop &
           "STOP: couplage entre atomes identiques Vpd_pi et Vdp_pi doivent etre egaux !"
        print "('    Vdd sigma = ',4F10.5)",SlaterKoster_par(n1,n2,12,:)
        print "('    Vdd pi =    ',4F10.5)",SlaterKoster_par(n1,n2,13,:)
        print "('    Vddd =      ',4F10.5)",SlaterKoster_par(n1,n2,14,:)
     enddo
  enddo




end SUBROUTINE Lec_SlaterKoster_par

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Intergrales de sauts de Slater & Koster (1954)
!
!     aout 2009 : AttenD
!
!     23/05/2008 : pour AttenC
!
!     15/05/2008 : Avec Attenuation du 19/05/2008
!                  pour orbital p uniquement
!
!     2005/6/8: pour s-s  s-p  s-d  p-p  p-d  (et symetriques)
!     
!     depend de la distance
!     
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      real(kind(1.d0)) function  &
           SlaterKoster_AttenD(orb1i,orb2i,x1,y1,z1,x2,y2,z2,dunite,v)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     orb* = 1 s,   
!            2 px,  3 py,  4 pz, 
!            5 dxy, 6 dyz, 7 dzx, 8 dx2-y2, 9 d3z2-r2 
!
!     orb1 : orbitale de l'atome en (x1,y1,z1) en dunite
!     orb2 : orbitale de l'atome en (x2,y2,z2) en dunite
!
!     dunite : unite des distances donnees dans x* y* et z* : en bohr (u.a.)
!

      use precision_reel 
      implicit none
!c     Variables d'entree
      integer orb1,orb2,orb1i,orb2i,l_orb1i,l_orb2i
      real(kind=prec) x1,x2,y1,y2,z1,z2,dunite

!c
      real(kind=prec),parameter :: ry_ev = 13.6058d0,ua = 0.529177d0
      real(kind=prec) Vsss,Vsps,Vpps,Vppp,Vsds,Vsdp,Vpds,Vpdp
      real(kind=prec) Vdds,Vddp,Vddd,SlaterKoster_param
      real(kind=prec) Vpss,Vdss,Vdps,Vdpp                 ! avril17

      real(kind=prec) l,m,n,d,q

      !real(kind=prec) apps,bpps,cpps,dpps,appp,bppp,cppp,dppp
      real(kind=prec) dpps,dppp,qpps,qppp,gpps,gppp
      real(kind=prec) r0_cut,l_cut,f_cut

      real(kind=prec),dimension(0:14,1:4) :: v   !avril17
      real(kind=prec) calcul_effet_dis           !8fev18

!	DONNEES :
!       p-p sigma :
!!        gpps = -0.4d0     ! couplage inter-plan (1er voisin)
        !qpps = 9.44625d0 ! pour  Atten
!!        qpps = 7.42807d0 ! pour  Atten2
!!	dpps = 1.36345d0 ! = dmin_pps = distance inter-plan

!	p-p pi :
!!        gppp = 3.0d0     ! couplage 1er voisin dans le plan
        !qppp = 4.0d0   ! pour Atten
!!        qppp = 3.1454d0   ! pour Atten2
!!	dppp = 0.57735d0 ! = dmin_ppp = distance 1er voisin dans le plan


!c      print*,'x1,y1,z1',x1,y1,z1
!c      print*,'x2,y2,z2',x2,y2,z2

      d = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1))
      if(d.lt.zero) then
      !stop '!!!! fonction: SlaterKoster: Atomes trop proche !!!'
       print*,'!!!! fonction: SlaterKoster: Atomes trop proche !!!'
      endif

      l = (x2-x1)/d
      m = (y2-y1)/d
      n = (z2-z1)/d


      d = d !* dunite  ! pour avoir d en u.a.

!!!     Function de cup of universel : d'apres Mehl PRB 50, R14695 (1998). ??? je ne la trouve pas !!
!!!     Voir aussi: Michael J. Mehl and Dimitrios A. Papaconstantopoulos, PRB 54 (1996) 4519-4530
!!!     Rmq.: pour la DOS cette fonction ne change presque rien
!!      r0_cut = 16.5d0 ! en u.a.
        ! Mehl96 : l_cut = 0.5 u.a.
        l_cut = 0.5 * ua / dunite ! pour avoir l_cut dans le meme unite de v et d (unite de a) 
                                  ! avec  dunite = a en Ang !!!!
!!      f_cut = 1.0d0 / (1 + exp((d-r0_cut)/l_cut))  

      !print*,'dunite, l_cut, v(1,4)=',dunite, l_cut, v(1,4)

      !Vpps = ( apps + bpps*d + cpps*d*d ) * exp(-dpps*dpps*d) *f_cut
      !Vppp = ( appp + bppp*d + cppp*d*d ) * exp(-dppp*dppp*d) *f_cut
!!       Vpps = gpps * exp(qpps*(1-d/dpps))  !!* f_cut
!!       Vppp = gppp * exp(qppp*(1-d/dppp))  !!* f_cut

      ! 8fev18
      !Vsss = v(1,1) * exp(v(1,2)*(1-d/v(1,3))) *1.0/(1.0+exp((d-v(1,4))/l_cut))
      !Vsps = v(2,1) * exp(v(2,2)*(1-d/v(2,3))) *1.0/(1.0+exp((d-v(2,4))/l_cut))
      !Vpss = v(3,1) * exp(v(3,2)*(1-d/v(3,3))) *1.0/(1.0+exp((d-v(3,4))/l_cut))  ! avril17
      !Vpps = v(4,1) * exp(v(4,2)*(1-d/v(4,3))) *1.0/(1.0+exp((d-v(4,4))/l_cut))
      !Vppp = v(5,1) * exp(v(5,2)*(1-d/v(5,3))) *1.0/(1.0+exp((d-v(5,4))/l_cut))
      !Vsds = v(6,1) * exp(v(6,2)*(1-d/v(6,3))) *1.0/(1.0+exp((d-v(6,4))/l_cut))
      !Vdss = v(7,1) * exp(v(7,2)*(1-d/v(7,3))) *1.0/(1.0+exp((d-v(7,4))/l_cut)) ! avril17
      !Vpds = v(8,1) * exp(v(8,2)*(1-d/v(8,3))) *1.0/(1.0+exp((d-v(8,4))/l_cut))
      !Vdps = v(9,1) * exp(v(9,2)*(1-d/v(9,3))) *1.0/(1.0+exp((d-v(9,4))/l_cut)) ! avril17
      !Vpdp = v(10,1) * exp(v(10,2)*(1-d/v(10,3))) *1.0/(1.0+exp((d-v(10,4))/l_cut))
      !Vdpp = v(11,1) * exp(v(11,2)*(1-d/v(11,3))) *1.0/(1.0+exp((d-v(11,4))/l_cut)) ! avril17
      !Vdds = v(12,1) * exp(v(12,2)*(1-d/v(12,3))) *1.0/(1.0+exp((d-v(12,4))/l_cut))
      !Vddp = v(13,1) * exp(v(13,2)*(1-d/v(13,3))) *1.0/(1.0+exp((d-v(13,4))/l_cut))
      !Vddd = v(14,1) * exp(v(14,2)*(1-d/v(14,3))) *1.0/(1.0+exp((d-v(14,4))/l_cut))
      Vsss = calcul_effet_dis(v(1,1),v(1,2),v(1,3),d) *1.0/(1.0+exp((d-v(1,4))/l_cut))
      Vsps = calcul_effet_dis(v(2,1),v(2,2),v(2,3),d) *1.0/(1.0+exp((d-v(1,4))/l_cut))
      Vpss = calcul_effet_dis(v(3,1),v(3,2),v(3,3),d) *1.0/(1.0+exp((d-v(1,4))/l_cut))
      Vpps = calcul_effet_dis(v(4,1),v(4,2),v(4,3),d) *1.0/(1.0+exp((d-v(1,4))/l_cut))
      Vppp = calcul_effet_dis(v(5,1),v(5,2),v(5,3),d) *1.0/(1.0+exp((d-v(1,4))/l_cut))
      Vsds = calcul_effet_dis(v(6,1),v(6,2),v(6,3),d) *1.0/(1.0+exp((d-v(1,4))/l_cut))
      Vdss = calcul_effet_dis(v(7,1),v(7,2),v(7,3),d) *1.0/(1.0+exp((d-v(1,4))/l_cut))
      Vpds = calcul_effet_dis(v(8,1),v(8,2),v(8,3),d) *1.0/(1.0+exp((d-v(1,4))/l_cut))
      Vdps = calcul_effet_dis(v(9,1),v(9,2),v(9,3),d) *1.0/(1.0+exp((d-v(1,4))/l_cut))
      Vpdp = calcul_effet_dis(v(10,1),v(10,2),v(10,3),d) *1.0/(1.0+exp((d-v(1,4))/l_cut))
      Vdpp = calcul_effet_dis(v(11,1),v(11,2),v(11,3),d) *1.0/(1.0+exp((d-v(1,4))/l_cut))
      Vdds = calcul_effet_dis(v(12,1),v(12,2),v(12,3),d) *1.0/(1.0+exp((d-v(1,4))/l_cut))
      Vddp = calcul_effet_dis(v(13,1),v(13,2),v(13,3),d) *1.0/(1.0+exp((d-v(1,4))/l_cut))
      Vddd = calcul_effet_dis(v(14,1),v(14,2),v(14,3),d) *1.0/(1.0+exp((d-v(1,4))/l_cut))



!avril17 corrige en janv18
!      ! Les calculs sont fait uniquement pour orb1 <= orb2, sinon "inversion" à la fin du calcul
!      if (orb1i.le.orb2i) then 
	 orb1 = orb1i
	 orb2 = orb2i
!      else 
!	 orb1 = orb2i
!	 orb2 = orb1i
!      endif

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     s-s
      if (orb1.eq.1.and.orb2.eq.1) then
         SlaterKoster_param = Vsss
         goto 100
      endif

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     s-p
      if (orb1.eq.1.and.orb2.eq.2) then
         SlaterKoster_param = l*Vsps
         goto 100
      endif
      if (orb1.eq.2.and.orb2.eq.1) then    ! avril17
         SlaterKoster_param = l*Vpss
         goto 100
      endif
      if (orb1.eq.1.and.orb2.eq.3) then
         SlaterKoster_param = m*Vsps
         goto 100
      endif
      if (orb1.eq.3.and.orb2.eq.1) then   ! avril17
         SlaterKoster_param = m*Vpss       
         goto 100
      endif
      if (orb1.eq.1.and.orb2.eq.4) then
         SlaterKoster_param = n*Vsps
         goto 100
      endif
      if (orb1.eq.4.and.orb2.eq.1) then     ! avril17
         SlaterKoster_param = n*Vpss
         goto 100
      endif

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     p-p
      if (orb1.eq.2.and.orb2.eq.2) then                  ! px - px 
         SlaterKoster_param = l*l*Vpps + (1.0 - l*l)*Vppp
         goto 100
      endif
      if (orb1.eq.3.and.orb2.eq.3) then                  ! py - py 
         SlaterKoster_param = m*m*Vpps + (1.0 - m*m)*Vppp
         goto 100
      endif
      if (orb1.eq.4.and.orb2.eq.4) then                  ! pz - pz 
         SlaterKoster_param = n*n*Vpps + (1.0 - n*n)*Vppp
         goto 100
      endif

      if ((orb1.eq.2.and.orb2.eq.3).or.(orb1.eq.3.and.orb2.eq.2)) then ! px - py  ! avril17
      !if (orb1.eq.2.and.orb2.eq.3) then
         SlaterKoster_param = l*m*(Vpps - Vppp)
         goto 100
      endif
      if ((orb1.eq.2.and.orb2.eq.4).or.(orb1.eq.4.and.orb2.eq.2)) then ! px - pz    ! avril17
      !if (orb1.eq.2.and.orb2.eq.4) then
         SlaterKoster_param = l*n*(Vpps - Vppp)
         goto 100
      endif
      if ((orb1.eq.3.and.orb2.eq.4).or.(orb1.eq.4.and.orb2.eq.3)) then  ! py - pz   ! avril17
      !if (orb1.eq.3.and.orb2.eq.4) then
         SlaterKoster_param = m*n*(Vpps - Vppp)
         goto 100
      endif

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     s-d
      if (orb1.eq.1.and.orb2.eq.5) then
         SlaterKoster_param = l*m*Vsds*sqrt(3.0)
         goto 100
      endif
      if (orb1.eq.5.and.orb2.eq.1) then                ! avril17
         SlaterKoster_param = l*m*Vdss*sqrt(3.0)  
         goto 100
      endif
      if (orb1.eq.1.and.orb2.eq.6) then
         SlaterKoster_param = m*n*Vsds*sqrt(3.0)
         goto 100
      endif
      if (orb1.eq.6.and.orb2.eq.1) then               ! avril17
         SlaterKoster_param = m*n*Vdss*sqrt(3.0)
         goto 100
      endif
      if (orb1.eq.1.and.orb2.eq.7) then
         SlaterKoster_param = n*l*Vsds*sqrt(3.0)
         goto 100
      endif
      if (orb1.eq.7.and.orb2.eq.1) then              ! avril17
         SlaterKoster_param = n*l*Vdss*sqrt(3.0)
         goto 100
      endif
      if (orb1.eq.1.and.orb2.eq.8) then
         SlaterKoster_param = (l*l - m*m)*Vsds*0.5*sqrt(3.0)
         goto 100
      endif
      if (orb1.eq.8.and.orb2.eq.1) then                        ! avril17
         SlaterKoster_param = (l*l - m*m)*Vdss*0.5*sqrt(3.0)
         goto 100
      endif
      if (orb1.eq.1.and.orb2.eq.9) then
         SlaterKoster_param = ( n*n - 0.5*(l*l + m*m) )*Vsds
         goto 100
      endif
      if (orb1.eq.9.and.orb2.eq.1) then                        ! avril17
         SlaterKoster_param = ( n*n - 0.5*(l*l + m*m) )*Vdss
         goto 100
      endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     p-d
       if (orb1.eq.2.and.orb2.eq.5) then  ! px - dxy
         SlaterKoster_param = l*l*m*Vpds*sqrt(3.0) + m*(1-2.0*l*l)*Vpdp
         goto 100
      endif
       if (orb1.eq.5.and.orb2.eq.2) then  ! dxy - px                     ! avril17
         SlaterKoster_param = l*l*m*Vdps*sqrt(3.0) + m*(1-2.0*l*l)*Vdpp
         goto 100
      endif
      if (orb1.eq.3.and.orb2.eq.6) then  ! py - dyz
         SlaterKoster_param = m*m*n*Vpds*sqrt(3.0) + n*(1-2.0*m*m)*Vpdp
         goto 100
      endif
      if (orb1.eq.6.and.orb2.eq.3) then  ! dyz - py                     ! avril17
         SlaterKoster_param = m*m*n*Vdps*sqrt(3.0) + n*(1-2.0*m*m)*Vdpp
         goto 100
      endif
      if (orb1.eq.4.and.orb2.eq.7) then  ! pz - dzx
         SlaterKoster_param = n*n*l*Vpds*sqrt(3.0) + l*(1-2.0*n*n)*Vpdp
         goto 100
      endif
      if (orb1.eq.7.and.orb2.eq.4) then  ! dzx - pz                     ! avril17
         SlaterKoster_param = n*n*l*Vdps*sqrt(3.0) + l*(1-2.0*n*n)*Vdpp
         goto 100
      endif
     
      if (orb1.eq.2.and.orb2.eq.6) then  ! px - dyz
         SlaterKoster_param = l*n*m*Vpds*sqrt(3.0) - 2.0*l*m*n*Vpdp
         goto 100
      endif
      if (orb1.eq.6.and.orb2.eq.2) then  ! dyz - px                    ! avril17
         SlaterKoster_param = l*n*m*Vdps*sqrt(3.0) - 2.0*l*m*n*Vdpp
         goto 100
      endif
      if (orb1.eq.3.and.orb2.eq.7) then  ! py - dzx
         SlaterKoster_param = l*n*m*Vpds*sqrt(3.0) - 2.0*l*m*n*Vpdp
         goto 100
      endif
      if (orb1.eq.7.and.orb2.eq.3) then  ! dzx - py                    ! avril17
         SlaterKoster_param = l*n*m*Vdps*sqrt(3.0) - 2.0*l*m*n*Vdpp
         goto 100
      endif
      if (orb1.eq.4.and.orb2.eq.5) then  ! pz - dxy
         SlaterKoster_param = l*n*m*Vpds*sqrt(3.0) - 2.0*l*m*n*Vpdp
         goto 100
      endif
      if (orb1.eq.5.and.orb2.eq.4) then  ! dxy - pz                    ! avril17
         SlaterKoster_param = l*n*m*Vdps*sqrt(3.0) - 2.0*l*m*n*Vdpp
         goto 100
      endif

      if (orb1.eq.2.and.orb2.eq.7) then  ! px - dzx
         SlaterKoster_param = l*l*n*Vpds*sqrt(3.0) + n*(1-2.0*l*l)*Vpdp
         goto 100
      endif
      if (orb1.eq.7.and.orb2.eq.2) then  ! dzx - px                         ! avril17
         SlaterKoster_param = l*l*n*Vdps*sqrt(3.0) + n*(1-2.0*l*l)*Vdpp
         goto 100
      endif
      if (orb1.eq.3.and.orb2.eq.5) then  ! py - dxy
         SlaterKoster_param = m*m*l*Vpds*sqrt(3.0) + l*(1-2.0*m*m)*Vpdp
         goto 100
      endif
      if (orb1.eq.5.and.orb2.eq.3) then  ! dxy - py                         ! avril17
         SlaterKoster_param = m*m*l*Vdps*sqrt(3.0) + l*(1-2.0*m*m)*Vdpp
         goto 100
      endif
      if (orb1.eq.4.and.orb2.eq.6) then  ! pz - dyz
         SlaterKoster_param = n*n*m*Vpds*sqrt(3.0) + m*(1-2.0*n*n)*Vpdp
         goto 100
      endif
      if (orb1.eq.6.and.orb2.eq.4) then  ! dyz - pz                        ! avril17
         SlaterKoster_param = n*n*m*Vdps*sqrt(3.0) + m*(1-2.0*n*n)*Vdpp
         goto 100
      endif
    
      if (orb1.eq.2.and.orb2.eq.8) then  ! px - dx2-y2
         SlaterKoster_param = l*(l*l-m*m)*Vpds*0.5*sqrt(3.0) + l*(1.0-l*l+m*m)*Vpdp
         goto 100
      endif
      if (orb1.eq.8.and.orb2.eq.2) then  ! dx2-y2 - px                            ! avril17
         SlaterKoster_param = l*(l*l-m*m)*Vdps*0.5*sqrt(3.0) + l*(1.0-l*l+m*m)*Vdpp
         goto 100
      endif
      if (orb1.eq.3.and.orb2.eq.8) then  ! py - dx2-y2
         SlaterKoster_param = m*(l*l-m*m)*Vpds*0.5*sqrt(3.0) - m*(1.0+l*l-m*m)*Vpdp  ! corrige le 32dec17
         goto 100
      endif
      if (orb1.eq.8.and.orb2.eq.3) then  ! dx2-y2 - py                            ! avril17
         SlaterKoster_param = m*(l*l-m*m)*Vdps*0.5*sqrt(3.0) - m*(1.0+l*l-m*m)*Vdpp   ! corrige le 32dec17
         goto 100
      endif
      if (orb1.eq.4.and.orb2.eq.8) then  ! pz - dx2-y2
         SlaterKoster_param = n*(l*l-m*m)*Vpds*0.5*sqrt(3.0) - n*(l*l-m*m)*Vpdp    ! corrige le 32dec17
         goto 100
      endif
      if (orb1.eq.8.and.orb2.eq.4) then  ! dx2-y2 - pz                            ! avril17
         SlaterKoster_param = n*(l*l-m*m)*Vdps*0.5*sqrt(3.0) - n*(l*l-m*m)*Vdpp    ! corrige le 32dec17
         goto 100
      endif

      if (orb1.eq.2.and.orb2.eq.9) then  ! px - d3z2-r2
         SlaterKoster_param = l*(n*n - 0.5*(l*l+m*m))*Vpds - l*n*n*Vpdp*sqrt(3.0)
         goto 100
      endif
      if (orb1.eq.9.and.orb2.eq.2) then  ! d3z2-r2 - px                           ! avril17
         SlaterKoster_param = l*(n*n - 0.5*(l*l+m*m))*Vdps - l*n*n*Vdpp*sqrt(3.0)
         goto 100
      endif
      if (orb1.eq.3.and.orb2.eq.9) then  ! py - d3z2-r2
         SlaterKoster_param = m*(n*n - 0.5*(l*l+m*m))*Vpds - m*n*n*Vpdp*sqrt(3.0)
         goto 100
      endif
      if (orb1.eq.9.and.orb2.eq.3) then  ! d3z2-r2 - py                           ! avril17        
         SlaterKoster_param = m*(n*n - 0.5*(l*l+m*m))*Vdps - m*n*n*Vdpp*sqrt(3.0)
         goto 100
      endif
      if (orb1.eq.4.and.orb2.eq.9) then  ! pz - d3z2-r2
         SlaterKoster_param = n*(n*n - 0.5*(l*l+m*m))*Vpds + n*(l*l+m*m)*Vpdp*sqrt(3.0)
         goto 100
      endif
      if (orb1.eq.9.and.orb2.eq.4) then  ! d3z2-r2  - pz                         ! avril17
         SlaterKoster_param = n*(n*n - 0.5*(l*l+m*m))*Vdps + n*(l*l+m*m)*Vdpp*sqrt(3.0)
         goto 100
      endif

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     d-d
       if (orb1.eq.5.and.orb2.eq.5) then  ! dxy - dxy
         SlaterKoster_param = 3.0*l*l*m*m*Vdds + (l*l+m*m-4.0*l*l*m*m)*Vddp +(n*n+l*l*m*m)*Vddd
         goto 100
       endif
       if (orb1.eq.6.and.orb2.eq.6) then  ! dyz - dyz
         SlaterKoster_param = 3.0*m*m*n*n*Vdds + (m*m+n*n-4.0*m*m*n*n)*Vddp +(l*l+m*m*n*n)*Vddd
         goto 100
       endif
       if (orb1.eq.7.and.orb2.eq.7) then  ! dzx - dzx   
         SlaterKoster_param = 3.0*n*n*l*l*Vdds + (n*n+l*l-4.0*n*n*l*l)*Vddp +(m*m+n*n*l*l)*Vddd
         goto 100
       endif


       if ((orb1.eq.5.and.orb2.eq.6).or.(orb1.eq.6.and.orb2.eq.5)) then  ! dxy - dyz           ! avril17
         SlaterKoster_param = 3.0*l*n*m*m*Vdds + l*n*(1.0-4.0*m*m)*Vddp + l*n*(m*m-1.0)*Vddd
         goto 100
       endif
       if ((orb1.eq.6.and.orb2.eq.7).or.(orb1.eq.7.and.orb2.eq.6)) then  ! dyz - dzx           ! avril17
         SlaterKoster_param = 3.0*l*m*n*n*Vdds + m*l*(1.0-4.0*n*n)*Vddp + m*l*(n*n-1.0)*Vddd
         goto 100
       endif

       if ((orb1.eq.5.and.orb2.eq.7).or.(orb1.eq.7.and.orb2.eq.5)) then  ! dxy - dzx           ! avril17
         SlaterKoster_param = 3.0*l*l*m*n*Vdds + m*n*(1.0-4.0*l*l)*Vddp + m*n*(l*l-1.0)*Vddd
         goto 100
       endif

       if ((orb1.eq.5.and.orb2.eq.8).or.(orb1.eq.8.and.orb2.eq.5)) then  ! dxy - dx2-y2        ! avril17
         SlaterKoster_param = 1.5*l*m*(l*l-m*m)*Vdds + 2.0*l*m*(m*m-l*l)*Vddp + 0.5*l*m*(l*l-m*m)*Vddd
         goto 100
       endif

       if ((orb1.eq.6.and.orb2.eq.8).or.(orb1.eq.8.and.orb2.eq.6)) then  ! dyz - dx2-y2        ! avril17
         SlaterKoster_param = 1.5*m*n*(l*l-m*m)*Vdds - m*n*(1.0+2.0*(l*l-m*m))*Vddp &          ! corrige le 31dec17
                            + m*n*(1.0+0.5*(l*l-m*m))*Vddd
         goto 100
       endif

       if ((orb1.eq.7.and.orb2.eq.8).or.(orb1.eq.8.and.orb2.eq.7)) then  ! dzx - dx2-y2       ! avril17
         SlaterKoster_param = 1.5*n*l*(l*l-m*m)*Vdds + n*l*(1.0-2.0*(l*l-m*m))*Vddp &
                            - n*l*(1.0+0.5*(l*l-m*m))*Vddd                                    ! corriger le 29dec17
         goto 100
       endif

       if ((orb1.eq.5.and.orb2.eq.9).or.(orb1.eq.9.and.orb2.eq.5)) then  ! dxy - d3z2-r2      ! avril17
         SlaterKoster_param = l*m*(n*n-0.5*(l*l+m*m))*Vdds*sqrt(3.0)    &
                            - 2.0*l*m*n*n*Vddp*sqrt(3.0)              &
                            + l*m*(1.0+n*n)*Vddd*0.5*sqrt(3.0)
         goto 100
       endif

       if ((orb1.eq.6.and.orb2.eq.9).or.(orb1.eq.9.and.orb2.eq.6)) then  ! dyz - d3z2-r2    ! avril17
         SlaterKoster_param = m*n*(n*n-0.5*(l*l+m*m))*Vdds*sqrt(3.0)    &
                            + m*n*(l*l+m*m-n*n)*Vddp*sqrt(3.0)              &              ! corrige le 29dec17
                            - m*n*(l*l+m*m)*Vddd*0.5*sqrt(3.0)
         goto 100
       endif
 
       if ((orb1.eq.7.and.orb2.eq.9).or.(orb1.eq.9.and.orb2.eq.7)) then  ! dzx - d3z2-r2    ! avril17
         SlaterKoster_param = l*n*(n*n-0.5*(l*l+m*m))*Vdds*sqrt(3.0)    &
                            + l*n*(l*l+m*m-n*n)*Vddp*sqrt(3.0)              &             ! corrige le 29dec17
                            - l*n*(l*l+m*m)*Vddd*0.5*sqrt(3.0)
         goto 100
       endif

       if (orb1.eq.8.and.orb2.eq.8) then  ! dx2-y2 - dx2-y2
         SlaterKoster_param = (l*l-m*m)*(l*l-m*m)*Vdds*0.75    &
                            + (l*l+m*m-(l*l-m*m)*(l*l-m*m))*Vddp              &
                            + (n*n+0.25*(l*l-m*m)*(l*l-m*m))*Vddd
         goto 100
       endif

       if ((orb1.eq.8.and.orb2.eq.9).or.(orb1.eq.9.and.orb2.eq.8)) then  ! dx2-y2 - d3z2-r2  ! avril17
         SlaterKoster_param = (l*l-m*m)*(n*n-0.5*(l*l+m*m))*Vdds*0.5*sqrt(3.0)    &
                            + n*n*(m*m-l*l)*Vddp*sqrt(3.0)              &
                            + (1.0+n*n)*(l*l-m*m)*Vddd*0.25*sqrt(3.0)
         goto 100
       endif
      
       if (orb1.eq.9.and.orb2.eq.9) then  ! d3z2-r2 - d3z2-r2
         SlaterKoster_param = (n*n-0.5*(l*l+m*m))*(n*n-0.5*(l*l+m*m))*Vdds    &
                            + 3.0*n*n*(l*l+m*m)*Vddp              &
                            + 0.75*(l*l+m*m)*(l*l+m*m)*Vddd
         goto 100
       endif
 
      stop'!!! Func.: SlaterKoster_AttenD: integrale de saut entre orbitales s, p et d uniquement !!!'
      return

 100  continue

  if (orb1i.gt.orb2i) then
	 if (orb1i.eq.1) then
	    l_orb1i = 0
	 else if (orb1i.le.4) then
            l_orb1i = 1
	 else if (orb1i.le.9) then
	    l_orb1i = 2
	 else 
	    stop 'Function: SlaterKoster_AttenD : orb1 > 9'
	 endif
	 if (orb2i.eq.1) then
	    l_orb2i = 0
	 else if (orb2i.le.4) then
            l_orb2i = 1
	 else if (orb2i.le.9) then
	    l_orb2i = 2
	 else 
	    stop 'Function: SlaterKoster_AttenD : orb2 > 9'
	 endif

     SlaterKoster_param = SlaterKoster_param * (-1.0)**(l_orb1i + l_orb2i)

	 !print*,'orb1i, orb2i =',orb1i, orb2i
	 !print*,'l_orb1i, l_orb2i,SlaterKoster_param =',l_orb1i, l_orb2i,SlaterKoster_param,(-1.0)**(l_orb1i + l_orb2i)
  endif

  SlaterKoster_AttenD  = SlaterKoster_param

  end function SlaterKoster_AttenD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 8fev18 Calcul de l'effet de la distance 
!
real(kind(1.d0)) function calcul_effet_dis(t,q,r0,d)
!
!    v1 = t(r=r0)
!    v2 = q
!    v3 = r0 = dmin
!    d = distance entre orbital
!    calcul_effet_dis = t(r=d)    output
!
!    if q=v2 > 0  : t(d) = t(r0) * exp(q(1 - d/rO))
!    if q=v2 < 0  : t(d) = t(r0) * (d/d0)^q 
!    if q=0 : t(d) = t(r0)
!

USE precision_reel

real(kind=prec) t,q,r0,d


if (q.ge.0) then  
  if (r0.eq.0) then
     !print*,'WARNING: in INPUT_SlaterKoster : dmin <> 0 ==> dmin = 1. !!!'
     r0=1.d0
  endif
  calcul_effet_dis = t * exp(q*(1-d/r0))
else if (q.lt.0) then
  if (r0.eq.0) then
    calcul_effet_dis = t
  else
    calcul_effet_dis = t * (d/r0)**q
  endif
endif

end function calcul_effet_dis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Ecriture_Hk(vec_k00,Hk)
!
!   2janv18 spin et SOC
!
!   Ecriture de Hk 
!
!      entree : vec_k
!      sortie : Hk = Hamiltonien H(k)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE precision_reel 
  USE constantes
  USE VarReseau
  USE VarHamil
  implicit none

  
                                           
  complex(kind=prec),dimension(1:ndim_max*nspin,1:ndim_max*nspin) :: Hk ! sortie

  integer i,i2,i_type,i2_type,j,orb1,orb2,nnspin,nnspin2
  real(kind=prec) R(1:dimension_pos_atom),vec_k00(1:dimension_reseau)
  real(kind=prec) vec_k(1:dimension_pos_atom),Rat_ij

  complex(kind=prec),parameter :: i_imag = (0.0d0,1.0d0)
  complex(kind=prec),parameter :: i_2pi =i_imag*6.28318530717958648d0
  integer ndim_res

vec_k(:) = 0.0
vec_k(1:dimension_reseau) = vec_k00(1:dimension_reseau)

! Ecriture du hamiltonien Hk:
Hk = (0.0d0,0.0d0)
! Element diagonaux initiaux (energie de site et terme intra-atome (SOC))
do i=1,natomes
   i_type = type_atome(i)
   do orb1=1,orbitales_lm(i_type,0)
     if (nspin.eq.1) then
	Hk(indice_orb(i,orb1),indice_orb(i,orb1)) & 
          = g(i,0,orb1,orb1,1,1)
     else                              ! 2janv18
      do orb2=1,orbitales_lm(i_type,0)
        !print*,"orb1,orb2,indice_orb(i,orb1),indice_orb(i,orb2)",orb1,orb2,indice_orb(i,orb1),indice_orb(i,orb2)
	do nnspin = 1,nspin                      
	do nnspin2 = 1,nspin 
	    Hk(indice_orb(i,orb1)+(nnspin-1)*ndim,indice_orb(i,orb2)+(nnspin2-1)*ndim) & 
		= g(i,0,orb1,orb2,nnspin,nnspin2)
	enddo
	enddo
      enddo
     endif
   enddo
enddo

do i=1,natomes
i_type = type_atome(i)
do j=1,nbvoisins(i)
   i2 = voisins(i,j,0)  ! numero du jeme voisin de i
   i2_type = type_atome(i2)
   ! R vecteur de la maille dans laquelle est l'atome i2 (voisin de i)
   R(:) = 0.0
   do ndim_res = 1,dimension_reseau
        R(:) = R(:) + a(ndim_res,:)*voisins(i,j,ndim_res) 
   enddo 
   do orb1=1,orbitales_lm(i_type,0)
   do orb2=1,orbitales_lm(i2_type,0)
      Hk(indice_orb(i,orb1),indice_orb(i2,orb2)) = & 
           Hk(indice_orb(i,orb1),indice_orb(i2,orb2)) + &
           (exp(i_2pi*DOT_PRODUCT(vec_k,R))) * g(i,j,orb1,orb2,1,1)

      if (nspin.eq.2) then  ! 2janv18
          !spin-orbite entre orbitale d'une meme atome uniquement  
	  !Hk(indice_orb(i,orb1),ndim+indice_orb(i2,orb2)) = & 
          !    Hk(indice_orb(i,orb1),ndim+indice_orb(i2,orb2)) + &
          !    (exp(i_2pi*DOT_PRODUCT(vec_k,R))) * g(i,j,orb1,orb2,1,nspin)
	  !
	  !Hk(ndim+indice_orb(i,orb1),indice_orb(i2,orb2)) = & 
          !    Hk(ndim+indice_orb(i,orb1),indice_orb(i2,orb2)) + &
          !    (exp(i_2pi*DOT_PRODUCT(vec_k,R))) * g(i,j,orb1,orb2,nspin,1)

	   Hk(ndim+indice_orb(i,orb1),ndim+indice_orb(i2,orb2)) = & 
               Hk(ndim+indice_orb(i,orb1),ndim+indice_orb(i2,orb2)) + &
               (exp(i_2pi*DOT_PRODUCT(vec_k,R))) * g(i,j,orb1,orb2,nspin,nspin)
      endif

   enddo
   enddo
enddo 
enddo

!print*,'Hk pour vec_k=', vec_k
!print*,'partie reelle de Hk ='
!do i=1,ndim*nspin
!   print*,real(real(Hk(i,1:ndim*nspin)))
!enddo
!print*,'partie imaginaire de Hk ='
!do i=1,ndim*nspin
!   print*,real(imag(Hk(i,1:ndim*nspin)))
!enddo
!stop

end SUBROUTINE Ecriture_Hk

