module mod_strat1

public :: repartition, branch_maker, down_grow_reser, vol_constr, mainte_cal, up_grow_mainte,&
		& upper_need_cal, birth_position, mech_constr_cal
!*****************************************************************************80  
!*****************************************************************************80  

!The following subroutines are the ones called within "flux_up"/"down" routines.
!
!Each of them takes as input the array of integer "schem". And they use one integer of that array 
!for a "select case". The routine will then behave according to the select case.
! 
!It means that "schem" is deciding the rules used to model the tree


!To understand their role, look at the flux_down/up routines

!===========================================================================================================================
! Pour le moment la plupart select case n'ont qu'un seul "case", mais ils se rempliront au fur et à mesure que l'on ajoute des modèles/stratégies.
! L'idée est que le module répertoriera toute ces différentes stratégies
! et que le tableau "schem" en input nous permet de selectionner celles que l'on veut avant chaque simulation
!===========================================================================================================================
CONTAINS 
!******************************************************************************80  
! schem(1) decides : -branch_maker (depends on vol_constr and birth position) *80
!                    -birth_position
!                    -score_position_selec
!	                                                                          *80
! schem(2) decides : -vol_constr								              *80
!					 -mainte_cal				
!                    -mech_constr_cal
!                    -weight_cal (not implemented yet)
!																	
! schem(3) decides : -down_grow_reser (depends on vol_constr and mech_constr_cal)					
!
! schem(4) decides : -up_grow_mainte (vol_constr and mainte_cal is inside too)
!
! schem(5) decides : -upper_need_cal
!                    -upper_need_vol
!
! schem(6) decides : -production
!                    -debut_printemps(not implemented yet)
!
! schem(7) decides : -repartition
!
! schem(8) decides : -cut_choice
!
!*****************************************************************************80



!====================================================================================
subroutine birth_position(a, list_point, space_above, desir2, schem)
!====================================================================================
  use precis_mod, only                      : prec => working_precis, precisint	
  use parameters_module
  use mod_tree
  implicit none
  type(branch), pointer, intent(in)                        :: a
  integer, intent(inout)                          :: desir2
  integer                                         :: ii,kk,temp
  logical, dimension(:,:), intent(inout)          :: space_above
  integer, dimension(:,:),intent(inout)           :: list_point
  integer, dimension(:), intent(in)             :: schem
  
  real(kind=prec)                               :: xg, yg, weight !coordinates of the barycenter
  integer, dimension(3)                         :: parent_coord
  type(branch), pointer                         :: father, bb

 !MEASURE OF THE BARYCENTER OF THE BRANCH, ITS BROTHERS AND NEPHEWS----------
  weight = 0.0_prec
  if(associated(a%parent%p))then
      parent_coord(1) = a%parent%p%int_x
      parent_coord(2) = a%parent%p%int_y
      parent_coord(3) = a%parent%p%int_z
      father => a%parent%p
      
      xg = 0.0_prec
      yg = 0.0_prec
      nullify(bb)
      block
      integer :: i, j
      real(kind=prec) :: ga, gb
      do i = 1, father%N_children
           bb => father%child(i)%p

           !weight = real(bb%volume + bb%upper_vol,kind=prec) + weight
           
           !gb = real(bb%volume + bb%upper_vol,kind=prec)
           !ga = weight - gb
           !xg = (xg*ga + bb%upbarycentre(1)*gb)/weight  !  NEW BARYCENTER
           !yg = (yg*ga + bb%upbarycentre(2)*gb)/weight  ! 
           
           weight = real(bb%volume,kind=prec) + weight
           
           gb = real(bb%volume)
           ga = weight - gb
           xg = (xg*ga + real(bb%int_x)*gb)/weight  !  NEW BARYCENTER
           yg = (yg*ga + real(bb%int_y)*gb)/weight  ! 
           
           if(bb%N_children>0)then
               do j = 1, bb%N_children
                    weight = weight + real(bb%child(j)%p%volume)
                    gb = real(bb%child(j)%p%volume)
                    ga = weight - gb
                    xg = (xg*ga + real(bb%child(j)%p%int_x)*gb)/weight  !  NEW BARYCENTER
                    yg = (yg*ga + real(bb%child(j)%p%int_y)*gb)/weight  ! 
               enddo
           endif
      enddo
      endblock

  else
      parent_coord(1) = a%int_x
      parent_coord(2) = a%int_y
      parent_coord(3) = a%int_z
      xg = real(parent_coord(1))
      yg = real(parent_coord(2))
  endif
 !MEASURE OF THE BARYCENTER OF THE BRANCH, ITS BROTHERS AND NEPHEWS----------
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !parent_coord(1) = 0          !
  !parent_coord(2) = 0          ! EXPERIMENTATION
  !parent_coord(3) = 0          !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !--------------------------------------------------------------------------

	!particularity of how the 26 direction around the branch are indexed:
	!these 26 adjacents cases are indexed by a 2-coordinates system instead of 3, with the index (5,2) being the location of the branch itself
	!Meaning of the indices/coordinates of "space_above" 
		! Plane above 
	    ! (1,3)  (2,3)  (3,3) 
	    ! (1,2)  (2,2)  (3,2)    
	    ! (1,1)  (2,1)  (3,1) 
	    !
	    ! Same plane:
	    ! (4,3)  (5,3)  (6,3)
	    ! (4,2)  (5,2)  (6,2)    
	    ! (4,1)  (5,1)  (6,1)  
	    !
	    !  Plane below:
	    ! (7,3)  (8,3)  (9,3)
	    ! (7,2)  (8,2)  (9,2)    
	    ! (7,1)  (8,1)  (9,1) 
	    
   temp = schem(1)
   if(schem(1)==4)then
       temp = 2
   endif
   if(schem(1)==6)then
       temp = 5
   endif
   
   select case(temp)
    case(2)  ! Position for creation of children is chosen randomly but with a bias for the level above
             ! case 4 is included in this case 2 : case 4 = apical lead
    	block
		real :: rundom1, rundom2, cal1, cal2, rundom3
		integer :: x_dum, y_dum, zz!, trunk_x, trunk_y, z_untiltrunk_debug
		logical :: apical_lead
		
		!=======if schem(1)==4 then we will make an apical lead==============
		if(schem(1)==4 .and. desir2>0)then
			
			!we verify if the branch is on the apical lead
			call apical_lead_test(a,apical_lead)
			
			!===============================
			!now we know whether the branch is part of the apical lead
			if(apical_lead)then
				if(.not. space_above(2,2))then
					list_point(1,desir2) = a%int_x
					list_point(2,desir2) = a%int_y
					list_point(3,desir2) = a%int_z + 1
					space_above(2, 2) = .true.
				
					desir2 = desir2 - 1
				endif
			endif
			
		endif
		!====================================================================
		
		
		if(ran_posit<0.0_prec) stop
		
		!cal1 = 1.0_prec + ran_posit + ran_posit*ran_posit
		!cal2 = (1.0_prec + ran_posit)/cal1
		cal1 = a%spacebabywish_z(1) + a%spacebabywish_z(2) + a%spacebabywish_z(3)
		cal2 = (a%spacebabywish_z(1) + a%spacebabywish_z(2))/cal1

		if(desir2<0)then
			print*, "bug_epsilon003332"
			stop
		endif
		loop_name3 : do 
			
			if(desir2==0) exit loop_name3
					
			call random_number(rundom1)
			call random_number(rundom2)
			call random_number(rundom3)
		
			x_dum = floor(rundom1*3.0)+1
			y_dum = floor(rundom2*3.0)+1
		
			if(rundom3<cal2)then
				x_dum = x_dum + 3
				if(rundom3 < a%spacebabywish_z(1)/cal1)then
					x_dum = x_dum + 3
				endif
			endif

			zz = (x_dum-1)/3
						
			
			!----DEBUG--------------------------------
			if(x_dum>9 .or. x_dum<1)then
				print*,"error231442"
				stop
			endif
			if(y_dum>3 .or. y_dum<1)then
				print*,"error233<33571442", y_dum
				stop
			endif
			!----DEBUG---------------------------------
			
			if(.not. space_above(x_dum, y_dum))then
				list_point(1,desir2) = a%int_x + x_dum - 2 - 3*zz
				list_point(2,desir2) = a%int_y + y_dum - 2
				list_point(3,desir2) = a%int_z + 1 - zz
				space_above(x_dum, y_dum) = .true.
				
				desir2 = desir2 - 1
			endif
			
		enddo loop_name3
		
			call random_number(rundom1) !
			call random_number(rundom2) ! je fais juste passer quelque numbre aleatoire
			call random_number(rundom3) !

		end block
		
	case(3) ! create children in a balanced way based on barycenter
		
		!Meaning of the indices of "space_above" 
		! Plane above 
	    ! (1,3)  (2,3)  (3,3) 
	    ! (1,2)  (2,2)  (3,2)    
	    ! (1,1)  (2,1)  (3,1) 
	    !
	    ! Same plane:
	    ! (4,3)  (5,3)  (6,3)
	    ! (4,2)  (5,2)  (6,2)    
	    ! (4,1)  (5,1)  (6,1)  
	    !
	    !  Plane below:
	    ! (7,3)  (8,3)  (9,3)
	    ! (7,2)  (8,2)  (9,2)    
	    ! (7,1)  (8,1)  (9,1) 
!------ If possible the new branch grows in the direction of its parent.---
		block
		integer :: x_dum, y_dum
		if(associated(a%parent%p))then
			x_dum = a%int_x - a%parent%p%int_x + 2 + 3*(a%parent%p%int_z - a%int_z + 1)
			y_dum = a%int_y - a%parent%p%int_y + 2
			
			!DEBUG--------------------------------
			if(x_dum>9 .or. x_dum<1)then
				print*,"error231442"
				stop
			endif
			if(y_dum>3 .or. y_dum<1)then
				print*,"error2333351442"
				stop
			endif
			if(x_dum - 3*(a%parent%p%int_z - a%int_z + 1)>3 .or. x_dum<1)then
				print*,"error231ze442"
				stop
			endif
			if(desir2<1)then
				print*,"erdez213"
				stop
			endif
			!DEBUG---------------------------------
			
			if(.not. space_above(x_dum, y_dum))then
				list_point(1,desir2) = a%int_x + x_dum-2 - 3*(a%parent%p%int_z - a%int_z + 1)
				list_point(2,desir2) = a%int_y + y_dum-2
				list_point(3,desir2) = a%int_z + a%int_z - a%parent%p%int_z
				space_above(x_dum, y_dum) = .true.
				
				weight = weight + 1.0_prec
				xg = (xg*(weight - 1.0) + real(list_point(1,desir2)))/weight  !  NEW BARYCENTER
				yg = (yg*(weight - 1.0) + real(list_point(2,desir2)))/weight  ! 
				
				desir2 = desir2 - 1
			endif
		endif
		end block
!--------------------------------------------------------------------------
		
		
		block
		integer :: x0, y0, dummyx, dummyy, inc, inc2
		integer,dimension(8,2) :: xylist = 0
		
		!From here on, we try to create the children such that the barycenter doesn't lean too far from the parent (try to keep some balance) 
		
		loop_name : do 
			
			if(desir2<1) exit loop_name
			
			dummyx = nint(10*(real(parent_coord(1)) - xg)) 
			dummyy = nint(10*(real(parent_coord(2)) - yg))
		
			if(abs(dummyx)+abs(dummyy)==0)then
				!print*,'A'
				dummyx = parent_coord(1) - a%int_x
				dummyy = parent_coord(2) - a%int_y
			endif
			
			if(dummyx/=0) dummyx = dummyx/abs(dummyx)
			if(dummyy/=0) dummyy = dummyy/abs(dummyy)			
			x0 = 2 + dummyx
			y0 = 2 + dummyy
			
			if(.not. space_above(x0, y0))then
				list_point(1,desir2) = a%int_x + x0 - 2 
				list_point(2,desir2) = a%int_y + y0 - 2
				list_point(3,desir2) = a%int_z + 1
				space_above(x0, y0) = .true.
				
				weight = weight + 1.0
				xg = (xg*(weight - 1.0) + real(list_point(1,desir2)))/weight  !  NEW BARYCENTER
				yg = (yg*(weight - 1.0) + real(list_point(2,desir2)))/weight  ! 		
	
				desir2 = desir2 - 1
				cycle loop_name
			endif
			
			!If space_above(x0,y0) is occupied : try to search an unoccupied space not too far from (x0,y0)
			!xylist(n,1) is the coord x while xylist(n,2) 
			
			if(x0+1>3)then
				xylist(1,1)=x0-1
				xylist(1,2)=y0
				if(y0-1<1)then
					xylist(2,2)=y0+1
					xylist(3,2)=y0+1
					xylist(4,2)=y0
					xylist(5,2)=y0+1
					xylist(6,2)=y0+2
					xylist(7,2)=y0+2
					xylist(8,2)=y0+2
				else
					xylist(2,2)=y0-1
					xylist(3,2)=y0-1
					xylist(4,2)=y0
					xylist(5,2)=y0-1
					xylist(6,2)=mod(y0-3+6,3)+1
					xylist(7,2)=mod(y0-3+6,3)+1
					xylist(8,2)=mod(y0-3+6,3)+1
				endif
				xylist(2,1)=x0
				xylist(3,1)=x0-1
				xylist(4,1)=x0-2 !mod(x0-3+6,3)+1
				xylist(5,1)=x0-2
				xylist(6,1)=x0-2
				xylist(7,1)=x0
				xylist(8,1)=x0-1
			else
				xylist(1,1)=x0+1
				xylist(1,2)=y0
				if(y0-1<1)then
					xylist(2,2)=y0+1
					xylist(3,2)=y0+1
					xylist(4,2)=y0
					xylist(5,2)=y0+1
					xylist(6,2)=y0+2
					xylist(7,2)=y0+2
					xylist(8,2)=y0+2
				else
					xylist(2,2)=y0-1
					xylist(3,2)=y0-1
					xylist(4,2)=y0
					xylist(5,2)=y0-1
					xylist(6,2)=mod(y0-3+6,3)+1
					xylist(7,2)=mod(y0-3+6,3)+1
					xylist(8,2)=mod(y0-3+6,3)+1
				endif
				xylist(2,1)=x0
				xylist(3,1)=x0+1
				xylist(4,1)=mod(x0+1+6,3)+1
				xylist(5,1)=mod(x0+1+6,3)+1
				xylist(6,1)=mod(x0+1+6,3)+1
				xylist(7,1)=x0
				xylist(8,1)=x0+1
			endif
			
			!DEBUGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
			do inc=1,8
				if(xylist(inc,1)<1 .or. xylist(inc,1)>3)then
					print*,"erffrf1820"
					stop
				endif
				if(xylist(inc,2)<1 .or. xylist(inc,2)>3)then
					print*,"erffrf182dezudue770"
					stop
				endif
			enddo
			!DEBUGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
			
			do inc2=0,2
				do inc=1,8
					if(.not. space_above(xylist(inc,1)+inc2*3, xylist(inc,2)))then
				
						list_point(1,desir2) = a%int_x + xylist(inc,1) - 2 
						list_point(2,desir2) = a%int_y + xylist(inc,2) - 2
						list_point(3,desir2) = a%int_z + 1 - inc2
						space_above(xylist(inc,1)+inc2*3, xylist(inc,2)) = .true.
				
						weight = weight + 1.0
						xg = (xg*(weight - 1.0) + real(list_point(1,desir2)))/weight  !  NEW BARYCENTER
						yg = (yg*(weight - 1.0) + real(list_point(2,desir2)))/weight  ! 		
	
						desir2 = desir2 - 1
						cycle loop_name
				 
					endif
				enddo
				if(inc2<2)then
					if(.not. space_above(x0+inc2*3+3, y0))then
						list_point(1,desir2) = a%int_x + x0 - 2 
						list_point(2,desir2) = a%int_y + y0 - 2
						list_point(3,desir2) = a%int_z - inc2
						space_above(x0+inc2*3+3, y0) = .true.
				
						weight = weight + 1.0
						xg = (xg*(weight - 1.0) + real(list_point(1,desir2)))/weight  !  NEW BARYCENTER
						yg = (yg*(weight - 1.0) + real(list_point(2,desir2)))/weight  ! 		
	
						desir2 = desir2 - 1
						cycle loop_name
					endif
				endif
			enddo
			
			
			if(desir2>0)then
				print*, "premier disease 404 in loop_name"
				stop
			else
				print*, "terminal disease 404"
				stop
			endif
			
		end do loop_name
		
		end block
			
	case(5) !Take into account the direction of the SUN and some balancing issue into a score

	!------ If possible the new branch grows in the direction of its parent.---
		block
		integer :: x_dum, y_dum
		if(associated(a%parent%p))then
			x_dum = a%int_x - a%parent%p%int_x + 2 + 3*(a%parent%p%int_z - a%int_z + 1)
			y_dum = a%int_y - a%parent%p%int_y + 2
			
			!DEBUG--------------------------------
			if(x_dum>9 .or. x_dum<1)then
				print*,"error231442"
				stop
			endif
			if(y_dum>3 .or. y_dum<1)then
				print*,"error2333351442"
				stop
			endif
			if(x_dum - 3*(a%parent%p%int_z - a%int_z + 1)>3 .or. x_dum<1)then
				print*,"error231ze442"
				stop
			endif
			if(desir2<1)then
				print*,"erdez213"
				stop
			endif
			!DEBUG---------------------------------
			
			if(.not. space_above(x_dum, y_dum))then
				list_point(1,desir2) = a%int_x + x_dum-2 - 3*(a%parent%p%int_z - a%int_z + 1)
				list_point(2,desir2) = a%int_y + y_dum-2
				list_point(3,desir2) = a%int_z + a%int_z - a%parent%p%int_z
				space_above(x_dum, y_dum) = .true.
				
				desir2 = desir2 - 1
			endif
		else
			x_dum = 2
			y_dum = 2
			if(.not. space_above(x_dum, y_dum))then
				list_point(1,desir2) = a%int_x 
				list_point(2,desir2) = a%int_y 
				list_point(3,desir2) = a%int_z + 1
				space_above(x_dum, y_dum) = .true.
				
				desir2 = desir2 - 1
			endif
		endif
		end block
	!--------------------------------------------------------------------------
	block
	integer, dimension(9,3) :: score_grid
	integer, dimension(27,2) :: indexx
	integer, dimension(2) ::indices
	integer :: zzz, ij, kl, maxx, n_v
	real :: shuffl
	
	call score_position_selec(a, space_above, score_grid, schem(1))
	
	
	do
		if(desir2==0) exit
		maxx = maxval(score_grid)
		
		if(maxx<-30)then
			print*,"errfferreur",maxx
			print*,score_grid
			print*, space_above
			stop
		endif
		n_v = 0
		indexx = 0
		indices = 0
		do ij=1,9
			do kl=1,3
				if(score_grid(ij,kl)==maxx)then
					n_v = n_v + 1
					indexx(n_v,1) = ij
					indexx(n_v,2) = kl
				endif
			enddo
		enddo
		
		if(n_v==0)then
			print*, "er042"
			stop
		endif
		call random_number(shuffl)
		shuffl = shuffl*real(n_v)
		ij = ceiling(shuffl)
		if(ij==0) ij = 1
		indices(1) = indexx(ij,1)
		indices(2) = indexx(ij,2)
		
		if(.not. space_above(indices(1), indices(2)))then
				
			zzz = (indices(1)-1)/3
			list_point(1,desir2) = a%int_x + indices(1) - 2 - 3*zzz
			list_point(2,desir2) = a%int_y + indices(2) - 2
			list_point(3,desir2) = a%int_z + 1 - zzz
			space_above(indices(1), indices(2)) = .true.
			
			desir2 = desir2 - 1
			
			!print*, score_grid, 1
			score_grid(indices(1),:) = score_grid(indices(1),:) - 3
			score_grid(indices(1),indices(2)) = -100000
			!print*, score_grid, indices
		else
			print*, "BUUUGGGG"
			print*, score_grid, desir2
			stop
		endif
		
	enddo
	
	end block
	case default
		print*, 'error brirth position'
		stop 
	end select
				
end subroutine birth_position
!==============================================================================

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
subroutine score_position_selec(a, space_above, score_grid, num)
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  use precis_mod, only                      : prec => working_precis, precisint	
  use parameters_module
  use mod_tree
implicit none
	type(branch), intent(in) :: a
	integer, intent(in) :: num
	logical, dimension(9,3), intent(in) :: space_above
	integer, dimension(9,3), intent(out) :: score_grid
	integer :: answer, iii, sampling, entier, jjj!, int_loop
	real(kind=prec) :: angle, radians, angle_moy, dummy_var!, photon
	integer, dimension(3) :: old_direc
	integer, dimension(0:6) :: convert
	integer :: int_angle, elimination
	integer, dimension(2) :: indices
	real(kind=prec) :: temporary
		!int_loop = desir2
		
			score_grid = 0
			!--Initialize score grid: already occupied space are excluded (less than 0 means exclusion)---
			block
			integer :: loopi, loopj
			do loopi=1,9
				do loopj=1,3
					if(space_above(loopi,loopj)) score_grid(loopi,loopj) = -100000
				enddo
			enddo
			end block
			!------------------------------------------------------------------------------------
			
			sampling = 10
			!sum_angle = 0.0_prec
			angle_moy = 0.0_prec
			!photon = 0.00000001_prec
			!------WE CALC THE DIRECTION FROM WHICH THE SUN COMES FROM (on average)
			do iii=0,2*sampling
				angle = real(iii,kind=prec)*90.0_prec/real(sampling,kind=prec)
				radians  = angle*atan(1.0_prec)/45.0_prec
				dummy_var = radians*sin(radians)
				
				if(angle>180.001_prec)then
					print*, "new er3r3"
					stop
				endif
				
				!---------------------------------------------------------------
				if(num==6) dummy_var = radians    ! HERE THE CASE 6 APPEARS
				!---------------------------------------------------------------
			
				!sum_angle = sum_angle + dummy_var
				call shading_at_equador(a%int_x, a%int_y, a%int_z, angle, a%owner%tree_point%the_space, answer)
				!---------------------! 
				if(answer<0)then      ! 
					print*,"er98eeeerf666"  ! DEBUGGMODE 
					stop              ! 
				endif                 ! 
				!---------------------! 
				if(answer>0) answer = 1
				angle_moy = angle_moy + (1.0_prec-real(answer,kind=prec))*dummy_var
			enddo
			angle_moy = angle_moy/real(2*sampling,kind=prec)
			!----------------------------------------------------------------
			
			!----------------------------------------------------------------------------------------------------------------
			!we determine the old direction (x,y,a) where a is the integer equivalent to the direction for "convert" function
			if(associated(a%parent%p))then
				old_direc(1) = a%int_x - a%parent%p%int_x + 2 + 3*(a%parent%p%int_z - a%int_z + 1)
				old_direc(2) = a%int_y - a%parent%p%int_y + 2
			else
				old_direc(1) = 2
				old_direc(2) = 2
			endif
			!----------------------------------------------------------------------------------------------------------------
			
			!----------NOW LET'S DETERMINE THE DIRECTION OF THE CHILD FROM A SCORE-SYSYTEM--------
			! (1,3)  (2,3)  (3,3)
			! (1,2)  (2,2)  (3,2)
			! (1,1)  (2,1)  (3,1)
			!
			! Same plane:
			! (4,3)  (5,3)  (6,3)
			! (4,2)  (5,2)  (6,2)
			! (4,1)  (5,1)  (6,1)
			!
			!  Plane below:
			! (7,3)  (8,3)  (9,3)
			! (7,2)  (8,2)  (9,2)
			! (7,1)  (8,1)  (9,1)
			convert(0) = 9
			convert(1) = 6
			convert(2) = 3
			convert(3) = 2
			convert(4) = 1
			convert(5) = 4
			convert(6) = 7
			
			old_direc(3) = -1
			do iii=0,6
				if(convert(iii)==old_direc(1))then
					old_direc(3) = iii
				endif
			enddo
			
			temporary = angle_moy*5.0_prec/(4.0_prec*atan(1.0_prec))
			if(temporary<=0.0_prec) temporary = temporary + 0.00001_prec
			int_angle = ceiling(temporary)
			!print*,int_angle, angle_moy
			!===DEBUGGGG==========
			if(int_angle < 1)then
				print*,"ferfrreffer"
				stop
			endif
			if(int_angle > 5)then
				print*,"ferdddedeedfrreffer"
				stop
			endif
			!=====================
			
			!~~~~~~~~~~~~~~THE CHILD WILL GROW TOWARD THE SUN~~~~~~~~~~~~~~~~~~~~~~~~~~~
			score_grid(convert(int_angle),:)=score_grid(convert(int_angle),:)+60
			!-------the adjacent locations are the second priority :
			score_grid(convert(int_angle+1),:)=score_grid(convert(int_angle+1),:)+40
			score_grid(convert(int_angle-1),:)=score_grid(convert(int_angle-1),:)+40
			!-------third priority is the rough direction :
			if(int_angle>1) score_grid(convert(int_angle-2),:)=score_grid(convert(int_angle-2),:)+20
			if(int_angle<5) score_grid(convert(int_angle+2),:)=score_grid(convert(int_angle+2),:)+20
			!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			
			!==============HERE WE GIVE POINTS TO ADJACENT DIRECTONS TO THE INITIAL ONE============================================
			if(old_direc(3)>-1)then
				if(old_direc(3)>0) score_grid(convert(old_direc(3)-1),old_direc(2)) = score_grid(convert(old_direc(3)-1),old_direc(2))+2
				if(old_direc(3)<6) score_grid(convert(old_direc(3)+1),old_direc(2)) = score_grid(convert(old_direc(3)+1),old_direc(2))+2
				if(old_direc(3)==0 .or. old_direc(3)==6) score_grid(8,old_direc(2)) = score_grid(8,old_direc(2))+2
			else
				if(old_direc(1)/=5 .and. old_direc(1)/=8)then
					print*,"frefererferf", old_direc(3), old_direc(1), old_direc(2)
					stop
				endif
				score_grid(old_direc(1)+1,old_direc(2)) = score_grid(old_direc(1)+1,old_direc(2))+2
				score_grid(old_direc(1)-1,old_direc(2)) = score_grid(old_direc(1)-1,old_direc(2))+2
			endif
			
			if(old_direc(2)>1) score_grid(old_direc(1),old_direc(2)-1) = score_grid(old_direc(1),old_direc(2)-1)+2
			if(old_direc(2)<3) score_grid(old_direc(1),old_direc(2)+1) = score_grid(old_direc(1),old_direc(2)+1)+2
			
			if(old_direc(2)==2)then
				if(old_direc(1)+3>9) score_grid(old_direc(1)+3,old_direc(2)) = score_grid(old_direc(1)+3,old_direc(2))+2
				if(old_direc(1)-3<1) score_grid(old_direc(1)-3,old_direc(2)) = score_grid(old_direc(1)-3,old_direc(2))+2
			endif
			!=======================================================================================================================
			
			block
			integer :: i0, j0, countting, x_ddum, y_ddum, z_ddum
			do i0 = 0,6
				countting = 0
				do j0=1,3
					if(score_grid(convert(i0),j0)>-1) countting = countting + 3
				enddo
				score_grid(convert(i0),:)=score_grid(convert(i0),:)+countting
				score_grid(convert(i0),2)=score_grid(convert(i0),2)+1
			enddo
			
!			do i0=1,9
!				do j0=1,3	
!					countting = 0
!				
!					z_ddum = (i0-1)/3
!					x_ddum = a%int_x + i0 - 2 - 3*z_ddum
!					y_ddum = a%int_y + j0 - 2
!					do iii=0,2*sampling
!						angle = real(iii,kind=prec)*90.0_prec/real(sampling,kind=prec)
!						radians  = angle*atan(1.0_prec)/45.0_prec			
!
!						call shading_at_equador(x_ddum, y_ddum, a%int_z-z_ddum+1, angle, a%owner%tree_point%the_space, answer)
!						!---------------------! 
!						if(answer<0)then      ! 
!							print*,"er98eeeerf666"  ! DEBUGGMODE 
!							stop              ! 
!						endif                 ! 
!						!---------------------! 
!						countting = countting + answer
!					
!					enddo
!					if(countting==2*sampling+1) score_grid(i0,j0) = -10000
!					
!				enddo
!			enddo
			end block
			!-------We only keep the best locations and eliminate the worse-----------
			!block
			!integer :: i1, i2
			!do i1 = 1,size(score_grid,1)
			!	do i2=1,size(score_grid,2)
			!		if(score_grid(i1,i2)<elimination) score_grid(i1,i2) = score_grid(i1,i2)-1000
			!	enddo
			!enddo
			!end block
			!------------------------------------------------------------------------
			
end subroutine score_position_selec

!======================================================================================
!------ Routine that decides when a branch creates a child(ren) / new leave(s) -------!
!======================================================================================
subroutine branch_maker(a,tree_object,prod,desir,generation,schem,volum_added)
!======================================================================================
  use precis_mod, only                      : prec => working_precis, precisint	
  use parameters_module
  use mod_tree
  implicit none
  type (branch), pointer, intent(inout)        :: a!, sommet
  type(arbre_pointer), pointer, intent(inout)  :: tree_object
  integer(kind=precisint), intent(inout)	   :: prod 
  integer, intent(inout)					   :: desir
  integer, intent(in)						   :: generation
  integer, dimension(:), intent(in)			   :: schem
  integer(kind=precisint),intent(inout)		   :: volum_added
  integer(kind=precisint)					   :: cost_cal, vol_init
  
  integer									   :: tempo
  
  cost_cal = extra_cost ! This sets the default cost of creation of 1 branch (can be changed in the select case, if you want something else)
  vol_init = 1          ! This sets the default volume of a newly born branch
  
  
  volum_added = 0  ! Counter that keeps track of the quantity of wood created (in volume unit)
  
  tempo = schem(1)  
  
  if(schem(1)==2)then  !
	tempo = 3          !
  endif                !
  
  if(schem(1)==4)then  !
	tempo = 3          !
  endif                !

  if(schem(1)==5)then  !
	tempo = 3          !
  endif                !
    
  if(schem(1)==6)then  !
	tempo = 3          !
  endif                !
	select case(tempo)
	
	case(0)		
		desir = min(prod/cost_cal, N_c_Max)
		if (a%N_children>0) then              
		
			if(schem(7)<5) then                           !
				if (desir/=0) then                        !
					print*, 'error at branch maker !!!'   ! DEBUG LINE : FOR THE CASES IN WHICH THE REPARTITION ROUTINE
					stop                                  !              DISTRIBUTES ALL THE RESSOURCE OF BRANCHES WITH CHILDREN
				endif 									  !
			endif
			
			desir = 0  ! Branch with children can't create more children
		endif
		
		if (desir>0) then
			block
			integer :: Nold_child
				Nold_child = a%N_children
				call new_branches(a,desir,generation,tree_object%peak_point)
				a%nb_leaves = 0  ! we remove the leaf when creating branches
				
				block												!
				integer :: i										!  By default when a branch is created its volume is 1
					do i=1,desir									!  This scheme allow you to give a different volume
						a%child(Nold_child+i)%p%volume = vol_init	!
						volum_added = volum_added + a%child(Nold_child+i)%p%volume
					enddo											!
				end block											!

			end block
			
			prod = prod - desir*cost_cal    ! Payment of the children created
			
		endif
		
	case(1) !------IDENTICAL TO CASE(0) EXCEPT THAT NOW LEAVES HAVE MAINTENANCE COSTS (LEAVES FALL IN WINTER AND HAS TO BE RENEWED IN SPRING)------
		desir = min(prod/cost_cal, N_c_Max)			
		
		if(a%N_children==0) then
			if (desir>0) then
				block
				integer :: Nold_child
					Nold_child = a%N_children
					call new_branches(a,desir,generation,tree_object%peak_point)
					a%nb_leaves = 0  ! we remove the leaf of the branch when it creates children
				
					block                                               !
					integer :: i                                        !  By default when a branch is created its volume is 1
						do i=1,desir                                    !  But this scheme allow you to give a different volume
							a%child(Nold_child+i)%p%volume = vol_init	!
							volum_added = volum_added + a%child(Nold_child+i)%p%volume
						enddo                                           !
					end block                                           !
				end block                                               !
			
				prod = prod - desir*cost_cal    ! Payment of the children created

			else
				if(a%nb_leaves/=1) then                        !
					print*, "error in case 1 of branch_maker"  ! DEBUG: IN THE CASE(6) (FOR REPARTITION) 
					stop                                       ! IT IS IMPOSSIBLE FOR AN EXTREMITY TO
				endif                                          ! NOT HAVE A LEAF AT THIS POINT
			    !----maintenance cost for the leaves!
				if(prod >= cost_leaf) then          !
					prod = prod - cost_leaf         !
					a%nb_leaves = 1                 !
				else                                !
					a%nb_leaves = 0                 !
				endif                               !
				!----maintenance cost for the leaves!
			endif
		endif
		
		
	case(3)  ! Adding space into the mix  
	! REMARQUE : LES CAS 2 ET 3 ONT LE MEME CODE (Mais spacial_recog se comportera differemment)
		block
		integer :: Nold_child, desir2, do_int, i
		logical, dimension(9,3)    :: space_around
		integer, dimension(:,:),allocatable     :: list_point3D
		type(great_space), pointer              :: space00
		
		desir = min(prod/cost_cal, N_c_Max)
		if (a%N_children>0) then
			do do_int=1,a%N_children
				if(.not. a%child(do_int)%p%death_flag)then
					desir = 0  ! Branch with living child(ren) can't create more children
					exit
				endif
			enddo
		endif
		                                         
		Nold_child = a%N_children
		
		space_around(:,:) = .false.
		desir2 = 0
		space00 => a%owner%tree_point%the_space
		
		if(desir<min_chil) desir = 0
		
		if(desir>0) call spacial_recog(a, desir2, space_around, space00) 
		!if(desir2<desir)then
		!	write(unit=20,fmt="(i7,i7,i7, i7)") a%int_x, a%int_y, a%int_z, generation
		!endif
		                                ! This routine determines which are the locations where our branch can grow to 
		                                ! desir2 = number of spaces/directions available
		desir = min(desir2,desir)       ! space_above = matrix that locates these direction
		desir2 = desir                  ! the different terms of space_above in fonction of the direction is
		                                ! (1,3)  (2,3)  (3,3)
		                                ! (1,2)  (2,2)  (3,2) 
		                                ! (1,1)  (2,1)  (3,1)
		                                ! the value 0 denotes an available poistion/direction and 1 an unavailable
		                            ! crea is a term used for debugging
		!print*,'hi'
        if(desir>0)then        
	                                                                   
	     ! Now after calling this subroutine: space_above is operational
	     ! and dummy1 points to the plan above
			call new_branches(a,desir,generation,tree_object%peak_point)
			a%nb_leaves = 0  ! we remove the leaf when creating branches
			!print*,'hoi'
			!************************************************************************************************************
			! We determine the position of the newly created branches
			allocate(list_point3D(3,desir)) ! list_point(1,:) = abscisse, list_point(2,:) = ordinate, list_point(3,:) = height
			                              ! we list the coordinates we will give to the 'desir' branches we created.
			call birth_position(a, list_point3D, space_around, desir2, schem)
 			! We finished determining the position of the new branches
			!*************************************************************************************************************
			if(desir2/=0)then      !
				print*,"gros bugg" ! DEBUG
				stop               !
			endif                  !
			

			!print*,"slh"
			call  make_spacial_branch(a, space00, list_point3D, Nold_child, desir, space_around)
		endif
		                                                !
		                                                !  By default when a branch is created its volume is 1
		do i=1,desir									!  This scheme allow you to give a different volume
			a%child(Nold_child+i)%p%volume = vol_init	!
			volum_added = volum_added + a%child(Nold_child+i)%p%volume
		enddo											!
		end block										!
		
		prod = prod - desir*cost_cal    ! Payment of the children created

	case default
		print*,"error in branch_maker"
		stop 
	end select 
	
 end subroutine branch_maker
!=============================================================================!
!-----------END BRANCH MAKER--------------------------------------------------!
!=============================================================================!



!****************************************************************************!80
function vol_constr(a,schem)												 !80
!****************************************************************************!80
  use precis_mod, only                                : prec => working_precis, precisint 
  use parameters_module
  use mod_tree
  implicit none
  type (branch), pointer, intent(in)        :: a
  integer, dimension(:), intent(in)			:: schem
  integer(kind=precisint)                   :: vol_constr
  integer(kind=precisint) 					:: x,y
  integer                                   :: tempo
  !***************************************************************************80
	vol_constr = 0
	
	tempo=schem(2)
	if(tempo==2) tempo=1
	
	select case(tempo)
	case(0) 
		vol_constr = a%upper_leaves !(Leonardo rule)  
		!x=nint(leaves_weight*real(a%upper_leaves,kind=prec)**exp_constr)  !Alternalive, volume depends on weight above it
		!y=nint(volum_weight*real(a%upper_vol,kind=prec)**exp_constr)
		!vol_constr=x+y
	case(1)
		if(schem(1)<2)then
			print*, "error, schem(1) must be superior to 2 (spacial startegies)"
			stop
		endif

		vol_constr = a%upper_leaves
	case(3)
		vol_constr = a%upper_leaves
	case default
		print*,"error vol_constr"
		stop 
	end select

end function vol_constr
!****************************************************************************!80

!***********************************************************************************
function mainte_cal(a,schem)   	                                                 !**
  use precis_mod, only                      : prec => working_precis, precisint  !**
  use parameters_module
  use mod_tree
  implicit none
  type (branch), pointer, intent(in)        :: a
  integer, dimension(:), intent(in)			:: schem
  integer(kind=precisint)                   :: mainte_cal
  real(kind=prec)                           :: true_volume
  integer                                   :: tempo
  
  mainte_cal = 0
  
	tempo = schem(2)
	if(tempo==2) tempo = 0
	
	select case(tempo)
	!case(0)
	!	mainte_cal = nint(real(a%volume,kind=prec)*rel_mainte_vol)

    !exp_mainte = the proportion of the branch that needs maintenance (1/2 if only the coss-section perimeter is alive)
    
	case(0)
		true_volume = real(a%volume,kind=prec)!*a%mech_constr
		mainte_cal = nint((true_volume**exp_mainte)*rel_mainte_vol)  
		
		                                                                          
	case(1)
		if(schem(1)<2)then
			print*, "error, schem(1) must be > 2 (spacial startegies)"
			
			stop
		endif

		true_volume = real(a%volume,kind=prec)*(1.0_prec+a%mech_constr)
		mainte_cal = nint((true_volume**exp_mainte)*rel_mainte_vol)  
	
	case(3)
		!true_volume = real(a%volume,kind=prec)!*a%mech_constr
		mainte_cal = nint(a%volume*rel_mainte_vol)  
	case default
		print*,"error mainte_cal"
		stop 
	end select
		

end function mainte_cal
!**********************************************************************************


!**********************************************************************************
subroutine mech_constr_cal(a,schem)
!**********************************************************************************
  use precis_mod, only                      : prec => working_precis, precisint 
  use parameters_module
  use mod_tree
  implicit none
  type (branch), pointer, intent(in)        :: a
  integer, dimension(:), intent(in)			:: schem
  integer                                   :: tempo
  real(kind=prec),dimension(2)              :: calc1
  real(kind=prec)                           :: calc2, weight
  real(kind=prec)                           :: xx, yy

	tempo = schem(2)
	if(tempo==2) tempo=1
	
	select case(tempo)
    
	case(0)
		
	
	case(1)
		if(schem(1)<2)then
			print*, "error, schem(1) must be > 2 (spacial startegies)"
			stop
		endif
		
		if(associated(a%parent%p))then
			xx = real(a%parent%p%int_x,kind=prec)
			yy = real(a%parent%p%int_y,kind=prec)
		else
			if(.not. associated(a,a%owner%trunk_point))then
				print*,"error unknown 4004"
				stop
			endif
			xx = 0.0_prec
			yy = 0.0_prec
		endif
		
		calc1(1) = a%upbarycentre(1) - xx
		calc1(2) = a%upbarycentre(2) - yy
		calc2 = sqrt(calc1(1)*calc1(1) + calc1(2)*calc1(2))
		weight = real(a%volume + a%nb_leaves*leaves_weight + a%upper_weight,kind=prec)
		
		a%mech_constr = weight*calc2*mech_coef !+ 1.0_prec
		
	case(3)
		
	case default
		print*,"error constr_cal"
		stop 
	end select

!************************************************************************************
end subroutine mech_constr_cal
!************************************************************************************


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Grow using energy and keep some reserve (the leftover will be given to the parent !
subroutine down_grow_reser(a,prod,schem)	                                     !!!!

! Routine called during flux_down: after receiving the flux down from the children the
! branch can grow using energy and keep some reserve (the leftover will be given to the parent)

  use precis_mod, only                      : prec => working_precis, precisint	 !!!!
  use parameters_module														     !!!!
  use mod_tree																     !!!!
  implicit none																     !!!!
  type (branch), pointer, intent(inout)        :: a							     !!!!
  integer, dimension(:), intent(in)			   :: schem						     !!!!
  integer(kind=precisint), intent(inout)	   :: prod						     !!!!
  integer(kind=precisint)					   :: desir1, desir2, temp!, des_cost1!!!
  real(kind=prec)                              :: cost_temp 
  integer                                      :: tempo
  !cost_temp incorporates the fact that a branch that grows on the side will need to create more wood to endure gravity.
  
  tempo = schem(3)
  if(schem(3)==2)tempo = 1
  
  select case(tempo)

  case(0)
    ! First we grow if we can and need to ----------
    desir1 = nint( real( vol_constr(a,schem) - a%volume, kind=prec) * frac_desir_down) 
    desir1 = max(0,desir1)
    if(cost_volum/=0)then
		cost_temp = real(cost_volum,kind=prec)!*a%mech_constr
		temp = int(real(prod,kind=prec)/cost_temp)
	else
		cost_temp = 0.0_prec
		temp = desir1
	endif
    desir1 = min(temp,desir1)
    a%volume = a%volume + desir1
    prod = prod - int(real(desir1,kind=prec)*cost_temp)
    if(prod<0) then 
		print*,"error nan",prod,temp,cost_volum,desir1
		stop 
	endif
	
	!reserve_______________________________________
    desir2=min(prod,a%reserve_wish)
    
    desir2=min(desir2,a%reserve_max)
    a%reserve=desir2
    prod=prod-desir2


  case(1)  ! Same as case 1 but we add mech failure at the fkux down
    ! First we grow if we can and need to ----------
    desir1 = nint( real( vol_constr(a,schem) - a%volume, kind=prec) * frac_desir_down) 
    desir1 = max(0,desir1)
    if(cost_volum/=0)then
		cost_temp = real(cost_volum,kind=prec)!*a%mech_constr
		temp = int(real(prod,kind=prec)/cost_temp)
	else
		cost_temp = 0.0_prec
		temp = desir1
	endif
    desir1 = min(temp,desir1)
    a%volume = a%volume + desir1
    prod = prod - int(real(desir1,kind=prec)*cost_temp) 
    if(prod<0) then 
		print*,"error nan",prod,temp,cost_volum,desir1
		stop 
	endif
	
	!reserve_______________________________________
    desir2=min(prod,a%reserve_wish)
    
    desir2=min(desir2,a%reserve_max)
    a%reserve=desir2
    prod=prod-desir2
    
    !----------Mechanic Failure ADDED----------------!
	block
	real(kind=prec) :: indicator
	indicator = base_rupture + rupture_coef*real(a%volume,kind=prec)**rupture_exp
		! 
		! (a%mech_constr - base_rupt) / rupt_coef = vol 
		!
		!
		
	if(a%mech_constr>indicator) then
		!if(.NOT. associated(a%parent%p%child(1)%p,a))then
			a%death_flag = .true.
		!endif
		a%upper_leaves = 1
		a%upper_vol = 0
		a%upper_weight = 0
		a%volume = 0
		a%nb_leaves = 0
		if(associated(a%parent%p))then
			a%upbarycentre(1) = a%parent%p%int_x 
			a%upbarycentre(2) = a%parent%p%int_y 
			a%upbarycentre(3) = a%parent%p%int_z
		else
			if(.not. associated(a%owner%trunk_point,a))then
				print*, "error non-trunk without children"
				stop
			endif
		endif
		if(schem(3)==2) then
			a%upper_weight = leaves_weight
		endif

	endif
	end block
	
	
    !-------------------------------------------!
!  case(2)  ! SAME AS 1 BUT RESERVE MAX DEPENDS ON VOL
    ! First we grow if we can and need to ----------
!    desir1 = nint( real( vol_constr(a,schem) - a%volume, kind=prec) * frac_desir_down) 
!    desir1 = max(0,desir1)
!    if(cost_volum/=0)then
!		temp = prod/cost_volum
!	else
!		temp = desir1
!	endif
!    desir1 = min(temp,desir1)
!    a%volume = a%volume + desir1
!    prod = prod - desir1*cost_volum
!    if(prod<0) then 
!		print*,"error nan"
!		stop 
!	endif
	
!	a%reserve_max = res_prod_rel*prod_leaf*a%volume
    
	!reserve_______________________________________
!    temp = nint(real(prod+a%reserve_wish)*0.5)  !replace 0.5 by a%wish_coeff later
!    desir2 = min(prod,temp)
    
!    desir2=min(desir2,a%reserve_max)
!    a%reserve = desir2
!    prod = prod - desir2
	
  case default
		print*,"error in down_grow_reser"
		stop 
  end select

end subroutine down_grow_reser
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!----------------------------------------------------------------------------------
subroutine up_grow_mainte(a,prod,schem)								!!!

!Routine called during flux_up: The branch may be forced to grow to the volume it has to grow
!                       The routine also decides its death (or not) by issuing a death_flag 

  use precis_mod, only                      : prec => working_precis, precisint !!!
  use parameters_module	                                                        !!!
  use mod_tree                                                                  !!!
  implicit none                                                                 !!!
  type (branch), pointer, intent(inout)        :: a
  integer, dimension(:), intent(in)			   :: schem
  !integer, intent(inout)			     	   :: death_flag
  integer(kind=precisint), intent(inout)	   :: prod
  integer(kind=precisint)					   :: desir, desir_cost
  real(kind=prec)                              :: cost_temp
  
  integer									   :: TEMPO	
  
	TEMPO =	schem(4)
  
  select case(TEMPO)

  case(0)		
		! We grow until we required volume to survive. DEATH if it can't grow to the correct size
		! GROWTH IN VOLUME : here, desir = volume augmentation it needs to survive !
		!************************************************************************  !
		desir = vol_constr(a,schem) - a%volume									   !	
		!************************************************************************  !
		desir = max(0,desir)
		a%volume = a%volume + desir
		cost_temp = real(cost_volum,kind=prec)!*a%mech_constr
		desir_cost = int(real(desir,kind=prec)*cost_temp) 
		prod = prod - desir_cost
		
		!+++++++++Maintenance payment (FULL)+++++++++++
		prod = prod - a%maintenance
		
  case(1)		
		! We grow until we required volume to survive. DEATH if it can't grow to the correct size or if the moments forces is too much
		! GROWTH IN VOLUME : here, desir = volume augmentation it needs to survive !
		!************************************************************************  !
		desir = vol_constr(a,schem) - a%volume									   !	
		!************************************************************************  !
		desir = max(0,desir)
		a%volume = a%volume + desir
		cost_temp = real(cost_volum,kind=prec)!*a%mech_constr
		desir_cost = int(real(desir,kind=prec)*cost_temp) 
		prod = prod - desir_cost
		
		!+++++++++Maintenance payment (FULL)+++++++++++
		block
		real(kind=prec) :: indicator
		indicator = base_rupture + rupture_coef*real(a%volume,kind=prec)**rupture_exp
		
		! 
		! (a%mech_constr - base_rupt) / rupt_coef = vol 
		!
		!
		
		if(a%mech_constr>indicator) then
			!if(.NOT. associated(a%parent%p%child(1)%p,a))then
				a%death_flag = .true.
			!endif
		endif
		end block
		
		prod = prod - a%maintenance
		
  case(2)
  		block
		integer :: temp
			temp=nint(real(a%flux_up)*0.5)
			if(desir>prod-temp) a%death_flag=.true. ! Using more than half the flux up for emergency growth means death in this model
		end block
		! GROWTH IN VOLUME : here, desir = volume augmentation it needs to survive !
		!************************************************************************  !
		desir = vol_constr(a,schem) - a%volume									   !	
		!************************************************************************  !
		desir = max(0,desir)
		a%volume = a%volume + desir
        desir_cost = desir*cost_volum
		prod = prod - desir_cost

		!+++++++++Maintenance payment (FULL)+++++++++++
		prod = prod - a%maintenance
		
  case default
		print*,"error in down_grow_reser"
		stop 
  end select
end subroutine up_grow_mainte
!----------------------------------------------------------------------------------



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function upper_need_cal(a,schem)                                               !%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  use precis_mod, only                      : prec => working_precis, precisint	
  use parameters_module
  use mod_tree
  implicit none
  type (branch), pointer, intent(in)        :: a
  integer, dimension(:), intent(in)	        :: schem
  integer(kind=precisint)                   :: upper_need_cal
  integer                                   :: i
  
  upper_need_cal = 0
  
  select case(schem(5))
	case(0)   ! CASE IN WHICH WE DON'T LIE
		upper_need_cal = max(0,vol_constr(a,schem) - a%volume)*cost_volum + a%maintenance
	case(1) 
		upper_need_cal = max(0,vol_constr(a,schem) - a%volume)*cost_volum + a%maintenance
		do i=1,a%N_children
			upper_need_cal = upper_need_cal + a%child(i)%p%upper_needs
		enddo
		stop
	case default
		print*,"error upper_need_cal "
		stop 
	end select

end function upper_need_cal
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!##########################################################################################
function production(a,tree_object,generation,schem)                                      !#
!##########################################################################################
  use precis_mod, only                                : prec => working_precis, precisint 
  use parameters_module
  use mod_tree	
	implicit none
  type (branch), pointer, intent(in)        :: a
  type(arbre_pointer), pointer, intent(in)  :: tree_object 
  integer, dimension(:), intent(in)			:: schem	
  integer, intent(in)						:: generation
  integer(kind=precisint)					:: production
  
  select case(schem(6))
  
  case(0)  !+++Production par défaut : SEULEMENT LES EXTREMITES PEUVENT AVOIR 1 FEUILLE (ET 1 SEULEMENT)+++ 
  
	!DEBUG CODE TO ENSURE ONLY EXTREMITIES HAVE A (UNIQUE) LEAF
	if (a%nb_leaves>0) then
		if(a%N_children>0) then
			print*, "production error 2: a non extremity has a leaf"
			stop
		endif
		if (a%nb_leaves>1) then
			print*, "production error 1: more than 1 leaf in generation", generation
			print*,a%nb_leaves
			stop
		endif
	endif
	if (schem(4)==1) then                                             !
		if(a%N_children<1) then                                       !
			if(a%nb_leaves/=1)then                                    !
				print*, "production error 3: an extremity has no leaf"! VERY SPECIAL DEBUG LINES
				stop                                                  !
			endif                                                     !
		endif                                                         !
	endif                                                             !
	!DEBUG END  -----------------------------------------------
	production = a%nb_leaves*Prod_leaf 
	
  case(1) !---------IDEM QUE CASE 0-------!MAIS AVEC UN ASPECT OMBRE/SOLEIL DIURNE
	if(schem(1)==0 .or. schem(1)==1)then
		print*, "erreor12444"
		stop
	endif
	
	if(a%nb_leaves>0)then
	
		block
		integer :: answer, iii, sampling
		real(kind=prec) :: angle, photon, sum_sinus, radians
	
		sampling = 10
		sum_sinus = 0.0_prec
		photon = 0.00000001_prec
		do iii=1,2*sampling-1
			angle = real(iii,kind=prec)*90.0_prec/real(sampling,kind=prec)
			radians  = angle*atan(1.0_prec)/45.0_prec
			sum_sinus = sum_sinus + sin(radians)
			call shading_at_equador(a%int_x, a%int_y, a%int_z, angle, a%owner%tree_point%the_space, answer)
			!---------------------! 
			if(answer<0)then      ! 
				print*,"er98666"  ! 
				stop              ! 
			endif                 ! 
			!---------------------! 
			if(answer>0) answer = 1
			photon = photon + (1.0_prec-real(answer,kind=prec))*sin(radians) 
		enddo
	
		production = a%nb_leaves*nint(photon/sum_sinus*real(Prod_leaf,kind=prec))
		if(production<0)then
			print*, photon, sum_sinus
			print*,"error prod negative"
			stop
		endif
		end block
	else
		production = 0 
	endif
  
  case(2) !---------IDEM QUE CASE 1-------!MAIS LA LUMINOSITE NE DEPEND PAS DE L'ANGLE (SKYLIGHT)
	if(schem(1)==0 .or. schem(1)==1)then
		print*, "erreor12444"
		stop
	endif
	
	if(a%nb_leaves>0)then
	
		block
		integer :: answer, iii, sampling
		real(kind=prec) :: angle, photon, sum_photon!, radians
	
		sampling = 20
		sum_photon = 0.0_prec
		photon = 0.00000001_prec
		do iii=1,2*sampling-1
			angle = real(iii,kind=prec)*90.0_prec/real(sampling,kind=prec)
			!radians  = angle*atan(1.0_prec)/45.0_prec
			sum_photon = sum_photon + 1.0_prec
			call shading_at_equador(a%int_x, a%int_y, a%int_z, angle, a%owner%tree_point%the_space, answer)
			!---------------------! 
			if(answer<0)then      ! 
				print*,"er98666"  ! 
				stop              ! 
			endif                 ! 
			!---------------------! 
			if(answer>0) answer = 1   ! a single branch will block all light
			photon = photon + (1.0_prec-real(answer,kind=prec))
		enddo
	
		production = a%nb_leaves*nint(photon/sum_photon*real(Prod_leaf,kind=prec))
		if(production<0)then
			print*, photon, sum_photon
			print*,"error prod negative"
			stop
		endif
		end block
	else
		production = 0 
	endif
  case(3)	!SHADOW/LIGHT VERTICAL
  
  	if(schem(1)==0 .or. schem(1)==1)then
		print*, "erreor1443R2444"
		stop
	endif
	
	if(a%nb_leaves>0)then
	
		block
		integer :: answer!, iii, sampling
		!real(kind=prec) :: angle, photon, sum_photon!, radians
		
		call shading_at_equador(a%int_x, a%int_y, a%int_z, 90.0_prec, a%owner%tree_point%the_space, answer)
		if(answer>0) answer = 1
		
		production = a%nb_leaves*Prod_leaf*(1-answer)
		if(production<0)then
			!print*, photon!, sum_photon
			print*,"error prod3 negative"
			stop
		endif
		end block
	else
		production = 0 
	endif
	
  case(4)   !SHADOW/LIGHT VERTICAL + CARDINAL DIRECTION
    if(schem(1)==0 .or. schem(1)==1)then
		print*, "erreor1443RRR2444"
		stop
	endif
	
	if(a%nb_leaves>0)then
	
		block
		integer :: answer!, iii, sampling
		real(kind=prec) :: photon, vertical_portion!, sum_photon!, radians
		vertical_portion = 0.2_prec
		
		photon = 0.0_prec
		!sum_photon = photon
		!------VERTICAL LIGHT-------------
		call shading_at_equador(a%int_x, a%int_y, a%int_z, 90.0_prec, a%owner%tree_point%the_space, answer)
		if(answer==0) then
			photon = vertical_portion
		endif
		!---------------------------------
		
		!-----EAST/WEST LIGHT-------------------
		call shading_at_equador(a%int_x, a%int_y, a%int_z, 0.0_prec, a%owner%tree_point%the_space, answer)
		if(answer==0) then
			photon = photon + (1.0_prec-vertical_portion)/5.0_prec
		endif
		call shading_at_equador(a%int_x, a%int_y, a%int_z, 180.0_prec, a%owner%tree_point%the_space, answer)
		if(answer==0) then
			photon = photon + (1.0_prec-vertical_portion)/5.0_prec
		endif
		!---NORTH/SOUTH LIGHT-------------------
		call shading_northsouth(a%int_x, a%int_y, a%int_z, .true., a%owner%tree_point%the_space, answer)
		if(answer==0) then
			photon = photon + (1.0_prec-vertical_portion)/5.0_prec
		endif
		call shading_northsouth(a%int_x, a%int_y, a%int_z, .false., a%owner%tree_point%the_space, answer)
		if(answer==0) then
			photon = photon + (1.0_prec-vertical_portion)/5.0_prec
		endif
		!--------------------------------------
		
		production = a%nb_leaves*nint(real(Prod_leaf,kind=prec)*photon)
		if(production<0)then
			print*, photon!, sum_photon
			print*,"error prod3 negative"
			stop
		endif
		end block
	else
		production = 0 
	endif
  case default
	print*,"error production", schem(6)
	stop 
	
  end select  
end function production
!##########################################################################################



!============================================================================80000 0
!---------  Repartition routine  --------------------------------------------80000 0
!============================================================================80000 0
subroutine repartition(a,prod,schem,generation) ! calcule le flux a remonter aux enfants de a 0
!============================================================================80000 0
  use precis_mod, only                          : prec => working_precis,precisint
  use parameters_module, only                   :  exp_egoisme,exp_reward,N_c_Max,random_flag,mech_coef
  use mod_tree
  implicit none
  type (branch), pointer, intent(inout)        :: a
  integer(kind=precisint),intent(inout)        :: prod
  integer, dimension(:), intent(in)			   :: schem	
  integer, intent(in)                          :: generation
  real(kind=prec)                              :: sum_flux_down
  integer(kind=precisint)                      :: prod_old, tempo, tempo2
  integer									   :: i, nchildren
  
  integer									   :: tempo_case
  !*****************************************************************************80
 
    sum_flux_down = 0.0_prec
    prod_old = prod
    
    if(prod<0) then
		print*,"error-1",prod
		stop 
    endif
    
    !=========================================================
    !REMARK : REPARTIRION IS ONLY CALLED IF a%N_children > 0 !
    !=========================================================
    
    tempo_case = schem(7)
	
	nchildren = 0
    
    !-------------------------------------------------------
    ! One can add new repartition scheme in the select case
    !-------------------------------------------------------
    select case (tempo_case)
    
    case (0) !--------------DEFAULT CASE-------NON DISSIPATIF _ PAS D'ALEATOIRE _ EXPOSANT EGO/ALTR----------
		
		do i=1,a%N_children
			if(.not. a%child(i)%p%death_flag)then
				sum_flux_down = sum_flux_down+ real(a%child(i)%p%flux_down,kind=prec)**exp_reward&
				&*real(a%child(i)%p%upper_needs,kind=prec)**exp_egoisme				
				nchildren = nchildren + 1
			endif
		enddo
		if (nint(sum_flux_down)==0) sum_flux_down=1.0
		do i=1,a%N_children
			if(.not. a%child(i)%p%death_flag)then
				a%child(i)%p%flux_up = int(floor((1.0*real(a%child(i)%p%flux_down,kind=prec)**exp_reward*&
				&real(a%child(i)%p%upper_needs,kind=prec)**exp_egoisme/sum_flux_down)*real(prod_old,kind=prec)))
				prod=prod-a%child(i)%p%flux_up
			endif
		enddo
		
! 		keep the excess/deficit of production due to the euclidean division
		if (prod.ne.0) then
			if(nchildren>0)then
				a%reserve = min(prod,a%reserve_max)		
			
				if(a%reserve<0) then											!
					print*,"repartition error", a%reserve, prod, prod_old, tempo! DEBUG LINES																							!
					stop 														!
				endif															!
			
				prod=0
			endif
		endif
! 		Redistribute the excess/deficit of production due to nint randomly
!		if (prod.ne.0) then
!			i=floor(a%N_children*rand())+1
!			a%child(i)%p%flux_up= a%child(i)%p%flux_up+prod
!			prod=0
!		endif
		
	case(1)    ! ------RANDOM REDISTRIBUTION OF RESSOURCE-----
	
	block
	real(kind=prec),dimension(a%N_children) :: rando
	real(kind=prec)  :: sum_rando
	rando = 0.0_prec
	sum_rando = 0.0_prec
	do i=1,a%N_children
		if(.not. a%child(i)%p%death_flag)then
			nchildren = nchildren + 1
			call random_number(rando(i))
			sum_rando = sum_rando + rando(i)
		endif
	enddo
	
	do i=1,a%N_children
		if(.not. a%child(i)%p%death_flag)then
			rando(i) = rando(i)/sum_rando
			a%child(i)%p%flux_up = int(floor(rando(i)*real(prod_old,kind=prec)))
			prod = prod - a%child(i)%p%flux_up
			!print*, prod_old
		endif
	enddo	
	
	!keep the excess/deficit of production due to the euclidean division
	if (prod .ne. 0) then
		if(nchildren>0)then
			a%reserve = min(prod,a%reserve_max)		
			
			if(a%reserve<0) then											!
				print*,"repartition error1", a%reserve, prod, prod_old! DEBUG LINES
				print*, rando			
				print*,"tttaaaaaaefefzfefz"
				print*, sum_rando
				do i=1,a%N_children
					print*,a%child(i)%p%flux_up
				enddo!
				stop 														!
			endif															!
		endif
		prod=0
	endif
	
	endblock
    case (2)  !-------CAS SYMETRIQUE NON DISSIPATIF--------!
		block
		integer(kind=precisint) :: debug_int
		do i=1,a%N_children
			if(.not. a%child(i)%p%death_flag) nchildren = nchildren + 1
		enddo
		tempo = 0
		tempo = prod/nchildren
		debug_int = mod(prod,nchildren)
		
		if (nint(sum_flux_down)==0) sum_flux_down=1.0d0
		do i=1,a%N_children
			if(.not. a%child(i)%p%death_flag)then
				a%child(i)%p%flux_up = tempo
				prod = prod - a%child(i)%p%flux_up
			endif
		enddo
		
		if(prod/=debug_int) then		!
			print*, "failure debug int"	! DEBUG LINES
			stop 	                    !
		endif							!
		
! 		keep the excess/deficit of production due to the euclidean division
		if (prod.ne.0) then
			if(nchildren>0)then
				a%reserve = min(prod,a%reserve_max)		
			
				if(a%reserve<0) then											!
					print*,"repartition error", a%reserve, prod, prod_old, tempo! DEBUG LINES																							!
					stop 														!
				endif															!
			
				prod=0
			endif
		endif
		end block

    case(3) !--------------ponctual anisotropic distribution--------------
	block
	real :: coinflip, isotropy
	integer :: dummy_var, dummy_flag, transfert
		
		isotropy = 0.99
		
		dummy_flag = 0
		call random_number(coinflip)
		if(generation/=20) coinflip = 0.1
		if(coinflip > isotropy) then
			if(a%N_children>=1) then
				if(a%child(1)%p%N_children==0)then
					if(random_flag<1) then 				! A modifier
						dummy_flag = 1					!
						random_flag = random_flag + 1	!
						print*,"hhahahahahahahahaaahahah&
						&hahhahhahhzhhehhhezhhzzhsazssazs"
						print*,"fhzfhzefezuhzeuifhzefiuh"
						print*,"fhzfhzefezuhzeuifhzefiuh"
						print*,"fhzfhzefezuhzeuifhzefiuh"
						print*,"fhzfhzefezuhzeuifhzefiuh"
						print*,"fhzfhzefezuhzeuifhzefiuh"
						print*,"fhzfhzefezuhzeuifhzefiuh"
					endif
				endif
			endif
		endif
		
		block
		real(kind=prec) :: flux_bas, proport
		
		do i=1,a%N_children
			if(a%child(i)%p%flux_down > 0) then
				flux_bas = real(a%child(i)%p%flux_down,kind=prec)
			else if (a%child(i)%p%flux_down==0) then
				flux_bas = 0.00001d0
			else
				stop
				print*,"error 301"
			endif
			sum_flux_down = sum_flux_down + flux_bas**exp_reward&
			&*real(a%child(i)%p%upper_needs,kind=prec)**exp_egoisme
		enddo
		
		if (nint(sum_flux_down)==0) sum_flux_down=1.0
		do i=1,a%N_children
		
			flux_bas = real(a%child(i)%p%flux_down,kind=prec)
			proport = (flux_bas**exp_reward)*&
			&(real(a%child(i)%p%upper_needs,kind=prec)**exp_egoisme)/sum_flux_down
			
			a%child(i)%p%flux_up = nint(proport*real(prod_old,kind=prec))
			prod=prod-a%child(i)%p%flux_up
		enddo
 		!Keep the excess/deficit of production due to nint randomly
		if (prod.ne.0) then
			!i=floor(a%N_children*rand())+1
			!a%child(i)%p%flux_up= a%child(i)%p%flux_up+prod
			!prod=0
			if(prod>0)then
				a%reserve = min(prod,a%reserve_max)		
				if(a%reserve<0) then											!
					print*,"repartition error", a%reserve, prod, prod_old, tempo! DEBUG LINES																							!
					stop 														!
				endif															!
				prod=0
			else
				prod=0
			endif
		endif
		
		end block
		
		if(dummy_flag == 1) then		
			a%child(1)%p%flux_up = a%child(1)%p%flux_up + (N_c_Max+1)*extra_cost	
		endif
		
	end block
    case(4)   ! Same as CASE 0 but also include %mech_constr
		do i=1,a%N_children
			if(.not. a%child(i)%p%death_flag)then
				sum_flux_down = sum_flux_down+ real(a%child(i)%p%flux_down,kind=prec)**exp_reward&
				&*real(a%child(i)%p%upper_needs,kind=prec)**exp_egoisme*(mech_coef+a%child(i)%p%mech_constr)
			endif
		enddo
		if (nint(sum_flux_down)==0) sum_flux_down=1.0
		do i=1,a%N_children
			if(.not. a%child(i)%p%death_flag)then
				a%child(i)%p%flux_up = int(floor((1.0*real(a%child(i)%p%flux_down,kind=prec)**exp_reward*&
				&real(a%child(i)%p%upper_needs,kind=prec)**exp_egoisme*(mech_coef+a%child(i)%p%mech_constr)&
				&/sum_flux_down)*real(prod_old,kind=prec)))
				prod=prod-a%child(i)%p%flux_up
			endif
		enddo
		
! 		keep the excess/deficit of production due to the euclidean division
		if (prod.ne.0) then
			a%reserve = min(prod,a%reserve_max)		
			
			if(a%reserve<0) then											!
				print*,"repartition error", a%reserve, prod, prod_old, tempo! DEBUG LINES																							!
				stop 														!
			endif															!
			
			prod=0
		endif
    case(5)  ! largely identical to case 0 but with some apical lead variable
    block
    logical  :: apical_lead, apical_child
    integer  :: leads, n_childlead
    !=========verify if we are on the apical lead====================
    call apical_lead_test(a,apical_lead)
	leads = 0
	if(apical_lead)then	
		do i=1,a%N_children
			call apical_lead_test(a%child(i)%p,apical_child)
			if(apical_child)then
				leads = leads + 1
				n_childlead = i
			endif
		enddo
    endif
    if(leads/=1)then
		apical_lead = .false.
    endif
    !---debug-------
	if(apical_lead)then
		block
		logical :: debuglog
		call apical_lead_test(a%child(n_childlead)%p,debuglog)
		if(.not. debuglog)then
			print*,"repatition bug apical lead"
			stop
		endif
		end block
	endif
    !---------------
    !if(generation>5)apical_lead = .false.
    !================================================================
    !we have now determined whether the branch is o the apical lead (apical_lead)
    !================================================================
    
		do i=1,a%N_children
			if(.not. a%child(i)%p%death_flag)then
				sum_flux_down = sum_flux_down+ real(a%child(i)%p%flux_down,kind=prec)**exp_reward&
				&*real(a%child(i)%p%upper_needs,kind=prec)**exp_egoisme				
				nchildren = nchildren + 1
			endif
		enddo
		! here if we are on the apical lead, the branch will be guaranteed to give 5% of all the flux to its successor
		if(apical_lead) sum_flux_down = 1.10_prec*sum_flux_down
		if (nint(sum_flux_down)==0) sum_flux_down=1.0
		do i=1,a%N_children
			if(.not. a%child(i)%p%death_flag)then
				a%child(i)%p%flux_up = int(floor((1.0*real(a%child(i)%p%flux_down,kind=prec)**exp_reward*&
				&real(a%child(i)%p%upper_needs,kind=prec)**exp_egoisme/sum_flux_down)*real(prod_old,kind=prec)))
				if(apical_lead)then
					if(i==n_childlead) a%child(i)%p%flux_up = a%child(i)%p%flux_up+(0.10_prec/1.10_prec)*real(prod_old,kind=prec)
				endif
				prod=prod-a%child(i)%p%flux_up
			endif
		enddo
		
! 		keep the excess/deficit of production due to the euclidean division
		if (prod.ne.0) then
			if(nchildren>0)then
				a%reserve = min(prod,a%reserve_max)		
			
				if(a%reserve<0) then											!
					print*,"repartition error", a%reserve, prod, prod_old, tempo! DEBUG LINES																							!
					stop 														!
				endif															!
			
				prod=0
			endif
		endif
	end block
    case(6) ! here apical lead = first born
    
		block
		logical :: apical_lead
		apical_lead = .false.
		do i=1,a%N_children
			if(.not. a%child(i)%p%death_flag)then
				sum_flux_down = sum_flux_down+ real(a%child(i)%p%flux_down,kind=prec)**exp_reward&
				&*real(a%child(i)%p%upper_needs,kind=prec)**exp_egoisme				
				nchildren = nchildren + 1
				if(i==1) apical_lead=.true.
			endif
		enddo
		! here if we are on the apical lead, the branch will be guaranteed to give 5% of all the flux to its successor
		if(apical_lead) sum_flux_down = 1.10_prec*sum_flux_down
		if (nint(sum_flux_down)==0) sum_flux_down=1.0
		do i=1,a%N_children
			if(.not. a%child(i)%p%death_flag)then
				a%child(i)%p%flux_up = int(floor((1.0*real(a%child(i)%p%flux_down,kind=prec)**exp_reward*&
				&real(a%child(i)%p%upper_needs,kind=prec)**exp_egoisme/sum_flux_down)*real(prod_old,kind=prec)))
				if(apical_lead)then
					if(i==1) a%child(i)%p%flux_up = a%child(i)%p%flux_up+(0.10_prec/1.10_prec)*real(prod_old,kind=prec)
				endif
				prod=prod-a%child(i)%p%flux_up
			endif
		enddo
		end block
		
! 		keep the excess/deficit of production due to the euclidean division
		if (prod.ne.0) then
			if(nchildren>0)then
				a%reserve = min(prod,a%reserve_max)		
			
				if(a%reserve<0) then											!
					print*,"repartition error", a%reserve, prod, prod_old, tempo! DEBUG LINES																							!
					stop 														!
				endif															!
			
				prod=0
			endif
		endif
		
	case(7)  ! symtric redistrib  but with some apical lead variable
    block
    logical  :: apical_lead, apical_child
    integer  :: leads, n_childlead
    
    !=========verify if we are on the apical lead====================
    call apical_lead_test(a,apical_lead)
	leads = 0
	if(apical_lead)then	
		do i=1,a%N_children
			call apical_lead_test(a%child(i)%p,apical_child)
			if(apical_child)then
				leads = leads + 1
				n_childlead = i
			endif
		enddo
    endif
    if(leads/=1)then
		apical_lead = .false.
    endif
    !---debug-------
	if(apical_lead)then
		block
		logical :: debuglog
		integer :: tempo2
		call apical_lead_test(a%child(n_childlead)%p,debuglog)
		if(.not. debuglog)then
			print*,"repatition bug apical lead"
			stop
		endif
		end block
	endif
    !---------------
    !if(generation>5)apical_lead = .false.
    !================================================================
    !we have now determined whether the branch is o the apical lead (apical_lead)
    !================================================================
    
		do i=1,a%N_children
			if(.not. a%child(i)%p%death_flag)then
				nchildren = nchildren + 1
			endif
		enddo
		! here if we are on the apical lead, the branch will be guaranteed to give 5% of all the flux to its successor
		if(apical_lead) tempo2 = floor(0.10_prec*real(prod,kind=prec))
		
		tempo = 0
		tempo = (prod-tempo2)/nchildren
		!debug_int = mod(prod,nchildren)
		
		if (nint(sum_flux_down)==0) sum_flux_down=1.0d0
		do i=1,a%N_children
			if(.not. a%child(i)%p%death_flag)then
				a%child(i)%p%flux_up = tempo
				if(apical_lead)then
					if(i==n_childlead)then
						a%child(i)%p%flux_up = a%child(i)%p%flux_up + tempo2
					endif
				endif
				prod = prod - a%child(i)%p%flux_up
			endif
		enddo
		
		
! 		keep the excess/deficit of production due to the euclidean division
		if (prod.ne.0) then
				a%reserve = min(prod,a%reserve_max)		
			
				if(a%reserve<0) then											!
					print*,"repartition error", a%reserve, prod, prod_old, tempo! DEBUG LINES																							!
					stop 														!
				endif															!
			
				prod=0
		endif
	end block
    case default
		block
		!real(kind=prec) :: sum_flux_down2, uno, due
		!integer :: prop_flux_up_ego
		!sum_flux_down2 = 0.0_prec
		do i=1,a%N_children
			if(.not. a%child(i)%p%death_flag)then
				sum_flux_down = sum_flux_down &
				&+ (rel_mainte_vol*real(a%child(i)%p%upper_leaves,kind=prec)**exp_mainte)**exp_egoisme
				
				!sum_flux_down2=sum_flux_down2+real(a%child(i)%p%flux_down,kind=prec)**(exp_egoisme*exp_mainte)
								
				nchildren = nchildren + 1
			endif
		enddo
		if (nint(sum_flux_down)==0) sum_flux_down=1.0
		do i=1,a%N_children
			if(.not. a%child(i)%p%death_flag)then
				a%child(i)%p%flux_up = int(floor(((rel_mainte_vol*real(a%child(i)%p%upper_leaves,kind=prec)&
				&**exp_mainte)**exp_egoisme/sum_flux_down)*real(prod_old,kind=prec)))
				prod=prod-a%child(i)%p%flux_up
				
				!prop_flux_up_ego = int(floor((real(a%child(i)%p%flux_down,kind=prec)**(exp_egoisme*exp_mainte)&
				!&/sum_flux_down2)*real(prod_old,kind=prec)))
				
				!uno = prop_flux_up_ego/prod_old
				!due = a%child(i)%p%flux_up/prod_old
				!if(abs(uno-due)>0.05)then
				!	print*,"hee-hooo", a%child(i)%p%flux_up, prop_flux_up_ego
				!	print*,"second", a%child(i)%p%upper_leaves, a%child(i)%p%volume
				!	print*,(rel_mainte_vol*real(a%child(i)%p%upper_leaves,kind=prec)&
				!	&**exp_mainte)**exp_egoisme/sum_flux_down, &
				!	&real(a%child(i)%p%flux_down,kind=prec)**(exp_egoisme*exp_mainte)/sum_flux_down2
				!	!stop
				!endif
			endif
		enddo
		
! 		keep the excess/deficit of production due to the euclidean division
		if (prod.ne.0) then
			if(nchildren>0)then
				a%reserve = min(prod,a%reserve_max)		
			
				if(a%reserve<0) then											!
					print*,"repartition error", a%reserve, prod, prod_old, tempo! DEBUG LINES																							!
					stop 														!
				endif															!
			
				prod=0
			endif
		endif
		end block
! 		Redistribute the excess/deficit of production due to nint randomly
!		if (prod.ne.0) then
!			i=floor(a%N_children*rand())+1
!			a%child(i)%p%flux_up= a%child(i)%p%flux_up+prod
!			prod=0
!		endif
		
		!print*,"ERROR repartition"
		!stop 
    end select


    if(prod<0)then
		print*,"negative prod at repartition"
		stop 
    endif
    
    return
end subroutine repartition
!============================================================================80000 0

!------------------------------------------------------------------------------------
subroutine cut_choice(a, sommet, prod, n_cut, continue_flag, flag, schem)
!------------------------------------------------------------------------------------
  use precis_mod, only                         : precisint
  use parameters_module
  use mod_tree
  implicit none
  type (branch), pointer, intent(inout)        :: a, sommet
  integer(kind=precisint), intent(in)		   :: prod
  integer(kind=precisint), intent(inout)	   :: n_cut
  integer, intent(inout)					   :: continue_flag,flag
  integer, dimension(:), intent(in)			   :: schem
  type (branch), pointer					   :: b => null(), parent_of_a

	select case(schem(8))

	!Here we decide if a branch dies or remain alive (this select case must be consistent with the select case of constr_grow)

	case(0)	! -------- Dead branches are not replaced by a bud--------------
		if (prod<0 .or. a%death_flag) then       
			if (associated(a%before%p)) then
				b => a%before%p
				parent_of_a => a%parent%p
				!print *, 'CALL CUT_BRANCH: ', a%petit_nom, &
				!& '  // PROD, MAINTENANCE, DESIR: ', prod, a%maintenance, desir
				call cut_branch(a,sommet,n_cut)
				a => b
				continue_flag=0
			else
				!print *, 'FLUX_UP_DITRIB cutting the trunk: ', &
				! '  // PROD, MAINTENANCE, DESIR: ', prod, a%maintenance
				!print *, 'Terminating program'
				continue_flag=0
				flag=1
			endif
		endif
    
    case(1)  !Case where the a bud come out at the location of the dead branch
		if (prod<0 .or. a%death_flag) then       
			if (associated(a%before%p)) then
				b => a%before%p
				parent_of_a => a%parent%p

				call cut_branch(a,sommet,n_cut)
				a => b
				continue_flag = 0
				
				if(parent_of_a%N_children<1) then                     ! If all the children are dead
					parent_of_a%nb_leaves = parent_of_a%nb_leaves + 1 ! then a leaf appear.
				endif                                                 ! I.E. No extremities can remain leaveless !!!
				
			else
				continue_flag=0
				flag=1
			endif
		endif

    case(2)  !Case where the a bud come out at the location of the dead branch by consuming a reserve of energy
		if (prod<0 .or. a%death_flag) then       
			if (associated(a%before%p)) then
				b => a%before%p
				parent_of_a => a%parent%p

				call cut_branch(a,sommet,n_cut)
				a => b
				continue_flag = 0
				
				if(parent_of_a%reserve >= cost_leaf .and. parent_of_a%N_children<1) then !
					parent_of_a%reserve = parent_of_a%reserve - cost_leaf                ! A leaf is created if there is enough reserve.
					parent_of_a%nb_leaves = parent_of_a%nb_leaves + 1                    !
				endif									                                 !
				
			else
				continue_flag=0
				flag=1
			endif
		endif
    
	case default
		print*,"error in cut choice", schem(8)
		stop 
	end select 

end subroutine cut_choice
!------------------------------------------------------------------------------------


end module mod_strat1
