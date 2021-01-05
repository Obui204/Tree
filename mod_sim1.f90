module mod_sim1

public :: flux_up_distrib, flux_down_distrib

contains

! Here, we have the flux_down and flux_up routines dictating the dynamics of the tree. 
! 
! The routines are themselves dependant on the subroutines stored in the module mod_strat1 which are
! the real ones deciding how the rules of the simulations are.
! 
! The code of the routines flux_up/down simply represents the skeleton. 
!
!
!The following comments is just a memo :
!----------------------------------------------------------------------------
! SUMMARY OF THE RELATIONS BWT SUBROUTINES									!
!																			!
!																			!
!  flux_down_distrib calls the following subroutines :						!
!			-production (mod_strat1)										!
! 			-down_grow_reser (mod_strat1)									!
!			-mainte_cal (mod_strat1)										!
!			-vol_constr (mod_strat1)										!
!           -mech_constr_cal (mod_strat1)                                   !
!																			!
!  flux_up_distrib calls the following subroutines :						!
!			-vol_constr through up_grow_mainte (mod_strat1)					!
!			-cut_branches through cut_choice (mod_tree1)					!
!			-new_branches through branch_maker (mod_tree1)					!
!			-up_grow_mainte (mod_strat1)									!
!			-cut_choice (mod_strat1)										!
!			-repartition (mod_strat1)										!
!			-branch_maker (mod_strat1)										!
!																			!
!----------------------------------------------------------------------------

!*****************************************************************************80
subroutine flux_down_distrib(schem,tree_object,hei,generation)               !80
!*****************************************************************************80
  use precis_mod, only                          : prec => working_precis, precisint
  use parameters_module
  use mod_tree
  use mod_strat1, only							: vol_constr, mainte_cal, down_grow_reser, upper_need_cal, production, mech_constr_cal
  implicit none
  integer, intent(in)                          :: hei, generation
  integer, dimension(:), intent(in)			   :: schem
  type(arbre_pointer), pointer, intent(inout)  :: tree_object 
  !type (branch), pointer, intent(inout)        :: sommet
!*****************************************************************************80
! Distributes the ressources created from the nb_leaves
  integer                                      :: desir, i, mass_weight
  integer(kind=precisint)					   :: prod!, upper_weight
  real(kind=prec),dimension(3)                 :: calc_upbarycentre
  type (branch), pointer                       :: a, bb
  
  !This is the "flux down" process: the leaves create the resource. And the branch send it down to its parent (but can keep some for itself to grow in volume)
  !then the parent do the same with thei own parent etc. The following for-loop codes this process.
  
  tree_object%tree_point%generat = generation
  a => tree_object%peak_point
 do
   	! prod = a%nb_leaves*Prod_leaf*a%distance/hei + a%reserve ! production function of the height
    
    
    prod = a%reserve
    if(prod<0)then
		print*,"Er010"
		stop
    endif
    
    !at this stage we update the knowledge of the branch ( number of leaves above it, quantity of volume above it, quantity of ressource its children are given him)
    !--------------------------------------------------------------------------------
    
    a%reserve = 0
    a%upper_leaves = a%nb_leaves
    a%upper_vol = 0
    a%upper_needs = 0
    a%upper_weight = 0
    
    calc_upbarycentre(:) = 0.0_prec  ! dummy var to calculate a%upper_barycentre

    do i=1,a%N_children
		prod = prod + a%child(i)%p%flux_down                       ! compute the ressouce given by the children
		a%upper_leaves = a%upper_leaves + a%child(i)%p%upper_leaves         ! compute desired diameter to sustain upper branches 		
		a%upper_vol = a%upper_vol + a%child(i)%p%volume + a%child(i)%p%upper_vol   ! compute the volume total above him
		
		mass_weight = a%child(i)%p%volume + a%child(i)%p%nb_leaves*leaves_weight + a%child(i)%p%upper_weight
		a%upper_weight = a%upper_weight +  mass_weight
		!a%upper_needs = a%upper_needs + a%child(i)%p%upper_needs     ! compute the total "needs" its descendant are asking for 
		
		
		calc_upbarycentre(:) = calc_upbarycentre(:) + a%child(i)%p%upbarycentre(:)*mass_weight

		if(prod<0) then																!
			print*, prod,a%child(i)%p%flux_down,i , "reserve error 2 at flux_down"	! DEBUG LINES
			stop   																	!
		endif																		!
    end do
    !-------Calculate the barycenter of the branch and its descendants=======!
    mass_weight = a%volume + a%nb_leaves*leaves_weight
    calc_upbarycentre(1) = calc_upbarycentre(1) + real(a%int_x*mass_weight,kind=prec)
    calc_upbarycentre(2) = calc_upbarycentre(2) + real(a%int_y*mass_weight,kind=prec)
    calc_upbarycentre(3) = calc_upbarycentre(3) + real(a%int_z*mass_weight,kind=prec)
    
    a%upbarycentre(:) = calc_upbarycentre(:)/real(a%upper_weight + mass_weight,kind=prec)
    !-------------------------------------===================================!
    
    call mech_constr_cal(a,schem)

!=========DIAGNOSTIC======================================
!    if(generation>38)then
!		print*, sqrt((a%upbarycentre(1)-a%int_x)**2+(a%upbarycentre(2)-a%int_y)**2), a%volume, a%distance, a%N_children
!    endif
!=========================================================
    prod = prod + production(a,tree_object,generation,schem)
    ! If the branch has some leaves the function "production" will decide how much ressource they are producing

    if(prod<0) then																			!
		print*, a%reserve,production(a,tree_object,generation,schem),a%nb_leaves, "reserve 1 error at flux_down"	! DEBUG LINES
		stop   																				!
    endif
    
    if(a%volume<1) then
		print*, "Error1"
		stop  
    endif
    
    call down_grow_reser(a,prod,schem)
    ! the routine grow_reserve decides how much the branch will use to grow (or keep in reserve) before giving the remainder to its parent
    
    if(a%volume<1 .and. .not. a%death_flag) then      !
		print*, "Error2"     ! DEBUG
		stop                 !
    endif                    !
    
    a%flux_down = prod    ! The remaining ressource will be given as flux down to the parent
    
    if(a%flux_down<0) then                      !
    	print*, a%flux_down, "flux down error"  ! DEBUG
    	stop                                    !
    endif                                       !
    
    

    !--------------------------------------------------------------------------------------

    !======Below we calculate the "upper_needs" which is the amount the branch will demand to its parents during the flux_up process=======
    
    a%maintenance = mainte_cal(a,schem)
    ! "mainte_cal" calculate the maintenance cost the branch will have to pay during winter (during the flux_up period)

    !call maintenance_down() !NOT IMPLEMENTED YET
    
    a%upper_needs = a%upper_needs + upper_need_cal(a,schem)
    !if(a%upper_needs/=a%maintenance)then
	!	print*,"huhhui",a%maintenance, a%upper_needs, vol_constr(a,schem)
	!	stop
    !endif
    !if(abs(a%flux_down/prod_leaf-nint((real(a%maintenance)/rel_mainte_vol)**2))>1)then
	!	print*,"STTTOOPPP", a%flux_down, a%maintenance, a%volume, a%upper_leaves
	!	print*, a%flux_down/prod_leaf-nint((real(a%maintenance)/rel_mainte_vol)**2), prod_leaf, rel_mainte_vol
	!	print*,a%flux_down/prod_leaf, nint(real(a%maintenance)/rel_mainte_vol)**2
	!	stop
    !endif
    ! upper need is what the branch will ask to the parent during winter (flux up process/routine)
    ! and "upper_need_cal" is the function that computes what the branch will ask its own need.
    ! Lies can be implemented in this function.
    
	if(a%death_flag)then   !
		a%maintenance = 0  ! If the branch is dead, there is no maintenance or needs
		a%upper_needs = 0  !
	endif                  !
    !=======================================================================================================================================


    !we now close the loop
    if (associated(a%before%p)) then
		a => a%before%p
    else
		exit
    end if
 end do

  return
end subroutine flux_down_distrib
!*****************************************************************************80


!*******************************************************************************************************!80
subroutine flux_up_distrib(schem,tree_object,generation,n_cre,n_cut1,n_cut2,flag,h,voltot,restot,childless)     !80
!*******************************************************************************************************!80
  use precis_mod, only                          : prec => working_precis, precisint 
  use parameters_module, only                   : extra_cost, N_c_Max, exp_egoisme, rel_mainte_vol,&
												 volum_weight, leaves_weight
  use mod_tree
  use mod_strat1, only							: branch_maker, repartition, up_grow_mainte, &
												cut_choice, mainte_cal, mech_constr_cal
  implicit none
  type(arbre_pointer), pointer, intent(inout)  :: tree_object 
  integer, intent(in)                          :: generation
  integer, dimension(:), intent(in)			   :: schem
  integer(kind=precisint) , intent(out)		   :: n_cre,n_cut1,n_cut2,childless,voltot,restot!,h_cut!,flag
  integer, intent(out)						   :: h,flag

  !type (branch), pointer				       :: tronc, sommet

! Distributes the ressources created from the nb_leaves
  integer(kind=precisint)               :: prod, temp, kk!, cost_cal
  integer                               :: counti1, h_cut!,prod_old,sum_flux_down,
  integer                               :: desir, continue_flag, counti2!,max_vol, i_maxvol,i
  !real                                 :: sum_real
  type (branch), pointer               	:: a, b, parent_of_a
  integer								:: is_before_the_parent
  is_before_the_parent = 0
  
 	a => tree_object%trunk_point
  
  counti1 = 0										!
  do												!
  	counti1 = 1 + counti1							!
    if (associated(a,tree_object%peak_point)) then	!
      exit											! DEBUG lines
    else											!
      a => a%after%p								!
    endif											!
  enddo                                             !
  counti2 = 0										!
  
  voltot = 0        ! The values of these variables are the data we
  restot = 0        ! will write in our output files :
  childless = 0	    !
  h = -1	            ! voltot = volume total de l'abre, restot = reserve totale, h = hauteur de l'arbre
  n_cre=0           ! , childless = number of extremities, n_cre = number of branches birthed this generation
  n_cut1=0          ! , n_cut1 = number of branches which died this generation
  n_cut2=0          ! , n_cut2 = nb of branches which died because of mech failure
  h_cut=0           ! will mark the average height of branches cut by starvation
  deallocate(tree_object%tree_point%branches_par_etage)                      
  allocate(tree_object%tree_point%branches_par_etage(tree_object%tree_point%hauteur+1))
  tree_object%tree_point%branches_par_etage = 0
  deallocate(tree_object%tree_point%feuilles_par_etage)                      
  allocate(tree_object%tree_point%feuilles_par_etage(tree_object%tree_point%hauteur+1))
  tree_object%tree_point%feuilles_par_etage = 0

  flag=0
  a => tree_object%trunk_point
  b => null()		! pointer b is used to skip newly created children (cf lines of code below)
  
  a%flux_up = a%flux_down
  !write(*,*) 'Flux rebound: ', a%flux_up, a%volume
 do
	
	!INITIALIZE DEFAULT FLAG STATUS
	is_before_the_parent = 0
	if(associated(a%before%p,a%parent%p)) is_before_the_parent = 1
	!death_flag = 0
	continue_flag = 1
	!_____________________________
	
	if(a%maintenance.ne.mainte_cal(a,schem)) then												        !-------------!
		print *, "pb up/down maintenance", a%maintenance, a%volume, mainte_cal(a,schem), counti2&	    !			  !
		&, a%N_children, a%distance, a%int_x, a%int_y, a%parent%p%N_children, generation, a%generation  ! DEBUG lines !
		stop  																					     	!			  !
	endif																							    !-------------!

	prod = a%flux_up + a%reserve
    a%reserve = 0	
	if(.not. a%death_flag) call up_grow_mainte(a,prod,schem)

    !if(generation>126)then
	!	if(associated(a%parent%p))then
	!		if(a%parent%p%N_children>1)then
	!			print*, a%distance, a%int_z, a%upbarycentre(1)-a%int_x, a%upbarycentre(2)-a%int_y	
	!		endif
	!	endif
    !endif
    
	block
	logical  :: temporary_flag
	integer  :: previous_ncut1!, namae, par
	integer  :: h_cut_dummy
	previous_ncut1 = 0
	temporary_flag = .false.
	if(a%death_flag) then
			temporary_flag = .true.
			previous_ncut1 = n_cut1
	endif
	
	!----scheme for calculating h_cut the average height of starved branches----!
	block
	type(branch), pointer :: ccb
	ccb => a
	h_cut_dummy = ccb%int_z
	do
		if (associated(ccb%after%p)) then
			ccb => ccb%after%p
		else
			exit
		end if
		if(ccb%distance>a%distance)then
			h_cut_dummy = h_cut_dummy + ccb%int_z
		else
			exit
		endif
	enddo
	end block
	!---------------------------------------------------------------------------!
	!if(prod<0)then
	!	print*,a%flux_up
	!	print*,a%maintenance
	!endif
 	call cut_choice(a,tree_object%peak_point,prod,n_cut1,continue_flag,flag,schem)
	! and gives as output continue_flag (indicating if we did kill the branch)
	! and flag (indicating whether it is the trunk we have cut)
	if(temporary_flag)then
		n_cut2 = n_cut1 - previous_ncut1 + n_cut2 
!	 print*, n_cut2, n_cut1
	else if(continue_flag==0)then
		h_cut = h_cut + h_cut_dummy
	endif
	end block
	if(flag==1) return
	
	b => a		! pointer b is used to skip any newly created branches from the flux_up process (cf end of the loop)
	!********************************************************************
	
	if(continue_flag==1) then
		!-------------------------------!
		counti2=counti2+1				! DEBUGG line
		!-------------------------------!
		!call reserve_bef_repart(   (to be implemented)
		tree_object%tree_point%branches_par_etage(a%distance) =&
		& tree_object%tree_point%branches_par_etage(a%distance) + 1
		! Now we distribute the flux up to children or/and keep some in reserve
		if (a%N_children>0) then
			call repartition(a,prod,schem,generation)
		endif
		! End of distribution
		desir=0
		! Creation of children : here desir = nb of new branches it wants to make
		!===================================================================
		call branch_maker(a,tree_object,prod,desir,generation,schem,temp)
		!Inside branch_maker we should have something like: prod = prod - desir*cost_of_a_single_branch 
		! call leaves_maker(a,tree_object, (to be implemented)
		
		a%reserve = a%reserve + prod

		a%reserve=min(a%reserve,a%reserve_max) 
		
		!---STORING INFO CONCERNING THE TREE (ITS TOTAL VOLUME, HEIGHT, NUMBER OF BRANCHES ETC.) so that we can write them on some output files---!
		voltot = voltot + temp				! temp is the volume of the newly created children !
		voltot = voltot + a%volume
		restot = restot + a%reserve			! REMARK : DO NOT GIVE ANY DEFAULT RESERVE TO THE CHILDREN
		childless = childless + desir
		if(a%N_children==0) then
			childless = childless + 1
			tree_object%tree_point%feuilles_par_etage(a%distance) =&
			& tree_object%tree_point%feuilles_par_etage(a%distance) + 1
		endif
		n_cre = desir + n_cre
		h = max(h,a%distance)
		
		!--------------------------------------------------------------------------------------------------------!
		
		if(desir>0) then
			!h = max(h,a%distance + 1)
			h = max(h,a%child(a%N_children)%p%distance)                                           !
			tree_object%tree_point%branches_par_etage(a%distance+1) = &       !
			& tree_object%tree_point%branches_par_etage(a%distance+1) + desir ! We determine the height of the tree
			                                                                  ! and the number of branches of a given height
			tree_object%tree_point%feuilles_par_etage(a%distance+1) = &       !
			& tree_object%tree_point%feuilles_par_etage(a%distance+1) + desir !
			temp = a%N_children		! With b we skip the children created                         
			b => a%child(temp)%p	! 
		endif
        !=====================================================================
	else if(continue_flag==0) then
		!call emergency_bud(a,schem  (To be implemented)
		if(a%N_children==0 .and. is_before_the_parent==1) childless = childless + 1
		
	endif

	a => b   ! With b we skip the children created from the for-loop
    !we now close the loop
    if (associated(a%after%p)) then
		a => a%after%p
    else
		exit
    end if
 end do
	
	if(counti1-n_cut1 /= counti2) then	    ! -----------------
		print*,"error at counti"            ! 
		print*,counti1,counti2, n_cut1       ! DEBUG
		stop                                !
	endif                                   ! -----------------
  
  if(n_cut1-n_cut2>0)then
	write(unit=20,fmt="(i7,f19.10)") generation, real(h_cut)/real(n_cut1-n_cut2)
  endif
!  a => tree_object%trunk_point                      !
!  counti1 = 0										!
!  do												!
!  	counti1 = 1 + counti1							!
!    if (associated(a,tree_object%peak_point)) then	!
!      exit											! DEBUG lines
!    else											!
!      a => a%after%p								!
!    endif											!
!  enddo                                             !
  
!  	if(counti1 /= counti2 + n_cre) then	    ! -----------------
!		print*,"error at counti2"           ! 
!		print*,counti1, counti2, n_cre      ! DEBUG
!		stop                                !
!	endif                                   ! -----------------
 	!-----------------------------------------------------
 	!NOW WE UPDATE THE DATA OF OUR TREE
 	tree_object%tree_point%hauteur = h
 	block
 	integer :: iii
 	integer, dimension(tree_object%tree_point%hauteur) :: array_copy
 	
 	do iii=h+1, size(tree_object%tree_point%branches_par_etage)    !
		if(tree_object%tree_point%branches_par_etage(iii)/=0)then  ! DEBUG
			print*,"error end flux_up"                             !
			stop                                                   !
		endif                                                      !
 	enddo
 	do iii=1,h
		array_copy(iii) = tree_object%tree_point%branches_par_etage(iii)
 	enddo
 	deallocate(tree_object%tree_point%branches_par_etage)
 	allocate(tree_object%tree_point%branches_par_etage(h))
 	tree_object%tree_point%branches_par_etage = array_copy
 
 	endblock
 	
  return
end subroutine flux_up_distrib
!*******************************************************************************************************!80


!*****************************************************************************80
subroutine print_terminal(h,generation,n,v)
  use precis_mod, only                          : prec => working_precis, precisint
  use mod_tree
  implicit none
  integer, intent(in)                          :: generation, h
  integer(kind=precisint),intent(in)           :: n, v

		write(*,*) "Generation: ", generation, " / height: ", h, " / # branches: ", n, " / total volume: ", v
		
end subroutine print_terminal
!*****************************************************************************80

end module mod_sim1
