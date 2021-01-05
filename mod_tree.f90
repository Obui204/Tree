!--------------------------------------------------------------------------!
module mod_tree    !  In this module all tree characteristics are defined  !
  !----------------!-------------------------------------------------------!
  !                                                                        !
  !  We try to make a tree with integers numbers characteristics           !
  !                                                                        !
  !------------------------------------------------------------------------!
  
      !       !==================================================================================================================!      
     ! !      !REMARK : A TREE WITH A TRUNK ALONE IS CONSIDERED TO BE OF HEIGHT 1 !!!                                            !     
    ! ! !     !                                                                                                                  !    
   !  !  !    !                                                                                                                  !   
  !   !   !   !SO THE TRUE " HEIGHT " OF THE TREE (DEPENDING ON HOW YOU DEFINE IT) MAY ACTUALLY BE  " true_height = hauteur - 1 "!   
 !!!!!!!!!!!  !==================================================================================================================!  
  
  use precis_mod, only                         : precisint, prec => working_precis
  use parameters_module
  
!  type seed
!	integer                                   :: x_seed
!	integer                                   :: y_seed
!  end type seed
  
!  type seed_pointer
!	type(seed), pointer                      :: p => null()
!  end type seed_pointer
  
  type great_space
	
	logical, dimension(:,:,:), allocatable       :: xyz_coord != 0
	
	logical, dimension(:,:), allocatable         :: seedbed
!	type(seed_pointer), dimension(:), allocatable   :: list_of_seed
	
  end type great_space

  !Gives main caracteristics for each tree
  type arbre
	integer(kind=precisint)	                  :: nbre_branches = 1
	integer(kind=precisint)	                  :: nbre_extremite = 1    
	integer(kind=precisint)                   :: nbre_feuilles = -1    ! yet to be implemented
	integer	                                  :: hauteur = 1
	integer                                   :: generat = 1              ! the current time 
	integer(kind=precisint),allocatable	      :: branches_par_etage(:) 
	integer(kind=precisint),allocatable       :: feuilles_par_etage(:) ! yet to be implemented inside the program
	type(great_space), pointer                :: the_space => null()
  end type arbre
  
! Arbre pointer: Pointe vers un arbre (utile si simulation de foret)
  type arbre_pointer
    type (arbre), pointer                     :: tree_point => null()   ! Pointer pointing the variable "arbre" corresponding to our tree
    type (branch), pointer                    :: trunk_point => null()  ! Pointer for the first item (the trunk) of the linked list that lists the branches of our tree
    type (branch), pointer                    :: peak_point => null()   ! Pointer for the last item of the linked list that lists the branches of our tree
  end type arbre_pointer 
  
  type branch_pointer
    type (branch), pointer                    :: p => null()
  end type branch_pointer
  
  type branch                             	  ! Definition of a new type: branch
    integer                                   :: generation = 0      			! generation is the time at which the branch was created
    integer                                   :: nb_leaves = -5				    ! We put some leaves at the end of a branch
    integer                                   :: distance = 1        			! Distance of the branch to the trunk (the distance of the trunk to itself is set to 1)
    !integer                                   :: extremal_distance = 0          ! Distance with his furthest descendant
    integer(kind=precisint)	                  :: volume = 1          		    ! Size of the branch
    integer(kind=precisint)                   :: upper_leaves = 1    			! number of leaves above this branch
    integer(kind=precisint)                   :: upper_needs = 0     			! Sum of the wishes of the descendence
    integer                                   :: reserve = 0         			!-----------------------! Reserve of "food" kept in the branch (set at 0 after discovering a race that left the tree alway growing)
    integer                                   :: reserve_max = -5     			!  MORE OR LESS UNUSED  ! Maximal size of the reserve
    integer                                   :: reserve_wish = -5   			!-----------------------! What we wish to achieve as a reserve
    integer(kind=precisint)                   :: maintenance = -10     			! Yearly cost to survive as is (set default as zero if you do not want to create a BUG)
    integer(kind=precisint)                   :: flux_down = 0       			! Flux coming from children (Spring period)
    integer(kind=precisint)                   :: flux_up = 0         			! Flux going to Children (Winter and maintenance period)
    !integer								  :: flux_max = -5					! Maximum flux it can handle
    integer                                   :: N_children = 0      			! Number of children
    integer                                   :: petit_nom
    !real(kind=prec)						  :: health=1.0_prec				! NOT USED YET
    logical                                   :: death_flag = .false.
    
    integer                                   :: int_x = 0        !Coordonnées sur l'espace discret
    integer                                   :: int_y = 0        !
    integer                                   :: int_z = 1        !integer that indicates whether the branch's position is ill-defined or not (0 means ill-defined, 1 means defined)
        
    real(kind=prec)                           :: mech_constr = 0.0_prec	
    integer(kind=precisint)                   :: upper_vol = 0					! quantity of volume of "wood" above the branch (NOT INCLUDING ITS OWN VOLUME)
    integer(kind=precisint)                   :: upper_weight = 0               ! weight above the branch (not including its own	
    real(kind=prec),dimension(3)              :: upbarycentre = 0.0_prec        ! barycenter of the branch and all its descendants
    
    real(kind=prec), dimension(3)             :: spacebabywish_z                ! Dertermine whether how much a branch wants to make babies above it or not
                                                                                ! by default spacebabywish_z is entirely determined by the parameter ran_posit
!--------MAY BE DELETED-------------------------------------------------------------------------------------------------------------------------------!
    real                                      :: x = 0.0             ! x position in the plane of the extremity of the branch (used for drawing tree) !
    real                                      :: y = 0.0             ! y position in the plane of the extremity of the branch (used for drawing tree) !
    real                                      :: z = 1.0             ! z position in the plane of the extremity of the branch (used for drawing tree) !
    real                                      :: theta = 0.0         ! inclination of the branch (between -pi/2 and pi/2)  (used for drawing tree)    !
    real                                      :: phi = 0.0           ! inclination of the branch (between -pi/2 and pi/2)  (used for drawing tree)    !
!------------------------------------------------------------------------------------------------------------------------------------------------------
    type (branch_pointer)                     :: parent              ! Single parent of the branch
    type (branch_pointer), allocatable        :: child(:)		     ! Children's address
    type (branch_pointer)                     :: before              ! Place in the tree chain
    type (branch_pointer)                     :: after               ! Place in the tree chain
    !type (point2D), pointer                   :: spacial_location => null()   ! pointer toward its location in the list of points
    type (arbre_pointer), pointer             :: owner => null()
  end type branch
  

public :: new_branches, cut_branch, branch, height, spacial_recog, opposite_side_of_3x3,&
		alloc_branch, deall, alloc_arbre, alloc_arbre_pointer, spacial_deallocate_branch,&
		make_spacial_branch, apical_lead_test, column_generator, shading_at_equador, shading_northsouth


!*******************************************************************************!
CONTAINS  !*********************************************************************!
!*******************************************************************************!

!*********************************************************************
!routine to allocate a branch
subroutine alloc_branch(point, tree_object)
  use parameters_module, only               : rw, init_nb_leaves, prod_leaf, res_prod_rel,&
											flux_res_rel, N_c_Max
  implicit none
	type(branch), intent(inout)				         :: point
	type(arbre_pointer), pointer, intent(inout)      :: tree_object
	!integer, intent(in)					:: init_vol
	allocate(point%child(2*N_c_Max))
	point%reserve_wish = rw
	point%nb_leaves = init_nb_leaves
	point%reserve_max = nint(res_prod_rel*real(prod_leaf))
	if(point%reserve_max<0) then
		print*,"error res_max negative"
		stop
	endif
	
	point%spacebabywish_z(1) = 1.0_prec
	point%spacebabywish_z(2) = ran_posit
	point%spacebabywish_z(3) = ran_posit*ran_posit
	
	point%owner => tree_object
	!point%flux_max = nint(flux_res_rel*real(point%reserve_max))
	!point%volume = init_vol
end subroutine alloc_branch
!*********************************************************************

!*********************************************************************
!routine to allocate an arbre(type)
subroutine alloc_arbre(point)
  ! Here we allocate the "arbre" variable corresponding to a tree composed of 
  ! a single trunk (default situation for a tree)
  use parameters_module, only               : init_nb_leaves
  implicit none
	type(arbre), intent(inout)				:: point
	
	point%nbre_feuilles = init_nb_leaves
	allocate(point%feuilles_par_etage(1))
	allocate(point%branches_par_etage(1))
	
end subroutine alloc_arbre
!*********************************************************************

!*********************************************************************
!routine to allocate/initiate a tree
subroutine alloc_arbre_pointer(point)
  implicit none
	type(arbre_pointer), pointer, intent(inout)  :: point
	
	allocate(point%trunk_point)
	call alloc_branch(point%trunk_point, point)
	point%trunk_point%reserve = rel_res_init*extra_cost
	point%trunk_point%volume = 1
	point%trunk_point%distance = 1
	!tronc%spatial_def = 1
	
	point%peak_point => point%trunk_point
	
	allocate(point%tree_point)
	call alloc_arbre(point%tree_point)
	
end subroutine alloc_arbre_pointer
!*********************************************************************

!*********************************************************************
function opposite_side_of_3x3(entier)  result (entier2)
!*********************************************************************
 implicit none
	integer, intent(in)   :: entier
	integer               :: entier2
    if(entier==1) entier2=3
    if(entier==3) entier2=1
end function opposite_side_of_3x3
!*********************************************************************

!*********************************************************************
subroutine spacial_recog(a,desir2, space_around, space0)
!********************************************************************
 implicit none
	!integer, intent(in)                       :: cas  ! Used in the select case(cas)
	type(great_space),          intent(in)    :: space0
	integer, intent(out)                      :: desir2
	logical, dimension(:,:),intent(out)       :: space_around
	type(branch), pointer, intent(in)         :: a
	integer                                   :: i, j, k, i1, j1, k1, dum_x, dum_y, dum_z, ent_dum1, ent_dum2
	
 ! integer :: crea indicates whether the plane above has is unoccupied (=0 for unocc, >0 for occ)
 ! desir2 is the number of location available
 ! space_around is a matix representation of the space above the branch we're dealing with
  
  	
	if(size(space_around,1)/=9)then
		print*,"er890"
		stop
		if(size(space_around,2)/=3)then
			print*,"errr1820"
			stop
		endif
	endif
	
	desir2 = 27   ! By default : a branch can create children in 26 directions. desir2 represet the nb of directions still available
	space_around(:,:) = .false. ! space_around is a spacialization of the available directions:
		! Indexation of the space_around array
	    ! Plane above:
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
	dum_z = a%int_z + 1
	dum_x = a%int_x - 1
	dum_y = a%int_y - 1
	
	!--------------------------------
	if(dum_z>N_z)then
		print*,"too low N_z"
		stop
	endif
	
	if(dum_y>N_xy-1)then
		print*,"too low N_y"
		stop
	endif
	
	if(dum_x>N_xy-1)then
		print*,"too low N_x"
		stop
	endif
	!--------------------------------
	
	do k = 0,2
		k1 = dum_z - k
		if(k1<1)then
			desir2 = desir2 - 9
			if(k/=2)then
				print*,"er9420"
				stop
			endif
			do ent_dum1 = 1, 3
				do ent_dum2 = 1, 3
					space_around(ent_dum1 + 3*k, ent_dum2) = .true.
				enddo
			enddo
			if(desir2 < 0)then
				print*, "ermodtree89"
				stop
			endif
			exit
		endif
		do i1 = dum_x, dum_x +2
			i = i1 - dum_x + 1 + 3*k
		
			do j1=dum_y, dum_y + 2
				j = j1 - dum_y + 1
			
				if(space0%xyz_coord(i1,j1,k1))then
					space_around(i,j) = .true.
					desir2 = desir2 - 1
				!else! if(cas==2) then
				!	if(k==2 .or. k==1)then
				!		space_around(i,j) = .true.
				!		desir2 = desir2 - 1
				!	endif
				endif
				
				
			enddo
		enddo
		
	enddo
	
	if(.NOT. space_around(5,2))then
		print*, "errrr092"
		stop
	endif
	
	
end subroutine spacial_recog
!*********************************************************************

!*********************************************************************
subroutine apical_lead_test(a,apical_lead)
!*********************************************************************
	implicit none
	type(branch), pointer, intent(in) :: a
	logical, intent(out)              :: apical_lead
	type(branch), pointer             :: bbb
	integer                           :: trunk_x, trunk_y, z_untiltrunk_debug
		
	!we need to verify the branch is on the apical lead => false or true
	
	!TWO CASES:  CASE 1 = Space exists : the apical lead is the vertical column springing from the root/trunk 
	!            CASE 2 = No space : the apical lead is the 1st child down to the root/trunk
	if(associated(a%owner%tree_point%the_space))then
		!----DEBUG--------------------
		!if(schem(1)<=1)then
		!	print*,"1dez3vv67hhbN"
		!	stop
		!endif
		!----DEBUG--------------------
		apical_lead = .false.
		trunk_x = a%owner%trunk_point%int_x
		trunk_y = a%owner%trunk_point%int_y
		z_untiltrunk_debug = a%int_z
			
		!============
		!we basically verify if the branch has the same xy coord as the trunk
		! as well as verify if the parent, grand parent etc. all have this same xy-coord
		! I.E. the apical lead must be a vertical column growing from the initial seed xy-coord
		if(abs(trunk_x-a%int_x)+abs(trunk_y-a%int_y)==0)then                              
			block
			type (branch), pointer :: cc
			cc => a                                                                           
			loop_name : do                                                                
				if(associated(cc%parent%p))then 
					cc => cc%parent%p
					if(abs(cc%int_x-trunk_x)+abs(cc%int_y-trunk_y)/=0) then 
						exit loop_name
					else 
						z_untiltrunk_debug = z_untiltrunk_debug - 1
					endif                                                                 
				else
					if(associated(a%owner%trunk_point,cc))then
						apical_lead = .true.
						exit loop_name
					else
						print*, "bug a branch without parent that isn't the trunk exist!!!"
						stop
					endif
				endif
			enddo loop_name
			end block
			!---------------DEBUGGING--------------------------
			if(apical_lead .and. z_untiltrunk_debug/=1)then
				print*, "bug in apical lead 1"
				stop
			endif
			if(.not. apical_lead)then
				if(z_untiltrunk_debug<2)then
					print*, "bug in apical lead 2"
					stop
				endif
			endif
			!-------------DEBUGGING--------------------------
		endif
	else
		!apical_lead = .true.
		bbb => a
		do
			if(associated(bbb%parent%p))then
				if(associated(bbb%parent%p%child(1)%p,bbb))then
					bbb => bbb%parent%p
				else
					apical_lead = .false.
					exit
				endif
			else
				apical_lead = .true.
				if(bbb%distance/=1)then
					print*, "error3411100033"
					stop
				endif
				exit
			endif
		enddo
	endif
endsubroutine apical_lead_test

!***********************************************************************************
subroutine make_spacial_branch(a, space00, list_point, Nold_child, desir, space_around)
!***********************************************************************************
  use precis_mod, only                      : prec => working_precis, precisint	
  use parameters_module
  		
	implicit none
	type (branch), pointer, intent(inout)        :: a
	type(great_space), pointer, intent(inout)    :: space00
	integer, intent(in)                          :: Nold_child, desir
	logical, dimension(:,:),intent(in)           :: space_around  

	integer, dimension(:,:), intent(in)          :: list_point
		
	! a : the branch
	! space00: the space
	! list_point : the coordinates of all points we want to create
	! 
	!
	if(desir<=0)then
		print*, "error 4041"
		stop
	endif
	

		!
		! LE CAS OU L ARBRE PEUT FAIRE POUSSER DES
		! ENFANTS SUR LE COTE OU EN DESSOUS
		!
		block
		integer :: dum_z, dum_x, dum_y, x0, y0, z0, i
		
		x0 = a%int_x
		y0 = a%int_y
		z0 = a%int_z
		!===========The first term of the do-loop is done manually===========!
		dum_x = list_point(1,desir+1-1)
		dum_y = list_point(2,desir+1-1) 
		dum_z = list_point(3,desir+1-1) 

		a%child(Nold_child+1)%p%int_x = dum_x   
		a%child(Nold_child+1)%p%int_y = dum_y  
		a%child(Nold_child+1)%p%int_z = dum_z

		if(space00%xyz_coord(dum_x, dum_y, dum_z))then  !
			print*, "eerrferccf"                           ! 
			print*,space00%xyz_coord(dum_x, dum_y, dum_z)  ! DEBUG
			stop                                           !
		endif                                              !

		space00%xyz_coord(dum_x, dum_y, dum_z) = .true.

		!==================================================================!
		
                                                   !  By default when a branch is created its volume is 1
		do i=2,desir
			dum_x = list_point(1,desir+1-i)
			dum_y = list_point(2,desir+1-i)
			dum_z = list_point(3,desir+1-i)
			a%child(Nold_child+i)%p%int_x = dum_x   ! Give the correct position
			a%child(Nold_child+i)%p%int_y = dum_y   !
			a%child(Nold_child+i)%p%int_z = dum_z   !

			if(space00%xyz_coord(dum_x, dum_y, dum_z))then   !
				print*, "eerrfterf"                             !
				print*,space00%xyz_coord(dum_x, dum_y, dum_z),i ! DEBUG
				stop                                            !
			endif                                               !

			space00%xyz_coord(dum_x, dum_y, dum_z) = .true.
			
			if(.not. space_around(dum_x-x0+2 + 3*(z0-dum_z+1), dum_y-y0+2))then
				print*, "immpossiddble!!!", dum_z
				print*, space_around
				stop
			endif
			
		enddo 
		end block 
!***********************************************************************************
end subroutine make_spacial_branch
!***********************************************************************************

!**********************************************************************************
subroutine shading_at_equador(a_x, a_y, a_z, angle, full_space, answer)
!***********************************************************************************
  use precis_mod, only                      : prec => working_precis, precisint	
  use parameters_module
 implicit none

	type(great_space),  intent(in)            :: full_space
	!type(branch), pointer, intent(in)         :: a
	integer, intent(in)                       :: a_x, a_y, a_z
	real(kind=prec), intent(in)               :: angle
	integer, intent(out)                      :: answer
	integer                                   :: i, j, k, run_x, run_z
	!real(kind=prec)                           :: radians
	
	answer = 0
	
	! Tuto :
	! 'angle' doesn't mean the 'incidence angle' of the sunrays
	!  It's the angle of the rays with the x-semi-axis (0;+infty) 
	! 
	! angle = 90 is when the incidence angle is 0 deg
	! angle = 45 and 135 are when the runrays are hitting perfectly diagonally
	! 
	! The other cases are a bit more complicated, and as such are treated in their own do-loop
	!
	if(abs(angle-90.0_prec)<0.001)then!-------------------------! Angle =90
		run_z = a_z + 1
		
		do
			if(run_z>=N_z-1)then  !
				exit              ! When we are out of bound of the 'full_space' we stop
			endif                 !
			
			if(full_space%xyz_coord(a_x, a_y, run_z))then
				answer = answer + 1	
			endif
			
			run_z = run_z + 1
		enddo
	
	else if(abs(abs(angle-90.0_prec)-45.0_prec)<0.001)then!-----! Angle = 45 or 135
		block
		integer :: increments
			if(angle < 90.0_prec)then    ! Angle = 45 deg
				increments = 1
			else !-----------------------! Angle = 135 deg
				increments = -1
			endif
			run_x = a_x + increments
			run_z = a_z + 1
			
			do
				if(run_z>=N_z-1)then
					exit
				endif
				if(abs(run_x)>=N_xy-1)then
					exit
				endif			
			
				if(full_space%xyz_coord(run_x, a_y, run_z))then
					answer = answer + 1	
				endif
			
				run_z = run_z + 1
				run_x = run_x + increments
			enddo
		end block
	
	else !----------------------------------------------------! Angle complicated
		block
		integer:: increments1, counting
		!vvlogical :: vertical_flag
		real(kind=prec)  :: radians, decimals_counting, decimals
			radians = angle*atan(1.0_prec)/45.0_prec
			if(angle < 90.0_prec)then    ! Angle < 90 deg
				increments1 = 1
			else !-----------------------! Angle > 90 deg
				increments1 = -1
				radians = 4.0_prec*atan(1.0_prec) - radians
			endif
			
			run_x = a_x
			run_z = a_z
			if(abs(angle-90.0_prec)<45.0_prec)then
				decimals = tan(radians)  
				counting = 0
				decimals_counting = 0.000001_prec
				decimals_counting = decimals_counting + decimals
				do	
					if(counting>floor(decimals_counting))then
						print*,"errorjo43178", counting, radians, angle
						stop
					endif
					
					if(counting==floor(decimals_counting))then
						decimals_counting = decimals_counting + decimals
						counting = counting + 1
						run_z = run_z + 1
						run_x = run_x + increments1
					else
						counting = counting + 1
						run_z = run_z + 1
					endif
					
					if(run_z>=N_z-1)then
						exit
					endif
					if(abs(run_x)>=N_xy-1)then
						exit
					endif	
					
					if(full_space%xyz_coord(run_x, a_y, run_z))then
						answer = answer + 1	
					endif
					
					
				enddo
				
			else
				counting = 0
				decimals_counting = 0.000001_prec
				if(abs(tan(radians))<decimals_counting) then
					decimals = 1.0_prec/decimals_counting
				else
					decimals = 1.0_prec/tan(radians)  !Here we take the cotan
				endif
				decimals_counting = decimals_counting + decimals
				do	
					if(counting>floor(decimals_counting))then
						print*,"errorde1234", run_x, run_z, counting, radians
						stop
					endif
					
					if(counting==floor(decimals_counting))then
						decimals_counting = decimals_counting + decimals
						counting = counting + 1
						run_z = run_z + 1
						run_x = run_x + increments1
					else
						counting = counting + 1
						run_x = run_x + increments1
					endif
					
					if(run_z>=N_z-1)then
						exit
					endif
					if(abs(run_x)>=N_xy-1)then
						exit
					endif	
					
					if(full_space%xyz_coord(run_x, a_y, run_z))then
						answer = answer + 1	
					endif
					
					
				enddo
				
				
			endif
			
		end block
	endif

!***********************************************************************************
end subroutine shading_at_equador
!**********************************************************************************

!**********************************************************************************
subroutine shading_northsouth(a_x, a_y, a_z, northsouth, full_space, answer)
!***********************************************************************************
  use precis_mod, only                      : prec => working_precis, precisint	
  use parameters_module
 implicit none

	type(great_space),  intent(in)            :: full_space
	!type(branch), pointer, intent(in)         :: a
	integer, intent(in)                       :: a_x, a_y, a_z
	logical, intent(in)                       :: northsouth
	integer, intent(out)                      :: answer
	integer                                   :: run_y, increm!, run_z
	!real(kind=prec)                           :: radians
	
	answer = 0
	
	! Tuto :
	! 'northsouth' = .true. means the rays are from the north
	!                .false. means from the south
	! 
	!
	if(northsouth)then    ! FROM NORTH
		increm = 1
	else !-----------------------! SOUTH
		increm = -1
	endif
		run_y = a_y + increm
		!run_z = a%int_z + 1
			
		do
		!	if(run_z>=N_z-1)then
		!		exit
		!	endif
			if(abs(run_y)>=N_xy-1)then
				exit
			endif			
			
			if(full_space%xyz_coord(a_x, run_y, a_z))then
				answer = answer + 1	
			endif
			
		!	run_z = run_z + 1
			run_y = run_y + increm
		enddo	

!***********************************************************************************
end subroutine shading_northsouth
!**********************************************************************************

!**********************************************************************************
subroutine column_generator(full_space, init_x, init_y, init_z, halfwidth, halflenght, height)
!***********************************************************************************
  use precis_mod, only                      : prec => working_precis, precisint	
  use parameters_module
 implicit none

	type(great_space),  intent(inout)           :: full_space
	integer,  intent(in)                      :: init_x, init_y, init_z, halfwidth, halflenght, height
	integer                                   :: i, j, k
	!real(kind=prec)                           :: radians
	
	!width = x-axis
	!lenght = y-axis
	
	do k = init_z, init_z + height
		if(k>N_z-2)then
			print*,"er12Z3"
			stop
		endif
		do i = init_x - halfwidth, init_x + halfwidth
			if(i>N_xy-2)then
				print*,"er13324R2Z3"
				stop
			endif
			do j = init_y - halflenght, init_y + halflenght
				if(j>N_xy-2)then
					print*,"er1333244RE23324R2Z3"
					stop
				endif
				
				full_space%xyz_coord(i,j,k) = .true.
				
			enddo
		enddo
	enddo
!***********************************************************************************
end subroutine column_generator
!**********************************************************************************


!******************************************************************************************!
!! BE CAREFUL THE BRANCH CREATED MAY BE IMCOMPLETE, IT IS LACKING THE FOOLOWING PROPERTY : !
!!																	                       !
subroutine new_branches(mother,Nb,generation,sommet) ! Creates Nb new branches in the tree !
!******************************************************************************************!
  use precis_mod, only                                : prec => working_precis
  use parameters_module, only                         : N_c_Max
!*******************************************************************************
  implicit none

  type (branch), pointer, intent(inout)              :: mother     !  parent branch
  integer, intent(in)                                :: generation
  type (branch), pointer, intent(inout)              :: sommet     ! is this part of  the top of the tree
  integer, intent(in)                                :: Nb         ! number of branches to create
  !integer, dimension(:), intent(in)					 :: schem
  !integer, intent(in)								 :: init_vol
!*****************************************************************************80
  !
  ! NEW_BRANCHES creates Nb new children branch with a leaf from a parent
  ! Nb needs to be equal r larger than 1
  !
  type (branch_pointer), dimension(Nb)               :: new
  type (branch_pointer), dimension(:), allocatable   :: old_child
  type (branch), pointer                             :: end_c => null()
  real(kind=prec), dimension(Nb)                     :: r
  real(kind=prec)                                    :: pi
  integer                                            :: N, i, Nf

  N = mother%N_children
  Nf = N + Nb
  pi=acos(-1.0_prec)
  
  ! DEBUG MODE
  ! test if it is feasible
  !
  if(Nb<=0)then
	print*, "error in branch maker: calling it when you don't make children"
	stop
  endif
  
  if ((N+Nb)>N_c_Max) then
    block
    integer :: counting
    counting = 0
    do i=1,N
       if(.not. mother%child(i)%p%death_flag) counting = counting + 1
    enddo
	if((counting+Nb)>N_c_max)then
	   print *, "Problem a branch is trying to have too many children"
       print *, "Terminating the tree program"
       stop
    endif
    endblock 
  endif
  !
  !
  
  !=============================================================================================================
  !Here we make so that mother%child(:) is now of size Nf !!!
  !PROBLEM : mother%child(:) may already contain branches !!!
  !Therefore we need to save its content in old_child(:) then deallocate then reallocate it then re-copy the data
  !
!  if (N > 0) then
!		if (.NOT. allocated(mother%child)) then                                         !
!			print*, "new branches: a branch with children has unallocated child-array"  ! DEBUG
!			stop                                                                        !
!		endif                                                                           !
!		
!		allocate(old_child(N))
!		do i=1,N
!			old_child(i)%p => mother%child(i)%p 
!			mother%child(i)%p  => null()
!		enddo
!		deallocate(mother%child)
!		allocate(mother%child(Nf))
!		do i=1,N
!			mother%child(i)%p => old_child(i)%p
!			old_child(i)%p  => null()
!		enddo
!		deallocate(old_child)
!  else
!		if(N<0)then                                     !
!			print*, "new branches: negative N_children" ! DEBUG
!			stop                                        !
!		endif                                           !
!		
!		allocate(mother%child(Nf))
!  endif
  !=============================================================================================================
  
  ! If we are here, we can create the children without worrying
  ! Now let's create the children, we wil put them after the mother in the chain
  
  ! Lets put some parameters
  call random_number(r)
  do i = 1,Nb
     !print *, ':0', N, Nb, i
		allocate(new(i)%p)
		call alloc_branch(new(i)%p, mother%owner)
    
		new(i)%p%parent%p => mother
		new(i)%p%generation = generation
		new(i)%p%distance = mother%distance + 1
		if (Nb==1) then
			new(i)%p%theta = mother%theta
			new(i)%p%phi = mother%phi
		else
        new(i)%p%theta = mother%theta+(real(i-1)/real(Nb-1)-0.5)*pi/2.0**(mother%distance-1)
 !       new(i)%p%theta = mother%theta+(real(i-1)/(Nb-1)-0.5)*2.0*pi/(sqrt(real(mother%distance)))
        new(i)%p%phi = abs(mother%phi+(real(i-1)/real(Nb-1)-0.5)*pi/(2.0*mother%distance))
 
		endif
		
		new(i)%p%x = mother%x+sin(new(i)%p%theta)*sin(new(i)%p%phi)
		new(i)%p%y = mother%y+cos(new(i)%p%theta)*sin(new(i)%p%phi)
		new(i)%p%z = mother%z+cos(new(i)%p%phi)    
		new(i)%p%petit_nom = floor(r(i)*1000)
		
		mother%child(N+i)%p  => new(i)%p
  end do

  ! Finally we need to reloop the chain
  new(Nb)%p%after%p => mother%after%p ! end of the chain relooped
  mother%after%p  => new(1)%p ! beginning of the chain relooped
  new(1)%p%before%p => mother ! some before stuff

  if (associated(mother,sommet)) then
    sommet => new(Nb)%p
  else
    new(Nb)%p%after%p%before%p => new(Nb)%p ! more before stuff
  end if

  ! now we need to link everything
  i=1
  do
    if (i>Nb-1) exit
    new(i)%p%after%p => new(i+1)%p
    new(i+1)%p%before%p => new(i)%p
    i=i+1
  end do

 ! mother%nb_leaves = 0  ! we remove the leaf when creating branches
  mother%N_children = Nf

end subroutine new_branches
!******************************************************************************************!

!*****************************************************************************80
subroutine cut_branch(start,sommet,n_cut)
!*****************************************************************************80
  use precis_mod, only                     : precisint
  implicit none
  type (branch), pointer, intent(inout)   :: start
  type (branch), pointer, intent(inout)   :: sommet
  integer(kind=precisint), intent(inout)  :: n_cut    ! Is the number of branches that died during all the current generation
  
!*****************************************************************************80
  !
  !! CUT_BRANCH does:
  !  - deallocate all branches after START
  !  - nullify mother%child for the sacrificed one
  !  - change mother%state, reorganize mother's children and number of children
  !  - Rebind chain.
  !
  type (branch), pointer                  :: mother, origin
  type (branch), pointer                  :: a, b
  type (branch_pointer), dimension(:), allocatable   :: temporaire
  integer                                 :: i, ii

  ! Are we trying to cut the trunk?
  if (associated(start%parent%p)) then
    mother => start%parent%p
  else
    print *, 'CUT_BRANCH cutting the trunk'
    print *, "Terminating program"
    call spacial_deallocate_branch(start)
    nullify(start%owner)
    deallocate(start)
    stop
  endif

  ! Which child is the branch we try to cut for its parent?
  ii=0
  do i=1,mother%N_children
    if (associated(mother%child(i)%p,start)) then
      ii=i
      exit
    endif
  enddo
  
  if (ii==0) then                                                                   !
    print *, "CUT_BRANCH cutting a desherited branch"                               ! DEBUG
    print *, "Terminating program"                                                  ! 
    print *, 'Cutting the branch: ', start%petit_nom, '  // Child #: ',ii           !
    print *, 'Cutting the child of branch: ', mother%petit_nom, '  // Child #: ',ii !
    stop                                                                            !
  endif

!  print *, 'Cutting the branch: ', start%petit_nom, '  // Child ',ii, 'out of ', mother%N_children


!====HERE WE CUT THE BRANCH, THEN GO THROUGH THE CHLDREN CUT THEM TOO, THEN REORGANIZE THE BEFORE-AFTER POINTERS====
  origin => start%before%p
  b => start
  if (associated(b,sommet)) then
    sommet => origin
    nullify(sommet%after%p)
    call spacial_deallocate_branch(b)
    nullify(b%owner)
    deallocate(b)
    n_cut = n_cut+1
  else
    do i=1,b%N_children
      if (associated(b%child(i)%p)) nullify(b%child(i)%p%parent%p)
    end do
    b => b%after%p
    call spacial_deallocate_branch(b%before%p)
    nullify(b%before%p%owner)
    deallocate(b%before%p)
    n_cut=n_cut+1
    nullify(b%before%p)
    nullify(start)
    do
      if (associated(b%parent%p)) then
        origin%after%p => b
        b%before%p => origin
        exit
      else
        if (associated(b,sommet)) then
          sommet => origin
          nullify(sommet%after%p)
          call spacial_deallocate_branch(b)
          nullify(b%owner)
          deallocate(b)
         n_cut=n_cut+1
          exit
        else
          do i=1,b%N_children
            if (associated(b%child(i)%p)) nullify(b%child(i)%p%parent%p)
          enddo
          b => b%after%p
          call spacial_deallocate_branch(b%before%p)
          nullify(b%before%p%owner)
          deallocate(b%before%p)
          n_cut=n_cut+1
          nullify(b%before%p)
        endif
      end if
    end do
  endif
!===========================================================================================================
  
  !---------HERE WE RESIZE AND REARRANGE MOTHER%CHILD(:) POINTERS----------------------------------------------------------!
  !----There are two cases : mother%child size is 1 (in which case we just deallocate it) or size > 1 (complicated case)---!
!  if(mother%N_children>1) then
!		allocate(temporaire(mother%N_children-1))
  
!		if (ii<mother%N_children) then        ! VERY IMPORTANT : THE CHILDREN MUST BE NUMBERED
!			do i=ii,mother%N_children-1			! ASSOCIATED CHILD(i)%P MUST BE THE FIRST N_CHILDREN
!				temporaire(i)%p => mother%child(i+1)%p
!				nullify(mother%child(i+1)%p)
!			end do
!		else
!			nullify(mother%child(ii)%p)
!		end if
!
!		do i=1,ii-1
!			temporaire(i)%p => mother%child(i)%p
!			nullify(mother%child(i)%p)
!		end do
		!Now temporaire(:) has become what mother%child(:) must become. 
!
!		!It leaves us with the task of resizing mother%child(:) then copy temporaire(:) back to mther%child(:)
!		deallocate(mother%child)
!		allocate(mother%child(mother%N_children-1))
!		do i=1,mother%N_children-1
!			mother%child(i)%p => temporaire(i)%p
!			nullify(temporaire(i)%p)
!		end do
!		deallocate(temporaire)
!	else
!		deallocate(mother%child)
!	end if
  !-------------------------------------------------------------------------------------
  
  if (ii<mother%N_children) then        ! VERY IMPORTANT : THE CHILDREN MUST BE NUMBERED
    do i=ii,mother%N_children-1			! ASSOCIATED CHILD(i)%P MUST BE THE FIRST N_CHILDREN
      mother%child(i)%p => mother%child(i+1)%p
      nullify(mother%child(i+1)%p)
    end do
  else
    nullify(mother%child(ii)%p)
  end if
  
  mother%N_children = mother%N_children - 1

end subroutine cut_branch
!*******************************************************************************



!§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
subroutine spacial_deallocate_branch(a)                                        !§
!§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
!Routine to use to delete the spacial component of a branch
!(must be called before deallocating any branches)
implicit none
	type(branch), pointer, intent(inout)        :: a
	type(great_space), pointer         :: dummy
	integer                            :: dum, xx, yy, zz, x0, y0, z0
	
	!====DEBUG DIAGNOSIS TO DELETE===============
	!logical,dimension(9,3)            :: space_around
	!integer                            :: desir2	
	!===========================================
	
	
	
	if(associated(a%owner%tree_point%the_space))then
		dummy => a%owner%tree_point%the_space
		!===DEBUG DIAGNOSIS TO DELETE==============
		!if(a%N_children==0)then
		!	call spacial_recog(a,desir2,space_around,dummy)
		!	write(unit=20,fmt="(5i7)") a%owner%tree_point%generat, desir2, a%int_x, a%int_y, a%int_z
		!endif
		!==========================================
		xx = a%int_x
		yy = a%int_y
		zz = a%int_z
		if(.not. dummy%xyz_coord(xx,yy,zz))then
			print*,"error(0) in spacial deall"
			print*, dummy%xyz_coord(xx,yy,zz)
			print*, xx, yy, zz
			stop
		endif
		
		dummy%xyz_coord(xx,yy,zz) = .false.
	endif
	
end subroutine spacial_deallocate_branch
!§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§



!*******************************************************************************!
! Routine used at the end of the program : We deallocate anything that wasn't   !
subroutine deall(a_tree_pointer, num)                                           !
	use precis_mod, only                    : prec => working_precis, precisint
	implicit none
	integer(kind=precisint), intent(inout)	:: num
	type(arbre_pointer), pointer, intent(inout)	:: a_tree_pointer
	type (branch), pointer 					:: start
	integer									:: chil, Nchiltronc
!	type (branch), pointer				    :: tronc
!	type (branch), pointer				    :: sommet
	
	Nchiltronc = a_tree_pointer%trunk_point%N_children
	do chil=1,Nchiltronc
		if(associated(a_tree_pointer%trunk_point%child(1)%p)) then
			start => a_tree_pointer%trunk_point%child(1)%p
			call cut_branch(start,a_tree_pointer%peak_point,num)
		else
			write(*,*) "BIG ERROR, you must rewrite the code !"
			stop
		endif
	enddo
	call spacial_deallocate_branch(a_tree_pointer%trunk_point)
	deallocate(a_tree_pointer%trunk_point)
	nullify(a_tree_pointer%trunk_point)
	nullify(a_tree_pointer%peak_point)
	
	if(associated(a_tree_pointer%tree_point%the_space))then
		deallocate(a_tree_pointer%tree_point%the_space%xyz_coord)
		deallocate(a_tree_pointer%tree_point%the_space)
		nullify(a_tree_pointer%tree_point%the_space)
	endif
	
	deallocate(a_tree_pointer%tree_point)
	nullify(a_tree_pointer%tree_point)
	
	deallocate(a_tree_pointer)
	
end subroutine deall
!*******************************************************************************!


!*****************************************************************************80
function height(sommet)                                                      !80
!*****************************************************************************80
  implicit none
  type (branch), pointer, intent(in)        :: sommet
  integer                                   :: height
  !*************************************************************************80
  type (branch), pointer                    :: a

  a => sommet
  height=0
  do
    height=max(a%distance,height)
    if (associated(a%before%p)) then
      a=> a%before%p
    else
      exit
    end if
  end do
end function height
!***************************************************************************80

!===============================================================================
!DEBUG_SUBrout

end module mod_tree
