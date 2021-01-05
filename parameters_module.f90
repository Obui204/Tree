!------------------------!-------------------------------------------------------------!
module parameters_module ! Modules storing all parameters and variables of the program !
!------------------------!-------------------------------------------------------------!

  use precis_mod, only                                : prec => working_precis

!------------------------!
! Branch specific stuff  !
!------------------------!

integer, save                  :: loop_choice = 0        !

real(kind=prec), save	       :: res_prod_rel = 2.0     ! relation btw prod_leaf and the default a%reserve_max of a branch (of volume 1) (reserve_max = prd_leaf * res_prod_rel)

real(kind=prec), save	       :: flux_res_rel = 2.0     ! a%flux_max = flux_res_rel * a%reserve_max . 

integer, save                  :: rel_res_init = 5       ! initial reserve of tronc = rel_res_init*cost 

integer, save                  :: N_generation = 20      ! Number of generations to be simulated

integer, save                  :: N_c_Max = 5            ! maximum number of children a branch can have

integer, save                  :: prod_leaf = 6          ! Default production of each leaf

integer, save                  :: init_nb_leaves = 1     ! Initial number of leaves

integer, save                  :: cost_volum = 10        ! Cost of creating 1 unit of volume

integer, save                  :: cost_leaf = 10         ! Cost of creating 1 operating leaf

integer, save                  :: extra_cost = 10        ! Cost minimal of a creating a new branch of volume 1 and with 1 leaf

real(kind=prec), save          :: leaves_weight = 1.0    ! self-explanatory

real(kind=prec),save           :: volum_weight = 0.0

integer, save                  :: rw = 0                 ! Reserve wish

real(kind=prec), save          :: mech_coef = 0.02       ! Coefficient describing the impact of bending moment => Calculate a%mech_constr

real(kind=prec), save          :: base_rupture = 1.2     ! Base a%mech_constr a branch can support 

real(kind=prec), save          :: rupture_coef = 0.02    ! How much more resistance to rupture the cross-section area of the branch is adding

real(kind=prec), save          :: rupture_exp = 1.5      ! Exponent for which power law follows the resistence to rupture in function of the branch cross-area 

real(kind=prec), save          :: exp_mainte = 1.3

integer, save                  :: min_chil

real(kind=prec), save          :: rel_mainte_vol= 1.0    ! Maintenance relative to volume (>=1)

real(kind=prec), save          :: frac_desir_down = 2.0  ! Fraction of volume growth to be satisfied on the way down (>=1)

real(kind=prec), save          :: exp_egoisme = 1.0      ! exponent for selfishness

real(kind=prec), save  	       :: exp_reward = 1.0       ! exponent for good production levels

!real(kind=prec), save         :: exp_trust = 1.0       ! TO BE IMPLEMENTED


real(kind=prec), save  	       :: ran_posit = 1.0 

integer, save 	               :: N_xy = 200

integer, save                  :: N_z = 400


integer, save       		   :: minN_c_Max = 5          ! These 3 varaibles are the increments and bounds for N_c_Max
integer, save       		   :: maxN_c_Max = 5		  ! when we will loop varying that parameter
integer, save       		   :: N_c_Max_inc = 1  		  !

integer, save       		   :: minprod_leaf = 6   	  !
integer, save       		   :: maxprod_leaf = 6        ! These 3 variables are the increments and bounds for prod_leaf
integer, save       		   :: prod_leaf_inc = 1       ! when we will do a loop varying the parameter prod_leaf

integer, save				   :: costextra_min = -2	  !
integer, save				   :: costextra_max = -2      ! Only useful when we do the case 4 of the tree main program. Determines the bounds for extra_cost of the simulation 
integer, save				   :: costextra_inc = 1

integer,save				   :: mincost_vol = -10		  ! Bounds of cost_vol
integer,save				   :: maxcost_vol = -10		  !
integer,save				   :: cost_vol_inc = -10	  ! And increment

real(kind=prec), save	       :: minrel_mainte_vol= 1.0  !
real(kind=prec), save	       :: maxrel_mainte_vol= 1.0  !
integer, save			       :: div_maint_vol = 1		  !

integer						   :: random_flag = 0		  ! do not count it in nb_param


!=========HAVE TO BE UPDATE EACH TIME A GLOBAL VAR IS ADDED=======================================================
integer, save				   :: nb_param=42   !number of global variables in that module (not including itself)

public                   :: N_generation, N_c_Max, prod_leaf, rw, rel_mainte_vol, exp_egoisme, exp_reward,&
							& init_nb_leaves, rel_res_init, frac_desir_down, minN_c_Max, maxN_c_Max, minprod_leaf, &
							& maxprod_leaf, minrel_mainte_vol, maxrel_mainte_vol, div_maint_vol,  N_c_Max_inc, &
							& prod_leaf_inc, loop_choice, exp_mainte, leaves_weight, volum_weight, res_prod_rel, flux_res_rel,&
							& mech_coef, nb_param, cost_vol_inc, mincost_vol, maxcost_vol, extra_cost, cost_volum, cost_leaf, &
							& costextra_min, costextra_max, costextra_inc, N_xy, N_z, base_rupture, rupture_coef, ran_posit, min_chil
							

end module parameters_module
