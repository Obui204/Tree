program main
  use precis_mod, only             :  prec => working_precis, precisint
  use parameters_module
  use init_module
  use mod_tree
  use mod_sim1  
  
  !*****************************************************************************80
  implicit none
  integer                           :: hmax, genmax
  integer(kind=precisint)           :: nmax

  type(arbre_pointer), pointer      :: tree_object
  integer, dimension(8)	            :: schem

  
	integer(kind=precisint)	:: childless,n_node,voltot,restot,n_cre,n_cut, n_cut2
	integer	                :: flag=0, inte,h,generation 
  
  
  !==================
  !We initialize the global variables stored in parameters_module
  ! using these 3 routines : they initialize their values by reading them 
  ! from text files named inputtree.basic, inputlog.basic and inputdata.basic
  call inputtree()
  call inputlog()
  call inputdata()
  !==================
  !Initialize the variable schem from the filetext inputstrat.basic
  call inputstrat(schem)
  
    allocate(tree_object)
  	call alloc_arbre_pointer(tree_object)  !initialize the tree

	  hmax=1
	  nmax=0_precisint
	  genmax=0
	  h=tree_object%trunk_point%distance
	  n_node=1
	  
	  
	do generation = 1, N_generation
		
		call flux_down_distrib(schem, tree_object, h, generation)
		!h is the height of our tree
		                     
		call flux_up_distrib(schem,tree_object,generation,n_cre,n_cut,n_cut2,flag,h,voltot,restot,childless)
		!this routine gives us as output:
		!   - n_cre : number of branches created during the process "flux_up"
		!   - n_cut : number of branches that died; n_cut2 : branches that died due to mechanical failure
		!   - flag : indicates whether the trunk is still alive or dead
		!   - voltot : volume total at the end of the flux_up process
		!   - restot : reserve total
		!   - childless : number of branches w/o any children (number of extremities)
		!  These output are useful if you want to see how the tree evolve

		n_node = n_node + n_cre - n_cut  ! variable indicating the number of branch the tree currently has
		
		if(flag==1) n_node = 0
		
		!======================================================
		!We print on the terminal the current height, number of branches and total volume of the tree
		call print_terminal(h,generation,n_node,voltot) 
		!======================================================
		
		nmax = max(nmax,n_node)  ! the maximal number of branches the tree has ever had during the simulation
		hmax = max(h,hmax)       ! the max height it has ever had

		if(flag==1) exit
		genmax = generation

	enddo
	
	print*, nmax, hmax, genmax
	
	!genmax = generation
  
	n_cut = 0
		
	call deall(tree_object, n_cut) !tree deallocation routine
  
  
end program main
