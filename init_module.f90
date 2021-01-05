module init_module

!=====================================================================!
! Module that regroups a collection of useful tools                   !
!                                                                     !
! Only the first 4 subroutines are essential for our tree simulation: !
! inputtree, inputlog, inputdata and inputstrat                       !
!                                                                     !
! The rest of the subroutines are tools for plotting graphs, writing  !
! files.                                                              !
!=====================================================================!


public :: inputtree, inputlog, inputdata, inputstrat, writeparam, write_gnuscript,&
			& filen ,bumpchar, gnuplotter, hist_int, int2char, move, int2charnext,&
			& write_octscript, octplotter

contains

!======================================================================================!
! The following 4 subroutines are routines that initialize the global variables that   !
! are stored within the module parameters_module.f90								   !
!																					   !
! These global variables are initialized from 4 text files called inputXXX.basic       !
! that must be stored inside a directory called input                             	   !
! the advantage with such subroutines is that the parameters becomes easily modifiable !
subroutine  inputtree()


  use precis_mod, only             :  prec => working_precis
  use parameters_module, only  : N_c_Max, init_nb_leaves, rw, frac_desir_down, &
								exp_egoisme, exp_reward, ran_posit, min_chil
  implicit none
  integer						:: nunit
  
  open(newunit=nunit, file = "input/inputtree.basic", status = "old",&
       action="read", form="formatted" , position="rewind")
       
     read(unit=nunit,fmt="(1i6)")    N_c_Max

     read(unit=nunit,fmt="(1i6)")    init_nb_leaves
     read(unit=nunit,fmt="(1i6)")    rw
     read(unit=nunit,fmt="(1i6)")    min_chil
     read(unit=nunit,fmt="(f8.0)")    frac_desir_down
     read(unit=nunit,fmt="(f8.0)")    ran_posit

     read(unit=nunit,fmt="(f8.0)") exp_egoisme
     read(unit=nunit,fmt="(f8.0)") exp_reward

  close(unit=nunit)

endsubroutine inputtree

subroutine  inputlog()

use precis_mod, only            :  prec => working_precis
use parameters_module, only : N_generation, prod_leaf, rel_mainte_vol, &
							  rel_res_init, exp_mainte, volum_weight, &
							  leaves_weight, res_prod_rel, flux_res_rel, mech_coef, &
							  cost_leaf, cost_volum, extra_cost, N_xy, N_z, rupture_coef, &
							  base_rupture, rupture_exp
  implicit none 
  integer					:: nunit
  
  open(newunit=nunit, file = "input/inputlog.basic", status = "old",&
       action="read", form="formatted" , position="rewind")

     read(unit=nunit,fmt="(1i6)")     N_generation
     read(unit=nunit,fmt="(1i6)")     prod_leaf
     read(unit=nunit,fmt="(1i6)")     extra_cost
     read(unit=nunit,fmt="(1i6)")     cost_leaf
     read(unit=nunit,fmt="(1i6)")     cost_volum
     read(unit=nunit,fmt="(f9.0)")	  leaves_weight
     read(unit=nunit,fmt="(f9.0)")	  volum_weight
     read(unit=nunit,fmt="(1i6)")     rel_res_init
     read(unit=nunit,fmt="(f9.0)")    mech_coef
     read(unit=nunit,fmt="(f9.0)")    base_rupture
     read(unit=nunit,fmt="(f9.0)")    rupture_coef
     read(unit=nunit,fmt="(f9.0)")    rupture_exp
     read(unit=nunit,fmt="(f9.0)")	  exp_mainte
     read(unit=nunit,fmt="(f9.0)")    rel_mainte_vol
     read(unit=nunit,fmt="(f9.0)")	  res_prod_rel
     read(unit=nunit,fmt="(f9.0)")	  flux_res_rel
     read(unit=nunit,fmt="(1i6)")     N_xy
     read(unit=nunit,fmt="(1i6)")     N_z

  close(unit=nunit)

  return 

endsubroutine inputlog

subroutine inputdata()

use precis_mod, only            :  prec => working_precis
use parameters_module, only : N_c_Max_inc, prod_leaf_inc, loop_choice,minprod_leaf,&
							maxprod_leaf, minrel_mainte_vol, maxrel_mainte_vol, div_maint_vol, &
							minN_c_Max, maxN_c_Max, mincost_vol, maxcost_vol, cost_vol_inc, costextra_min,&
							costextra_max, costextra_inc
  implicit none 

  integer					:: nunit
  
  open(newunit=nunit, file = "input/inputdata.basic", status = "old",&
       action="read", form="formatted" , position="rewind")
     read(unit=nunit,fmt="(1i6)")     loop_choice
     read(unit=nunit,fmt="(1i6)")     N_c_Max_inc
	 read(unit=nunit,fmt="(1i6)")     minN_c_Max
     read(unit=nunit,fmt="(1i6)")     maxN_c_Max
     read(unit=nunit,fmt="(1i6)")     prod_leaf_inc
     read(unit=nunit,fmt="(1i6)")     minprod_leaf
     read(unit=nunit,fmt="(1i6)")     maxprod_leaf
     read(unit=nunit,fmt="(f9.0)")    minrel_mainte_vol
     read(unit=nunit,fmt="(f9.0)")    maxrel_mainte_vol
     read(unit=nunit,fmt="(1i5)")     div_maint_vol
     read(unit=nunit,fmt="(i5)")      mincost_vol
     read(unit=nunit,fmt="(i5)")      maxcost_vol
     read(unit=nunit,fmt="(1i5)")     cost_vol_inc
     read(unit=nunit,fmt="(i5)")      costextra_min
     read(unit=nunit,fmt="(i5)")      costextra_max
     read(unit=nunit,fmt="(i5)")      costextra_inc
  close(unit=nunit)
  
end subroutine inputdata

subroutine inputstrat(schem)
use precis_mod, only            :  prec => working_precis
  implicit none 
  integer, dimension(:),intent(inout) 			:: schem
  integer										:: i, nunit
  open(newunit=nunit, file = "input/inputstrat.basic", status = "old",&
       action="read", form="formatted" , position="rewind")
  do i=1,size(schem)
     read(unit=nunit,fmt="(1i6)")    schem(i)
  enddo
  close(unit=nunit)
end subroutine inputstrat
!======================================================================================!


!-------------------------------!------------------------------------------------!----------------------------!
function int2char(p) result (c) ! returns the character equivalent of an integer (no integer higher than 9999)!
!-------------------------------!------------------------------------------------!----------------------------!
!   THE RETURNED INTEGER IS OF THE FORM 0001, 0002 etc
!
! Only checked with  positive integers
!
  integer,intent(in)                        :: p
!  character(len=int(log10(abs(real(p))))+1) :: c
  character(len=4) :: c
  character(len=*), parameter               :: numbers = "0123456789"
  integer                                   :: i,r,n,j

  n=p
  i=4
  do j=1,4
     r=n-10*(n/10)+1
     c(i:i)=numbers(r:r)
     n=n/10
     i=i-1
  end do

end function int2char

!-------------------------------!------------------------------------------------!----------------------------!
function int2charnext(p,i) result (c) ! returns the character equivalent of an integer (other subroutine) !
!-------------------------------!------------------------------------------------!----------------------------!
! ! THE RETURNED INTEGER IS OF THE FORM 01, 02, 03 if i = 2 etc
! Only checked with  positive integers
!
implicit none
  integer, intent(in)						:: i
  integer,intent(in)                        :: p
!  character(len=int(log10(abs(real(p))))+1) :: c
  character(len=i) :: c
  character(len=6) :: temp
  !character(len=*), parameter               :: numbers = "0123456789"
  integer                                   :: r,n,j

	if(p>=10**i)then
		print*,"invalid int2charnext argument"
		stop
	endif
	if(i>6)then
		print*,"invalid int2charnext argument"
		stop
	endif
	
  write(temp,'(i6)') p
  r = i - len(trim(ADJUSTL(temp)))
  do j=1,r
	c(j:j)="0"
  enddo
  c(r+1:i)=trim(ADJUSTL(temp))

end function int2charnext


!************************************************************
! The 2 following subroutines are simply there so that we can create text files that will have different names :
! For example , using filen and bumpchar allow us to successively create files with the names :
! name00data then name01data then name02data etc, for an arbitrary number of files (between 1 and 100) 
subroutine bumpchar(string,n)

!  Change the string "00" by "01", etc.
	  implicit none
      integer,intent(in)             :: n
      character(len=*),intent(inout) :: string
      integer                        :: m

      m = n
      do m = n,1,-1
         if (string(m:m) /= "9") then
            string(m:m) = char(ichar(string(m:m))+1)
            return
         endif
         string(m:m) = "0"
      enddo

endsubroutine bumpchar

!Write new file names     
subroutine filen(fname,naame)

!  File name must be of the form fname="name00.dat", and
!  the number of characters name=(number of characters in fname)
      implicit none
      integer, intent(in)             :: naame
      character(len=*), intent(inout) :: fname
      integer                         :: i, nn
      logical                         :: fexist
     
      nn = naame - 4
      do i = 1, 100
        inquire(file = fname, exist = fexist)
        if (fexist) then
           call bumpchar(fname, nn)
        else
          exit
        endif
      enddo

endsubroutine filen
!**********************************************************

!=======================================================================
! Subroutine that write the values we used for all the global variables
! during a given run of a program  in the directory 'fname'
subroutine writeparam(fname,option,schem)
 use precis_mod, only            :  prec => working_precis
 use parameters_module
 implicit none
 character(len=*),intent(in)			:: fname
 integer,intent(in)						:: option
 integer,dimension(:),intent(in)		:: schem
 integer								:: nunit
 
 open(newunit=nunit,status="new",action="write",file=fname)
  
     write(unit=nunit,fmt="(a,1i9)") "N_generation :",N_generation
     write(unit=nunit,fmt=*) " "
     write(unit=nunit,fmt="(a,1i9)") "Default prod_leaf :",prod_leaf
     write(unit=nunit,fmt=*) " "

     write(unit=nunit,fmt="(a,1i9)") "cost of 1 leaf :",cost_leaf
     write(unit=nunit,fmt="(a,1i9)") "cost of creating volume :",cost_volum
     write(unit=nunit,fmt="(a,1i9)") "cost of 1 full branch :",extra_cost
     write(unit=nunit,fmt=*) " "

     write(unit=nunit,fmt="(a,f19.7)") "leaves weight", leaves_weight
     write(unit=nunit,fmt="(a,f19.7)") "weight per volume", volum_weight
     write(unit=nunit,fmt=*) " "

     write(unit=nunit,fmt="(a,f19.7)") "relation bwt reserve max and prod per leaf (res_prod_rel)",res_prod_rel
     write(unit=nunit,fmt="(a,f19.7)") "relat bwt flux max and its max reserve (flux_res_rel)",flux_res_rel
     write(unit=nunit,fmt="(a,1i9)")  "rel_res_init :",rel_res_init
     write(unit=nunit,fmt=*) " "
     
     write(unit=nunit,fmt="(a,f19.7)") "Default rel_mainte_vol:",rel_mainte_vol    
     write(unit=nunit,fmt="(a,f19.7)") "exp_mainte, how the maintenance grows w/ volume", exp_mainte
     write(unit=nunit,fmt="(a,f19.7)")"mech_coeff, how the much constr",mech_coef
     write(unit=nunit,fmt="(a,f19.7)")"base_rupture :", base_rupture
     write(unit=nunit,fmt="(a,f19.7)")"rupture_coef, how resistance to rupture grows w/ volum", rupture_coef
      write(unit=nunit,fmt="(a,f19.7)")"rupture_exp, how resistance to rupture grows w/ volum", rupture_exp
     write(unit=nunit,fmt=*) " "

	 write(unit=nunit,fmt="(a,1i9)") "Default N_c_Max :",N_c_Max
     write(unit=nunit,fmt="(a,1i9)") "init_nb_leaves :",init_nb_leaves
     write(unit=nunit,fmt="(a,1i9)")  "rw :",rw
      write(unit=nunit,fmt="(a,1i9)")  "min_chil :",min_chil
     write(unit=nunit,fmt=*) " "
     
     write(unit=nunit,fmt="(a,f19.7)")  "frac_desir_down :",frac_desir_down
     write(unit=nunit,fmt="(a,f19.7)")  "ran_posit:", ran_posit
     write(unit=nunit,fmt="(a,f19.7)") "exp_egoisme :",exp_egoisme
     write(unit=nunit,fmt="(a,f19.7)") "exp_reward :",exp_reward
     write(unit=nunit,fmt=*) " "
     
     write(unit=nunit,fmt="(a,1i9)") "size of the xy surface, N_xy:",N_xy
     write(unit=nunit,fmt="(a,1i9)") "height N_z:", N_z

     write(unit=nunit,fmt=*) " "
     
     block
	 integer :: i
		do i=1,size(schem)
			write(unit=nunit,fmt="(a,a,a,1i4)") "strat ", int2charnext(i,2)," :",schem(i)
		enddo
	 end block
     write(unit=nunit,fmt=*) " "
     
  if(option==1) then
     write(unit=nunit,fmt="(a,1i9)") "min prod_leaf :",minprod_leaf
     write(unit=nunit,fmt="(a,1i9)") "max prod_leaf :",maxprod_leaf
     write(unit=nunit,fmt="(a,1i9)")   "Increment of prod_leaf", prod_leaf_inc
     write(unit=nunit,fmt="(a,1i9)") "min N_c_Max :",minN_c_Max
     write(unit=nunit,fmt="(a,1i9)") "max N_c_Max :",maxN_c_Max
     write(unit=nunit,fmt="(a,1i9)")   "Increment of N_c_Max :", N_c_Max_inc
	 write(unit=nunit,fmt="(a,f19.7)") "min rel_mainte_vol:",minrel_mainte_vol
     write(unit=nunit,fmt="(a,f19.7)") "max rel_mainte_vol:",maxrel_mainte_vol
     write(unit=nunit,fmt="(a,1i9)") "Div maintenance:",div_maint_vol
	 write(unit=nunit,fmt="(a,i9)") "min volume creation cost:",mincost_vol
     write(unit=nunit,fmt="(a,i9)") "max volume creation cost:",maxcost_vol
     write(unit=nunit,fmt="(a,1i9)") "Div on res max:",cost_vol_inc
	 write(unit=nunit,fmt="(a,i9)") "min creation cost (case 4 only):",costextra_min
     write(unit=nunit,fmt="(a,i9)") "max creation cost (case 4 only):",costextra_max
     write(unit=nunit,fmt="(a,i9)") "increment in costextra loop:",costextra_inc
     
     write(unit=nunit,fmt=*) " "
     write(unit=nunit,fmt=*) " "

     write(unit=nunit,fmt="(a,1i9)")  "loop_choice :", loop_choice
  else if(option==0) then
	!print*, "New simulation"
  else
	print*, "ILLICIT USE OF writeparam"
  endif
  
  if(nb_param/=42) then
	print*,"writeparam is not up to date please modify its code"
	stop
  endif
  
  close(unit=nunit)

end subroutine writeparam
!==========================================================

!######################################################################################
!The next routines create gnuplot script corresponding to a certain types of data files
!
subroutine write_gnuscript(fi,gscript,logoption,fileoption)
 use precis_mod, only            :  prec => working_precis
 implicit none
 character(len=*),dimension(:),intent(in)			:: fi
 character(len=*),intent(in)						:: gscript
 integer,intent(in)									:: logoption,fileoption
 integer											:: nunit, inte
 character(len=LEN(fi)+4),dimension(size(fi))		:: logfi
 
 open(newunit=nunit,status="new",action="write",file=gscript)
 write(unit=nunit,fmt="(a)") "set terminal png"
 
 select case(fileoption)
 
 case(1) ! Case 1: we want to plot the files called 'fi' in the main program
 	if (logoption==1) then
		do inte=1,size(fi)
			logfi(inte)="log_"//fi(inte)
		enddo
	endif
	!subdir//"/gscript.gp"
	if (logoption==1) 	write(unit=nunit,fmt="(a)") "set logscale y"
	do inte=1,size(fi)
		if (logoption==0) then
			write(unit=nunit,fmt="(a,a,a)") "set output '",fi(inte),".png'"
		else if (logoption==1) then
			write(unit=nunit,fmt="(a,a,a)") "set output '",logfi(inte),".png'"
		else
			print*,"error100"
			stop
		endif
		write(unit=nunit,fmt="(a,a,a)") "plot '",fi(inte),"' using 1:2 notitle with points pointtype 7"
		write(unit=nunit,fmt="(a)") "unset output"
	enddo
	
 
 case(2) ! make a script for the file named 'accross' in the main program
	
	do inte=1,size(fi) 
			write(unit=nunit,fmt="(a,a,a)") "set output '",fi(inte)//"_gen-max",".png'"
			write(unit=nunit,fmt="(a,a,a)") "plot '",fi(inte),"' using 1:2 notitle with boxes"
			write(unit=nunit,fmt="(a)") "unset output"
			
			write(unit=nunit,fmt="(a,a,a)") "set output '",fi(inte)//"_height-max",".png'"
			write(unit=nunit,fmt="(a,a,a)") "plot '",fi(inte),"' using 1:3 notitle with boxes"
			write(unit=nunit,fmt="(a)") "unset output"
			
			write(unit=nunit,fmt="(a,a,a)") "set output '",fi(inte)//"_node-max",".png'"
			write(unit=nunit,fmt="(a,a,a)") "plot '",fi(inte),"' using 1:4 notitle with boxes"
			write(unit=nunit,fmt="(a)") "unset output"
	enddo
 case(3)  ! For the 'cor' file
  !outdir//"/gscript.gp"
	write(unit=nunit,fmt="(a)") "set pointsize 5" 
	write(unit=nunit,fmt="(a,a,a)") "set output '",fi(1)//"_gen-max",".png'"
	write(unit=nunit,fmt="(a,a,a)") "plot '",fi(1),"' using 1:2:3 notitle with points pointtype 7 palette"
	write(unit=nunit,fmt="(a)") "unset output"
			
	write(unit=nunit,fmt="(a,a,a)") "set output '",fi(1)//"_height-max",".png'"
	write(unit=nunit,fmt="(a,a,a)") "plot '",fi(1),"' using 1:2:4 notitle with points pointtype 7 palette"
	write(unit=nunit,fmt="(a)") "unset output"
			
	write(unit=nunit,fmt="(a,a,a)") "set output '",fi(1)//"_node-max",".png'"
	write(unit=nunit,fmt="(a,a,a)") "plot '",fi(1),"' using 1:2:5 notitle with points pointtype 7 palette"
	write(unit=nunit,fmt="(a)") "unset output"
 case default
	print*, "error 101"
	stop
 end select
 close(unit=nunit)
end subroutine write_gnuscript
!######################################################################################

!--------------------------------------------------------------
!Simple routine that activates a gnuplot script at the directory you want
subroutine gnuplotter(subdir,gpfile)
  implicit none
  character(len=*),intent(in)		:: subdir, gpfile	
  integer							:: nunit
  open(newunit=nunit, action="write",status="new",form="FORMATTED", file="test1.sh")
	write(unit=nunit,fmt="(a)") "cd "//subdir
	write(unit=nunit,fmt="(a)") "gnuplot "//gpfile
  close(unit=nunit)
  call system ("chmod a+x test1.sh")
  call system ("./test1.sh")
  call system ("rm test1.sh")
end subroutine gnuplotter
!------------------------------------------------------------
!------------------------------------------------------------
subroutine move(subdir,newdir,files)
  implicit none
  character(len=*),intent(in)		:: subdir, files,newdir	
  integer							:: nunit
  open(newunit=nunit, action="write",status="new",form="FORMATTED", file="test2.sh")
	write(unit=nunit,fmt="(a)") "cd "//subdir
	write(unit=nunit,fmt="(a)") "mv "//files//" "//newdir
  close(unit=nunit)
  call system ("chmod a+x test2.sh")
  call system ("./test2.sh")
  call system ("rm test2.sh")
end subroutine move
!-----------------------------------------------------------

!888888888888888888888888888888888888888888888888888888888
subroutine hist_int(v_mi,v_ma,v_moy,bins,fil,filout)
!The subroutine make an histogramme from a list of INTEGER
!stored in a text file
!888888888888888888888888888888888888888888888888888888888
 use precis_mod, only            :  prec => working_precis
 implicit none
 integer, intent(out)						::v_mi, v_ma!, v_moy
 real, intent(out)							:: v_moy
 integer, intent(in)						:: bins
 character(len=*), intent(in)				::fil, filout
 
 integer									:: x,comp=0,io,i,wid,entier,nunit,nuna
 integer, dimension(0:bins-1)				::histo
 real(kind=prec)							:: mil, absc

 histo(:) = 0
 v_moy = 0.0
 !*************************************************************
! Determine size of the text file and its max/min value
! v_mi = min value, v_ma= max value of the data in the text file
 open(newunit=nunit, file = fil, status = "old",&
       action="read", form="formatted" , position="rewind")
	  Do
        read(unit=nunit,fmt="(1i8)",iostat=io) x
        if(comp==0) then
			v_mi=x
			v_ma=x
        endif
        If( io < 0 ) exit
        comp = comp + 1
        v_moy = v_moy + real(x)
        v_ma=max(x,v_ma)
        v_mi=min(x,v_mi)
      Enddo
  close(unit=nunit)
!**************************************************************
  v_moy = v_moy/real(comp)
!-------Decide on the binning of the histogram-----------------------------
!wid = width of the bins ; bins = number of bins maximum the histogram may have
  if(v_ma-v_mi>bins) then
	wid = ceiling(real(v_ma-v_mi,kind=prec)/real(bins,kind=prec))
  else
	wid=1
  endif
!By default we should obtain a histogram with number of bins equal to the integer "bins"
! But if the maximum value and the minimum value are too close then the else-statement
! reduce the numbers of bins.
!-----------------------------------------------------------
 mil=real(wid,kind=prec)/2.0_prec
 
  open(newunit=nunit, file = filout, status = "new",&
       action="write", form="formatted" , position="rewind")
	  
  open(newunit=nuna, file = fil, status = "old",&
       action="read", form="formatted" , position="rewind")
  histo(:)=0
										!
      do i=1,comp						! CREATE THE HISTOGRAM HERE
        read(unit=nuna,fmt="(1i8)") x		!
		if(x>v_ma) print*, "ERROR1"		!
		if(x<v_mi) print*, "ERROR2"		!
		entier=(x-v_mi)/wid				!
		histo(entier)=histo(entier)+1	!
      enddo
 
 !---WRITE THE HISTOGRAMME HERE------------------------
	 do i=0,bins-1
		absc=real(v_mi,kind=prec)+real(wid*i,kind=prec)+mil
		write(unit=nunit,fmt="(f12.5,1i8)") absc, histo(i)
	 enddo
      
 close(unit=nunit)
 close(unit=nuna)
end subroutine hist_int
!888888888888888888888888888888888888888888888888888888888


!subroutine readme1(fname)
! implicit none
! character(len=*),intent(in)			:: fname
!
!endsubroutine readme1

!===========================================================================
subroutine write_octscript(fi,gscript,fileoption)
!==========================================================================
 implicit none
 character(len=*),intent(in)						:: fi
 character(len=*),intent(in)						:: gscript
 integer,intent(in)									:: fileoption
 integer											:: nunit, inte
 character(len=LEN(fi))								:: tfi1,tfi2,tfi3,tfi4
 !character(len=LEN(fi)+4),dimension(size(fi))		:: logfi
 
 !print*,"hi"
 open(newunit=nunit,status="new",action="write",file=gscript)
 write(unit=nunit,fmt="(a)") "f=figure('Visible','off');"
 select case(fileoption)
 case(0)
 
	write(unit=nunit,fmt="(a)") "load "//fi
	write(unit=nunit,fmt="(a)") "plot("//fi//"(:,1),"//fi//"(:,2),'r.')"
	write(unit=nunit,fmt="(a)") "print -dpng "//fi//".png"
	
! case(1)
!	tfi1 = fi
!	tfi2 = fi
!	tfi3 = fi
!	tfi1(13:15) = int2charnext(1,3)
!	tfi2(13:15) = int2charnext(10,3)
!	tfi3(13:15) = int2charnext(30,3)
!	
!	write(unit=nunit,fmt="(a)") "hold on"
!		write(unit=nunit,fmt="(a)") "load "//tfi1
!		write(unit=nunit,fmt="(a)") "load "//tfi2
!		write(unit=nunit,fmt="(a)") "load "//tfi3	
!	write(unit=nunit,fmt="(a)") "plot("//tfi1//"(:,1),"//tfi1//"(:,2),'r.')"
!	write(unit=nunit,fmt="(a)") "plot("//tfi2//"(:,1),"//tfi2//"(:,2),'bo')"
!	write(unit=nunit,fmt="(a)") "plot("//tfi3//"(:,1),"//tfi3//"(:,2),'k+')"
!	write(unit=nunit,fmt="(a)") "print('MyPNG.png', '-dpng')"
!	write(unit=nunit,fmt="(a)") "exit"
 case default
	print*, "error 102"
	stop
 end select
 close(unit=nunit)
 
end subroutine write_octscript
!=======================================================================

!=======================================================================
subroutine octplotter(subdir,gpfile)
!=======================================================================
  implicit none
  character(len=*),intent(in)		:: subdir, gpfile	
  integer							:: nunit
  open(newunit=nunit, action="write",status="new",form="FORMATTED", file="test3.sh")
	write(unit=nunit,fmt="(a)") "cd "//subdir
	write(unit=nunit,fmt="(a)") "octave "//gpfile
  close(unit=nunit)
  call system ("chmod a+x test3.sh")
  call system ("./test3.sh")
  call system ("rm test3.sh")
end subroutine octplotter
!=======================================================================

!subroutine logy_convert(unite,

!end subroutine logy_convert

endmodule init_module
