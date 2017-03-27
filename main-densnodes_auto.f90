!Created by Stuart Davie
!********************************************************************************************************
! Version 1   : Takes sets of microsolvated systems represented in spherical coordinates and outputs nodes
!	        based on regions of high density. Outputed nodes are in cartesians, ready for xyz file usage
! Version 1.2 : Adapted by James L. McDonagh for inclusion into Hermes. September 2016.  
!
! CHANGE LOG
! Version 1.1: Linear overlap function added for averaging 
! Version 1.2: Minor alteration to modularize the program 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module den_nodes
  
  implicit none

contains

  subroutine density_nodes(TSname, AtomsBefore, NoWaters, Res, o_lap)

    character 			  ::  	TSName*50,DummyChar*50
    integer 			  ::  	AtomsBefore, NoWaters, configNo, i,j,config,cur_node, NodeNo
    double precision, allocatable ::    raw(:),waterfeatures(:,:),cartesians(:,:,:),min_list(:)
    double precision 		  ::	res,b2a,length,pi
    integer, allocatable 	  ::    density_grid(:,:,:),node_list(:,:),node_carts(:,:)
    integer 			  ::	x_min,x_max,y_min,y_max,z_min,z_max,x_atom,y_atom,z_atom
    integer 			  ::	O_lap,o_lapx,o_lapy,o_lapz,lap_val(3),ios
    parameter                           (b2a=0.529177249,pi=3.14159265359) ! b2a Bohr to Angstroms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!call system('ls') 									!
!print*,'Training Set Name?'								!
!read*,TSNAME										!
!print*,'No. of atoms before waters?'							!Read in the name of the training set file, the number of atoms 
!read*,AtomsBefore									!before the solvation shell (e.g. 13 for an alanine), the number 
!print*,'No. of waters?'								!of waters in the solvation shell, the resolution of the grid, and
!read*,NoWaters										!the amount of overlap between grid squares.
!print*,'Resolution?'									!
!read*,res										!
!print*,'Overlap?'									!
!read*,o_lap										!
!print*,'No. of Nodes?'									!
!read*,NodeNo										!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    open(unit=11,file='solv_input.txt',status='old',iostat=ios)
!    if(ios.eq.0) then
!       read(11,*)DummyChar,TSname
!       read(11,*)DummyChar,AtomsBefore
!       read(11,*)DummyChar,NoWaters
!       read(11,*)DummyChar,Res
!       read(11,*)DummyChar,o_lap
!    else
!       print*,'error reading solv_input.txt'
!    endif
!    close(11)
    NodeNo=NoWaters
    

    allocate(raw((AtomsBefore+3*NoWaters)*3-6+27))
    allocate(waterfeatures(NoWaters,9))
    allocate(node_list(NoWaters,3))
    allocate(node_carts(NoWaters,3))
    allocate(min_list(NoWaters))
    node_list=0;node_carts=0;configNo=0	

    open(unit=1, file=trim(TSNAME), status='old')				
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do while (.true.)							!
       read(1,*,err=100,end=100)raw					!
       configNo=configNo+1						!This loop reads through the training set file
    enddo								!to determine the number of configurations
100 continue								!
    if(configNo.eq.0) then						!
       print*,'No configurations counted. Check TS format and file name'!
       stop                                                             !
    endif								!
    rewind(1)								!
    print*,'total configurations=',configNo				!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(cartesians(configNo,NoWaters,3))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do config=1,configNo						     				        !
       read(1,*)raw												!
       do i=1,NoWaters												!
          do j=1,9												!
             waterfeatures(i,j)=raw((atomsBefore*3-6)+(9*(i-1)+j))						!
          enddo													!This loop reads in the training data and turns the
          cartesians(config,i,1)=waterfeatures(i,1)*sin(waterfeatures(i,2))*cos(waterfeatures(i,3))*b2a/res	!water features into cartesians for use later on.
          cartesians(config,i,2)=waterfeatures(i,1)*sin(waterfeatures(i,2))*sin(waterfeatures(i,3))*b2a/res	!Cartesians are scaled by b2a(bohr2angstrom) and by 
          cartesians(config,i,3)=waterfeatures(i,1)*cos(waterfeatures(i,2))*b2a/res				!the inverse resolution for grid. Grid is always 
       enddo													!1x1x1Ang so res will stretch or compress the system
    enddo													!as desired.
    close(1)													!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    x_min=floor(minval(cartesians(:,:,1)));x_max=ceiling(maxval(cartesians(:,:,1)))					!
    y_min=floor(minval(cartesians(:,:,2)));y_max=ceiling(maxval(cartesians(:,:,2)))					!
    z_min=floor(minval(cartesians(:,:,3)));z_max=ceiling(maxval(cartesians(:,:,3)))					!
														        !Finds the minimum and maximum x/y/z values required
    allocate(density_grid((x_max-x_min+1),(y_max-y_min+1),(z_max-z_min+1)))						!for volume to encompass full system. Rounds to lowest
													        	!and highest whole number so grid can be set up for
    print*,'x_min=',x_min,'x_max=',x_max										!density calc. Allocates density_grid array.
    print*,'y_min=',y_min,'y_max=',y_max										!
    print*,'z_min=',z_min,'z_max=',z_max										!
    print*,size(density_grid(:,1,1)),size(density_grid(1,:,1)),size(density_grid(1,1,:)),size(density_grid)		!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    open(unit=2,file='node_carts.txt', status='unknown')
    write(2,*)NodeNo;write(2,*)'Node list cartesians'
    open(unit=3,file='node_feats.txt', status='unknown')
    write(3,*)NodeNo;write(3,*)'Node list features'
!!**************************MAIN NODE ALLOCATION LOOP*******************************************!!!!!!!!!!!!!!!!!For each node,!!!!!!!!!!!!!!!!!!!!!
    do cur_node=1,NodeNo												                           !
       density_grid=0												                                   !Run through all configurations, find
       do config=1,configNo												                           !grid cell for each water, and add to array
          do j=1,NoWaters												                           !
             x_atom=nint(cartesians(config,j,1))-x_min+1	!Because even though range for dimension might be x_min	                           !	
             y_atom=nint(cartesians(config,j,2))-y_min+1	!to x_max, density_grid array has positions from	                           !		
             z_atom=nint(cartesians(config,j,3))-z_min+1	!1 to x_max-x_min+1					                           !
			                                                                                                                           !
             do o_lapx=-o_lap,o_lap                                                                                                                !
                do o_lapy=-o_lap,o_lap                                                                                                             !
                   do o_lapz=-o_lap,o_lap					                                                                   ! 
                      lap_val(1)=abs(o_lapx);lap_val(2)=abs(o_lapy);lap_val(3)=abs(o_lapz)				                           !	
                      if(x_atom+o_lapx.gt.0 .and. y_atom+o_lapy.gt.0 .and. z_atom+o_lapz.gt.0) then                                                !
                         if(x_atom+o_lapx.lt.(x_max-x_min+2) .and. y_atom+o_lapy.lt.(y_max-y_min+2) .and. &
				z_atom+o_lapz.lt.(z_max-z_min+2)) then   !	
                            density_grid(x_atom+o_lapx,y_atom+o_lapy,z_atom+o_lapz)= &                                                             !
                                 density_grid(x_atom+o_lapx,y_atom+o_lapy,z_atom+o_lapz)+ (o_lap-maxval(lap_val)+1)				   !	
                         endif                                                                                                                     !
                      endif												                           !
                   enddo                                                                                                                           !
                enddo                                                                                                                              !
             enddo											                                           !
          enddo													                                   !		
       enddo													                                   !			
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       print*,'total points=',sum(density_grid),NoWaters*ConfigNo							!
       print*,'max count=',maxval(density_grid)									        !
       Node_list(cur_node,:)=maxloc(density_grid)									!Find cell with maximum count (corresponding to
       node_carts(cur_node,1)=(node_list(cur_node,1)+(x_min-1))							        !location of maximum density). Write location.
       node_carts(cur_node,2)=(node_list(cur_node,2)+(y_min-1))							        !
       node_carts(cur_node,3)=(node_list(cur_node,3)+(z_min-1))							        !	
       write(2,*)' He',node_carts(cur_node,:)*res			!Nodes labelled He for ease of use with Jmol	!
       length=sqrt((node_carts(cur_node,1)*res/b2a)**2+(node_carts(cur_node,2)*res/b2a)**2   	&		        !
            +(node_carts(cur_node,3)*res/b2a)**2)			                                                !
       if(node_carts(cur_node,1).gt.0) then										!
          write(3,*)length,acos((node_carts(cur_node,3)*res/b2a)/length),				&		!
               atan(dble(node_carts(cur_node,2))/dble(node_carts(cur_node,1)))	!for use with other programs	        !
       elseif(node_carts(cur_node,2).gt.0) then                                                                         !
          write(3,*)length,acos((node_carts(cur_node,3)*res/b2a)/length),				&		!
               pi+atan(dble(node_carts(cur_node,2))/dble(node_carts(cur_node,1)))	!for use with other programs	!
       else                                                                                                             !
          write(3,*)length,acos((node_carts(cur_node,3)*res/b2a)/length),				&		!
               -pi+atan(dble(node_carts(cur_node,2))/dble(node_carts(cur_node,1)))                                      !
       endif                                                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       do config=1,configNo												!
          do j=1,NoWaters												!
             min_list(j)=sqrt((node_carts(cur_node,1)-cartesians(config,j,1))**2 &					!For each configuration, find distance between	
                  +(node_carts(cur_node,2)-cartesians(config,j,2))**2 &					                !each water and current node. Move the nearest 	
                  +(node_carts(cur_node,3)-cartesians(config,j,3))**2)					                !water outside of grid range so it does not	
          enddo													        !influence next node placement. This stops two nodes
          cartesians(config,minloc(min_list),:)=cartesians(config,minloc(min_list),:)-2*(x_max-x_min)		        !being next to each other due to high density, 	
       enddo														!when one node would capture the corresponding
    enddo														!structure alone.	
!!!!************************END MAIN NODE ALLOCATION LOOP*****************************************************!!!!!!!!!!!
    close(2)

  end subroutine density_nodes

end module den_nodes
