!Created by Stuart Davie
!********************************************************************************************************
! Version 1  : Takes sets of microsolvated systems represented in spherical coordinates and reorders water molecules
!               based on proximity to nodes. If no nodes are input, orders by distance to central molecule. Can also strip water
!		molecules from training set if desired.
! Version 1.1: Notifies user when no node file is present and creates example file. Gives user option to strip waters from set.
! Version 1.2 : Modified into modular format for incorporation into Hermes by J. L. Mcdonagh September 2016
! 
! CHANGE LOG 
! Version 1.1 Additional error checks
! Version 1.2 Modularized for incorporation into Hermes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module wat_order
  implicit none
  contains
  subroutine waterorder(TSname, AtomsBefore, NoWaters)
character 		::	TSname*50,DummyChar*50
integer 		::	AtomsBefore, NoWaters, configNo, i,j,water_min,NoNodes,NoStrip,ios
real, allocatable 	:: 	raw(:),rawsort(:),waterfeatures(:,:),waterfeatsort(:,:),waterfeatsort2(:,:)
real, allocatable 	:: 	cartesians(:,:),NodeCarts(:,:),NodeFeats(:,:),nodeorder(:)
integer, allocatable 	:: 	nodelist(:),orderlist(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!call system('ls')									!Read input parameters
!print*,'Training Set Name?'								!
!read*,TSNAME										!
!print*,'No. of atoms before waters?'							!		
!read*,AtomsBefore									!	
!print*,'No. of waters?'									!	
!read*,NoWaters
!print*,'Strip no. of waters? ("0" to leave all waters in)'
!read*,NoStrip									!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!open(unit=11,file='solv_input.txt',status='old',iostat=ios)
!if(ios.eq.0) then
!   read(11,*)DummyChar,TSname
!   read(11,*)DummyChar,AtomsBefore
!   read(11,*)DummyChar,NoWaters
!else
!   print*,'error reading solv_input.txt'
!endif
!close(11)
NoStrip=0
print*,'Number of waters to strip',NoStrip

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(raw((AtomsBefore+3*NoWaters)*3-6+27))						!
allocate(rawsort((AtomsBefore+(3*(NoWaters-NoStrip)))*3-6+27))					!
allocate(waterfeatures(NoWaters,9))							!
allocate(waterfeatsort(NoWaters,9))							!Allocate arrays and initialize variables
allocate(waterfeatsort2(NoWaters,9))							!
allocate(cartesians(NoWaters,9))							!
allocate(NodeCarts(NoWaters,3))								!
allocate(NodeFeats(NoWaters,3))								!	
allocate(nodelist(NoWaters))								!	
allocate(nodeorder(NoWaters))								!	
allocate(orderlist(NoWaters))								!		
configNo   =	0									!
NodeFeats  = 	0									!
rawsort    = 	0									!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(unit=4,file='node_feats.txt',status='old',err=100,iostat=ios)			!
read(4,*)NoNodes									!
read(4,*)										!Read input nodes if they exist
do i=1,NoNodes										!(otherwise NodeFeats=0 from above)
	read(4,*,err=100,end=100)NodeFeats(i,:)						!	
enddo											!
close(4)										!
100 continue										!
if(ios.gt.0)then
NoNodes=NoWaters									!
print*,'Something wrong with node_feats.txt. Positioning nodes at the origin'		!
open(unit=4,file='node_feats.txt_example',status='unknown')				!
write(4,*)NoNodes;write(4,*);								!
do i=1,NoNodes										!
	write(4,*)NodeFeats(i,:)							!	
enddo											!
close(4)										!
endif											!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


open(unit=3, file='Node_carts_ord.txt', status='unknown')
open(unit=1, file=trim(TSNAME), status='old')
open(unit=2, file=trim(TSNAME)//'_new', status='unknown')
!!*********************MAIN FEATURE REDEFINITION LOOP*************************************************!!!!!!!!!!!
do while (.true.)												!For every configuration:
   read(1,*,err=101,end=101)raw											!read in features
   configNo=configNo+1												!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do i=1,NoWaters												!
      do j=1,9								!					!	
         waterfeatures(i,j)=raw((atomsBefore*3-6)+(9*(i-1)+j))		!get the features of the waters		!convert water features to cartesians
      enddo								!					!
      cartesians(i,1)=waterfeatures(i,1)*sin(waterfeatures(i,2))*cos(waterfeatures(i,3))!convert to		!
      cartesians(i,2)=waterfeatures(i,1)*sin(waterfeatures(i,2))*sin(waterfeatures(i,3))!cartesians		!
      cartesians(i,3)=waterfeatures(i,1)*cos(waterfeatures(i,2))			!			!	
      NodeCarts(i,1)=NodeFeats(i,1)*sin(NodeFeats(i,2))*cos(NodeFeats(i,3))	!also convert node		!convert node features to cartesians		
      NodeCarts(i,2)=NodeFeats(i,1)*sin(NodeFeats(i,2))*sin(NodeFeats(i,3))	!coords to cartesians		!
      NodeCarts(i,3)=NodeFeats(i,1)*cos(NodeFeats(i,2))				!				!
      if(configNo.eq.1)   write(3,*)i,NodeCarts(i,1),NodeCarts(i,2),NodeCarts(i,3)!output node carts		!
   enddo													!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
   do i=1,NoWaters												!			
      do j=1,NoWaters												!for each node  {	
         nodeorder(j)=sqrt((cartesians(j,1)-NodeCarts(i,1))**2+(cartesians(j,2)-NodeCarts(i,2))**2+&		!find distance between waters and node		
         &(cartesians(j,3)-NodeCarts(i,3))**2)		!Make array of water to current node distances		!		
      enddo						!							!		
      do j=1,i-1				      !orderlist() is the list of already allocated water mols	!determine closest water not already taken		
         nodeorder(orderlist(j))=2*maxval(nodeorder)  !This loop sets such h20 mols to high value so they are 	!by other node
      enddo					      !not selected as closest to current node			!		
      orderlist(i)=minloc(nodeorder,1) 		!Finds closest h2o to node					!list water appropriately   }			
      waterfeatsort(i,:)=waterfeatures(orderlist(i),:)!Lists waters according to order				!	
   enddo													!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
   do i=1,NoWaters												!		
      do j=1,3							!						!		
         waterfeatsort2(i,j)=waterfeatsort(i,j)			!						!		
      enddo							!						!Reorder hydrogens based on distance	
      if(waterfeatsort(i,4).lt.waterfeatsort(i,7)) then		!IF first hydrogen is nearer than second	!to atom of interest
         waterfeatsort2(i,4:9)=waterfeatsort(i,4:9)		!list hydrogens in current order		!
      else							!ELSE						!
         waterfeatsort2(i,3:6)=waterfeatsort(i,6:9)		!swap the hydrogens in the list			!			
         waterfeatsort2(i,6:9)=waterfeatsort(i,3:6)		!						!	
      endif							!						!		
   enddo													!		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   rawsort(1:atomsbefore*3-6)=raw(1:atomsbefore*3-6)		!rawsort will be the node-defined		!	
   rawsort((atomsbefore*3-6+1):((atomsbefore+3*(NoWaters-NoStrip))*3-6))=0!feature set.Set h2o=0 so errors seeable!	
   rawsort(size(rawsort)-27+1:size(rawsort))=raw(size(raw)-27+1:size(raw))		
   do i=1,NoWaters-NoStrip											!		
      do j=1,9													!Get all new ordering
         rawsort((atomsBefore*3-6)+(9*(i-1)+j))=waterfeatsort2(i,j)	!New order from waterfeatsort2		!	
      enddo													!		
   enddo													!	
   write(2,*)rawsort							!Output new training set		!and output to appropriate file		
enddo														!		
101 continue													!	
!!*********************END MAIN FEATURE REDEFINITION LOOP*********************************************!!!!!!!!!!!
print*,'Total configs =',configNo						
call system('mv '//trim(TSNAME)//' '//trim(TSNAME)//'_original')		!Rename new and old training sets so they
call system('mv '//trim(TSNAME)//'_new '//trim(TSNAME))				!are immediately ferebus compatable

end subroutine waterorder

end Module wat_order
