!Created by Stuart Davie
!********************************************************************************************************
! Version 1  : Takes sets of microsolvated systems represented in spherical coordinates and reorders water molecules
!               based on proximity to nodes. If no nodes are input, orders by distance to central molecule. Can also strip water
!		molecules from training set if desired.
! Version 1.1: Modified into modular format for incorporation into Hermes by J. L. Mcdonagh September 2016
! 
! CHANGE LOG 
! Version 1.1  Modularized for incorporation into Hermes.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module nodorder
implicit none

contains

subroutine node_order(natoms, AtomsBefore, NoWaters)
Character		:: DummyChar*50
integer			:: KeepFeats,NoFeats,FeatNo,AtomsBefore,NoWaters,Counter,ios,curnode,featno2
integer 		:: NoNodes,NodeNo, natoms
integer			:: check
integer, allocatable	:: keep_list(:),node_list(:)
double precision, allocatable :: features(:,:)


!open(unit=11,file='solv_input.txt',status='old',iostat=ios)
!if(ios.eq.0) then
!   read(11,*)DummyChar
!   read(11,*)DummyChar,AtomsBefore
!   read(11,*)DummyChar,NoWaters
!   read(11,*)DummyChar
!   read(11,*)DummyChar
!   read(11,*)DummyChar
!   read(11,*)DummyChar
!else
!   print*,'error reading solv_input.txt'
!endif
! close(11)


open(unit=1,file='ordered_features.txt',status='old')
   read(1,*)DummyChar,DummyChar,DummyChar,DummyChar,DummyChar,KeepFeats
print*,KeepFeats,'important features'
   read(1,*)
   allocate(keep_list(KeepFeats))
   allocate(node_list(KeepFeats))
   do FeatNo=1,KeepFeats
      read(1,*)keep_list(FeatNo)
   enddo
close(1)

open(unit=12,file='node_feats.txt',status='old')
   read(12,*)NoNodes
print*,'Number of nodes =',NoNodes
allocate(features(NoNodes,3))
   read(12,*)
   do NodeNo=1,NoNodes
      read(12,*)features(NodeNo,:)
    !  print*,features(nodeno,:)
   enddo
close(12)

counter=0
do FeatNo=1,KeepFeats
   if(keep_list(FeatNo).gt.(AtomsBefore*3-6)) then
      curnode=1+(keep_list(FeatNo)-(AtomsBefore*3-6)-1)/9
      check=0
      do featno2=1,Counter
         if(curnode.eq.node_list(FeatNo2)) check=1 
      enddo
      if(check.eq.0) then
         counter=counter+1
         node_list(counter)=curnode
         print*,'important water',Curnode
      endif
   endif
enddo



open(unit=13,file='node_feats_cut.txt',status='unknown')
   write(13,*)Counter
   write(13,*)' Node list features'
   do NodeNo=1,Counter
      write(13,*)features(node_list(NodeNo),:)
   enddo
close(13)

open(unit=10,file='Change_natoms_in_FINPUT_to_value_in_this_file.txt',status='unknown')
print*,'Change natoms in FINPUT.txt to',AtomsBefore+Counter*3
write(10,*)AtomsBefore+Counter*3
natoms = AtomsBefore+Counter*3
close(10)
!open(unit=13,file='FINPUT.txt',status='old',iostat=ios)
!if(ios.ne.0) then
!   print*,'No FINPUT.txt file to edit'
!else
!   read(13,*)
!   write(13,*)'natoms',AtomsBefore+Counter*3
!endif
close(13)


end subroutine node_order

end module nodorder

