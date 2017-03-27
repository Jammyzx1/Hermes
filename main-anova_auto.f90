!Created by Stuart Davie
!********************************************************************************************************
! Version 1  :  Takes sets of microsolvated systems represented in spherical coordinates and reorders water molecules
!               based on proximity to nodes. If no nodes are input, orders by distance to central molecule. Can also strip water
!		molecules from training set if desired.
! Version 1.1:  Modified into modular format for incorporation into Hermes by J. L. Mcdonagh September 2016
! Version 1.2:  Allows a user defined minimum number of solvent molecules to surround the solute
! 
! CHANGE LOG 
! Version 1.1 Modularized for incorporation into Hermes
! Version 1.2 Allows user to define a minimum number of solvent molecules around a solute 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    module analysis_of_varience 
    implicit none
    contains 

    subroutine anova(Filename, AtomsBefore, NoWaters, Thresh, NoBins)

    character                     :: filename*30, AorF*1,DummyChar*50
    integer                       :: NoProps,NoFeats,ios,configNo,i,j,NoBins,BinNo,PropNo,FeatNo
    integer                       :: NoConfigs,BinAlo,curBin,FeatNo2,CurBin2,BinNo2,cross_switch
    double precision, allocatable :: featuresarray(:,:), propertiesarray(:,:),raw(:),TotalVar(:)
    double precision, allocatable :: maxminlist(:,:),meanlist(:,:,:),propmean(:),bvar(:,:),wvar(:,:,:)
    double precision, allocatable :: inner_means(:,:),inner_row_mean(:,:),max_array(:)
    integer, allocatable          :: binlocation(:,:),bincount(:,:),inner_count(:,:),Sel_Waters(:)
!multi_anova

    double precision, allocatable :: dbvar(:,:),dwvar(:)
    double precision              :: inner_total_mean,inner_bvar(2),inner_wvar,inner_total_var,max_val
    double precision		  :: thresh2,thresh
    integer			  :: atomsbefore,Nowaters,bmaxloc,counter,MinWaters,ccwat

!    call system('ls')
!    print*,'file?'
!    read*,filename
!    print*,'properties? (1-25)'
!    read*,NoProps
!    print*,NoProps,'properties'
!    print*,'No of atoms or no of features (a or f)?'
!    read*,AorF
!    if(AorF.eq.'a') then
!       print*,'No of Atoms?'
!       read*,NoFeats
!       NoFeats=3*NoFeats-6
!    else
!       print*,'No of Atoms?'
!       read*,NoFeats
!    endif
!    print*,NoFeats,'features'
!    print*,'No of bins'
!    read*,NoBins
!    print*,NoBins,'bins'
open(unit=11,file='solv_input.txt',status='old',iostat=ios)
if(ios.eq.0) then
   read(11,*)DummyChar,FileName
   print*,FileName
   read(11,*)DummyChar,AtomsBefore
   read(11,*)DummyChar,NoWaters
   read(11,*)DummyChar
   read(11,*)DummyChar
   read(11,*)DummyChar,Thresh
   print*,Thresh
   read(11,*)DummyChar,NoBins
   print*,NoBins
   read(11,*)DummyChar,MinWaters
else
   print*,'error reading solv_input.txt'
endif
close(11)    
    NoFeats=(AtomsBefore+3*NoWaters)*3-6
print*,NoFeats
NoProps=1
    cross_switch=0

    allocate(raw(NoFeats+NoProps))
    NoConfigs=0


    OPEN(UNIT=1, FILE=trim(FileName),STATUS='old',iostat=ios) !Open RAW file
    if(ios.ne.0) then
       print*,'error opening',trim(filename);stop
    endif
    do while(.true.)
       read(1,*,err=100,end=100)raw
       NoConfigs=NoConfigs+1
    enddo
100 continue
    print*,NoConfigs,' configurations'
    deallocate(raw)
    
    allocate(Sel_Waters(NoWaters))
    allocate(featuresarray(NoFeats,NoConfigs))
    allocate(propertiesarray(NoProps,NoConfigs))
    allocate(maxminlist(NoFeats,3))
    allocate(binlocation(NoFeats,NoConfigs))
    allocate(meanlist(NoBins,NoFeats,NoProps))
    allocate(bincount(NoBins,NoFeats))
    allocate(propmean(NoProps))
    allocate(bvar(NoFeats,NoProps))
    allocate(wvar(NoBins,NoFeats,NoProps))
    allocate(TotalVar(NoProps))
    allocate(max_array(NoFeats))
    meanlist=0;bincount=0;propmean=0;bvar=0;wvar=0;TotalVar=0

    rewind(1)
    do ConfigNo=1,NoConfigs
       read(1,*)featuresarray(:,ConfigNo),propertiesarray(:,ConfigNo)
       do PropNo=1,NoProps
          propmean(PropNo)=propmean(PropNo)+propertiesarray(PropNo,ConfigNo)
       enddo
    enddo
    close(unit=1)
    propmean=propmean/NoConfigs

    do PropNo=1,NoProps
       do ConfigNo=1,NoConfigs
          TotalVar(PropNo)=TotalVar(PropNo)+(propertiesarray(PropNo,ConfigNo)-propmean(PropNo))**2
       enddo
    enddo
    print*,'Total Sum of Squares deviation',TotalVar

    do FeatNo=1,NoFeats
       maxminlist(FeatNo,1)=maxval(featuresarray(FeatNo,:))
       maxminlist(FeatNo,2)=minval(featuresarray(FeatNo,:))
       maxminlist(FeatNo,3)=maxminlist(FeatNo,1)-maxminlist(FeatNo,2)
       do ConfigNo=1,NoConfigs
          BinAlo=ceiling((featuresarray(FeatNo,ConfigNo)-maxminlist(FeatNo,2))/(maxminlist(FeatNo,3))*NoBins)
          if(BinAlo.eq.0) BinAlo=1
          binlocation(FeatNo,ConfigNo)=BinAlo
          bincount(BinAlo,FeatNo)=bincount(BinAlo,FeatNo)+1
          do PropNo=1,NoProps
             meanlist(BinAlo,FeatNo,PropNo)=meanlist(BinAlo,FeatNo,PropNo)+(propertiesarray(PropNo,ConfigNo))
          enddo
       enddo
       do BinNo=1,NoBins
          if(bincount(BinNo,FeatNo).gt.0) then 
             do PropNo=1,NoProps
                meanlist(BinNo,FeatNo,PropNo)=meanlist(BinNo,FeatNo,PropNo)/bincount(BinNo,FeatNo)
                bvar(FeatNo,PropNo)=bvar(FeatNo,PropNo)+(meanlist(BinNo,FeatNo,PropNo)-propmean(PropNo))**2&
                &*(bincount(BinNo,FeatNo))           
             enddo
          endif
       enddo
    enddo

    do PropNo=1,NoProps
       do FeatNo=1,NoFeats
          do Configno=1,NoConfigs
             curBin=binlocation(FeatNo,ConfigNo)
             wvar(curBin,FeatNo,PropNo)=wvar(curBin,FeatNo,PropNo)&
                  &+(meanlist(curBin,FeatNo,PropNo)-(propertiesarray(PropNo,ConfigNo)))**2
          enddo  
       enddo
    enddo


    if(cross_switch.eq.1) then
    allocate(dbvar(NoFeats,NoProps));dbvar=0
    allocate(dwvar(NoProps));dwvar=0
    allocate(inner_means(NoBins,NoBins));inner_means=0
    allocate(inner_count(NoBins,NoBins));inner_count=0
    allocate(inner_row_mean(NoBins,2));inner_row_mean=0
    inner_bvar=0;inner_wvar=0;inner_total_var=0

!    do PropNo=1,NoProps
       !do FeatNo=1,NoFeats
          propno=1
          do Configno=1,NoConfigs
		featno=1
                featno2=2
                curbin=binlocation(FeatNo,ConfigNo)
                curbin2=binlocation(FeatNo2,ConfigNo)
              inner_means(curbin,curbin2)=inner_means(curbin,curbin2)+propertiesarray(PropNo,ConfigNo)
              inner_count(curbin,curbin2)= inner_count(curbin,curbin2)+1
          enddo  
          do BinNo=1,NoBins
             do BinNo2=1,NoBins
                if(inner_count(BinNo,BinNo2).gt.0) then
                   inner_means(BinNo,BinNo2)=inner_means(BinNo,BinNo2)/inner_count(BinNo,BinNo2)
                endif
             enddo
          enddo
          do BinNo=1,NoBins
             inner_row_mean(BinNo,1)=sum(inner_means(BinNo,:))/NoBins
             inner_row_mean(BinNo,2)=sum(inner_means(:,BinNo))/NoBins
          enddo        
          inner_total_mean=sum(inner_means(:,:))/NoBins/NoBins
         print*, inner_row_mean(:,1)
         print*, inner_row_mean(:,2)
         print*,inner_total_mean
         do BinNo=1,NoBins
            inner_bvar(1)=inner_bvar(1)+NoBins*(inner_row_mean(BinNo,1)-inner_total_mean)**2
            inner_bvar(2)=inner_bvar(2)+NoBins*(inner_row_mean(BinNo,2)-inner_total_mean)**2
            do BinNo2=1,NoBins
               inner_wvar=inner_wvar+(inner_means(BinNo,BinNo2)-inner_row_mean(BinNo,1)&
                          &-inner_row_mean(BinNo2,2)+inner_total_mean)**2
               inner_total_var=inner_total_var+(inner_means(BinNo,BinNo2)-inner_total_mean)**2
            enddo 
         enddo
 !        print*,inner_total_var
 !        print*,inner_bvar(:),inner_wvar
    endif




    open(unit=2,file='bin_means.txt',status='unknown')
    do PropNo=1,NoProps
       write(2,*)'  Property',PropNo;write(2,'(3A13)',advance='no')'Bin_Number','Bin_count','Global_mean'
       do FeatNo=1,NoFeats-1
          write(2,'(i4)',advance='no')FeatNo
       enddo
       write(2,'(i4)')NoFeats
       do BinNo=1,NoBins
          write(2,*)BinNo,propmean(PropNo),meanlist(BinNo,:,PropNo)
       enddo
       write(2,*)
    enddo
    close(2)

    open(unit=4,file='bin_counts.txt',status='unknown')
    write(4,'(2A13)',advance='no')'Bin_Number'
    do FeatNo=1,NoFeats-1
       write(4,'(i4)',advance='no')FeatNo
    enddo
    write(4,'(i4)')NoFeats
    do BinNo=1,NoBins
       write(4,*)BinNo,bincount(BinNo,:)
    enddo
    write(4,*)
    close(4)


    open(unit=3,file='ANOVA_1way.txt',status='unknown')
    do PropNo=1,NoProps
       write(3,*)'  Property',PropNo,' Total_SS ',TotalVar(PropNo)
       write(3,*)' Feature      BetweenSS      WithinSS      MSbetween      MSwithin   F' 
       do FeatNo=1,NoFeats
          write(3,*)FeatNo,bvar(FeatNo,PropNo),sum(wvar(:,FeatNo,PropNo)),bvar(FeatNo,PropNo)/(NoBins-1),&
                    &sum(wvar(:,FeatNo,PropNo))/(NoConfigs-NoBins),bvar(FeatNo,PropNo)/(NoBins-1)/&
                    &(sum(wvar(:,FeatNo,PropNo))/(NoConfigs-NoBins))
       enddo
       write(3,*)
    enddo
    close(3)

    Sel_Waters=0
    open(unit=5,file='ordered_features.txt',status='unknown')
    do PropNo=1,NoProps
       max_array(:)=bvar(:,PropNo)
       max_val=maxval(max_array)
       counter=0
       do FeatNo=1,NoFeats
          if(Maxloc(max_array,1).le.(AtomsBefore*3-6)) then
             if(bvar(Maxloc(max_array,1),PropNo)/max_val.gt.thresh) then
                max_array(Maxloc(max_array))=minval(max_array)-1
                counter=counter+1
	!	print*,Maxloc(max_array,1)
             endif
          else
	     if (sum(Sel_Waters).lt.MinWaters) then
		ccwat=ceiling(real(Maxloc(max_array,1)-(AtomsBefore*3-6))/9)
		if(Sel_Waters(ccwat).lt.1) Sel_Waters(ccwat)=1		
                max_array(Maxloc(max_array))=minval(max_array)-1
                counter=counter+1
             endif
	!     print*,Maxloc(max_array,1),sum(Sel_waters),ccwat
          endif
       enddo
       write(5,*)'  Property',PropNo,' Total_SS ',TotalVar(PropNo),'Count ',Counter,'threshold ',thresh
       write(5,*)' Feature      BetweenSS  %_of_max %_or_total' 
       max_array(:)=bvar(:,PropNo)
       max_val=maxval(max_array)
       sel_waters=0
       counter=0
       do FeatNo=1,NoFeats
          if(Maxloc(max_array,1).le.(AtomsBefore*3-6)) then
             if(bvar(Maxloc(max_array,1),PropNo)/max_val.gt.thresh) then
             write(5,*)Maxloc(max_array),bvar(Maxloc(max_array),PropNo),bvar(Maxloc(max_array),PropNo)/max_val*100,&
                       &bvar(Maxloc(max_array),PropNo)/TotalVar(PropNo)*100
                max_array(Maxloc(max_array))=minval(max_array)-1
                counter=counter+1
	!	print*,Maxloc(max_array,1)
             endif
          else
	     if (sum(Sel_Waters).lt.MinWaters) then
		ccwat=ceiling(real(Maxloc(max_array,1)-(AtomsBefore*3-6))/9)
		if(Sel_Waters(ccwat).lt.1) Sel_Waters(ccwat)=1		
             write(5,*)Maxloc(max_array),bvar(Maxloc(max_array),PropNo),bvar(Maxloc(max_array),PropNo)/max_val*100,&
                          &bvar(Maxloc(max_array),PropNo)/TotalVar(PropNo)*100
                max_array(Maxloc(max_array))=minval(max_array)-1
                counter=counter+1
             endif
	!     print*,Maxloc(max_array,1),sum(Sel_waters),ccwat
          endif
       enddo
       write(5,*)
    enddo
    close(5)
!
!    open(unit=4,file='within_bin_SS.txt',status='unknown')
!    do BinNo=1,NoBins
!       write(4,*)wvar(BinNo,:,1)
!    enddo
!    close(4)

  end subroutine anova

end module analysis_of_varience














