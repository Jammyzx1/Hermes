!Created by Stuart Davie
!********************************************************************************************************
! Version 1  : Takes set of true values and predictions of some property and calculates various prediction
!	       metrics. Searches for FINPUT to find output of program FEREBUS, to adequately account for 
!	       formatting. If FINPUT.txt not present, will ask user for file instead
! CHANGE LOG 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    program GetStats

    implicit none
    integer 		::	i,j,NoPreds,linecount,ios,NoProps,PropCount,StartProp,F_Switch,DummyInt
    character*100 	::	Title,DummyChar,outnames(25),InName,var_name,pass_name
    character 		::	curprop*2,atomname*5
    double precision  	::	meanabserror,meansquarederror,rmse,meantrue,meanpredicted,dbdum
    double precision 	::	relabserror,relsquarederror,relrmse,qsquared,trueXpred,true2,pred2   
    double precision 	::	k,kprime,r2,r02,r02prime,simpabserror,simpsquarederror,simprmse
    double precision 	::	r2num,r2den1,r2den2,r02num,r02primenum,r02den,r02primeden,rm2,meantrueerr
    double precision 	::	rm2prime,rm2mean,rm2delta,stdevtrue,stdevpred,errorstdev,MaxT,MinT,RangeT
    double precision, allocatable :: maearray(:),ordmaearray(:),true(:), predicted(:)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    F_Switch=0       									!As program will mostly be used for Popelier group
    open(unit=3,file='FINPUT.txt',status='old',iostat=ios)				!Ferebus output analysis, this checks for the Ferebus
    if(ios.ne.0) then									!input file
       print*,"Can't find FINPUT.txt";F_Switch=1;					!
    endif										!
    close(3)										!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(F_Switch.eq.0)then					!IF Ferebus input file found:
!---------------------------------------------------------------!
        open(unit=3,file='FINPUT.txt',status='old')		!Get prediction file name
        read(3,*)title						!
        close(3)						! 
!---------------------------------------------------------------! 
        var_name='predictions'              			!Get number of predictions
        call FINPUT_read(var_name,NoPreds,DummyChar,1)		!
        print*,var_name,NoPreds					!
!---------------------------------------------------------------!
        var_name='nproperties'              			!Get number of properties
        call FINPUT_read(var_name,NoProps,DummyChar,1)		!
        print*,var_name,NoProps					!
!---------------------------------------------------------------!
        var_name='starting_properties'              		!Get starting property
        call FINPUT_read(var_name,StartProp,DummyChar,1)	!
        print*,var_name,StartProp				!
!---------------------------------------------------------------!
        var_name='atoms'              				!Get atom name (to read correct prediction file)
        call FINPUT_read(var_name,DummyInt,pass_name,2)		!
        AtomName=trim(Pass_Name)				!
        print*,var_name,trim(AtomName)				!
!---------------------------------------------------------------!
    else							!ELSE
!---------------------------------------------------------------!
        call system('ls')					!
        print*,'Input filename?'				!Ask for prediction file name
        read*,InName						!
        open(unit=3,file=trim(InName),status='old',iostat=ios)	!Open it
        if(ios.ne.0) then					!
           print*,"Can't find "//trim(InName);			!
        endif							!
        linecount=0						!
        do while(.true.)					!
           read(3,*,iostat=ios)dbdum				!Count number of predictions
           if(ios.ne.0) goto 59					!
              linecount=linecount+1				!
        enddo							!
59      continue						!
        close(3)						!
        StartProp=1;NoProps=1;nopreds=linecount			!
	print*,'No of predictions =',nopreds			!Set required parameters
        if(linecount.eq.0) then					!
           print*,'reading error';stop				!
        endif     						!
    endif							!ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(maearray(nopreds));allocate(ordmaearray(nopreds));	!
    allocate(true(nopreds));allocate(predicted(nopreds));	!Allocate arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(F_Switch.eq.0) then										!IF Using Ferebus output:
	    outnames(1)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q00_PREDICTIONS.txt'		!
	    outnames(2)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q10_PREDICTIONS.txt'		!Prepare array of output file names
	    outnames(3)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q11c_PREDICTIONS.txt'		!
	    outnames(4)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q11s_PREDICTIONS.txt'		!
	    outnames(5)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q20_PREDICTIONS.txt'		!
	    outnames(6)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q21c_PREDICTIONS.txt'		!
	    outnames(7)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q21s_PREDICTIONS.txt'		!
	    outnames(8)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q22c_PREDICTIONS.txt'		!
	    outnames(9)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q22s_PREDICTIONS.txt'		!
	    outnames(10)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q30_PREDICTIONS.txt'		!
	    outnames(11)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q31c_PREDICTIONS.txt'	!
	    outnames(12)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q31s_PREDICTIONS.txt'	!
	    outnames(13)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q32c_PREDICTIONS.txt'	!
	    outnames(14)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q32s_PREDICTIONS.txt'	!
	    outnames(15)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q33c_PREDICTIONS.txt'	!
	    outnames(16)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q33s_PREDICTIONS.txt'	!
	    outnames(17)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q40_PREDICTIONS.txt'		!
	    outnames(18)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q41c_PREDICTIONS.txt'	!
	    outnames(19)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q41s_PREDICTIONS.txt'	!
	    outnames(20)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q42c_PREDICTIONS.txt'	!
	    outnames(21)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q42s_PREDICTIONS.txt'	!
	    outnames(22)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q43c_PREDICTIONS.txt'	!
	    outnames(23)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q43s_PREDICTIONS.txt'	!
	    outnames(24)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q44c_PREDICTIONS.txt'	!
	    outnames(25)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q44s_PREDICTIONS.txt'	!
	endif												!ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	do PropCount=StartProp,NoProps+StartProp-1						!Do for each property of interest:
!!!******************************MAIN LOOP******************************************************!
		if (F_Switch.eq.0) then								!IF using Ferebus output:
		   if(PropCount.lt.10) then							!\
		      write(curprop,"(I1)")PropCount						! \
		   else										!  Prepare string format for output
		      write(curprop,"(I2)")PropCount						! /
		   endif									!/
		   open(unit=2,file=trim(outnames(PropCount)),status='old',iostat=ios)		!\
                   if(ios.ne.0) then								! \
                      print*,"Can't find "//outnames(PropCount);cycle				!  |
                   endif									!  |Open predictions file
		   do j=1,5 									!  |and move to start of predictions
                      read(2,*)									!  |
                   enddo									!  |
		else										!  |
		   write(curprop,"(I1)")PropCount						!  |
		   open(unit=2,file=trim(InName),status='old')					! /
		endif										!/
!-----------------------------------------------------------------------------------------------!
												!
                meanabserror=0;meansquarederror=0;rmse=0;qsquared=0;trueXpred=0;true2=0  	!
                meantrue=0;meanpredicted=0;relabserror=0;relsquarederror=0;relrmse=0  		!Initialise metric variables
                pred2=0;r2=0;r02=0;r02prime=0;simpabserror=0;simpsquarederror=0;simprmse=0	!
                r2num=0;r2den1=0;r2den2=0;r02num=0;r02primenum=0;r02den=0;r02primeden=0		!
                rm2=0;meantrueerr=0;rm2prime=0;rm2mean=0;rm2delta=0;stdevtrue=0;stdevpred=0	!
                errorstdev=0									!
!-----------------------------------------------------------------------------------------------!
                do j=1,NoPreds									!
                    read(2,*)true(j),predicted(j)						!Read predictions
                enddo										!
                close(2)									!
!-----------------------------------------------------------------------------------------------!Start calculating properties:
		MaxT=maxval(true);MinT=minval(true)						!
		do j=1,NoPreds									!
                   meanabserror=meanabserror+abs(true(j)-predicted(j))/NoPreds			!Mean absolute error
                   meantrueerr=meantrueerr+(true(j)-predicted(j))/NoPreds			!Mean true error
                   meansquarederror=meansquarederror+(true(j)-predicted(j))**2/nopreds		!Mean squared error
                   meantrue=meantrue+true(j)/NoPreds						!Mean true value
                   meanpredicted=meanpredicted+predicted(j)/NoPreds				!Mean predicted value
                   true2=true2+true(j)*true(j)							!Sum of true squared
                   pred2=pred2+predicted(j)*predicted(j)					!Sum of predicted squared
                   trueXpred=trueXpred+true(j)*predicted(j)					!Sum of true*predicted
                enddo										!
                k=trueXpred/pred2								!rm2 parameter k
                kprime=trueXpred/true2								!rm2 parameter kprime
!-----------------------------------------------------------------------------------------------!
                do j=1,NoPreds									!Calculate properties that required precomputed means and sums
                    simpabserror=simpabserror+abs(meantrue-predicted(j))/NoPreds		!Absolute error between mean and predicted
                    simpsquarederror=simpsquarederror+(meantrue-true(j))**2/nopreds		!Absolute error between mean and true
                    stdevpred=stdevpred+(meanpredicted-predicted(j))**2/Nopreds			!Variance of predicted
                    errorstdev=errorstdev+(true(j)-predicted(j)-meantrueerr)**2/nopreds		!Variance of true errors
                    r2num=r2num+((true(j)-meantrue)*(predicted(j)-meanpredicted))		!Sum of (true-mean_true) * (pred - mean_pred) for Pearson coefficient
                    r2den1=r2den1+(predicted(j)-meanpredicted)**2				!Sum of predicted-mean_predicted squared
                    r2den2=r2den2+(true(j)-meantrue)**2						!Sum of true-mean_true squared
                    r02num=r02num+(true(j)-k*predicted(j))**2					!rm2 parameter
                    r02primenum=r02primenum+(predicted(j)-kprime*true(j))**2			!rm2 parameter
                enddo										!
!-----------------------------------------------------------------------------------------------!
                RangeT=MaxT-MinT                                                                !
                r2num=r2num**2									!
                r2=r2num/r2den1/r2den2								!Coefficient of determination (r^2)
                errorstdev=sqrt(errorstdev)							!Std dev of error (sqrt of error variance)
                rmse=sqrt(meansquarederror)							!Root mean square error
                stdevpred=sqrt(stdevpred)							!Std dev of pred (sqrt of pred variance)
                stdevtrue=sqrt(simpsquarederror)						!Std dev of true (sqrt of true variance)
                qsquared=1-meansquarederror/simpsquarederror					!q^2 coefficient
                r02=1-(r02num/r2den2)								!rm2 parameter
                r02prime=1-(r02primenum/r2den1)							!rm2 parameter
                rm2=r2*(1-sqrt((r2-r02)))							!rm2 parameter
                rm2prime=r2*(1-sqrt((r2-r02prime)))						!rm2 parameter
                rm2mean=(rm2+rm2prime)/2							!rm2 parameter
                rm2delta=abs(rm2-rm2prime)							!rm2 parameter
!-----------------------------------------------------------------------------------------------!
		do j=1,NoPreds									!Order arrays so percentiles can be calculated:
		   maearray(j)=abs(true(j)-predicted(j))					!Mean absolute error array
		enddo										!
		do j=1,Nopreds									!
		   ordmaearray(j)=minval(maearray,1)						!
		   maearray(minloc(maearray,1))=maearray(minloc(maearray,1))+maxval(maearray,1)	!
		enddo										!
!-----------------------------------------------------------------------------------------------!
		open(unit=1,file='prediction_evaluation_metrics_'//trim(curprop)//'.txt')	!Open output file	
                write(1,*)'true_Range                   =',RangeT				!Write output for each property	
                write(1,*)'mean_true_value              =',meantrue				!
                write(1,*)'true_value_stdev             =',stdevtrue				!
                write(1,*)'mean_predicted_value         =',meanpredicted			!
                write(1,*)'pred_value_stdev             =',stdevpred				!
                write(1,*)'mean_err                     =',meantrueerr				!
                write(1,*)'stdev_of_error               =',errorstdev                		!
                write(1,*)'************'							!
                write(1,*)'mean_absolute_error          =',meanabserror				!
                write(1,*)'mean_squared_error           =',meansquarederror			!
                write(1,*)'root_mean_squared_error      =',rmse					!
                write(1,*)'mae_%_Range		      	=',meanabserror/RangeT*100	        !
                write(1,*)'rmse_%_Range      		=',rmse/RangeT*100			!
                write(1,*)'************'							!
                write(1,*)'MAX_absolute_error      =',ordmaearray(nopreds)			!
                write(1,*)'99%_absolute_error      =',ordmaearray((99*nopreds/100))		!
                write(1,*)'95%_absolute_error      =',ordmaearray((95*nopreds/100))		!
                write(1,*)'90%_absolute_error      =',ordmaearray((90*nopreds/100))		!
                write(1,*)'50%_absolute_error      =',ordmaearray((50*nopreds/100))		!
                write(1,*)'************'							!
                write(1,*)'q^2__(target>0.5)            =',qsquared				!
                write(1,*)'r^2__(target>0.5)            =',r2					!
                write(1,*)'************'							!
                write(1,*)'average_rm^2__(target>0.5)   =',rm2mean       			!
                write(1,*)'rm^2_delta__(target<0.2)     =',rm2delta       			!
                write(1,*)'rm^2                         =',rm2					!
                write(1,*)"reverse_rm^2                 =",rm2prime				!
                write(1,*)''									!
                write(1,*)''									!
                write(1,*)''									!	
                write(1,*)'For more discussion see'						!	
                write(1,*)'Further exploring rm2 metrics for validation of QSPR models'		!	
                write(1,*)'Ojha et al. 2011. doi:10.1016/j.chemolab.2011.03.011'		!	
                write(1,*)"http://www.sciencedirect.com/science/article/pii/S016974391100061X"	!	
            close(unit=1)									!Close output file
!!!******************************MAIN LOOP END**************************************************!
        enddo
        
    end program GetStats
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine FINPUT_read(input_name,input_value,input_char,switch)	!
    character*50 input_name,dummy,input_char,char_hold			!
    integer input_value,linecount,ios,switch				!
									!
    open(unit=3,file='FINPUT.txt',status='old')				!
    linecount=0								!
    do while (.true.)							!
       read (3, *,iostat=ios) dummy					!
       linecount=linecount+1						!
       if(trim(dummy).eq.(trim(input_name)))   goto 100			!
       if(ios.ne.0) exit						!
    enddo								!
100 continue   								!
    rewind(3)								!
    if(switch.eq.1) then						!
    do i=1,linecount-1							!
       read(3,*)							!
    enddo								!
    read(3,*)dummy,input_value						!
    elseif(switch.eq.2) then						!
    do i=1,linecount							!
       read(3,*)							!
    enddo								!
       read(3,*)dummy,char_hold						!
       input_char=trim(char_hold)					!
    endif								!
    close(3)								!
    return								!
    end subroutine							!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
