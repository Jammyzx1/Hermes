! Created by From Stuart Davie 2015 Univeristy of Manchester
!****************************************************************************
!
! Originally created by Stuart.
! Version 1.1
! CHANGE LOG 
! Version 1  : To calculate statistics from krigining multipoles. 
! Version 1.1: Modular format and incorporation in the Hermes program.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    module krigingcoefficient

      implicit none

    contains

      subroutine krigcoefficent
        integer i,j,NoPreds,linecount,ios,ios2,ios3,nprops,propcount
        character*100 Title,dummychar,dum,dum2,atomname,filename,outnames(25),curprop
        double precision summax, summax2, summax3, sumc1, sumc2, sumc3, stdev1, stdev2
        double precision true, predicted, meanabserror,squarederror,rmse,meantrue,meanpredicted
        double precision relabserror,relsquarederror,relrmse,qsquared,trueXpred,true2,pred2   
        double precision k,kprime,r2,r02,r02prime,simpabserror,simpsquarederror,simprmse
        double precision r2num,r2den1,r2den2,r02num,r02primenum,r02den,r02primeden,rm2,meantrueerr
        double precision rm2prime,rm2mean,rm2delta,stdevtrue,stdevpred,errorstdev,maxp,minp,rangep

                
        open(unit=3,file='FINPUT.txt',status='old',iostat=ios2)
        if(ios2.ne.0) then
           print*,"Can't find FINPUT.txt";stop
        endif
    
        read(3,*)TITLE
        rewind(3)
    
        dum='predictions'
        linecount=0
        do while (.true.)
           read (3, *,iostat=ios) dum2
           linecount=linecount+1
           if(trim(dum2).eq.(trim(dum)))   goto 55
           if(ios.ne.0) exit
        enddo
55      continue   
        rewind(3)
        do i=1,linecount-1
           read(3,*)
        enddo
        read(3,*)dummychar,nopreds
        rewind(3)
        dum='atoms'
        linecount=0
        do while (.true.)
           read (3, *,iostat=ios) dum2
           linecount=linecount+1
           if(trim(dum2).eq.(trim(dum)))   goto 56
           if(ios.ne.0) exit
        enddo
56      continue   
        rewind(3)
        do i=1,linecount
           read(3,*)
        enddo
        read(3,*)dummychar,atomname
        rewind(3)
        dum='nproperties'
        linecount=0
        do while (.true.)
           read (3, *,iostat=ios) dum2
           linecount=linecount+1
           if(trim(dum2).eq.(trim(dum)))   goto 57
           if(ios.ne.0) exit
        enddo
57      continue   
        rewind(3)
        do i=1,linecount-1
           read(3,*)
        enddo
        read(3,*)dummychar,nprops
        close(3)


        outnames(1)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q00_PREDICTIONS.txt'
        outnames(2)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q10_PREDICTIONS.txt'
        outnames(3)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q11c_PREDICTIONS.txt'
        outnames(4)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q11s_PREDICTIONS.txt'
        outnames(5)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q20_PREDICTIONS.txt'
        outnames(6)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q21c_PREDICTIONS.txt'
        outnames(7)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q21s_PREDICTIONS.txt'
        outnames(8)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q22c_PREDICTIONS.txt'
        outnames(9)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q22s_PREDICTIONS.txt'
        outnames(10)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q30_PREDICTIONS.txt'
        outnames(11)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q31c_PREDICTIONS.txt'
        outnames(12)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q31s_PREDICTIONS.txt'
        outnames(13)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q32c_PREDICTIONS.txt'
        outnames(14)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q32s_PREDICTIONS.txt'
        outnames(15)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q33c_PREDICTIONS.txt'
        outnames(16)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q33s_PREDICTIONS.txt'
        outnames(17)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q40_PREDICTIONS.txt'
        outnames(18)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q41c_PREDICTIONS.txt'
        outnames(19)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q41s_PREDICTIONS.txt'
        outnames(20)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q42c_PREDICTIONS.txt'
        outnames(21)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q42s_PREDICTIONS.txt'
        outnames(22)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q43c_PREDICTIONS.txt'
        outnames(23)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q43s_PREDICTIONS.txt'
        outnames(24)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q44c_PREDICTIONS.txt'
        outnames(25)='./'//trim(title)//'_kriging_'//trim(atomname)//'_q44s_PREDICTIONS.txt'

        do propcount=1,nprops
           
           open(unit=10,file='curpropno',status='unknown')
           write(10,*)propcount;close(10)
           open(unit=10,file='curpropno',status='old')
           read(10,*)curprop;close(10)
           
           open(unit=2,file=outnames(propcount),status='old',iostat=ios3)
           if(ios3.ne.0) then
              print*,"Can't find",filename ;stop
           endif
       
           open(unit=1,file='prediction_evaluation_metrics_'//trim(curprop)//'.txt')
           do j=1,5 
              read(2,*)
           enddo
                 
           meanabserror=0;squarederror=0;rmse=0;qsquared=0;trueXpred=0;true2=0  
           meantrue=0;meanpredicted=0;relabserror=0;relsquarederror=0;relrmse=0  
           pred2=0;r2=0;r02=0;r02prime=0;simpabserror=0;simpsquarederror=0;simprmse=0
           r2num=0;r2den1=0;r2den2=0;r02num=0;r02primenum=0;r02den=0;r02primeden=0
           rm2=0;meantrueerr=0;rm2prime=0;rm2mean=0;rm2delta=0;stdevtrue=0;stdevpred=0
           errorstdev=0
           do j=1,NoPreds
              read(2,*)true,predicted
              if(j.eq.1) then 
                 maxp=true;minp=true; endif
                 if(true.gt.maxp)maxp=true
                 if(true.lt.minp)minp=true
                 meanabserror=meanabserror+abs(true-predicted)/NoPreds
                 meantrueerr=meantrueerr+(true-predicted)/NoPreds
                 squarederror=squarederror+(true-predicted)**2/nopreds
                 meantrue=meantrue+true/NoPreds
                 meanpredicted=meanpredicted+predicted/NoPreds
                 true2=true2+true*true
                 pred2=pred2+predicted*predicted
                 trueXpred=trueXpred+true*predicted
              enddo
              rmse=sqrt(squarederror)
              k=trueXpred/pred2
              kprime=trueXpred/true2
              rewind(2)
              do j=1,5
                 read(2,*)
              enddo
              do j=1,NoPreds
                 read(2,*)true,predicted
                 simpabserror=simpabserror+abs(meantrue-predicted)/NoPreds
                 simpsquarederror=simpsquarederror+(meantrue-true)**2/nopreds
                 stdevpred=stdevpred+(meanpredicted-predicted)**2/Nopreds
                 errorstdev=errorstdev+(true-predicted-meantrueerr)**2/nopreds
                 r2num=r2num+((true-meantrue)*(predicted-meanpredicted))
                 r2den1=r2den1+(predicted-meanpredicted)**2
                 r2den2=r2den2+(true-meantrue)**2
                 r02num=r02num+(true-k*predicted)**2
                 r02primenum=r02primenum+(predicted-kprime*true)**2
              enddo
              r2num=r2num**2
              r2=r2num/r2den1/r2den2
              r02=1-(r02num/r2den2)
              r02prime=1-(r02primenum/r2den1)
              errorstdev=sqrt(errorstdev)
              stdevpred=sqrt(stdevpred)
              stdevtrue=sqrt(simpsquarederror)
              qsquared=1-squarederror/simpsquarederror
              rm2=r2*(1-sqrt((r2-r02)))
              rm2prime=r2*(1-sqrt((r2-r02prime)))
              rm2mean=(rm2+rm2prime)/2
              rm2delta=abs(rm2-rm2prime)
              
          
          
              write(1,*)'true_range                   =',maxp-minp
              write(1,*)'mean_true_value              =',meantrue
              write(1,*)'true_value_stdev             =',stdevtrue
              write(1,*)'mean_predicted_value         =',meanpredicted
              write(1,*)'pred_value_stdev             =',stdevpred
              write(1,*)'mean_err                     =',meantrueerr
              write(1,*)'stdev_of_error               =',errorstdev                
              write(1,*)'************'
              write(1,*)'mean_absolute_error          =',meanabserror
              write(1,*)'mean_squared_error           =',squarederror
              write(1,*)'root_mean_absolute_error     =',rmse
              write(1,*)'************'
              write(1,*)'q^2__(target>0.5)            =',qsquared
              write(1,*)'************'
              write(1,*)'average_rm^2__(target>0.5)   =',rm2mean       
              write(1,*)'rm^2_delta__(target<0.2)     =',rm2delta       
              write(1,*)'r^2                          =',r2
              write(1,*)'rm^2                         =',rm2
              write(1,*)"reverse_rm^2                 =",rm2prime
              write(1,*)''
              write(1,*)''
              write(1,*)''
              write(1,*)'For more discussion see'
              write(1,*)'Further exploring rm2 metrics for validation of QSPR models'
              write(1,*)'Ojha et al. 2011. doi:10.1016/j.chemolab.2011.03.011'
              write(1,*)"http://www.sciencedirect.com/science/article/pii/S016974391100061X"
          
              close(unit=2) 
              close(unit=1)
           enddo
        
      end subroutine krigcoefficent
      
    end module krigingcoefficient
