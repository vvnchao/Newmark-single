c    This program is for newmark analysis
c    with:
c
c    strike and dipping angle of slope
c    and frictional parameters
c
c    Input file is CWB txt format: PA03_20210418141140.00_CWB.txt
c
c
c    user should prepare the parameter list file:
c
c    xxxx.yyyy.mmdd.dop.txt
c
c    run the main program:
c
c    > ./Newmark-single PA03_20210418141140.00_CWB.txt
c   
c-----------------------------------------------------------------------------------
c     September 01 2023 created by vvn weian chao, associate professor
c     Email: vvnchao@gmail.com; vvnchao@nycu.edu.tw
c     Department of Civil Engineering, National Yang Ming Chiao Tung University
c-----------------------------------------------------------------------------------
      program newmark
        parameter(npmax=327680)
c--------------------------
c   For input txt:
c--------------------------        
        real t(npmax),acv(npmax),acn(npmax),ace(npmax)
        real dt
        character(len=100)  fn   !--  PA03_20210418141140.00_CWB.txt    
      
c-------------------------------
c   For projected accelerations:
c-------------------------------
        real aN(npmax)  !--Acceleration normal to sliding plane (positive upward)
        real aH(npmax)  !--Horizontal acceleration along slope dipping direction
        real aD(npmax)  !--Acceleration in the down-dip direction (positive downward) 
      
c------------------------
c   For Newmark analysis:
c------------------------
        real S(npmax)    !-- acceleration of the sliding block
        real vel(npmax)  !-- cumulative velocity
        real disp(npmax) !-- cumulative displacement
        real mud(npmax)
        real v0,d0,mu
      
      
c-------------------------
c   For Option Parameters:
c-------------------------
        CHARACTER(len=100) par,par1      
        REAL str  !-- strike of slope
        REAL dip  !-- dipping angle of slope
        REAL mur  !-- residual friction coef.
        REAL mup  !-- peak friction coef.
        REAL mus  !-- steady-state friction coef.
        REAL Dc   !-- critical distance [m]
        REAL Vc   !-- critical velocity [m/s]
        REAL muss
        LOGICAL lstr,ldip,lmur,lmup,lmus,lDc,lVc
      
c--------------------
c   For Command Mode:
c--------------------
        INTEGER*4 bat,system
        CHARACTER(len=300) cmd
        
c--------------------------
c   For call time function:
c--------------------------
        CHARACTER(len=40) datex,timex           
      
c--------------------------
c   parameters:
c--------------------------
        double precision pi,rad
        double precision dipc,dips,strc,strs
        
        pi=atan(1.0)*4.0
	rad=pi/180.
     
      
         call date_and_time(datex,timex)
         write(*,'(/a,a4,a,a2,a,a2,a,a2,a,a2,a,a4,/)')
     1      '  Time progression start. Current date and time:   ',
     2      datex(1:4),'-',datex(5:6),'-',datex(7:8),'.  ',
     3      timex(1:2),':',timex(3:4),':',timex(5:8)      
      
      
c   INPUTS:
c---------------------------
c   Get number of arguments:
       narg=iargc()
       if(narg.lt.1) THEN
        call usage
       endif      
      
       call getarg(1,par)

       open(1,file=par,status='old')
       do i=1,11
       read(1,*)  !-- header
       enddo
       i=1
       do
       read(1,*,iostat=istat) t(i),acv(i),acn(i),ace(i)  !-- vertical north east
       if(istat.ne.0) exit
       i=i+1
       enddo
       close(1)
       np=i-1
       print*,'Total data points: ',np
       if(np.gt.npmax) stop 'ERROR: np > npmax!! '
       dt=t(2)-t(1)  !-- sampling interval
      
c   Baseline Correction: high-pass filter with corner frequency 0.075 Hz
       call IIRFILT(acv(1:np),np,'BUTTER',3,'HP',0.075,0.0,dt,1)       
       call IIRFILT(acn(1:np),np,'BUTTER',3,'HP',0.075,0.0,dt,1)
       call IIRFILT(ace(1:np),np,'BUTTER',3,'HP',0.075,0.0,dt,1)      
       print*,' filtering...'
      



c   DEFAULT VALUE:
c--------------------
       lstr=.false.
       ldip=.false.
       lmur=.false.
       lmup=.false.
       lmus=.false.
       lDc=.false.
       lVc=.false.
      
      
c   Get arguments (optional):
c----------------------------
       DO iarg=2,narg
        call getarg(iarg,par)
        l=len(par)
        IF(par(1:4).eq.'str=') THEN 
        read(par(5:l),*) str
        lstr=.true.
        ELSE IF(par(1:4).eq.'dip=') THEN
        read(par(5:l),*) dip
        ldip=.true.
        ELSE IF(par(1:4).eq.'mur=') THEN
        read(par(5:l),*) mur
        lmur=.true.
        ELSE IF(par(1:4).eq.'mup=') THEN
        read(par(5:l),*) mup
        lmup=.true.
        ELSE IF(par(1:4).eq.'mus=') THEN
        read(par(5:l),*) mus
        lmus=.true.
        ELSE IF(par(1:3).eq.'Dc=') THEN
        read(par(4:l),*) Dc
        lDc=.true.
        ELSE IF(par(1:3).eq.'Vc=') THEN
        read(par(4:l),*) Vc
        lVc=.true.        
        ELSE
        call usage
        stop
        ENDIF
       ENDDO      
      
       if(.not.lstr) stop 'ERROR: check str='
       if(.not.ldip) stop 'ERROR: check dip='      
       if(.not.lmur) stop 'ERROR: check mur='         
       if(.not.lmup) stop 'ERROR: check mup='       
       if(.not.lmus) stop 'ERROR: check mus='       
       if(.not.lDc) stop 'ERROR: check Dc='       
       if(.not.lVc) stop 'ERROR: check Vc='      
      
       print*,'used str=',str
       print*,'used dip=',dip
       print*,'used mur=',mur
       print*,'used mup=',mup
       print*,'used mus=',mus
       print*,'used Dc=',Dc
       print*,'used Vc=',Vc
      
	dipc=dcos(dip*rad)
        dips=dsin(dip*rad)
	strc=dcos(str*rad)
	strs=dsin(str*rad)       
      
c   Projected accelerations:
c---aH: Horizontal acceleration along slope dipping direction
c---aN: Acceleration normal to sliding plane (positive upward)
c---aD: Acceleration in the down-dip direction (positive downward) 
       aH=0.0
       aN=0.0
       aD=0.0
       do i=1,np
         aD(i)=ace(i)*dipc*strc-acn(i)*dipc*strs-acv(i)*dips
         aN(i)=ace(i)*dips*strc-acn(i)*dips*strs+acv(i)*dipc
         aH(i)=ace(i)*strc-acn(i)*strs
       enddo
      
      
c   output
       open(1,file='newmark_oup.txt',status='unknown')
       write(1,*) 'time,Vertical,North,East,aH,aD,aN,S,vel,disp,mu'

      
c   Start to Newmark analysis
       v0=0.       !-- initial velocity [m/s]
       d0=0.       !-- initial displacement [m]
       S=0.
       vel=0.
       disp=0.
       mud=0.
       mu=mur
      
       do i=1,np
         S(i)=(980.0*dips-aD(i))-mu*(980.0*dipc+aN(i))
         if(S(i).gt.0.0.or.v0.gt.0.0) then
           v0=v0+S(i)*dt/100.0
           if(v0.lt.0.0) v0=0.0
           d0=d0+v0*dt+S(i)*dt*dt/2.0
c           print*,i*dt,S(i),v0,d0
           if(v0.gt.Vc) then
             muss=mus+(mup-mus)*exp(-(v0-Vc)/Vc)
             mu=muss+(mup-muss)*exp(log(0.05)*(d0/Dc))
c             print*,v0,Vc,muss,mu,d0,Dc
           else
             mu=mur
           endif
         endif
         vel(i)=v0
         disp(i)=d0
         mud(i)=mu
         write(1,44) t(i),acv(i),acn(i),ace(i),aH(i),aD(i),aN(i),S(i),ve
     &l(i),disp(i),mud(i)

       enddo
      
      
   44 FORMAT(f8.4,7(1x,f10.3),2f15.3,f8.4)   
      
       close(1)
       print*,' '
       print*,'OUTPUT: newmark_oup.txt'
       print*,' '
      
      end program newmark
      
      
c---------------------------------------------------------------------------       
       SUBROUTINE usage
        write(*,*)' '
        write(*,*)'--------------------------------------------------'
        write(*,*)' '
        write(*,*)'PLOT DOMAIN: '
        write(*,*)' '
        write(*,*)' '        
        write(*,*)'USAGE:./Newmark-single PA03_20210418141140.00_CWB.txt
     & parameter_list '
        write(*,*)' '
        write(*,*)'# parameter_list: '
        write(*,*)' '
        write(*,*)'str=  : strike of slope'
        write(*,*)'dip=  : dipping angle of slope'
        write(*,*)'mur=  : residual friction coef.'
        write(*,*)'mup=  : peak friction coef.'
        write(*,*)'mus=  : steady-state friction coef.'
        write(*,*)'Dc=   : critical distance [m]'
        write(*,*)'Vc=   : critical velocity [m/s]'
        write(*,*)' '
        write(*,*)'Example: '
        write(*,*)'./Newmark-single PA03_20210418141140.00_CWB.txt str=5
     &4.0 dip=24.0 mur=0.44 mup=0.55 mus=0.15 Dc=11.0 Vc=0.396 '
        write(*,*)' '
        write(*,*)'AUTHOR: vvn Weian Chao, Sep. 01 2023'
        write(*,*)' '
        write(*,*)'NOTES: Please, send bugs, comments, improvements'
        write(*,*)' .... to vvnchao@gmail.com;vvnchao@nycu.edu.tw'
        write(*,*)' '
        write(*,*)'WEBSITE:  '
        write(*,*)' '
         STOP
        RETURN
        
       END SUBROUTINE usage      
