!----------------------------------------------------------------------         
      subroutine rdr2daux(io,dout,dinp,mxnx,mxny,nx,ny,lcmprs,clab12)           
!----------------------------------------------------------------------         
!                                                                               
! --- CALPUFF    Version: TNG-7.0.0    Level: 140521          RDR2DAUX          
!                D. Strimaitis                                                  
!                                                                               
! --- PURPOSE:  Read NX*NY words of 2-D real array                              
!               (possibly with compression)                                     
!                                                                               
! --- INPUTS:                                                                   
!               IO - integer     - Unit number of input file                    
!      DINP(nx,ny) - real array  - Data array from file                         
!        MXNX,MXNY - integers    - Dimensions of output data array              
!            NX,NY - integers    - Dimensions of input data array               
!           LCMPRS - logocal     - Compression flag                             
!                                                                               
! --- OUTPUT:                                                                   
!  DOUT(mxnx,mxny) - real array  - Data array to calling routine                
!           CLAB12 - character*12- Data label                                   
!                                                                               
! --- RDR2DAUX called by:  AUX1, RDAUX                                          
! --- RDR2DAUX calls:      UNCOMPRS                                             
!----------------------------------------------------------------------         
      real dinp(nx,ny)                                                          
      real dout(mxnx,mxny)                                                      
      character*12 clab12                                                       
      character*15 clab15                                                       
      logical lcmprs                                                            
                                                                                
      nxy=nx*ny                                                                 
      mxnxy=mxnx*mxny                                                           
      nchar=12                                                                  
                                                                                
      if(mxnxy.LT.nxy) then                                                     
         write(*,*)'ERROR in RDR2DAUX:  actual input array will not ',         &
     &             'fit in output array!'                                       
         write(*,*)'Input  array NX, NY:  ',nx,ny                               
         write(*,*)'Output array NX, NY:  ',mxnx,mxny                           
         stop                                                                   
                                                                                
      elseif(nx.EQ.mxnx .AND. ny.EQ.mxny) then                                  
! ---    Input and output arrays have same shape so simple read if not          
! ---    compressed.   Use DINP as the work array if compressed.                
         if(lcmprs) then                                                        
            read(io) n                                                          
            call UNCOMPRS(dinp,n,io,nxy,nchar,clab12,clab15,dout)               
         else                                                                   
            read(io) clab12,dout                                                
         endif                                                                  
                                                                                
      else                                                                      
! ---    Input array smaller than output array, so read and then                
! ---    transfer by element.  Use output array as the work array               
! ---    if compressed.                                                         
         if(lcmprs) then                                                        
            read(io) n                                                          
            call UNCOMPRS(dout,n,io,nxy,nchar,clab12,clab15,dinp)               
         else                                                                   
            read(io) clab12,dinp                                                
         endif                                                                  
                                                                                
         do j=1,ny                                                              
            do i=1,nx                                                           
               dout(i,j)=dinp(i,j)                                              
            enddo                                                               
         enddo                                                                  
                                                                                
      endif                                                                     
                                                                                
      return                                                                    
      end                                                                       
                                                                                
!----------------------------------------------------------------------         
      subroutine rdi2daux(io,iout,iinp,mxnx,mxny,nx,ny,lcmprs,clab12)           
!----------------------------------------------------------------------         
!                                                                               
! --- CALPUFF    Version: TNG-7.0.0    Level: 140521          RDI2DAUX          
!                D. Strimaitis                                                  
!                                                                               
! --- PURPOSE:  Read NX*NY words of 2-D integer array                           
!               (no compression allowed)                                        
!                                                                               
! --- INPUTS:                                                                   
!               IO - integer     - Unit number of input file                    
!      IINP(nx,ny) - int. array  - Data array from file                         
!        MXNX,MXNY - integers    - Dimensions of output data array              
!            NX,NY - integers    - Dimensions of input data array               
!           LCMPRS - logocal     - Compression flag                             
!                                                                               
! --- OUTPUT:                                                                   
!  IOUT(mxnx,mxny) - int. array  - Data array to calling routine                
!           CLAB12 - character*12- Data label                                   
!                                                                               
! --- RDI2DAUX called by:  AUX1, RDAUX                                          
! --- RDI2DAUX calls:      none                                                 
!----------------------------------------------------------------------         
      integer iinp(nx,ny)                                                       
      integer iout(mxnx,mxny)                                                   
      character*12 clab12                                                       
      character*15 clab15                                                       
      logical lcmprs                                                            
                                                                                
      nxy=nx*ny                                                                 
      mxnxy=mxnx*mxny                                                           
      nchar=12                                                                  
                                                                                
      if(lcmprs) then                                                           
         write(*,*)'ERROR in RDI2DAUX: Data compression option ',              &
     &             'is NOT implemented'                                         
         stop                                                                   
      endif                                                                     
                                                                                
      if(mxnxy.LT.nxy) then                                                     
         write(*,*)'ERROR in RDI2DAUX:  actual input array will not ',         &
     &             'fit in output array!'                                       
         write(*,*)'Input  array NX, NY:  ',nx,ny                               
         write(*,*)'Output array NX, NY:  ',mxnx,mxny                           
         stop                                                                   
                                                                                
      elseif(nx.EQ.mxnx .AND. ny.EQ.mxny) then                                  
! ---    Input and output arrays have same shape so simple read                 
         read(io) clab12,iout                                                   
                                                                                
      else                                                                      
! ---    Input array smaller than output array, so read and then                
! ---    transfer by element.                                                   
         read(io) clab12,iinp                                                   
         do j=1,ny                                                              
            do i=1,nx                                                           
               iout(i,j)=iinp(i,j)                                              
            enddo                                                               
         enddo                                                                  
                                                                                
      endif                                                                     
                                                                                
      return                                                                    
      end                                                                       
                                                                                
!----------------------------------------------------------------------         
      subroutine avgcldmr(qcz,tkz,patmz,qcup,zupbot,zuptop,zface,nzp1,         &
     &                    zbot,ztop,ldb,cldamr,fzcld,cldt,cldp)                 
!----------------------------------------------------------------------         
!                                                                               
! --- CALPUFF    Version: TNG-7.0.0    Level: 140521          AVGCLDMR          
!                D. Strimaitis                                                  
!                                                                               
! --- PURPOSE:  Obtain average LWC from vertical profile data                   
!               Average is "in-cloud", so only LWC>0 contributes                
!               T,P "in-cloud" average is calculated the same way               
!               and assumes that T,P are constant in each layer                 
!                                                                               
! --- INPUTS:                                                                   
!       QCZ(mxnz) - real array - Cloud water mixing ratio (g/kg) profile        
!       TKZ(mxnz) - real array - Temperature (K) profile @ layer heights        
!     PATMZ(mxnz) - real array - Pressure (atm) profile @ layer heights         
!            QCUP - real       - Cloud water mixing ratio (g/kg) aloft          
!         (ZUPBOT - real       - Bottom (mAGL) of cloud layers aloft)           
!         (ZUPTOP - real       - Top (mAGL) of cloud layers aloft)              
!                                ZUPBOT/ZUPTOP not currently used               
!   ZFACE(mxnzp1) - real array - Cell face heights (m) for each layer           
!            NZP1 - integer    - Number of cell face heights (NZ + 1)           
!            ZBOT - real       - Bottom (mAGL) of layer to be averaged          
!            ZTOP - real       - Top (mAGL) of layer to be averaged             
!             LDB - logical    - Debug output flag                              
!                                                                               
!     Parameters:                                                               
!           MXNZ, MXNZP1, IO6                                                   
!                                                                               
! --- OUTPUT:                                                                   
!          CLDAMR - real       - Average in-cloud liquid water mixing           
!                                ratio (g/kg) for layer                         
!           FZCLD - real       - Fraction of interval ZTOP-ZBOT with            
!                                LWC>0                                          
!            CLDT - real       - Associated temperature (K)                     
!            CLDP - real       - Associated pressure (atm)                      
!                                                                               
! --- AVGCLDMR called by:  CHEM                                                 
! --- AVGCLDMR calls:      ZFIND                                                
!----------------------------------------------------------------------         
! --- Include parameters                                                        
      include 'params.puf'                                                      
                                                                                
      real qcz(mxnz),tkz(mxnz)                                                  
      real zface(mxnzp1),patmz(mxnz)                                            
      logical ldb
! --- wangzhm declared and save
      real,intent(in) :: qcup,zupbot,zuptop,zbot,ztop
      integer,intent(in) :: nzp1
      real,intent(out) :: cldamr,fzcld,cldt,cldp
      real,save :: zabot,zatop,zbar,sumqc,sumdz,sumtq,sumpq,&
             & dz,dzq
      integer,save :: ibot,ibar,itop,i1,i2,i

! --- wangzhm omp: set threadprivate                                           
      !$OMP THREADPRIVATE(zabot,zatop,zbar,sumqc,sumdz,sumtq,sumpq,&
      !$OMP& dz,dzq,ibot,ibar,itop,i1,i2,i) 
                                                                                
! --- Averaging height range                                                    
! --- Set averaging limits to model domain faces (keep orig range)              
      zabot=zbot                                                                
      zatop=ztop                                                                
      zabot=MAX(zabot,zface(1))                                                 
      zabot=MIN(zabot,zface(nzp1))                                              
      zatop=MAX(zatop,zface(1))                                                 
      zatop=MIN(zatop,zface(nzp1))                                              
                                                                                
! --- Find the grid layers containing the bottom, average, and top              
      zbar=0.5*(zabot+zatop)                                                    
      call ZFIND(zabot,zface,nzp1,ibot)                                         
      call ZFIND(zbar, zface,nzp1,ibar)                                         
      call ZFIND(zatop,zface,nzp1,itop)                                         
                                                                                
! --- Initially no cloud water                                                  
      cldamr=0.0                                                                
      fzcld=0.0                                                                 
      cldt=tkz(ibar)                                                            
      cldp=patmz(ibar)                                                          
                                                                                
      if(ibot.EQ.itop) then                                                     
! ---    Just 1 layer for average                                               
         cldamr=qcz(ibot)                                                       
         if(cldamr.GT.0.0) fzcld=1.0                                            
         cldt=tkz(ibot)                                                         
         cldp=patmz(ibot)                                                       
                                                                                
      else                                                                      
! ---    2 or more layers                                                       
! ---    Get average in-cloud mixing ratio (do not average zeroes)              
         sumqc=0.0                                                              
! ---    Cloud thickness                                                        
         sumdz=0.0                                                              
! ---    In-cloud T,P weighted by LWC                                           
         sumtq=0.0                                                              
         sumpq=0.0                                                              
                                                                                
! ---    Partial layer from ZABOT to ZFACE(ibot+1)                              
         if(qcz(ibot).GT.0.0) then                                              
            dz=zface(ibot+1)-zabot                                              
            sumdz=sumdz+dz                                                      
            dzq=dz*qcz(ibot)                                                    
            sumqc=sumqc+dzq                                                     
            sumtq=sumtq+dzq*tkz(ibot)                                           
            sumpq=sumpq+dzq*patmz(ibot)                                         
         endif                                                                  
                                                                                
! ---    Partial layer from ZFACE(itop) to ZATOP                                
         if(qcz(itop).GT.0.0) then                                              
            dz=zatop-zface(itop)                                                
            sumdz=sumdz+dz                                                      
            dzq=dz*qcz(itop)                                                    
            sumqc=sumqc+dzq                                                     
            sumtq=sumtq+dzq*tkz(itop)                                           
            sumpq=sumpq+dzq*patmz(itop)                                         
         endif                                                                  
                                                                                
! ---    Remaining are full layers between IBOT+1 and ITOP-1                    
         i1=ibot+1                                                              
         i2=itop-1                                                              
         if(i2.GE.i1) then                                                      
            do i=i1,i2                                                          
               if(qcz(i).GT.0.0) then                                           
                  dz=zface(i+1)-zface(i)                                        
                  sumdz=sumdz+dz                                                
                  dzq=dz*qcz(i)                                                 
                  sumqc=sumqc+dzq                                               
                  sumtq=sumtq+dzq*tkz(i)                                        
                  sumpq=sumpq+dzq*patmz(i)                                      
               endif                                                            
            enddo                                                               
         endif                                                                  
                                                                                
! ---    Compute average                                                        
         if(sumdz.GT.0.0) then                                                  
            cldamr=sumqc/sumdz                                                  
            cldt=sumtq/sumqc                                                    
            cldp=sumpq/sumqc                                                    
            fzcld=sumdz/(ztop-zbot)                                             
         endif                                                                  
      endif                                                                     
                                                                                
! *** Not Active ***                                                            
!c --- Treatment for interaction with clouds at top of model domain             
!c --- Averaging layer must extend into the top model layer                     
!c --- (use linear weight to soften transition)                                 
      weight=0.0                                                                
!      if(qcup.GT.0.0) then                                                     
!         zlo=zface(nzp1-1)                                                     
!         if(ztop.GT.zlo) then                                                  
!            weight=(ztop-zlo)/(zface(nzp1)-zlo)                                
!            weight=AMIN1(weight,1.0)                                           
!         endif                                                                 
!c ---    Bottom of any cloud layer aloft must touch the model-top              
!         if(zupbot.LE.zface(nzp1)) then                                        
!c ---       Use the maximum of the average computed above and the              
!c ---       weighted average of the cloud water aloft                          
!            cldamr=AMAX1(cldamr,weight*qcup)                                   
!         endif                                                                 
!      endif                                                                    
! *** Not Active ***                                                            
                                                                                
      if(ldb) then                                                              
         write(io6,*) 'AVGCLDMR: zbot,ztop,cldamr = ',zbot,ztop,cldamr          
         write(io6,*) '          weight,qcup      = ',weight,qcup               
         write(io6,*) '          T(K),P(atm)      = ',cldt,cldp                 
         write(io6,*) '          zabot,zatop      = ',zabot,zatop               
         write(io6,*) '          ibot,itop        = ',ibot,itop                 
         write(io6,*) '          sumqc,sumdz,fzcld= ',sumqc,sumdz,fzcld         
      endif                                                                     
                                                                                
      return                                                                    
      end                                                                       
! --- wangzhm declared
!     real,intent(in) :: qcup,zupbot,zuptop,zbot,ztop
!     integer,intent(in) :: nzp1
!     real,intent(out) :: cldamr,fzcld,cldt,cldp
!     real :: zabot,zatop,zbar,sumqc,sumdz,sumtq,sumpq,&
!            & dz,dzq
!     integer :: ibot,ibar,itop,i1,i2,i                                                                                
!----------------------------------------------------------------------         
      subroutine makep3d                                                        
!----------------------------------------------------------------------         
!                                                                               
! --- CALPUFF    Version: TNG-7.0.0    Level: 140521           MAKEP3D          
!                D. Strimaitis                                                  
!                                                                               
! --- PURPOSE:  Estimate the 3D pressure (atm) at the face heights              
!                                                                               
! --- INPUTS:                                                                   
!     Common block /GRIDNEST/ variables:                                        
!         ngrid                                                                 
!     Common block /METHD/ variables:                                           
!         nxm(mxmetdom),nym(mxmetdom),nzm,zfacem(mxnzp1)                        
!     Common block /METHR/ variables:                                           
!         temp2d(mxnx,mxny,mxmetdom),rho2d(mxnx,mxny,mxmetdom),                 
!         tmet(mxnx,mxny,mxnz,mxmetdom)                                         
!                                                                               
!     Parameters:                                                               
!           MXNZ, MXNZP1, MXNX, MXNY, MXMETDOM                                  
!                                                                               
! --- OUTPUT:                                                                   
!     Common block /METHR/ variables:                                           
!         pmet(mxnx,mxny,mxnzp1,mxmetdom)                                       
!                                                                               
! --- MAKEP3D called by:  COMP                                                  
! --- MAKEP3D calls:      none                                                  
!----------------------------------------------------------------------         
! --- Include parameters                                                        
      include 'params.puf'                                                      
                                                                                
      include 'gridnest.puf'                                                    
      include 'methd.puf'                                                       
      include 'methr.puf'                                                       
                                                                                
! --- P2=P1* EXP(-(g/R)(z2-z1)/Tv)                                              
! --- g=9.81 m2/s2   R=287.0 J/(kg K)                                           
! --- Tv is average virtual temperature for layer z2-z1 (use layer T)           
      data gbyr/0.0341812/                                                      
                                                                                
! --- Loop over met domains                                                     
      do im=1,ngrid                                                             
                                                                                
! ---    Set surface pressure                                                   
         kz=1                                                                   
         do jy=1,nym(im)                                                        
            do ix=1,nxm(im)                                                     
! ---          kg-Molar volume(m3) at ambient T,P computed from density         
               vkgmol = 28.97/rho2d(ix,jy,im)                                   
               pmet(ix,jy,kz,im)=(22.4141/vkgmol)*                             &
     &                           (temp2d(ix,jy,im)/273.15)                      
            enddo                                                               
         enddo                                                                  
                                                                                
! ---    Set pressure profiles from temperature profiles                        
         do kz=2,nzm+1                                                          
            kzm1=kz-1                                                           
            do jy=1,nym(im)                                                     
               do ix=1,nxm(im)                                                  
                  f=-gbyr*(zfacem(kz)-zfacem(kzm1))                             
                  pmet(ix,jy,kz,im)=pmet(ix,jy,kzm1,im)*                       &
     &                              EXP(f/tmet(ix,jy,kzm1,im))                  
               enddo                                                            
            enddo                                                               
         enddo                                                                  
                                                                                
      enddo                                                                     
                                                                                
      return                                                                    
      end                                                                       
