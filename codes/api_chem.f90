!----------------------------------------------------------------------         
! --- API_CHEM -- CALPUFF subroutines added for API (MCHEM=6,7)                 
!----------------------------------------------------------------------         
!                                                                               
! --- CALPUFF    Version: TNG-7.1.0    Level: 140521          API_CHEM          
!                                                                               
! --- PURPOSE: This collection of routines was prepared for API by              
!              AER to add chemical transformation and aerosol options           
!              to CALPUFF.  They are bundled as an include-file to              
!              facilitate future modular upgrades.  Note that another           
!              set of subroutines (isorropia.for) has not been placed           
!              in api_chem.for and must be 'included' to complete               
!              the code for the related options.                                
!                                                                               
! --- IMPLEMENTATION:  Several modifications made to these new options          
!              have been assigned either new user-selected modeling             
!              option variables (changed via the control file) or               
!              local logical operators that may be switched in an               
!              individual subroutine to restore the original logic.             
!              CALPUFF will need to be recompiled if a local logical            
!              variable is changed.                                             
!                                                                               
!              The following variable settings will restore the original        
!              features and logic of the API code (for testing):                
!                                                                               
!              Control File                                                     
!                     MLWC = 0 Do not read gridded cloud water from file        
!                     MNH3 = 1 Read monthly ammonia vertical profiles           
!                  MAVGNH3 = 0 Do not average ammonia vertical profiles         
!                              across puff                                      
!                   RNITE1 = 0.0 for no heterogeneous SO2 transformation        
!              Local Logical in Subroutine                                      
!                     L_KGPCM = .FALSE. in CHEMRIV6, CHEMRIV7                   
!                 L_TNO3FLOOR = .FALSE. in CHEMRIV6, CHEMRIV7                   
!                 L_RAINCLOUD = .FALSE. in CHEMRIV6, CHEMRIV7                   
!                     L_SCAV6 = .FALSE. in WET                                  
!                                                                               
! -----------------------------                                                 
! --- CONTENT:                                                                  
! -----------------------------                                                 
! --- AER routines based on existing CALPUFF routines                           
!      subroutine chemriv6                                                      
!      subroutine chemriv7                                                      
!      subroutine chmriv6                                                       
!      subroutine chmriv7                                                       
!      subroutine setbckoc                                                      
!                                                                               
! --- AER ISORROPIA interface routine                                           
!     subroutine isodriver                                                      
!     (include 'isorropia.for')                                                 
!                                                                               
! --- AER CALTECH SOA routines                                                  
!     subroutine soadriver                                                      
!     subroutine caltech_soa                                                    
!     subroutine soasub                                                         
!                                                                               
! --- AQUEOUS-CHEMISTRY ROUTINES BASED ON RADM/CMAQ                             
!     subroutine aqradm                                                         
!     function hlconst                                                          
!     function index1                                                           
! -----------------------------                                                 
!                                                                               
!                                                                               
! --- UPDATE                                                                    
!                                                                               
! **********************************************************************        
! --- Exponent, Inc. Updates:                                                   
! **********************************************************************        
! --- V6.41-V6.42_x1.1 140521  : Call ISORROPIA with the METASTABLE             
!                                control mode ON                                
!                                Modified: ISODRIVER                            
!                      140521  : Use iterative procedure in calling             
!                                CALTECH_SOA to revise absorbing                
!                                organic mass based on current aerosol          
!                                mass as equilibrium is sought                  
!                                Modified: CHEMRIV7, SOADRIVER                  
!                      140521  : Use cloud fraction to set puff mass            
!                                altered when using AUX-file LWC                
!                                Use local temperature, pressure for AQ         
!                                Adjust mass transformed by fraction of         
!                                puff that overlays cloud water layers          
!                                Modified: CHEMRIV6, CHEMRIV7                   
!                      140521  : Enforce NO3<=TNO3 to avoid negative            
!                                HNO3 concentrations due to precision in        
!                                HNO3=TNO3-NO3 operation                        
!                                Also guard against NO3<0                       
!                                Modified: CHEMRIV6, CHEMRIV7                   
!                      140521  : Add minimum RH and SO4 for ISORROPIA           
!                                Modified: CHEMRIV6, CHEMRIV7                   
! **********************************************************************        
!                                                                               
! --- V6.4-V6.41   110301 (DGS): Treat aqueous-phase transformation             
!                                case of precip without cloud water             
!                                (from AUX file) by using a default             
!                                lwc=0.5g/m3 with at least 10% cloud            
!                                cover                                          
!                                Modified: CHEMRIV6, CHEMRIV7                   
!                                                                               
! --- V6.302-V6.4  101025 (DGS): Update ISORROPIA to V2.1                       
!                                Replaced: (include file isorropia.for)         
!                                          (include file isrpia.inc)            
!                                Modified: ISODRIVER                            
!                                                                               
! --- V5.8-V6.302  100917 (DGS): Restructure value**-n to value**(-n)           
!                                to satisfy LF95 compiler                       
!                                Modified: CHMRIV6, CHMRIV7                     
!                  100917 (DGS): Place NH3 profiles into /CHEMDAT/              
!                                and add chembkz to chembk logic                
!                  100917 (DGS): pass RSHET(fraction/hr) for setting            
!                                heterogeneous SO2 reaction rate                
!                                Modified: CHEMRIV6, CHEMRIV7                   
!                                          CHMRIV6, CHMRIV7                     
!                  100917 (DGS): Pass local cloud water mixing ratio            
!                                when available from CALMET 3D file             
!                                Modified: CHEMRIV6, CHEMRIV7                   
!                  100917 (DGS): Add local logical to convert LWC               
!                                passed to AQRADM from g/m3 to kg/m3            
!                                Modified: CHEMRIV6, CHEMRIV7                   
!                  100917 (DGS): Add local logical to use the minimum           
!                                TNO3 as a floor rather than as a               
!                                breakpoint for NO3=0.0                         
!                                Modified: CHEMRIV6, CHEMRIV7                   
!                  100917 (DGS): Add local logical to set a minimum             
!                                cloud cover (10%) whenever there is            
!                                liquid precipitation, and                      
!                                add but do not activate code to                
!                                change effective scavenging                    
!                                coefficients to limit their action             
!                                in any timestep to just the cloud              
!                                fraction affecting the puff mass.              
!                                Modified: CHEMRIV6, CHEMRIV7                   
!----------------------------------------------------------------------         
      subroutine chemriv6(delt,qin,coz,ctnh3,ch2o2,maqchem,rshet,temp,         &
     &                    rhum,rhoair,pivol,zlen,zpuf,cldamr,fzcld,            &
     &                    cldt,cldp,cloud,prate,zcoef,nspec,ldb1,io6,          &
     &                    rh_isrp,so4_isrp,                                    &
     &                    rate,scav)                                            
!----------------------------------------------------------------------         
!                                                                               
! --- CALPUFF    Version: TNG-7.1.0    Level: 140521          CHEMRIV6          
!                P. Karamchandani, AER (Adapted from chemriv)                   
!                                                                               
! --- PURPOSE:  This routine sets up call to RIVAD chemical                     
!               transformation subroutines, and passes results to               
!               calling program                                                 
!                                                                               
! --- UPDATES:                                                                  
!                                                                               
! --- V6.41-V6.42_x1.1 140521  : Use cloud fraction to set puff mass            
!                                altered when using AUX-file LWC                
!                                Use local temperature, pressure for AQ         
!                                Adjust mass transformed by fraction of         
!                                puff that overlays cloud water layers          
!                                Condition NO3 to be no larger than TNO3        
!                                and no smaller than zero                       
!                                Add minimum RH and SO4 for ISORROPIA           
!                                                                               
! --- V6.302-V6.41 110301 (DGS): Treat aqueous-phase transformation             
!                                case of precip without cloud water             
!                                (from AUX file) by using a default             
!                                lwc=0.5g/m3 with at least 10% cloud            
!                                cover                                          
! --- V5.8-V6.302  100917 (DGS): Pass RSHET(fraction/hr) for setting            
!                                heterogeneous SO2 reaction rate                
!                         (DGS): Add local cloud water mixing ratio             
!                                to input argument list (MLWC option)           
!                         (DGS): Add local logical to convert LWC               
!                                passed to AQRADM from g/m3 to kg/m3            
!                         (DGS): Add local logical to use the minimum           
!                                TNO3 as a floor rather than as a               
!                                breakpoint for NO3=0.0                         
!                         (DGS): Add local logical to set a minimum             
!                                cloud cover (10%) whenever there is            
!                                liquid precipitation for use in the            
!                                aqueous phase reaction/scavenging,             
!                                and add but do not activate code to            
!                                change the effective scavenging                
!                                coefficients by limiting their action          
!                                in any timestep to just the cloud              
!                                fraction affecting the puff mass.              
!                                                                               
! --- INPUTS:                                                                   
!         DELT - real    - integration time interval (hours)                    
!  QIN(mxspec) - real    - Pollutant mass (g) in the puff                       
!                            QIN(1) = SO2                                       
!                            QIN(2) = SO4                                       
!                            QIN(3) = NO                                        
!                            QIN(4) = NO2                                       
!                            QIN(5) = TNO3 (HNO3 + NO3)                         
!                            QIN(6) = NH4NO3                                    
!                            QIN(+) = Not Used Here                             
!                                                                               
!          COZ - real    - puff ozone concentration (ppb)                       
!        CTNH3 - real    - background ammonia concentration (ppb)               
!        CH2O2 - real    - puff H2O2 concentration (ppb)                        
!      MAQCHEM - integer - Aqueous phase transformation flag                    
!                            0 = aqueous phase transformation                   
!                                not modeled                                    
!                            1 = transformation rates adjusted                  
!                                for aqueous phase reactions                    
!        RSHET - real    - SO2 heterogeneous loss rate (fraction/hr)            
!         TEMP - real    - temperature (deg. K)                                 
!         RHUM - real    - relative humidity (percent)                          
!       RHOAIR - real    - surface air density (kg/m**3)                        
!        PIVOL - real    - Reciprocal of puff volume (1/m**3)                   
!         ZLEN - real    - Puff vertical length scale (m)                       
!         ZPUF - real    - Puff/Slug elevation (m MSL)                          
!       CLDAMR - real    - Average cloud water mixing ratio (g/kg)              
!                                                                               
! --- 6.42_x1.1                                                                 
!        FZCLD - real    - Fraction of puff layer in cloud                      
!         CLDT - real    - Average cloud temperature (K)                        
!         CLDP - real    - Average cloud pressure (atm)                         
!      RH_ISRP - real    - Minimum relative humidity (%)                        
!                          for ISORROPIA                                        
!     SO4_ISRP - real    - Minimum SO4 (g/m3) for ISORROPIA                     
!                                                                               
!        CLOUD - real    - Cloud cover (tenths)                                 
!        PRATE - real    - Precip. rate (mm/hr)                                 
!        ZCOEF - real    - Cosine of solar zenith angle                         
!        NSPEC - real    - number of species                                    
!         LDB1 - logical - Control variable for printing of debug               
!                          information                                          
!          IO6 - integer - Fortran unit number of printed output                
!                                                                               
!                                                                               
! --- OUTPUT:                                                                   
!  QIN(mxspec) - real    - Pollutant mass (g) in the puff                       
!                            QIN(1) = SO2                                       
!                            QIN(2) = SO4                                       
!                            QIN(3) = NO                                        
!                            QIN(4) = NO2                                       
!                            QIN(5) = TNO3 (HNO3 + NO3)                         
!                            QIN(6) = NH4NO3                                    
!                            QIN(+) = Not Used Here                             
!      RATE(4) - real    - Transformation rates (percent/hour)                  
!                            R(1) -- SO2 loss rate                              
!                            R(2) -- NOX loss rate                              
!                            R(3) -- TNO3 formation rate                        
!                            R(4) -- NO  loss rate                              
!     SCAV(6) - real    -  Scavenging coefficients (1/s)                        
!                            SCAV(1) -- SO2                                     
!                            SCAV(2) -- SO4                                     
!                            SCAV(3) -- NO                                      
!                            SCAV(4) -- NO2                                     
!                            SCAV(5) -- HNO3                                    
!                            SCAV(6) -- NO3                                     
!                                                                               
!                                                                               
! --- CHEMRIV6 called by: CHEM                                                  
! --- CHEMRIV6 calls:     PHOT, CHMRIV6, ISODRIVER, AQRADM                      
!----------------------------------------------------------------------         
!                                                                               
      implicit none                                                             
                                                                                
! --- Arguments 
! --- wangzhm save                                                                
      integer nspec
      integer maqchem, io6                                               
      real delt, coz, ctnh3, ch2o2, temp, rhum, rhoair, pivol, zlen             
      real zpuf, cloud, prate, zcoef, rshet                                     
      real qin(nspec), scav(nspec)                                              
      real rate(4)                                                              
      real cldamr                                                               
                                                                                
! --- 6.42_x1.1                                                                 
      real cldt,cldp,fzcld
      real,save :: fcloud2                                              
      real rh_isrp,so4_isrp        
! --- wangzhm save                                                   
      real,save :: cldta,cldpa                                                          
      real,save :: confcta,rhoaira                                                      
      real,save :: test                                                                 
                                                                                
      logical ldb1                                                              

! --- wangzhm save                                                                                
! --- Locals                                                                    
      real,save :: ppb(6),ppbi(6),q(6)                                                  
      real rmwt(6)                                                              
      real,save :: lwc    ! Liquid water content, g/m3                                  
      real,save :: cozm, ch2o2m, ctnh3m  ! Concentrations in mols/mols air units        
      real,save :: tno3floor                                                            
! --- 6.42_x1.1                                                                 
      real,save :: tso4min                                                              
! --- Local controls                                                            
      logical l_kgpcm, l_tno3floor, l_raincloud                                 
                                                                                
! --- Concentration array (moles/mole of air) for AQRADM                        
      real,save :: con(6)                                                               
                                                                                
! --- Note: TNO3 is weighted as NO3                                             
      data rmwt/64.,96.,30.,46.,62.,62./                                        
                                                                                
      real,save :: dt, vkgmol, confct, rk1, zcoefb, o3ppm                               
      real,save :: presur, f, ppbix,  ppbi4, ppbx                                       
      real,save :: delno, tso4, tno3, tnh3, rhfrac, pno3                                
      real,save :: patm, taucld, rhoairm                                                
      integer,save :: is, iss, j, ii, i                                                 
                                                                                
      real,save :: pcloud, fcloud                                                       
                                                                                
! ----------------------                                                        
! --- Set local controls                                                        
! ----------------------                                                        
! --- These enable 3 changes to the code to be reverted to their                
! --- original condition.                                                       
! --- 1.  Selecting l_kgpcm=.TRUE. converts the LWC passed to AQRADM            
! ---     from g/m3 to kg/m3, since AQRADM assumes it to be kg/m3.              
! ---     Selecting l_kgpcm=.FALSE. retains the original code and               
! ---     passes LWC to AQRADM in g/m3.                                         
      data l_kgpcm/.TRUE./                                                      
! --- 2.  Selecting l_tno3floor=.TRUE. uses the cut-off TNO3                    
! ---     concentration as a floor, and computes NO3 corresponding              
! ---     to this floor for all TNO3 less than this floor.  The                 
! ---     resulting ratio of NO3/TNO3(floor) is then multiplied by              
! ---     the actual TNO3 to obtain the final NO3 concentration.                
! ---     Selecting l_tno3floor=.FALSE. retains the original code               
! ---     and NO3=0.0 for all TNO3<TNO3(floor).                                 
      data l_tno3floor/.TRUE./                                                  
! --- 3.  Selecting l_raincloud=.TRUE. sets a minimum cloud cover               
! ---     of 10% whenever liquid precipitation is non-zero.  This               
! ---     only applies to the aqueous phase option and forces wet               
! ---     removal whenever there is rain.  When cloud cover is zero,            
! ---     aqueous conversion and wet removal are also zero.                     
! ---     (This also changes the effective scavenging coefficients              
! ---     by limiting their action in any timestep to just the                  
! ---     cloud fraction affecting the puff mass.)--NA                          
      data l_raincloud/.TRUE./                                                  
! ----------------------                                                        

! --- wangzhm omp: set threadsprivate                                           
      !$OMP THREADPRIVATE(cldta,cldpa,confcta,rhoaira,test,ppb,ppbi,q,&
      !$OMP& lwc,cozm,ch2o2m,ctnh3m,tno3floor,tso4min,&
      !$OMP& con,dt,vkgmol,confct,rk1,zcoefb,o3ppm,presur,f,ppbix,ppbi4,ppbx,delno,tso4,tno3,&
      !$OMP& tnh3,rhfrac,pno3,patm,taucld,rhoairm,is,iss,j,ii,i,pcloud,fcloud,fcloud2) 
                                                                                
! --- Transfer mass data to local array                                         
      iss = MIN(nspec,6)                                                        
      do is = 1, iss                                                            
         q(is) = qin(is)                                                        
      end do                                                                    
      do is = (iss+1), 6                                                        
         q(is) = 0.0                                                            
      end do                                                                    
                                                                                
! --- Compute initial concentrations as PPB                                     
! --- kg-Molar volume (m^3) at ambient T,P computed from air density            
      vkgmol = 28.97/rhoair                                                     
! --- Conversion factor to ppm in RIVAD                                         
      confct = (1.0e-3)*vkgmol*pivol                                            
      do is = 1, iss                                                            
         ppbi(is) = (1.0e09)*confct*q(is)/rmwt(is)                              
      end do                                                                    
                                                                                
! --- Get NO2 photolysis rate                                                   
      call PHOT(zpuf,cloud,zcoef,rk1,zcoefb)                                    
                                                                                
! --- Do conversions needed for RIVAD input arguments                           
! --- RIVAD weights TNO3 as HNO3: scale mass by 63/62=1.016129                  
      q(5) = q(5)*1.016129                                                      
! --- ozone in PPM rather than PPB                                              
      o3ppm = 1.0e-03*coz                                                       
! --- Implied pressure (mb) -- use T0=273. P0=1013. as in RIVAD                 
      patm = (22.4/vkgmol)*(temp/273.)                                          
      presur = 1013.*patm                                                       
                                                                                
      if(ldb1) then                                                             
         write(io6,*)                                                           
         write(io6,*)'CHEMRIV6:'                                                
         write(io6,*)'    rk1,rhoair,O3ppb = ',rk1,rhoair,coz                   
         write(io6,*)'     TNH3ppb,H2O2ppb = ',ctnh3, ch2o2                     
         write(io6,*)'    presur,temp,rhum = ',presur,temp,rhum                 
         write(io6,*)'  cloud,prate,cldamr = ',cloud,prate,cldamr               
         write(io6,*)'     zlen,zmsl,zcoef = ',zlen,zpuf,zcoef                  
         write(io6,*)' vkgmol,confct,pivol = ',vkgmol,confct,pivol              
                                                                                
! --- 6.42_x1.1                                                                 
         write(io6,*)'     fzcld,cldt,cldp = ',fzcld,cldt,cldp                  
         write(io6,*)'    rh_isrp,so4_isrp = ',rh_isrp,so4_isrp                 
                                                                                
         write(io6,*)                                                           
!         write(io6,101)                                                        
101      format('ppb =  SO2',7x,'SO4',7x,'NO',7x,'NO2',7x,'TNO3',              &
     &              6x,'NO3')                                                   
         write(io6,*)'Starting Concs -----'                                     
         write(io6,101)                                                         
         write(io6,'(3x,6e10.2)') (ppbi(j),j=1,6)                               
         write(io6,*)'Starting Puff Mass (g) -----'                             
         write(io6,'(3x,6f10.2)') (qin(j),j=1,6)                                
      endif                                                                     
                                                                                
! --- Call RIVAD module                                                         
      call CHMRIV6(q,temp,presur,rhum,o3ppm,rk1,zcoef,delt,confct,             &
     &             rshet)                                                       
! --- ozone in PPB                                                              
      coz = 1.0e03*o3ppm                                                        
                                                                                
! --- Re-weight TNO3 as NO3                                                     
      q(5) = q(5)*0.984127                                                      
                                                                                
      if(ldb1) then                                                             
         write(io6,*)'After CHMRIV6, Puff Mass (g) -----'                       
         write(io6,'(3x,6f10.2)') (q(j),j=1,6)                                  
      endif                                                                     
                                                                                
      if (maqchem.EQ.1) then                                                    
!                                                                               
! --- initialize scavenging coefficients                                        
         do ii = 1, 6                                                           
            scav(ii) = 0.                                                       
         end do                                                                 
                                                                                
! --- Use local cloud water mixing ratio if valid, or                           
! --- Assign liquid water content if there is cloud cover and the               
! --- temperature is above freezing                                             
         fcloud = cloud * 0.1    ! fractional cloud cover                       
         pcloud = 0.0                                                           
         if(l_raincloud .AND. temp.GT.273.15) then                              
! ---       Force at least 10% cloud cover when there is rain                   
            if(prate.GT.0.0) fcloud=MAX(fcloud,0.1)                             
         endif                                                                  
                                                                                
         if(fcloud.GT.0.0) pcloud = prate/fcloud                                
                                                                                
! --- 6.42_x1.1                                                                 
         rhoaira=rhoair                                                         
         confcta=confct                                                         
         cldta=temp                                                             
         cldpa=patm                                                             
! ---    Fraction of puff within cloud horizontally and vertically              
         fcloud2=fcloud                                                         
         if(cldamr.GE.0.0) then                                                 
! ---       Valid cloud water provided (use it)                                 
! ---       Use average cloud temperature and pressure                          
            cldta=cldt                                                          
            cldpa=cldp                                                          
            rhoaira=rhoair*(cldpa/patm)*(temp/cldta)                            
            confcta=(1.0e-3)*(28.97/rhoaira)*pivol                              
! ---       Convert: g/kg(air) * rhoair(kg/m3) ==> g/m3(air)                    
            lwc=cldamr*rhoaira                                                  
! ---       Set minimum cloud fraction to 10%                                   
            fcloud=MAX(fcloud,0.1)                                              
            pcloud=prate/fcloud                                                 
! ---       Fraction of puff within cloud horizontally and vertically           
            fcloud2=fcloud*fzcld                                                
                                                                                
         elseif (fcloud.gt.0. .and. temp.gt.273.15) then                        
            if (prate .gt. 0.) then                                             
               lwc = 0.5                                                        
            else                                                                
               lwc = 0.1                                                        
            end if                                                              
         else                                                                   
            lwc = 0.                                                            
         end if                                                                 
                                                                                
         if(l_kgpcm) then                                                       
! ---       Units for lwc in AQRADM are kg/m3                                   
            lwc=0.001*lwc                                                       
         endif                                                                  
                                                                                
         if (lwc .gt. 0.) then                                                  
!                                                                               
! --- Calculate concs in moles/mole air units                                   
            do ii = 1, 6                                                        
                                                                                
! --- 6.42_x1.1                                                                 
               con(ii) = MAX(confcta*q(ii)/rmwt(ii),0.)                         
                                                                                
            end do                                                              
!                                                                               
! --- Get HNO3 conc from total nitrate and PM nitrate                           
            con(5) = MAX(con(5) - con(6),0.)                                    
!                                                                               
            cozm = coz * 1.E-9                                                  
            ch2o2m = ch2o2 * 1.E-9                                              
            ctnh3m = ctnh3 * 1.E-9                                              
!                                                                               
! --- Call aqueous-phase chemistry module                                       
            taucld = delt * 3600.   ! timestep in seconds                       
                                                                                
! --- 6.42_x1.1                                                                 
            rhoairm = rhoaira * 1.E3/28.97   ! air density in moles/m3          
            if(ldb1) then                                                       
               write(io6,*)'AQRADM:cldta,cldpa,pcloud,lwc,cozm,'//             &
     &                     'ch2o2m,ctnh3m,con6= '                               
               write(io6,*)'called:   ',cldta,cldpa,pcloud,lwc,cozm,           &
     &                      ch2o2m,ctnh3m,con                                   
            endif                                                               
            call AQRADM(cldta,cldpa,taucld,pcloud,lwc,cozm,ch2o2m,             &
     &                  ctnh3m,con,rhoairm,zlen,scav)                           
            if(ldb1) then                                                       
               write(io6,*)'returned: ',cldta,cldpa,pcloud,lwc,cozm,           &
     &                      ch2o2m,ctnh3m,con                                   
               write(io6,*)'    scav: ',scav                                    
            endif                                                               
                                                                                
! --- Assign adjusted SO2 and SO4 concs; adjust for cloud cover                 
            q(1) = q(1)*(1. - fcloud2) + fcloud2*con(1)*rmwt(1)/confcta         
            q(2) = q(2)*(1. - fcloud2) + fcloud2*con(2)*rmwt(2)/confcta         
!                                                                               
! --- Adjusted oxidant concs                                                    
            coz = coz*(1. - fcloud2) + fcloud2*cozm*1.E9                        
            ch2o2 = ch2o2*(1. - fcloud2) + fcloud2*ch2o2m*1.E9                  
                                                                                
!                                                                               
! --- Adjusted scavenging rates                                                 
            if (prate .gt. 0.) then                                             
! ---          Possible alternate method (NA)                                   
! ---           if(l_raincloud .AND. fcloud.LT.0.99) then                       
! ---              do ii = 1, 6                                                 
! ---                 scav(ii)=-ALOG(1.0 - fcloud*                              
! --- &                        (1.0 - EXP(-scav(ii)*taucld) ))/taucld           
! ---              end do                                                       
! ---           else                                                            
                  do ii = 1, 6                                                  
                     scav(ii) = scav(ii) * fcloud                               
                  end do                                                        
! ---           endif                                                           
            end if                                                              
         end if                                                                 
      end if                                                                    
                                                                                
! --- Condition mass results                                                    
      do is = 1,iss                                                             
         q(is) = MAX(q(is),0.0)                                                 
      end do                                                                    
                                                                                
! --- Compute ending concentrations as PPB                                      
      do is = 1,iss                                                             
         ppb(is) = (1.0e09)*confct*q(is)/rmwt(is)                               
      end do                                                                    
                                                                                
! --- Compute conversion rates (%/hr) for QA review                             
      do i = 1,4                                                                
         rate(i) = 0.0                                                          
      end do                                                                    
                                                                                
! --- Compute these logs only if debug output is ON                             
      if(ldb1) then                                                             
         dt = .01*delt                                                          
! ---    SOX:                                                                   
         if(ppbi(1).GT.0.0 .AND. ppb(1).GT.0.0) rate(1) = -ALOG(ppb(1)/        &
     &                                                   ppbi(1))/dt            
! ---    NOX:                                                                   
         ppbix = ppbi(4) + ppbi(3)                                              
         ppbx = ppb(4) + ppb(3)                                                 
         if(ppbix.GT.0.0 .AND. ppbx.GT.0.0) rate(2) = -ALOG(ppbx/              &
     &                                                   ppbix)/dt              
! ---    TNO3:  Conversion of NO2 after NO:NO2 process                          
         ppbi4 = ppbi(4) + ppbi(3) - ppb(3)                                     
         if(ppbi4.GT.0.0 .AND. ppb(4).GT.0.0) rate(3) = -ALOG(ppb(4)/          &
     &                                                   ppbi4)/dt              
! ---    NO:                                                                    
         delno = ppb(3) - ppbi(3)                                               
         if(delno.LT.0.0) then                                                  
            if(ppbi(3).GT.0.0 .AND. ppb(3).GT.0.0) rate(4) =                   &
     &                          ALOG(ppb(3)/ppbi(3))/dt                         
         elseif(delno.GT.0.0) then                                              
            if(ppbi4.GT.0.0 .AND. ppb(4).GT.0.0) rate(4) = -ALOG(ppbi4/        &
     &                                                    ppbi(4))/dt           
         endif                                                                  
      endif                                                                     
                                                                                
! --- Inorganic aerosol equilibrium with ISORROPIA                              
      if(q(5).gt.0.0)then                                                       
!                                                                               
! --- concs in mols/m3                                                          
                                                                                
! --- 6.42_x1.1                                                                 
! ---   Apply so4(g/m3) concentration constraint from control file              
        tso4=q(2)*pivol                                                         
        tso4=MAX(tso4,so4_isrp)/rmwt(2)                                         
                                                                                
        tso4 = MAX(tso4,1.E-12)                                                 
        tno3 = q(5) * pivol / rmwt(6)                                           
        tnh3 = ctnh3 * rhoair / 28.97E6                                         
                                                                                
! --- 6.42_x1.1                                                                 
! ---   Apply RH(%) constraint from control file to compute RH fraction         
        rhfrac=0.01*MAX(rhum,rh_isrp)                                           
                                                                                
        if(tno3.gt.1.E-10)then                                                  
          if(ldb1) then                                                         
                                                                                
! --- 6.42_x1.1                                                                 
             write(io6,*)'ISODRIVER:tso4,tno3,tnh3,rhfrac,temp     = '          
             write(io6,*)'called:   ',tso4,tno3,tnh3,rhfrac,temp                
                                                                                
          endif                                                                 
          call isodriver(tso4,tno3,tnh3,rhfrac,temp,pno3)                       
          if(ldb1) then                                                         
             write(io6,*)'returned: ',tso4,tno3,tnh3,rhfrac,temp,pno3           
          endif                                                                 
                                                                                
        elseif(l_tno3floor) then                                                
          tno3floor=1.E-10                                                      
          if(ldb1) then                                                         
             write(io6,*)                                                      &
     &             'ISODRIVER:tso4,tno3floor,tnh3,rhfrac,temp,pno3= '           
                                                                                
! --- 6.42_x1.1                                                                 
             write(io6,*)'called:   ',                                         &
     &                        tso4,tno3floor,tnh3,rhfrac,temp                   
                                                                                
          endif                                                                 
          call isodriver(tso4,tno3floor,tnh3,rhfrac,temp,pno3)                  
          if(ldb1) then                                                         
             write(io6,*)'returned: ',                                         &
     &                        tso4,tno3floor,tnh3,rhfrac,temp,pno3              
          endif                                                                 
          pno3=pno3*(tno3/tno3floor)                                            
                                                                                
        else                                                                    
          pno3 = 0.                                                             
        endif                                                                   
        q(6) = pno3 * rmwt(6) / pivol                                           
                                                                                
! --- 6.42_x1.1                                                                 
! ---   Condition NO3 (q(6)) to be no smaller than zero                         
        if(q(6).LT.0.0) then                                                    
           write(io6,*)'CHEMRIV6 Warning: reset NO3(g) from ',                 &
     &                 q(6),' to ZERO'                                          
           q(6)=0.0                                                             
        endif                                                                   
! ---   Condition NO3 to be no larger than TNO3 (q(5))                          
! ---   (q(5)>0.0 in this block)                                                
        test=q(6)/q(5)-1.0                                                      
        if(test.GT.1.0e-06) then                                                
           write(io6,*)'CHEMRIV6 Warning: reset NO3(g) from ',                 &
     &                 q(6),' to ',q(5)                                         
        endif                                                                   
        q(6)=MIN(q(5),q(6))                                                     
                                                                                
      endif                                                                     
                                                                                
! --- Transfer revised mass data to original array                              
      do is = 1,iss                                                             
         qin(is) = q(is)                                                        
      end do                                                                    
                                                                                
      if(ldb1) then                                                             
         write(io6,*)'After ISODRIVER, Puff Mass (g) -----'                     
         write(io6,'(3x,6f10.2)') (qin(j),j=1,6)                                
         write(io6,*)                                                           
      endif                                                                     
                                                                                
      return                                                                    
      end                                                                       
                                                                                
!----------------------------------------------------------------------         
      subroutine chemriv7(delt,qin,coz,ctnh3,ch2o2,maqchem,rshet,temp,         &
     &                    rhum,rhoair,pivol,zlen,zpuf,cldamr,fzcld,            &
     &                    cldt,cldp,cloud,prate,zcoef,nspec,ldb1,io6,          &
     &                    rh_isrp,so4_isrp,                                    &
     &                    rate,scav)                                            
!----------------------------------------------------------------------         
!                                                                               
! --- CALPUFF    Version: TNG-7.1.0    Level: 140521          CHEMRIV7          
!                P. Karamchandani, AER (Adapted from chemriv)                   
!                                                                               
! --- PURPOSE:  This routine sets up call to RIVAD chemical                     
!               transformation subroutines, and passes results to               
!               calling program                                                 
!                                                                               
! --- UPDATES:                                                                  
!                                                                               
! --- V6.41-V6.42_x1.1 140521  : Use cloud fraction to set puff mass            
!                                altered when using AUX-file LWC                
!                                Use local temperature, pressure for AQ         
!                                Adjust mass transformed by fraction of         
!                                puff that overlays cloud water layers          
!                                Condition NO3 to be no larger than TNO3        
!                                and no smaller than zero                       
!                                Add list file unit to SOADRIVER                
!                                Add minimum RH and SO4 for ISORROPIA           
!                                                                               
! --- V6.302-V6.41 110301 (DGS): Treat aqueous-phase transformation             
!                                case of precip without cloud water             
!                                (from AUX file) by using a default             
!                                lwc=0.5g/m3 with at least 10% cloud            
!                                cover                                          
! --- V5.8-V6.302  100917 (DGS): Pass RSHET(fraction/hr) for setting            
!                                heterogeneous SO2 reaction rate                
!                         (DGS): Add local cloud water mixing ratio             
!                                to input argument list (MLWC option)           
!                         (DGS): Add local logical to convert LWC               
!                                passed to AQRADM from g/m3 to kg/m3            
!                         (DGS): Add local logical to use the minimum           
!                                TNO3 as a floor rather than as a               
!                                breakpoint for NO3=0.0                         
!                         (DGS): Add local logical to set a minimum             
!                                cloud cover (10%) whenever there is            
!                                liquid precipitation for use in the            
!                                aqueous phase reaction/scavenging,             
!                                and add but do not activate code to            
!                                change the effective scavenging                
!                                coefficients by limiting their action          
!                                in any timestep to just the cloud              
!                                fraction affecting the puff mass.              
!                                                                               
! --- INPUTS:                                                                   
!         DELT - real    - integration time interval (hours)                    
!  QIN(mxspec) - real    - Pollutant mass (g) in the puff                       
!                            QIN(1) = SO2                                       
!                            QIN(2) = SO4                                       
!                            QIN(3) = NO                                        
!                            QIN(4) = NO2                                       
!                            QIN(5) = TNO3 (HNO3 + NO3)                         
!                            QIN(6) = NH4NO3                                    
!                            QIN(7)  = Primary OC (POC)                         
!                            QIN(8)  = TOL (Toluene)                            
!                            QIN(9)  = TOLAER1 (Condensable product)            
!                            QIN(10) = TOLAER2 (Condensable product)            
!                            QIN(11) = ATOLA1 (SOA 1 from TOL)                  
!                            QIN(12) = ATOLA2 (SOA 2 from TOL)                  
!                            QIN(13) = XYL (Xylene)                             
!                            QIN(14) = XYLAER1 (Condensable product)            
!                            QIN(15) = XYLAER2 (Condensable product)            
!                            QIN(16) = AXYLA1 (SOA 1 from XYL)                  
!                            QIN(17) = AXYLA2 (SOA 2 from XYL)                  
!                            QIN(18) = ALKH (Higher alkanes)                    
!                            QIN(19) = ALKHAER (Condensable product)            
!                            QIN(20) = ALKHA (SOA 1 from ALKH)                  
!                            QIN(21) = PAH                                      
!                            QIN(22) = PAHAER1 (Condensable product)            
!                            QIN(23) = PAHAER2 (Condensable product)            
!                            QIN(24) = APAHA1 (SOA 1 from PAH)                  
!                            QIN(25) = APAHA2 (SOA 2 from PAH)                  
!                                                                               
!          COZ - real    - puff ozone concentration (ppb)                       
!        CTNH3 - real    - background ammonia concentration (ppb)               
!        CH2O2 - real    - puff H2O2 concentration (ppb)                        
!      MAQCHEM - integer - Aqueous phase transformation flag                    
!                            0 = aqueous phase transformation                   
!                                not modeled                                    
!                            1 = transformation rates adjusted                  
!                                for aqueous phase reactions                    
!        RSHET - real    - SO2 heterogeneous loss rate (fraction/hr)            
!         TEMP - real    - temperature (deg. K)                                 
!         RHUM - real    - relative humidity (percent)                          
!       RHOAIR - real    - surface air density (kg/m**3)                        
!        PIVOL - real    - Reciprocal of puff volume (1/m**3)                   
!         ZLEN - real    - Puff vertical length scale (m)                       
!         ZPUF - real    - Puff/Slug elevation (m MSL)                          
!       CLDAMR - real    - Average cloud water mixing ratio (g/kg)              
!                                                                               
! --- 6.42_x1.1                                                                 
!        FZCLD - real    - Fraction of puff layer in cloud                      
!         CLDT - real    - Average cloud temperature (K)                        
!         CLDP - real    - Average cloud pressure (atm)                         
!      RH_ISRP - real    - Minimum relative humidity (%)                        
!                          for ISORROPIA                                        
!     SO4_ISRP - real    - Minimum SO4 (g/m3) for ISORROPIA                     
!                                                                               
!        CLOUD - real    - Cloud cover (tenths)                                 
!        PRATE - real    - Precip. rate (mm/hr)                                 
!        ZCOEF - real    - Cosine of solar zenith angle                         
!        NSPEC - real    - number of species                                    
!         LDB1 - logical - Control variable for printing of debug               
!                          information                                          
!          IO6 - integer - Fortran unit number of printed output                
!                                                                               
!                                                                               
! --- OUTPUT:                                                                   
!  QIN(mxspec) - real    - Pollutant mass (g) in the puff                       
!                            QIN(1)  = SO2                                      
!                            QIN(2)  = SO4                                      
!                            QIN(3)  = NO                                       
!                            QIN(4)  = NO2                                      
!                            QIN(5)  = TNO3 (HNO3 + NO3)                        
!                            QIN(6)  = NH4NO3                                   
!                            QIN(7)  = Primary OC (POC)                         
!                            QIN(8)  = TOL (Toluene)                            
!                            QIN(9)  = TOLAER1 (Condensable product)            
!                            QIN(10) = TOLAER2 (Condensable product)            
!                            QIN(11) = ATOLA1 (SOA 1 from TOL)                  
!                            QIN(12) = ATOLA2 (SOA 2 from TOL)                  
!                            QIN(13) = XYL (Xylene)                             
!                            QIN(14) = XYLAER1 (Condensable product)            
!                            QIN(15) = XYLAER2 (Condensable product)            
!                            QIN(16) = AXYLA1 (SOA 1 from XYL)                  
!                            QIN(17) = AXYLA2 (SOA 2 from XYL)                  
!                            QIN(18) = ALKH (Higher alkanes)                    
!                            QIN(19) = ALKHAER (Condensable product)            
!                            QIN(20) = ALKHA (SOA 1 from ALKH)                  
!                            QIN(21) = PAH                                      
!                            QIN(22) = PAHAER1 (Condensable product)            
!                            QIN(23) = PAHAER2 (Condensable product)            
!                            QIN(24) = APAHA1 (SOA 1 from PAH)                  
!                            QIN(25) = APAHA2 (SOA 2 from PAH)                  
!      RATE(8) - real    - Transformation rates (percent/hour)                  
!                            R(1) -- SO2 loss rate                              
!                            R(2) -- NOX loss rate                              
!                            R(3) -- TNO3 formation rate                        
!                            R(4) -- NO  loss rate                              
!                            R(5) -- TOL loss rate                              
!                            R(6) -- XYL loss rate                              
!                            R(7) -- ALKH loss rate                             
!                            R(8) -- PAH loss rate                              
!     SCAV(6) - real    -  Scavenging coefficients (1/s)                        
!                            SCAV(1) -- SO2                                     
!                            SCAV(2) -- SO4                                     
!                            SCAV(3) -- NO                                      
!                            SCAV(4) -- NO2                                     
!                            SCAV(5) -- HNO3                                    
!                            SCAV(6) -- NO3                                     
!                                                                               
! --- CHEMRIV7 called by: CHEM                                                  
! --- CHEMRIV7 calls:     PHOT, CHMRIV7, ISODRIVER, SOADRIVER, AQRADM           
!----------------------------------------------------------------------         
!                                                                               
      implicit none                                                             

! --- wangzhm save
                                                                                
! --- Arguments                                                                 
      integer nspec
      integer maqchem, io6                                               
      real delt, coz, ctnh3, ch2o2, temp, rhum, rhoair, pivol, zlen             
      real zpuf, cloud, prate, zcoef, rshet                                     
      real qin(nspec), scav(nspec)                                              
      real rate(8)                                                              
      real cldamr                                                               
                                                                                
! --- 6.42_x1.1                                                                 
      real cldt,cldp,fzcld
      real,save :: fcloud2                                              
      real rh_isrp,so4_isrp 
                                                          
      real,save :: cldta,cldpa                                                          
      real,save :: confcta,rhoaira                                                      
      real,save :: test                                                                 
                                                                                
      logical ldb1                                                              
                                                                                
! --- Locals                                                                    
      real,save :: ppb(25),ppbi(25),q(25)                                               
      real rmwt(25)                                                             
      real,save :: lwc    ! Liquid water content, g/m3                                  
      real,save :: cozm, ch2o2m, ctnh3m  ! Concentrations in mols/mols air units        
      real,save :: tno3floor                                                            
! --- 6.42_x1.1                                                                 
      real,save :: tso4min                                                              
! --- Local controls                                                            
      logical l_kgpcm, l_tno3floor, l_raincloud                                 
                                                                                
! --- Concentration array (moles/mole of air) for AQRADM                        
      real,save :: con(6)                                                               
                                                                                
! --- Note: TNO3 is weighted as NO3                                             
      data rmwt/64.,96.,30.,46.,62.,62.,180.,                                  &
     &          5*92.,5*106.,3*226.,5*156./                                     
                                                                                
      real,save :: dt, vkgmol, confct, rk1, zcoefb, o3ppm                               
      real,save :: presur, f, ppbix,  ppbi4, ppbx                                       
      real,save :: delno, tso4, tno3, tnh3, rhfrac, pno3                                
      real,save :: patm, taucld, rhoairm                                                
      integer,save :: is, iss, j, ii, i                                                 
!                                                                               
      real,save :: pcloud, fcloud                                                       
                                                                                
! ----------------------                                                        
! --- Set local controls                                                        
! ----------------------                                                        
! --- These enable 2 changes to the code to be reverted to their                
! --- original condition.                                                       
! --- 1.  Selecting l_kgpcm=.TRUE. converts the LWC passed to AQRADM            
! ---     from g/m3 to kg/m3, since AQRADM assumes it to be kg/m3.              
! ---     Selecting l_kgpcm=.FALSE. retains the original code and               
! ---     passes LWC to AQRADM in g/m3.                                         
      data l_kgpcm/.TRUE./                                                      
! --- 2.  Selecting l_tno3floor=.TRUE. uses the cut-off TNO3                    
! ---     concentration as a floor, and computes NO3 corresponding              
! ---     to this floor for all TNO3 less than this floor.  The                 
! ---     resulting ratio of NO3/TNO3(floor) is then multiplied by              
! ---     the actual TNO3 to obtain the final NO3 concentration.                
! ---     Selecting l_tno3floor=.FALSE. retains the original code               
! ---     and NO3=0.0 for all TNO3<TNO3(floor).                                 
      data l_tno3floor/.TRUE./                                                  
! --- 3.  Selecting l_raincloud=.TRUE. sets a minimum cloud cover               
! ---     of 10% whenever liquid precipitation is non-zero.  This               
! ---     only applies to the aqueous phase option and forces wet               
! ---     removal whenever there is rain.  When cloud cover is zero,            
! ---     aqueous conversion and wet removal are also zero.                     
! ---     (This also changes the effective scavenging coefficients              
! ---     by limiting their action in any timestep to just the                  
! ---     cloud fraction affecting the puff mass.)--NA                          
      data l_raincloud/.TRUE./                                                  
! ----------------------                                                        
! --- wangzhm omp: set threadsprivate                                           
      !$OMP THREADPRIVATE(cldta,cldpa,confcta,rhoaira,test,ppb,ppbi,q,&
      !$OMP& lwc,cozm,ch2o2m,ctnh3m,tno3floor,tso4min,&
      !$OMP& con,dt,vkgmol,confct,rk1,zcoefb,o3ppm,presur,f,ppbix,ppbi4,ppbx,delno,&
      !$OMP& tso4,tno3,tnh3,rhfrac,pno3,patm,taucld,rhoairm,is,iss,j,ii,i,pcloud,fcloud,fcloud2)                  
                                                                                       
! --- Transfer mass data to local array                                         
      iss = MIN(nspec,25)                                                       
      do is = 1, iss                                                            
         q(is) = qin(is)                                                        
      end do                                                                    
      do is = (iss+1), 25                                                       
         q(is) = 0.0                                                            
      end do                                                                    
                                                                                
! --- Compute initial concentrations as PPB                                     
! --- kg-Molar volume (m^3) at ambient T,P computed from air density            
      vkgmol = 28.97/rhoair                                                     
! --- Conversion factor to ppm in RIVAD                                         
      confct = (1.0e-3)*vkgmol*pivol                                            
      do is = 1, iss                                                            
         ppbi(is) = (1.0e09)*confct*q(is)/rmwt(is)                              
      end do                                                                    
                                                                                
! --- Get NO2 photolysis rate                                                   
      call PHOT(zpuf,cloud,zcoef,rk1,zcoefb)                                    
                                                                                
! --- Do conversions needed for RIVAD input arguments                           
! --- RIVAD weights TNO3 as HNO3: scale mass by 63/62=1.016129                  
      q(5) = q(5)*1.016129                                                      
! --- ozone in PPM rather than PPB                                              
      o3ppm = 1.0e-03*coz                                                       
! --- Implied pressure (mb) -- use T0=273. P0=1013. as in RIVAD                 
      patm = (22.4/vkgmol)*(temp/273.)                                          
      presur = 1013.*patm                                                       
                                                                                
      if(ldb1) then                                                             
         write(io6,*)'CHEMRIV7:'                                                
         write(io6,*)'    rk1,rhoair,O3ppb = ',rk1,rhoair,coz                   
         write(io6,*)'     TNH3ppb,H2O2ppb = ',ctnh3, ch2o2                     
         write(io6,*)'    presur,temp,rhum = ',presur,temp,rhum                 
         write(io6,*)'  cloud,prate,cldamr = ',cloud,prate,cldamr               
         write(io6,*)'     zlen,zmsl,zcoef = ',zlen,zpuf,zcoef                  
         write(io6,*)' vkgmol,confct,pivol = ',vkgmol,confct,pivol              
                                                                                
! --- 6.42_x1.1                                                                 
         write(io6,*)'     fzcld,cldt,cldp = ',fzcld,cldt,cldp                  
         write(io6,*)'    rh_isrp,so4_isrp = ',rh_isrp,so4_isrp                 
                                                                                
         write(io6,*)                                                           
!         write(io6,101)                                                        
101      format('ppb =  SO2',7x,'SO4',7x,'NO',7x,'NO2',7x,'TNO3',              &
     &              6x,'NO3')                                                   
         write(io6,*)'Starting Concs -----'                                     
         write(io6,101)                                                         
         write(io6,'(3x,6e10.2)') (ppbi(j),j=1,6)                               
      endif                                                                     
                                                                                
! --- Call RIVAD module                                                         
      call CHMRIV7(q,temp,presur,rhum,o3ppm,rk1,zcoef,delt,confct,             &
     &             rshet)                                                       
! --- ozone in PPB                                                              
      coz = 1.0e03*o3ppm                                                        
                                                                                
! --- Re-weight TNO3 as NO3                                                     
      q(5) = q(5)*0.984127                                                      
                                                                                
      if(ldb1) then                                                             
         write(io6,*)'After CHMRIV6, Puff Mass (g) -----'                       
         write(io6,'(3x,6f10.2)') (q(j),j=1,6)                                  
      endif                                                                     
                                                                                
      if (maqchem.EQ.1) then                                                    
!                                                                               
! --- initialize scavenging coefficients                                        
         do ii = 1, 6                                                           
            scav(ii) = 0.                                                       
         end do                                                                 
                                                                                
! --- Use local cloud water mixing ratio if valid, or                           
! --- Assign liquid water content if there is cloud cover and the               
! --- temperature is above freezing                                             
         fcloud = cloud * 0.1    ! fractional cloud cover                       
         pcloud = 0.0                                                           
         if(l_raincloud .AND. temp.GT.273.15) then                              
! ---       Force at least 10% cloud cover when there is rain                   
            if(prate.GT.0.0) fcloud=MAX(fcloud,0.1)                             
         endif                                                                  
                                                                                
         if(fcloud.GT.0.0) pcloud = prate/fcloud                                
                                                                                
! --- 6.42_x1.1                                                                 
         rhoaira=rhoair                                                         
         confcta=confct                                                         
         cldta=temp                                                             
         cldpa=patm                                                             
! ---    Fraction of puff within cloud horizontally and vertically              
         fcloud2=fcloud                                                         
         if(cldamr.GE.0.0) then                                                 
! ---       Valid cloud water provided (use it)                                 
! ---       Use average cloud temperature and pressure                          
            cldta=cldt                                                          
            cldpa=cldp                                                          
            rhoaira=rhoair*(cldpa/patm)*(temp/cldta)                            
            confcta=(1.0e-3)*(28.97/rhoaira)*pivol                              
! ---       Convert: g/kg(air) * rhoair(kg/m3) ==> g/m3(air)                    
            lwc=cldamr*rhoaira                                                  
! ---       Set minimum cloud fraction to 10%                                   
            fcloud=MAX(fcloud,0.1)                                              
            pcloud=prate/fcloud                                                 
! ---       Fraction of puff within cloud horizontally and vertically           
            fcloud2=fcloud*fzcld                                                
                                                                                
         elseif (fcloud.gt.0. .and. temp.gt.273.15) then                        
            if (prate .gt. 0.) then                                             
               lwc = 0.5                                                        
            else                                                                
               lwc = 0.1                                                        
            end if                                                              
         else                                                                   
            lwc = 0.                                                            
         end if                                                                 
                                                                                
         if(l_kgpcm) then                                                       
! ---       Units for lwc in AQRADM are kg/m3                                   
            lwc=0.001*lwc                                                       
         endif                                                                  
                                                                                
         if (lwc .gt. 0.) then                                                  
!                                                                               
! --- Calculate concs in moles/mole air units                                   
            do ii = 1, 6                                                        
                                                                                
! --- 6.42_x1.1                                                                 
               con(ii) = MAX(confcta*q(ii)/rmwt(ii),0.)                         
                                                                                
            end do                                                              
!                                                                               
! --- Get HNO3 conc from total nitrate and PM nitrate                           
            con(5) = MAX(con(5) - con(6),0.)                                    
!                                                                               
            cozm = coz * 1.E-9                                                  
            ch2o2m = ch2o2 * 1.E-9                                              
            ctnh3m = ctnh3 * 1.E-9                                              
!                                                                               
! --- Call aqueous-phase chemistry module                                       
            taucld = delt * 3600.   ! timestep in seconds                       
                                                                                
! --- 6.42_x1.1                                                                 
            rhoairm = rhoaira * 1.E3/28.97   ! air density in moles/m3          
            if(ldb1) then                                                       
               write(io6,*)'AQRADM:cldta,cldpa,pcloud,lwc,cozm,'//             &
     &                     'ch2o2m,ctnh3m,con6= '                               
               write(io6,*)'called:   ',cldta,cldpa,pcloud,lwc,cozm,           &
     &                      ch2o2m,ctnh3m,con                                   
            endif                                                               
            call AQRADM(cldta,cldpa,taucld,pcloud,lwc,cozm,ch2o2m,             &
     &                  ctnh3m,con,rhoairm,zlen,scav)                           
            if(ldb1) then                                                       
               write(io6,*)'returned: ',cldta,cldpa,pcloud,lwc,cozm,           &
     &                      ch2o2m,ctnh3m,con                                   
               write(io6,*)'    scav: ',scav                                    
            endif                                                               
                                                                                
! --- Assign adjusted SO2 and SO4 concs; adjust for cloud cover                 
            q(1) = q(1)*(1. - fcloud2) + fcloud2*con(1)*rmwt(1)/confcta         
            q(2) = q(2)*(1. - fcloud2) + fcloud2*con(2)*rmwt(2)/confcta         
!                                                                               
! --- Adjusted oxidant concs                                                    
            coz = coz*(1. - fcloud2) + fcloud2*cozm*1.E9                        
            ch2o2 = ch2o2*(1. - fcloud2) + fcloud2*ch2o2m*1.E9                  
                                                                                
!                                                                               
! --- Adjusted scavenging rates                                                 
            if (prate .gt. 0.) then                                             
! ---          Possible alternate method (NA)                                   
! ---           if(l_raincloud .AND. fcloud.LT.0.99) then                       
! ---              do ii = 1, 6                                                 
! ---                 scav(ii)=-ALOG(1.0 - fcloud*                              
! --- &                        (1.0 - EXP(-scav(ii)*taucld) ))/taucld           
! ---              end do                                                       
! ---           else                                                            
                  do ii = 1, 6                                                  
                     scav(ii) = scav(ii) * fcloud                               
                  end do                                                        
! ---           endif                                                           
            end if                                                              
         end if                                                                 
      end if                                                                    
                                                                                
! --- Condition mass results                                                    
      do is = 1,iss                                                             
         q(is) = MAX(q(is),0.0)                                                 
      end do                                                                    
                                                                                
! --- Compute ending concentrations as PPB                                      
      do is = 1,iss                                                             
         ppb(is) = (1.0e09)*confct*q(is)/rmwt(is)                               
      end do                                                                    
                                                                                
! --- Compute conversion rates (%/hr) for QA review                             
      do i = 1,8                                                                
         rate(i) = 0.0                                                          
      end do                                                                    
                                                                                
! --- Compute these logs only if debug output is ON                             
      if(ldb1) then                                                             
         dt = .01*delt                                                          
! ---    SOX:                                                                   
         if(ppbi(1).GT.0.0 .AND. ppb(1).GT.0.0) rate(1) = -ALOG(ppb(1)/        &
     &                                                   ppbi(1))/dt            
! ---    NOX:                                                                   
         ppbix = ppbi(4) + ppbi(3)                                              
         ppbx = ppb(4) + ppb(3)                                                 
         if(ppbix.GT.0.0 .AND. ppbx.GT.0.0) rate(2) = -ALOG(ppbx/              &
     &                                                   ppbix)/dt              
! ---    TNO3:  Conversion of NO2 after NO:NO2 process                          
         ppbi4 = ppbi(4) + ppbi(3) - ppb(3)                                     
         if(ppbi4.GT.0.0 .AND. ppb(4).GT.0.0) rate(3) = -ALOG(ppb(4)/          &
     &                                                   ppbi4)/dt              
! ---    NO:                                                                    
         delno = ppb(3) - ppbi(3)                                               
         if(delno.LT.0.0) then                                                  
            if(ppbi(3).GT.0.0 .AND. ppb(3).GT.0.0) rate(4) =                   &
     &                          ALOG(ppb(3)/ppbi(3))/dt                         
         elseif(delno.GT.0.0) then                                              
            if(ppbi4.GT.0.0 .AND. ppb(4).GT.0.0) rate(4) = -ALOG(ppbi4/        &
     &                                                    ppbi(4))/dt           
! ---    TOL:                                                                   
         if(ppbi(8).GT.0.0 .AND. ppb(8).GT.0.0) rate(5) = -ALOG(ppb(8)/        &
     &                                                    ppbi(8))/dt           
! ---    XYL:                                                                   
         if(ppbi(13).GT.0.0 .AND. ppb(13).GT.0.0) rate(6) =                    &
     &                                  -ALOG(ppb(13)/ppbi(13))/dt              
! ---    ALKH:                                                                  
         if(ppbi(18).GT.0.0 .AND. ppb(18).GT.0.0) rate(7) =                    &
     &                                  -ALOG(ppb(18)/ppbi(18))/dt              
! ---    PAH:                                                                   
         if(ppbi(21).GT.0.0 .AND. ppb(21).GT.0.0) rate(8) =                    &
     &                                  -ALOG(ppb(21)/ppbi(21))/dt              
         endif                                                                  
      endif                                                                     
                                                                                
! --- Inorganic aerosol equilibrium with ISORROPIA                              
      if(q(5).gt.0.0)then                                                       
!                                                                               
! --- concs in mols/m3                                                          
                                                                                
! --- 6.42_x1.1                                                                 
! ---   Apply so4(g/m3) concentration constraint from control file              
        tso4=q(2)*pivol                                                         
        tso4=MAX(tso4,so4_isrp)/rmwt(2)                                         
                                                                                
        tso4 = MAX(tso4,1.E-12)                                                 
        tno3 = q(5) * pivol / rmwt(6)                                           
        tnh3 = ctnh3 * rhoair / 28.97E6                                         
                                                                                
! --- 6.42_x1.1                                                                 
! ---   Apply RH(%) constraint from control file to compute RH fraction         
        rhfrac=0.01*MAX(rhum,rh_isrp)                                           
                                                                                
        if(tno3.gt.1.E-10)then                                                  
          if(ldb1) then                                                         
                                                                                
! --- 6.42_x1.1                                                                 
             write(io6,*)'ISODRIVER:tso4,tno3,tnh3,rhfrac,temp     = '          
             write(io6,*)'called:   ',tso4,tno3,tnh3,rhfrac,temp                
                                                                                
          endif                                                                 
          call isodriver(tso4,tno3,tnh3,rhfrac,temp,pno3)                       
          if(ldb1) then                                                         
             write(io6,*)'returned: ',tso4,tno3,tnh3,rhfrac,temp,pno3           
          endif                                                                 
                                                                                
        elseif(l_tno3floor) then                                                
          tno3floor=1.E-10                                                      
          if(ldb1) then                                                         
             write(io6,*)                                                      &
     &             'ISODRIVER:tso4,tno3floor,tnh3,rhfrac,temp,pno3= '           
                                                                                
! --- 6.42_x1.1                                                                 
             write(io6,*)'called:   ',                                         &
     &                        tso4,tno3floor,tnh3,rhfrac,temp                   
                                                                                
          endif                                                                 
          call isodriver(tso4,tno3floor,tnh3,rhfrac,temp,pno3)                  
          if(ldb1) then                                                         
             write(io6,*)'returned: ',                                         &
     &                        tso4,tno3floor,tnh3,rhfrac,temp,pno3              
          endif                                                                 
          pno3=pno3*(tno3/tno3floor)                                            
                                                                                
        else                                                                    
          pno3 = 0.                                                             
        endif                                                                   
        q(6) = pno3 * rmwt(6) / pivol                                           
                                                                                
! --- 6.42_x1.1                                                                 
! ---   Condition NO3 (q(6)) to be no smaller than zero                         
        if(q(6).LT.0.0) then                                                    
           write(io6,*)'CHEMRIV6 Warning: reset NO3(g) from ',                 &
     &                 q(6),' to ZERO'                                          
           q(6)=0.0                                                             
        endif                                                                   
! ---   Condition NO3 to be no larger than TNO3 (q(5))                          
! ---   (q(5)>0.0 in this block)                                                
        test=q(6)/q(5)-1.0                                                      
        if(test.GT.1.0e-06) then                                                
           write(io6,*)'CHEMRIV6 Warning: reset NO3(g) from ',                 &
     &                 q(6),' to ',q(5)                                         
        endif                                                                   
        q(6)=MIN(q(5),q(6))                                                     
                                                                                
      endif                                                                     
!                                                                               
! --- SOA equilibrium                                                           
! --- 6.42_x1.1                                                                 
      call soadriver(io6,q,temp,presur,pivol)                                   
!                                                                               
! --- Transfer revised mass data to original array                              
      do is = 1,iss                                                             
         qin(is) = q(is)                                                        
      end do                                                                    
                                                                                
      return                                                                    
      end                                                                       
                                                                                
!----------------------------------------------------------------------         
      subroutine chmriv6(PM,TEMPER,PRESUR,RH,O3PUFF,RK1,ZCOEF,                 &
     &                   TSTP,CONFCT,RSHET)                                     
!----------------------------------------------------------------------         
!                                                                               
! --- CALPUFF    Version: TNG-7.1.0    Level: 100917           CHMRIV6          
!                P. Karamchandani, AER (Adapted from chemriv)                   
!                                                                               
! --- ADOPTED FROM ARM3 (see banner below)                                      
!                Sulfate & nitrate conversion uses exponential fcn              
!                Unused variables are deactivated                               
!                                                                               
! --- UPDATE                                                                    
! --- V5.8-V6.302  100917 (DGS): Restructure value**-n to value**(-n)           
!                                to satisfy LF95 compiler                       
!                  100917 (DGS): Add heterogeneous reaction rate for            
!                                SO4 and pass RSHET(fraction/hr)                
! --- V5.8         071025  (PK): Original                                       
!                                                                               
!----------------------------------------------------------------------         
      IMPLICIT NONE                                                             
      SAVE                                                                      
!                                                                               
!DECK.CHMRIV                                                                    
!                                                                               
!  DATE: NOVEMBER 1987                                                          
!  VERSION: ARM3-1.0                                                            
!                                                                               
!     CALCULATE TRANSFORMATION OF SO2 TO SO4 AND NO2 TO HNO3-NO3                
!     FOR ONE TIME STEP, TSTP                                                   
!     BASED ON RIVAD CHEMICAL MECHANISM                                         
!     ADPATED FROM THE RIVAD MODEL                                              
!                                                                               
!                                                                               
!  INPUT ARGUMENTS:                                                             
!                   PM     R  MASS OF SPECIES (I) IN PUFF (G)                   
!                             SO2, SO4, NO, NO2, HNO3-NO3, TSP                  
!                   TEMPER R  TEMPERATURE AT PLUME HEIGHT (K)                   
!                   PRESUR R  PRESSURE AT PLUME HEIGHT (MB)                     
!                   RH     R  RELATIVE HUMIDITY AT PLUME (%)                    
!                   O3PUFF R  PUFF OZONE CONCENTRATION (PPM)                    
!                   RK1    R  NO2 PHOTOLYSIS RATE CONSTANT (PPM/MIN)            
!                   ZCOEF  R  COSINE OF SOLAR ZENITH ANGLE                      
!                   TSTP   R  TIME STEP (HOURS)                                 
!                   CONFCT R  CONVERSION FACTOR FOR UG/M3 TO PPM-MWT            
!                   RSHET  R  SO2 heterogeneous loss rate (fraction/hr)         
!                                                                               
!  OUTPUT ARGUMENTS:                                                            
!                   PM     R  NEW MASS OF PUFF DUE TO CHANGES BASED ON          
!                                                                               
!                                                                               
!  SUBROUTINES CALLED:                                                          
!                                                                               
!  CALLED BY:   CHEMRIV6                                                        
!                                                                               
                                                                                
! --- Constants                                                                 
      REAL COEF1                  ! Molec/cc to ppm conv factor coefficient     
      PARAMETER ( COEF1 = 7.33981E+15 )                                         
                                                                                
      REAL CONSTC                 ! Constant for falloff type reaction          
      PARAMETER ( CONSTC = 0.6 )                                                
                                                                                
      REAL TI300                  ! 1.0 / 300.                                  
      PARAMETER ( TI300 = 1.0 / 300.0 )                                         
                                                                                
! --- Arguments                                                                 
      REAL PM(6)                                                                
      REAL TEMPER, PRESUR, RH, O3PUFF, RK1, ZCOEF, TSTP, CONFCT, RSHET          
 
! --- wangzhm save 
                                                                                
! --- Locals                                                                    
      REAL,save :: TFACT, R26, ROHM, QJ, PHIKK, RCONST                                  
      INTEGER,save :: L                                                                 
      REAL,save :: CNO, CNO2, CNOX, TAMB, PAMB, RH1, H2O, SUM                           
      REAL,save :: XNO, XNOX, XO3, XOX, XNO2, RNO2X, XNO3, XN2O5                        
      REAL,save :: RSULF, RNITR, RKX, RNO3, RNITRN, RNITRD, SULFN                       
      REAL,save :: XOHMAX, XOH                                                          
      REAL,save :: A0, A, B, C                                                          
                                                                                
      REAL,save :: RK0               ! K0 in falloff rate expressions                   
      REAL,save :: RKINF             ! KINF in falloff rate expressions                 
      REAL,save :: XEND              ! Exponent in falloff rate expressions             
!                                                                            
      REAL,save :: KSO2OH            ! SO2 + OH rate constant                           
      REAL,save :: KNO2OH            ! NO2 + OH rate constant                           
!                                                                               
      REAL*8,save :: TINV, CFACT, RFACT                                                 
!                                                                               
      REAL,save :: PPM(5),SMWT(5)                                                       
      LOGICAL,save :: LFIRST                                                            
      DATA LFIRST/.TRUE./                                                       
      DATA SMWT/64.0,96.0,30.0,46.0,63.0/                                       
!                                                                               
!     FIRST TIME THROUGH SET UP SOME GLOBAL RATE CONSTANTS                      
!
! --- wangzhm omp: set threadsprivate                                           
      !$OMP THREADPRIVATE(TFACT,R26,ROHM,QJ,PHIKK,RCONST,L,&
      !$OMP& CNO,CNO2,CNOX,TAMB,PAMB,RH1,H2O,SUM,&
      !$OMP& XNO,XNOX,XO3,XOX,XNO2,RNO2X,XNO3,XN2O5,&
      !$OMP& RSULF,RNITR,RKX,RNO3,RNITRN,RNITRD,SULFN,&
      !$OMP& XOHMAX,XOH,A0,A,B,C,RK0,RKINF,XEND,KSO2OH,KNO2OH,&
      !$OMP& TINV,CFACT,RFACT,ppm,smwt,lfirst)
                                                                               
      IF (LFIRST) THEN                                                          
         TFACT = 1./273.                                                        
         R26 = 1./26.4                                                          
         ROHM = 4.87E-7/1.32E-3                                                 
         LFIRST = .FALSE.                                                       
      END IF                                                                    
!                                                                               
      DO L = 1,5                                                                
         PPM(L) = 1.E6*CONFCT*PM(L)/SMWT(L)                                     
      END  DO                                                                   
!                                                                               
      CNO = PPM(3)                                                              
      CNO2 = PPM(4)                                                             
      CNOX = CNO + CNO2                                                         
!                                                                               
      TAMB = TEMPER                                                             
      PAMB = PRESUR/1013.0                                                      
      TINV = 1. / TAMB                                                          
      CFACT = COEF1 * PAMB * TINV                                               
!                                                                               
      IF (RK1.LE.0.0) THEN                                                      
         QJ = 0.0                                                               
         PHIKK = 0.0                                                            
      ELSE                                                                      
         PHIKK = RK1*R26                                                        
         QJ = (1.338E-3)*ZCOEF**2.74                                            
      ENDIF                                                                     
      IF (TAMB.GE.273.) THEN                                                    
         RCONST = 18.02*(597.3-.566*(TAMB-273.))/1.9869                         
      ELSE                                                                      
         RCONST = 6133.17                                                       
      END IF                                                                    
      RH1 = MIN(RH,95.)                                                         
      RH1 = MAX(RH1,0.)                                                         
      H2O = (6030.*.01*RH1/PAMB)*EXP(RCONST*(TFACT-TINV))                       
!                                                                               
!  DO SIMPLE CHEMISTRY                                                          
!                                                                               
      XNOX = CNOX                                                               
      XOX = O3PUFF + CNO2                                                       
      SUM = XOX + XNOX + PHIKK                                                  
      XNO2 = 0.5*(SUM-SQRT(ABS(SUM*SUM-4.*XNOX*XOX)))                           
      XNO2 = MIN(XNO2,XNOX)                                                     
      XNO2 = MAX(XNO2,0.)                                                       
! dgs      XNO=XNOX-XNO2                                                        
! dgs      IF (XNO.LT.0.0) XNO=0.0                                              
      RNO2X = 0.                                                                
      IF (XNOX.GT.0.) RNO2X = XNO2/XNOX                                         
                                                                                
! dgs Use NO2/NOX ratio to test for zero NO (single precision)                  
      xno = xnox - xno2                                                         
      if(rno2x.GE.0.999999) XNO=0.0                                             
                                                                                
      XO3 = XOX - XNO2                                                          
      XO3 = MAX(XO3,0.)                                                         
                                                                                
      IF (ZCOEF.LT.0.06975) THEN                                                
!                                                                               
!  NIGHTTIME CHEMISTRY                                                          
!                                                                               
         RSULF = 0.                                                             
!                                                                               
!  IMPLEMENT NEW NIGHTTIME CHEMISTRY - 1/86                                     
!                                                                               
         RNITR = 0.0                                                            
         XO3 = MAX(XO3,0.)                                                      
         IF (XO3.GT.0.) THEN                                                    
            RKX = 1780./(1.9E-6*H2O + 3.12)                                     
            RNO3 = 0.086                                                        
            A0 = 0.59 + 1780. -3.12*RKX                                         
            A = 2.*RKX*RNO3-A0                                                  
            B = RNO3 + 0.0474*XO3 +XNO2*A0                                      
            C = -0.0474*XNO2*XO3                                                
!                                                                               
            XNO3 = -2.*C/(B+SQRT(B*B-4.*A*C))                                   
            XN2O5 = XNO3*XNO2*RKX                                               
            RNITRN = 2.*1.9E-6*H2O*XN2O5*60.*TSTP                               
            RNITRN = MIN(XNO2,RNITRN)                                           
            PPM(5) = PPM(5) + RNITRN                                            
            XNO2 = XNO2 - RNITRN                                                
         END IF                                                                 
      ELSE                                                                      
!                                                                               
!  DO DAYTIME CHEMISTRY                                                         
!                                                                               
         XOHMAX = ROHM*QJ                                                       
         XOH = XOHMAX                                                           
         IF (PPM(1).NE.0..OR.XNO2.NE.0.)                                       &
     &   XOH = 2.*QJ*3.4E5*H2O*XO3/((4.45E10+3.4E5*H2O)*                       &
     &       (2000.*PPM(1) + 14000.*XNO2))                                      
         XOH = MIN(XOH,XOHMAX)                                                  
                                                                                
! --- SO2 + OH rate constant (falloff expression)                               
         RK0 = 1.0E+06 * CFACT * 3.0E-31 * ( TAMB * TI300 )**(-3.3)             
         RKINF = 1.5E-12                                                        
         XEND = 1.0 / ( 1.0 + ( LOG10( RK0 / RKINF ) )**2 )                     
         KSO2OH = ( RK0 / ( 1.0 + RK0 / RKINF ) ) * CONSTC**XEND                
                                                                                
! --- NO2 + OH rate constant (falloff expression)                               
         RK0 = 1.0E+06 * CFACT * 2.6E-30 * ( TAMB * TI300 )**(-3.2)             
         RKINF = 2.4E-11 * ( TAMB * TI300 )**(-1.3)                             
         XEND = 1.0 / ( 1.0 + ( LOG10( RK0 / RKINF ) )**2 )                     
         KNO2OH = ( RK0 / ( 1.0 + RK0 / RKINF ) ) * CONSTC**XEND                
!                                                                               
! --- Convert rate constants from molec-cc-1 s-1 to ppm-1 hr-1 units            
         rfact = 2.64e19 * (pamb * tinv)                                        
         kso2oh = kso2oh * rfact                                                
         kno2oh = kno2oh * rfact                                                
         rsulf = 1. - EXP(-kso2oh*xoh*tstp)                                     
         rnitr = 1. - EXP(-kno2oh*xoh*tstp)                                     
         RNITRD = RNITR*XNO2                                                    
         RNITRD = MIN(XNOX,RNITRD)                                              
         PPM(5) = PPM(5)+RNITRD                                                 
         XNOX = XNOX-RNITRD                                                     
         XNOX =  MAX(XNOX,0.)                                                   
         XNO = XNOX*(1.-RNO2X)                                                  
         XNO2 = XNOX*RNO2X                                                      
                                                                                
      END IF                                                                    
!                                                                               
! --- Add heterogeneous rate wso4=rshet*ppm(1)*tstp                             
      SULFN = PPM(1)*RSULF + rshet*ppm(1)*tstp                                  
      SULFN = MIN(PPM(1),SULFN)                                                 
      PPM(2) = PPM(2) + SULFN                                                   
      PPM(1) = PPM(1) - SULFN                                                   
!                                                                               
!     UPDATE NEW STEADY-STATE NO AND NO2 VALUES                                 
!                                                                               
      PPM(3) = XNO                                                              
      PPM(4) = XNO2                                                             
!                                                                               
!     UPDATE NEW O3 VALUE                                                       
      O3PUFF = XO3                                                              
!                                                                               
!     CONVERT BACK TO GRAMS OF SPECIES IN PUFF                                  
!                                                                               
      DO L = 1,5                                                                
         PM(L) = 1.E-6*PPM(L)*SMWT(L)/CONFCT                                    
      END DO                                                                    
      RETURN                                                                    
      END                                                                       
                                                                                
!----------------------------------------------------------------------         
      subroutine chmriv7(PM,TEMPER,PRESUR,RH,O3PUFF,RK1,ZCOEF,                 &
     &                   TSTP,CONFCT,RSHET)                                     
!----------------------------------------------------------------------         
!                                                                               
! --- CALPUFF    Version: TNG-7.1.0    Level: 100917           CHMRIV7          
!                P. Karamchandani, AER (Adapted from chemriv)                   
!                                                                               
! --- ADOPTED FROM ARM3 (see banner below)                                      
!                Sulfate & nitrate conversion uses exponential fcn              
!                Unused variables are deactivated                               
!                                                                               
! --- UPDATE                                                                    
! --- V5.8-V6.302  100917 (DGS): Restructure value**-n to value**(-n)           
!                                to satisfy LF95 compiler                       
!                         (DGS): Add heterogeneous reaction rate for            
!                                SO4 and pass RSHET(fraction/hr)                
! --- V5.8         071025  (PK): Original                                       
!                                                                               
!----------------------------------------------------------------------         
      IMPLICIT NONE                                                             
      SAVE                                                                      
!                                                                               
!DECK.CHMRIV                                                                    
!                                                                               
!  DATE: NOVEMBER 1987                                                          
!  VERSION: ARM3-1.0                                                            
!                                                                               
!     CALCULATE TRANSFORMATION OF SO2 TO SO4 AND NO2 TO HNO3-NO3                
!     FOR ONE TIME STEP, TSTP                                                   
!     BASED ON RIVAD CHEMICAL MECHANISM                                         
!     ADPATED FROM THE RIVAD MODEL                                              
!                                                                               
!                                                                               
!  INPUT ARGUMENTS:                                                             
!                   PM     R  MASS OF SPECIES (I) IN PUFF (G)                   
!                             SO2, SO4, NO, NO2, HNO3-NO3, TSP                  
!                   TEMPER R  TEMPERATURE AT PLUME HEIGHT (K)                   
!                   PRESUR R  PRESSURE AT PLUME HEIGHT (MB)                     
!                   RH     R  RELATIVE HUMIDITY AT PLUME (%)                    
!                   O3PUFF R  PUFF OZONE CONCENTRATION (PPM)                    
!                   RK1    R  NO2 PHOTOLYSIS RATE CONSTANT (PPM/MIN)            
!                   ZCOEF  R  COSINE OF SOLAR ZENITH ANGLE                      
!                   TSTP   R  TIME STEP (HOURS)                                 
!                   CONFCT R  CONVERSION FACTOR FOR UG/M3 TO PPM-MWT            
!                   RSHET  R  SO2 heterogeneous loss rate (fraction/hr)         
!                                                                               
!  OUTPUT ARGUMENTS:                                                            
!                   PM     R  NEW MASS OF PUFF DUE TO CHANGES BASED ON          
!                                                                               
!                                                                               
!  SUBROUTINES CALLED:                                                          
!                                                                               
!  CALLED BY:   CHEMRIV7                                                        
!                                                                               
                                                                                
! --- Constants                                                                 
      REAL COEF1                  ! Molec/cc to ppm conv factor coefficient     
      PARAMETER ( COEF1 = 7.33981E+15 )                                         
                                                                                
      REAL CONSTC                 ! Constant for falloff type reaction          
      PARAMETER ( CONSTC = 0.6 )                                                
                                                                                
      REAL TI300                  ! 1.0 / 300.                                  
      PARAMETER ( TI300 = 1.0 / 300.0 )                                         
                                                                                
! --- Arguments                                                                 
      REAL PM(25)                                                               
      REAL TEMPER, PRESUR, RH, O3PUFF, RK1, ZCOEF, TSTP, CONFCT, RSHET          
                                                                                
! --- Locals 

! --- wangzhm save                                                                   
      REAL,save :: TFACT, R26, ROHM, QJ, PHIKK, RCONST                                  
      INTEGER,save :: L                                                                 
      REAL,save :: CNO, CNO2, CNOX, TAMB, PAMB, RH1, H2O, SUM                           
      REAL,save :: XNO, XNOX, XO3, XOX, XNO2, RNO2X, XNO3, XN2O5                        
      REAL,save :: RSULF, RNITR, RKX, RNO3, RNITRN, RNITRD, SULFN                       
      REAL,save :: XOHMAX, XOH                                                          
      REAL,save :: A0, A, B, C                                                          
                                                                                
      REAL,save :: RK0               ! K0 in falloff rate expressions                   
      REAL,save :: RKINF             ! KINF in falloff rate expressions                 
      REAL,save :: XEND              ! Exponent in falloff rate expressions             
!                                                                           
      REAL,save :: KSO2OH            ! SO2 + OH rate constant                           
      REAL,save :: KNO2OH            ! NO2 + OH rate constant                           
      REAL,save :: KTOLOH            ! TOL + OH rate constant                           
      REAL,save :: KXYLOH            ! XYL + OH rate constant                           
      REAL,save :: KALKHOH           ! ALKH + OH rate constant                          
      REAL,save :: KPAHOH            ! PAH + OH rate constant                           
!                                                                               
      REAL,save :: DTOL, DXYL, DALKH, DPAH                                              
!                                                                               
      REAL*8,save :: TINV, CFACT, RFACT                                                 
!                                                                               
      REAL,save :: PPM(25),SMWT(25)                                                     
      LOGICAL,save :: LFIRST                                                            
      DATA LFIRST/.TRUE./                                                       
      DATA SMWT/64.0,96.0,30.0,46.0,63.0,62.0,180.0,                           &
     &          5*92.0,5*106.0,3*226.0,5*156.0/                                 

! --- wangzhm omp: set threadsprivate                                           
      !$OMP THREADPRIVATE(TFACT,R26,ROHM,QJ,PHIKK,RCONST,L,&
      !$OMP& CNO,CNO2,CNOX,TAMB,PAMB,RH1,H2O,SUM,XNO,XNOX,XO3,XOX,XNO2,RNO2X,XNO3,XN2O5,&
      !$OMP& RSULF,RNITR,RKX,RNO3,RNITRN,RNITRD,SULFN,XOHMAX,XOH,A0,A,B,C,RK0,RKINF,XEND,&
      !$OMP& KSO2OH,KNO2OH,KTOLOH,KXYLOH,KALKHOH,KPAHOH,DTOL,DXYL,DALKH,DPAH,TINV,CFACT,RFACT,&
      !$OMP& PPM,SMWT,LFIRST)                                                                                 
!                                                                               
!     FIRST TIME THROUGH SET UP SOME GLOBAL RATE CONSTANTS                      
!                                                                               
      IF (LFIRST) THEN                                                          
         TFACT = 1./273.                                                        
         R26 = 1./26.4                                                          
         ROHM = 4.87E-7/1.32E-3                                                 
         LFIRST = .FALSE.                                                       
      END IF                                                                    
!                                                                               
      DO L = 1,5                                                                
         PPM(L) = 1.E6*CONFCT*PM(L)/SMWT(L)                                     
      END  DO                                                                   
!                                                                               
      DO L = 8,25                                                               
         PPM(L) = 1.E6*CONFCT*PM(L)/SMWT(L)                                     
      END  DO                                                                   
!                                                                               
      CNO = PPM(3)                                                              
      CNO2 = PPM(4)                                                             
      CNOX = CNO + CNO2                                                         
!                                                                               
      TAMB = TEMPER                                                             
      PAMB = PRESUR/1013.0                                                      
      TINV = 1. / TAMB                                                          
      CFACT = COEF1 * PAMB * TINV                                               
!                                                                               
      IF (RK1.LE.0.0) THEN                                                      
         QJ = 0.0                                                               
         PHIKK = 0.0                                                            
      ELSE                                                                      
         PHIKK = RK1*R26                                                        
         QJ = (1.338E-3)*ZCOEF**2.74                                            
      ENDIF                                                                     
      IF (TAMB.GE.273.) THEN                                                    
         RCONST = 18.02*(597.3-.566*(TAMB-273.))/1.9869                         
      ELSE                                                                      
         RCONST = 6133.17                                                       
      END IF                                                                    
      RH1 = MIN(RH,95.)                                                         
      RH1 = MAX(RH1,0.)                                                         
      H2O = (6030.*.01*RH1/PAMB)*EXP(RCONST*(TFACT-TINV))                       
!                                                                               
!  DO SIMPLE CHEMISTRY                                                          
!                                                                               
      XNOX = CNOX                                                               
      XOX = O3PUFF + CNO2                                                       
      SUM = XOX + XNOX + PHIKK                                                  
      XNO2 = 0.5*(SUM-SQRT(ABS(SUM*SUM-4.*XNOX*XOX)))                           
      XNO2 = MIN(XNO2,XNOX)                                                     
      XNO2 = MAX(XNO2,0.)                                                       
! dgs      XNO=XNOX-XNO2                                                        
! dgs      IF (XNO.LT.0.0) XNO=0.0                                              
      RNO2X = 0.                                                                
      IF (XNOX.GT.0.) RNO2X = XNO2/XNOX                                         
                                                                                
! dgs Use NO2/NOX ratio to test for zero NO (single precision)                  
      xno = xnox - xno2                                                         
      if(rno2x.GE.0.999999) XNO = 0.0                                           
                                                                                
      XO3 = XOX - XNO2                                                          
      XO3 = MAX(XO3,0.)                                                         
      IF (ZCOEF.LT.0.06975) THEN                                                
!                                                                               
!  NIGHTTIME CHEMISTRY                                                          
!                                                                               
         RSULF = 0.                                                             
!                                                                               
!  IMPLEMENT NEW NIGHTTIME CHEMISTRY - 1/86                                     
!                                                                               
         RNITR = 0.0                                                            
         XO3 = MAX(XO3,0.)                                                      
         IF (XO3.GT.0.) THEN                                                    
            RKX = 1780./(1.9E-6*H2O + 3.12)                                     
            RNO3 = 0.086                                                        
            A0 = 0.59 + 1780. -3.12*RKX                                         
            A = 2.*RKX*RNO3-A0                                                  
            B = RNO3 + 0.0474*XO3 +XNO2*A0                                      
            C = -0.0474*XNO2*XO3                                                
!                                                                               
            XNO3 = -2.*C/(B+SQRT(B*B-4.*A*C))                                   
            XN2O5 = XNO3*XNO2*RKX                                               
            RNITRN = 2.*1.9E-6*H2O*XN2O5*60.*TSTP                               
            RNITRN = MIN(XNO2,RNITRN)                                           
            PPM(5) = PPM(5) + RNITRN                                            
            XNO2 = XNO2 - RNITRN                                                
         END IF                                                                 
      ELSE                                                                      
!                                                                               
!  DO DAYTIME CHEMISTRY                                                         
!                                                                               
         XOHMAX = ROHM*QJ                                                       
         XOH = XOHMAX                                                           
         IF (PPM(1).NE.0..OR.XNO2.NE.0.)                                       &
     &   XOH = 2.*QJ*3.4E5*H2O*XO3/((4.45E10+3.4E5*H2O)*                       &
     &     (2000.*PPM(1) + 14000.*XNO2))                                        
         XOH = MIN(XOH,XOHMAX)                                                  
                                                                                
! --- SO2 + OH rate constant (falloff expression)                               
         RK0 = 1.0E+06 * CFACT * 3.0E-31 * ( TAMB * TI300 )**(-3.3)             
         RKINF = 1.5E-12                                                        
         XEND = 1.0 / ( 1.0 + ( LOG10( RK0 / RKINF ) )**2 )                     
         KSO2OH = ( RK0 / ( 1.0 + RK0 / RKINF ) ) * CONSTC**XEND                
!                                                                               
! --- NO2 + OH rate constant (falloff expression)                               
         RK0 = 1.0E+06 * CFACT * 2.6E-30 * ( TAMB * TI300 )**(-3.2)             
         RKINF = 2.4E-11 * ( TAMB * TI300 )**(-1.3)                             
         XEND = 1.0 / ( 1.0 + ( LOG10( RK0 / RKINF ) )**2 )                     
         KNO2OH = ( RK0 / ( 1.0 + RK0 / RKINF ) ) * CONSTC**XEND                
!                                                                               
! --- TOL + OH rate constant                                                    
         ktoloh = 2.1e-12 * EXP(322.*tinv)                                      
!                                                                               
! --- XYL + OH rate constant                                                    
         kxyloh = 1.7e-11 * EXP(116.*tinv)                                      
!                                                                               
! --- ALKH + OH rate constant                                                   
         kalkhoh = 1.97e-11                                                     
!                                                                               
! --- PAH + OH rate constant                                                    
         kpahoh = 7.7e-11                                                       
!                                                                               
! --- Convert rate constants from molec-cc-1 s-1 to ppm-1 hr-1 units            
         rfact = 2.64e19 * (pamb * tinv)                                        
         kso2oh = kso2oh * rfact                                                
         kno2oh = kno2oh * rfact                                                
         rsulf = 1. - EXP(-kso2oh*xoh*tstp)                                     
         rnitr = 1. - EXP(-kno2oh*xoh*tstp)                                     
         RNITRD = RNITR*XNO2                                                    
         RNITRD = MIN(XNOX,RNITRD)                                              
         PPM(5) = PPM(5) + RNITRD                                               
         XNOX = XNOX - RNITRD                                                   
         XNOX = MAX(XNOX,0.)                                                    
         XNO = XNOX*(1.-RNO2X)                                                  
         XNO2 = XNOX*RNO2X                                                      
                                                                                
         ktoloh = ktoloh * rfact                                                
         kxyloh = kxyloh * rfact                                                
         kalkhoh = kalkhoh * rfact                                              
         kpahoh = kpahoh * rfact                                                
!                                                                               
! --- Change in SOA precursor concentrations                                    
         dtol = (1. - EXP(-ktoloh*xoh*tstp)) * ppm(8)                           
         dxyl = (1. - EXP(-kxyloh*xoh*tstp)) * ppm(13)                          
         dalkh = (1. - EXP(-kalkhoh*xoh*tstp)) * ppm(18)                        
         dpah = (1. - EXP(-kpahoh*xoh*tstp)) * ppm(21)                          
                                                                                
         ppm(8) = ppm(8) - dtol          ! TOL                                  
         ppm(9) = ppm(9) + 0.071*dtol    ! TOLAER1                              
         ppm(10) = ppm(10) + 0.138*dtol  ! TOLAER2                              
                                                                                
         ppm(13) = ppm(13) - dxyl        ! XYL                                  
         ppm(14) = ppm(14) + 0.038*dxyl  ! XYLAER1                              
         ppm(15) = ppm(15) + 0.167*dxyl  ! XYLAER2                              
                                                                                
         ppm(18) = ppm(18) - dalkh       ! ALKH                                 
         ppm(19) = ppm(19) + 1.173*dalkh ! ALKHAER                              
                                                                                
         ppm(21) = ppm(21) - dpah        ! PAH                                  
         ppm(22) = ppm(22) + 0.156*dpah  ! PAHAER1                              
         ppm(23) = ppm(23) + 0.777*dpah  ! PAHAER2                              
      END IF                                                                    
                                                                                
! --- Add heterogeneous rate wso4=rshet*ppm(1)*tstp                             
      SULFN = PPM(1)*RSULF + rshet*ppm(1)*tstp                                  
      SULFN = MIN(PPM(1),SULFN)                                                 
      PPM(2) = PPM(2) + SULFN                                                   
      PPM(1) = PPM(1) - SULFN                                                   
!                                                                               
!     UPDATE NEW STEADY-STATE NO AND NO2 VALUES                                 
!                                                                               
      PPM(3) = XNO                                                              
      PPM(4) = XNO2                                                             
!                                                                               
!     UPDATE NEW O3 VALUE                                                       
      O3PUFF = XO3                                                              
!                                                                               
!     CONVERT BACK TO GRAMS OF SPECIES IN PUFF                                  
!                                                                               
      DO L = 1,5                                                                
         PM(L) = 1.E-6*PPM(L)*SMWT(L)/CONFCT                                    
      END DO                                                                    
!                                                                               
      DO L = 8,25                                                               
         PM(L) = 1.E-6*PPM(L)*SMWT(L)/CONFCT                                    
      END DO                                                                    
!                                                                               
      RETURN                                                                    
      END                                                                       
                                                                                
!---------------------------------------------------------------------          
      subroutine setbckoc(ndathr)                                               
!---------------------------------------------------------------------          
!                                                                               
! --- CALPUFF    Version: TNG-7.1.0    Level: 071025          SETBCKOC          
!                P. Karamchandani, AER (Adapted from setsoa)                    
!                                                                               
! --- PURPOSE:  Background OC concentration (ug/m3) for current hour            
!                                                                               
! --- INPUTS:                                                                   
!            NDATHR - integer    - YYJJJHH date-time for hour                   
!                                                                               
!     Common Block /CHEMDAT/ variables:                                         
!        BCKPMF(12),OFRAC(12)                                                   
!                                                                               
! --- OUTPUT:                                                                   
!                                                                               
!     Common Block /NEWSOA/ variables: bckoc                                    
!                                                                               
! --- SETBCKOC called by:  COMP                                                 
! --- SETBCKOC calls:      GRDAY                                                
!----------------------------------------------------------------------         
      include 'params.puf'                                                      
      include 'chemdat.puf'                                                     
      include 'newsoa.puf'                                                      
                                                                                
! --- Extract month from date-time                                              
      iyr=ndathr/100000                                                         
      ijul=ndathr/100 - 1000*iyr                                                
      call GRDAY(io6,iyr,ijul,imo,iday)                                         
                                                                                
! --- Set background OC data for this month                                     
      bckoc=bckpmf(imo)*ofrac(imo)                                              
                                                                                
      return                                                                    
      end                                                                       
                                                                                
!----------------------------------------------------------------------         
      subroutine isodriver(so4,no3,nh3,nrh,ntempk,pno3c)                        
!----------------------------------------------------------------------         
!                                                                               
! --- CALPUFF    Version: TNG-7.1.0    Level: 140521         ISODRIVER          
!                P. Karamchandani, AER                                          
!                                                                               
! --- PURPOSE:  Driver routine to calculate gas/particle equilibrium            
!               using ISORROPIA aerosol equilibrium module                      
!               This version uses double precision                              
!                                                                               
! --- UPDATES:                                                                  
!                                                                               
! --- V6.4-V6.42_x1.1  140521  : Call ISOROPIA with METASTABLE ON               
!                                                                               
! --- V6.302-V6.4  101025 (DGS): Revise ISOROPIA arguments for V2.1             
!                                                                               
! --- INPUTS:                                                                   
!            so4    - real - Total sulfate concentration (mole/m3)              
!            no3    - real - Total nitrate concentration (mole/m3)              
!            nh3    - real - Total ammonia concentration (mole/m3)              
!            nrh    - real - Fractional relative humidity                       
!            ntempk - real - Temperature (K)                                    
!                                                                               
! --- OUTPUT:                                                                   
!            pno3c  - real - Particle-phase equilibrium NO3                     
!                            concentration (mole/m3)                            
!                                                                               
! --- ISODRIVER called by:  CHEMRIV6, CHEMRIV7                                  
! --- ISODRIVER calls:      ISOROPIA                                            
!----------------------------------------------------------------------         
                                                                                
!***********************************************************************        
!                                                                      *        
!  REVISION HISTORY:                                                   *        
!     Version 1.0 developed November 2006 by PK, AER, Inc. for         *        
!     CALPUFF for API Contract No. 2006-102376                         *        
!                                                                      *        
!***********************************************************************        
!                                                                               
!...........  INCLUDES                                                          
                                                                                
      INCLUDE 'isrpia.inc'                                                      
                                                                                
!...........  ARGUMENTS and their descriptions                                  
                                                                                
      REAL SO4             ! Total sulfate (moles/m3)                           
      REAL NO3             ! Total nitrate (moles/m3)                           
      REAL NH3             ! Total ammonia (moles/m3)                           
      REAL NRH             ! Relative humidity as fraction                      
      REAL NTEMPK          ! Temperature in Kelvin                              
                                                                                
      REAL PNO3C           ! Particle-phase NO3 concentration (mole/m**3)       
                                                                                
! --- Locals  

! --- wangzhm save                                                                  
      REAL,save :: GNO3C           ! Gas-phase concentration of HNO3 in mole/m**3 air   
                                                                                
      INTEGER NCTRL, NOTHER                                                     
!V1.7 PARAMETER(NCTRL = 2,NOTHER = 6)                                           
!V2.1 PARAMETER(NCTRL = 2,NOTHER = 9)                                           
      PARAMETER(NCTRL = 2,NOTHER = 9)                                           
                                                                                
! --- Gas-phase concentration array in moles/m**3 air                           
      REAL*8,save ::  GAS(NGASAQ)                                                       
                                                                                
! --- Aqueous-phase concentration array in moles/m**3 air                       
      REAL*8,save ::  AERLIQ(NIONS+NGASAQ+2)                                            
                                                                                
! --- Solid-phase concentration array in moles/m**3 air                         
      REAL*8,save ::  AERSLD(NSLDS)                                                     
                                                                                
! --- Flag for different types of problems solved and                           
! --- different state of aerosols (deliquescent or metastable)                  
      REAL*8,save ::  CNTRL(NCTRL)                                                      
                                                                               &
! Sol&tion information array (see ISOCOM.f for details)                         
      REAL*8,save ::  OTHER(NOTHER)                                                     
                                                                               &
! Tot&l species concentrations in moles/m**3 air                                
      REAL*8,save ::  WI(NCOMP)                                                         
                                                                               &
! Tem&erature and humidity                                                      
      REAL*8,save ::  TEMPI, RHI
      
! --- wangzhm omp: set threadsprivate                                           
      !$OMP THREADPRIVATE(GNO3C,GAS,AERLIQ,AERSLD,CNTRL,OTHER,WI,TEMPI,RHI)    
      
                                                                 
!***********************************************************************        
!  begin body of subroutine                                                     
!                                                                               
! *** Assign input total concentrations to WI                                   
                                                                                
      WI(1) = 0.    ! Na concentration                                          
      WI(2) = so4                                                               
      WI(3) = nh3                                                               
      WI(4) = no3                                                               
      WI(5) = 0.    ! Cl concentration                                          
!V2.1 Additional 3 are not used here                                            
      WI(6) = 0.    ! total calcium concentration                               
      WI(7) = 0.    ! total potassium concentration                             
      WI(8) = 0.    ! total magnesium concentration                             
                                                                                
      RHI   = NRH                                                               
      TEMPI = NTEMPK                                                            
!                                                                               
! *** CALL ISORROPIA                                                            
!                                                                               
      CNTRL(1) = 0.     ! 0 = FORWARD PROBLEM,  1 = REVERSE PROBLEM             
                                                                                
! --- 6.42_x1.1                                                                 
      CNTRL(2) = 1.     ! 0 = SOLID + LIQUID AEROSOL,  1 = METASTABLE           
! --- wangzhm omp: critical
      !$OMP CRITICAL (call_isoropia)                                                                               
      CALL ISOROPIA (WI,  RHI,  TEMPI,   CNTRL,                                &
     &               W,   GAS,  AERLIQ,  AERSLD,  SCASE,  OTHER)  
      !$OMP END CRITICAL (call_isoropia)
                   
!                                                                               
! *** SAVE RESULTS                                                              
                                                                               &
! Gas&phase HNO3 (moles/m3)                                                     
      GNO3C  = GAS(2)                                                           
                                                                               &
! Par&icle-phase NO3 (moles/m3)                                                 
      PNO3C  = MAX(0.,NO3 - GNO3C)                                              
!                                                                               
      RETURN                                                                    
      END                                                                       
                                                                                
                                                                                
!----------------------------------------------------------------------         
      subroutine soadriver(ilst,q,tempk,press,pivol)                            
!----------------------------------------------------------------------         
                                                                                
!***********************************************************************        
!  FUNCTION: Driver routine to calculate SOA gas-particle equilibrium  *        
!                                                                      *        
!                                                                               
! --- UPDATES:                                                                  
!                                                                               
! --- V6.4-V6.42_x1.1  140521  : Iterate on CALTECH_SOA call to update          
!                                absorbing organic mass                         
!                                                                               
!  INPUTS:                                                             *        
!     ilst   - List-file unit                                          *        
!     q      - Mass of species in puff                                 *        
!     tempk  - Temperature (K)                                         *        
!     press  - Pressure at plume height (mb)                           *        
!     pivol  - Reciprocal of puff volume (1/m**3)                      *        
!                                                                      *        
!  RETURN VALUES:                                                      *        
!     q      - Adjusted mass of species in puff                        *        
!                                                                      *        
!  REVISION HISTORY:                                                   *        
!     Version 1.0 developed November 2006 by PK, AER, Inc. for         *        
!     CALPUFF for API Contract No. 2006-102376                         *        
!***********************************************************************        
      IMPLICIT NONE                                                             
!                                                                               
!...........  INCLUDES                                                          
                                                                                
      INCLUDE 'newsoa.puf'                                                      
      INCLUDE 'soadat.puf'                                                      
                                                                                
!...........  ARGUMENTS and their descriptions                                  
                                                                                
! --- 6.42_x1.1                                                                 
      integer ilst         ! List-file unit                                     
                                                                                
      REAL q(25)           ! Puff masses (g)                                    
      REAL tempk           ! Temperature in Kelvin                              
      REAL press           ! Pressure in mb                                     
      REAL pivol           ! reciprocal of puff volume (1/m3)                   

! --- wangzhm save
                                                                               &
! Loc&l variables                                                               
      integer,save :: iorg ! loop index                                                &
! Sec&ndary organic aerosol concentrations (ug/m3)                              
! --- Gas-phase compounds                                                       
      REAL,save :: gasorg(NORG)                                                         
! --- Particle-phase compounds                                                  
      REAL,save :: partorg(NORG)                                                        
! --- Total (gas-phase + particle-phase)                                        
      REAL,save :: worg(NORG+1)  ! additional species at end for primary OC             
                                                                                
! --- Conversion factor                                                         
      REAL,save :: cfact                                                                
                                                                                
! --- 6.42_x1.1                                                                 
! --- OA mass iteration                                                         
      integer,save :: iter, niter                                                       
      real,save :: psum1, psum2, pdiff                                                  

! --- wangzhm omp: set threadsprivate                                           
      !$OMP THREADPRIVATE(iorg,gasorg,partorg,worg,cfact,iter,niter,psum1,psum2,pdiff) 
                                                                                
!***********************************************************************        
!  begin body of subroutine                                                     
!                                                                               
! *** Convert input puff masses (g) to puff concentrations (ug/m3)              
                                                                                
      cfact = pivol * 1.E6                                                      
                                                                               &
! Gas&s                                                                         
      gasorg(1) = q(9) * cfact    ! TOLAER1                                     
      gasorg(2) = q(10) * cfact   ! TOLAER2                                     
      gasorg(3) = q(14) * cfact   ! XYLAER1                                     
      gasorg(4) = q(15) * cfact   ! XYLAER2                                     
      gasorg(5) = q(19) * cfact   ! ALKHAER                                     
      gasorg(6) = q(22) * cfact   ! PAHAER1                                     
      gasorg(7) = q(23) * cfact   ! PAHAER2                                     
                                                                               &
! Par&icles                                                                     
      partorg(1) = q(11) * cfact   ! ATOLA1                                     
      partorg(2) = q(12) * cfact   ! ATOLA2                                     
      partorg(3) = q(16) * cfact   ! AXYLA1                                     
      partorg(4) = q(17) * cfact   ! AXYLA2                                     
      partorg(5) = q(20) * cfact   ! AALKHA                                     
      partorg(6) = q(24) * cfact   ! APAHA1                                     
      partorg(7) = q(25) * cfact   ! APAHA2                                     
                                                                               &
! Tot&l                                                                         
      do iorg = 1, NORG                                                         
        worg(iorg) = gasorg(iorg) + partorg(iorg)                               
      end do                                                                    
                                                                               &
! Org&nic absorbing mass                                                        
      worg( NORG + 1 ) = MAX(q(7)*cfact + bckoc, 0.)                            
                                                                               &
! Cal&ulate organic aerosol formation using Caltech SOA Module                  
                                                                                
! --- 6.42_x1.1                                                                 
! --- Iterate on this call until difference in particulate OC total <1%         
! --- PDIFF is 1% criterion for difference                                      
! --- NITER is set large but finite to trap endless loop                        
! --- ITER is current loop counter                                              
      pdiff=0.01                                                                
      niter=10000                                                               
      iter=0                                                                    
10    psum1=0.0                                                                 
      do iorg = 1, NORG                                                         
        psum1=psum1+partorg(iorg)                                               
      end do                                                                    
                                                                                
      call caltech_soa(worg, gasorg, partorg, tempk)                            
                                                                                
! --- 6.42_x1.1                                                                 
      iter=iter+1                                                               
      psum2=-psum1                                                              
      do iorg = 1, NORG                                                         
        psum2=psum2+partorg(iorg)                                               
      end do                                                                    
      if(psum1.GT.0.0) then                                                     
         if(iter.LE.niter) then                                                 
            if(ABS(psum2/psum1) .GT. pdiff) goto 10                             
         else                                                                   
! ---       Report too many iterations warning and stop iterating               
            write(ilst,*)'SOADRIVER:  Warning ...'                              
            write(ilst,*)'  Iterations exceed NITER = ',niter,                 &
     &                   ' when computing particulate SOA '                     
            write(ilst,*)'  Target fractional SOA change = ',pdiff              
            write(ilst,*)'  Last fractional SOA change = ',psum2/psum1          
         endif                                                                  
      elseif(worg(norg+1).GT.0.0 .AND. iter.EQ.1) then                          
         goto 10                                                                
      endif                                                                     
                                                                                
! --- Assign adjusted concentrations back to puff masses                        
                                                                             
      cfact = 1./cfact                                                       
! Gases                                                                      
      q(9) = gasorg(1) * cfact    ! TOLAER1                                  
      q(10) = gasorg(2) * cfact   ! TOLAER2                                  
      q(14) = gasorg(3) * cfact   ! XYLAER1                                  
      q(15) = gasorg(4) * cfact   ! XYLAER2                                  
      q(19) = gasorg(5) * cfact   ! ALKHAER                                  
      q(22) = gasorg(6) * cfact   ! PAHAER1                                  
      q(23) = gasorg(7) * cfact   ! PAHAER2                                  
                                                                             
! Particles                                                                  
      q(11) = partorg(1) * cfact   ! ATOLA1                                  
      q(12) = partorg(2) * cfact   ! ATOLA2                                  
      q(16) = partorg(3) * cfact   ! AXYLA1                                     
      q(17) = partorg(4) * cfact   ! AXYLA2                                     
      q(20) = partorg(5) * cfact   ! AALKHA                                     
      q(24) = partorg(6) * cfact   ! APAHA1                                     
      q(25) = partorg(7) * cfact   ! APAHA2                                     
!                                                                               
      RETURN                                                                    
      END                                                                       
                                                                                
!----------------------------------------------------------------------         
      SUBROUTINE CALTECH_SOA (WORG, GASORG, PARTORG, CURTEMP)                   
!----------------------------------------------------------------------         
                                                                                
!***********************************************************************        
!  FUNCTION: Program to simulate Secondary Organic Aerosols            *        
!     Total no. of species = 4 + 34 = 38                               *        
!     Partition equation                                               *        
!     Kom, i = (Ai/Mo)/Gi                                              *        
!     Ai = particle-phase concentration (ug/m**3 air)                  *        
!     Gi = gas-phase concentration (ug/m**3 air)                       *        
!     MSUM = sum Ai + primary organics                                 *        
!     Kom, i has the units of m3/ug                                    *        
!  PRECONDITION REQUIRED: called from subr. AEROEQ                     *        
!  RETURN VALUES:                                                      *        
!     PARTORG(I) - Particulate concentration of organic species i      *        
!     GASORG(I)  - Gas-phase concentration of organic species i        *        
!  REVISION HISTORY:                                                   *        
!     Written by Betty K. Pun of AER, Inc. for EPRI's Aerosol          *        
!         Module Implementation Project, May, 1999                     *        
!     Revised by Yang Zhang of AER, Inc. for EPRI's Aerosol            *        
!         Module Implementation Project based on Models3's             *        
!         coding standard July, 1999                                   *        
!     Revised by Betty K. Pun of AER, Inc. November, 99                *        
!         for incorporation into 3-D model under EPRI.                 *        
!         To reduce computational requirement and to                   *        
!         facilitate numerical solution, simultaneous equations        *        
!         are not solved in this version.  Instead, the partition of   *        
!         each organic aerosol at each time depends on the amount of   *        
!         material in the organic phase at the step before the partition*       
!     Modified by Betty K. Pun, AER to include temperature dependence  *        
!         based on generic Hvap, February/March 2002                   *        
!     Revised by BKP March 2006, change reference temperature from     *        
!             310 to 298K                                              *        
!     Updated by PK December 2006, for implementation in CALPUFF       *        
!             for API Contract No. 2006-102376                         *        
!  REFERENCES:                                                         *        
!    1. Odum et al.,  97.  EST 31:1890                                 *        
!    2. Griffin et al.,  99.  JGR 104:3555                             *        
!                                                                      *        
!***********************************************************************        
                                                                                
      IMPLICIT NONE                                                             
                                                                                
!...........  INCLUDES                                                          
      INCLUDE 'soadat.puf'                                                      
                                                                                
!...........  ARGUMENTS and their descriptions and some other variables         
                                                                                
      REAL WORG(NORG+1)   ! Total organic species concentration                 
      REAL PARTORG(NORG)  ! Particulate concentration of organic species i      
      REAL GASORG(NORG)   ! Gas-phase concentration of organic species i        
      REAL CURTEMP        ! Puff temperature used in modifying parititon        
                          ! constants                                           
                                                                                
!..........  Local Variables                                                    
                                                                                
      REAL MP             ! Total organic species concentration                 
      REAL MSUM           ! Sum of particle-phase concentration (ug/m**3 air)   
      INTEGER I           ! Organic species index                               
                                                                                
!***********************************************************************        
!  begin body of subroutine                                                     
                                                                                
      MP = WORG(NORG + 1)              ! OC (primary+background)                
      MSUM = MP                                                                 
                                                                                
      DO I = 1,  NORG                                                           
         IF (WORG(I) .LT. 0.0) THEN                                           
            WRITE (*,  95000) I                                                 
95000       FORMAT ('ERROR: Negative conc. read for species',  I2)              
            WRITE (*,  *) 'Negative Conc. = ',WORG(I)                           
            call flush(6)                                                       
            STOP                                                                
         END IF                                                                 
         MSUM = MSUM + PARTORG(I)      !add starting SOA to MSUM                
      END DO                                                                    
                                                                                
      IF (MSUM .GT. 0.0) THEN                                                 
         call SOASUB (MSUM, WORG, PARTORG, GASORG, CURTEMP)                     
      ENDIF                                                                     
                                                                                
      DO I = 1,  NORG                                                           
         IF (PARTORG(I) .LT. 0.0) THEN                                        
            PARTORG(I) = 0.0                                                    
            GASORG(I) = WORG(I)                                                 
         END IF                                                                 
      END DO                                                                    
                                                                                
      RETURN                                                                    
      END                                                                       
                                                                                
!----------------------------------------------------------------------         
      SUBROUTINE SOASUB (MSUM, WORG, PARTORG, GASORG, CURTEMP)                  
!----------------------------------------------------------------------         
                                                                                
!***********************************************************************        
!  FUNCTION: Program to calculate the equilibria                       *        
!  PRECONDITION REQUIRED: called from subr. CALTECHSOA                 *        
!  RETURN VALUES:                                                      *        
!     PARTORG(J) - Particule-phase concentrations (microgram/m**3)     *        
!     GASORG(J) - Gas-phase concentrations (microgram/m**3)            *        
!  REVISION HISTORY:                                                   *        
!     Written by Betty K. Pun of AER, Inc. for EPRI's Aerosol          *        
!         Module Implementation Project, May, 1999                     *        
!     Revised by Yang Zhang of AER, Inc. for EPRI's Aerosol            *        
!         Module Implementation Project based on Models3's             *        
!         coding standard July, 1999                                   *        
!     Modified by Betty K. Pun of AER, Inc, for Models-3 application   *        
!         Analytical solution calculated November, 99                  *        
!     Modified by Betty K. Pun, AER to include temperature dependence  *        
!         based on generic Hvap, February/March 2002                   *        
!***********************************************************************        
                                                                                
      IMPLICIT NONE                                                             
                                                                                
!...........  INCLUDES                                                          
                                                                                
      INCLUDE 'soadat.puf'                                                      
                                                                                
!...........  ARGUMENTS and their descriptions                                  
                                                                                
      REAL MSUM           ! Sum of particle-phase concentration (ug/m**3 air)   
      REAL WORG(NORG+1)   ! Total organic species concentration                 
      REAL PARTORG(NORG)  ! Particulate concentration of organic species i      
      REAL GASORG(NORG)   ! Gas-phase concentration of organic species i        
      REAL CURTEMP        ! Grid cell temperature used in modifying parititon   
                          ! constants                                           
!...........  Local variables                                                   
      INTEGER J           ! Organic species index                               
      REAL KCORR(NORG)    ! Partition coefficients, corrected for temperature   
                          ! [m**3/ug]                                           
      REAL TEXPT          ! temperature where KOM (experimental values)         
                          ! are obtained                                        
      PARAMETER (TEXPT = 298.0)                                                 
                                                                                
! ..........  Other Variables                                                   
                                                                                
!***********************************************************************        
!  begin body of subroutine                                                     
                                                                                
      IF (MSUM .EQ. 0.0) THEN                                                 
         WRITE(*,  *)    'Error in SOA: MSUM = 0.0'                             
         STOP                                                                   
      END IF                                                                    
!                                                                               
! *** Calculate the equilibrium concentration                                   
!     GASORG(J) = WORG(J) - PARTORG(J)                                          
!     G         = C      -  A                                                   
!     (A/M)/G = K is equal to A/(C-A) = KM is equal to A = CKM - AKM            
!     Therefore A = CKM / (1 + KM)                                              
                                                                                
! add code to correct for CURTEMPerature, assuming KOM is experimental          
! values obtained at 310 K (changed to a ref temp of 298K)                      
      DO J = 1, NORG                                                            
         KCORR(J) = EXP(HVAP(J)*(1/CURTEMP - 1/TEXPT)/RKJMOLK)                 &
     &                 * KOM(J) * CURTEMP/TEXPT                                 
         PARTORG(J) = WORG(J) * KCORR(J)* MSUM /                               &
     &                    (1. + MSUM * KCORR(J))                                
         GASORG(J) = WORG (J) - PARTORG(J)                                      
      END DO                                                                    
                                                                                
      RETURN                                                                    
      END                                                                       
                                                                                
                                                                                
!----------------------------------------------------------------------         
      subroutine aqradm ( temp, pres_atm, taucld, prcrate, wcavg, coz,         &
     &                    ch2o2, ctnh3, conc, rhoair, len, scav )               
!----------------------------------------------------------------------         
!                                                                               
!  DESCRIPTION:                                                                 
!    Compute concentration changes in cloud due to aqueous chemistry            
!    Adapted from RADM Cloud implementation in CMAQ/SCICHEM for CALPUFF         
!    by PK, AER, January 2007 for API Contract No. 2006-102376                  
!                                                                               
!  Reference:                                                                   
!     Walcek & Taylor, 1986, A theoretical Method for computing                 
!      vertical distributions of acidity and sulfate within cumulus             
!      clouds, J. Atmos Sci.,  Vol. 43, no. 4 pp 339 - 355                      
!                                                                               
!  Called by:  CHEMRIV6, CHEMRIV7                                               
!                                                                               
!  Calls the following functions:  HLCONST                                      
!                                                                               
!  ARGUMENTS     TYPE      I/O       DESCRIPTION                                
!  ---------     ----  ------------  --------------------------------           
!   CONC(6)      real  input&output  Concentration for species i=1,11           
!                                    (1) = SO2   conc (mol/mol of SO2)          
!                                    (2) = SO4   conc (mol/mol of SO4)          
!                                    (3) = NO    conc (mol/mol of NO)           
!                                    (4) = NO2   conc (mol/mol of NO2)          
!                                    (5) = HNO3  conc (mol/mol of HNO3)         
!                                    (6) = NO3   conc (mol/mol of NO3)          
!-----------------------------------------------------------------------        
                                                                                
      IMPLICIT NONE                                                             
                                                                                
!...........PARAMETERS and their descriptions:                                  
!      include 'params.puf'                                                    &
! tem&orarily set io6=2 here since params.puf variables are not declared,      &
! res&lting in compile errors.                                                  
      INTEGER IO6                                                               
      PARAMETER (IO6 = 2) ! from params.puf                                     
                                                                               &
! num&er of oxidizing reactions                                                 
      INTEGER NUMOX                                                             
      PARAMETER (NUMOX =  3)   ! H2O2, O3, and Fe-Mn catalyzed                  
                                                                               &
! min&mum and maximum pH                                                        
      REAL PHMIN, PHMAX                                                         
      PARAMETER (PHMIN = 0.0001, PHMAX = 10.0)                                  
                                                                               &
! min&mum concentration                                                         
      REAL CONCMIN                                                              
      PARAMETER (CONCMIN = 1.0E-30)                                             
                                                                               &
! con&ert seconds to hours                                                      
      REAL SEC2HR                                                               
      PARAMETER (SEC2HR = 1.0 / 3600.0)                                         
!                                                                              &
! Mol&r volume at STP [ L/mol ] Non MKS units                                   
      REAL MOLVOL                                                               
      PARAMETER (MOLVOL = 22.41410)                                             
                                                                               &
! Sta&dard Temperature [ K ]                                                    
      REAL STDTEMP                                                              
      PARAMETER (STDTEMP = 273.15)                                              
                                                                               &
! den&ity of water at 20 C and 1 ATM (kg/m3)                                    
      REAL H2ODENS                                                              
      PARAMETER (H2ODENS = 1000.0)                                              
                                                                               &
! Mol&cular weights of Fe and Mn                                                
      REAL MW_FE, MW_MN                                                         
      PARAMETER (MW_FE = 55.8, MW_MN = 54.9)                                    
                                                                                
! CO2 background concentration (moles/mole air)                                 
      REAL CCO2                                                                 
      PARAMETER (CCO2 = 3.4E-04)   ! 340 ppm                                    
                                                                                
!...........ARGUMENTS and their descriptions                                    
                                                                                
      REAL         TEMP            ! temperature (K)                            
      REAL         PRES_ATM        ! pressure (atm)                             
      REAL         TAUCLD          ! timestep for cloud (s)                     
      REAL         PRCRATE         ! precip rate (mm/hr)                        
      REAL         WCAVG           ! liquid water content (kg/m3)               
      REAL         COZ             ! ozone concentration (mol/mol)              
      REAL         CH2O2           ! H2O2 concentration (mol/mol)               
      REAL         CTNH3           ! total ammonium  concentration (mol/mol)    
      REAL         CONC   ( 6 )    ! species concentrations (mol/molV)          
      real         rhoair          ! air density, moles/m3                      
                                                                                
! --- puff length scale (for scavenging coefficient calculations) (m)           
      real len                                                                  
                                                                                
! --- Scavenging coefficients (1/s)                                             
      real scav( 6 )                                                            

! --- wangzhm save
                                                                                
!...........LOCAL VARIABLES (scalars) and their descriptions:                   
                                                                                
      REAL,save ::       RT              ! gas const * temperature (liter atm/mol)    
      REAL,save ::       RECIPAP1        ! one over pressure (/atm)                   
      REAL,save ::       ONE_OVER_TEMP   ! 1.0 / TEMP                                 
                                                                                
! --- Gas liquid equilibria                                                     
                                                                                
      REAL,save ::         PH2O20          ! total H2O2 partial pressure (atm)          
      REAL,save ::         PH2O2F          ! gas only H2O2 partial pressure (atm)       
      REAL,save ::         H2O2H           ! Henry's Law Constant for H2O2              
      REAL,save ::         H2O2L           ! H2O2 conc in cloudwater (mol/liter)        
                                                                                
      REAL,save ::         PHNO30          ! total HNO3 partial pressure (atm)          
      REAL,save ::         PHNO3F          ! gas only HNO3 partial pressure (atm)       
      REAL,save ::         HNO3H           ! Henry's Law Constant for HNO3              
      REAL,save ::         HNO31           ! First dissociation constant for HNO3       
      REAL,save ::         HNO31H          ! HNO31*HNO3H                                
      REAL,save ::         HNO3L           ! HNO3 conc in cloudwater (mol/liter)        
                                                                                
      REAL,save ::         PNH30           ! total NH3 partial pressure (atm)           
      REAL,save ::         PNH3F           ! gas only NH3 partial pressure (atm)        
      REAL,save ::         NH3H            ! Henry's Law Constant for NH3               
      REAL,save ::         NH31            ! First dissociation constant for NH3        
      REAL,save ::         NH3DH2O         !                                            
      REAL,save ::         NH31HDH         !                                            
      REAL,save ::         NH3L            ! NH3 conc in cloudwater (mol/liter)         
                                                                                
      REAL,save ::         PO30            ! total O3 partial pressure (atm)            
      REAL,save ::         PO3F            ! gas only O3 partial pressure (atm)         
      REAL,save ::         O3H             ! Henry's Law Constant for O3                
      REAL,save ::         O3L             ! O3 conc in cloudwater (mol/liter)          
                                                                                
      REAL,save ::         PSO20           ! total SO2 partial pressure (atm)           
      REAL,save ::         PSO2F           ! gas only SO2 partial pressure (atm)        
      REAL,save ::         SO2H            ! Henry's Law Constant for SO2               
      REAL,save ::         SO21            ! First dissociation constant for SO2        
      REAL,save ::         SO22            ! Second dissociation constant for SO2       
      REAL,save ::         SO212           ! SO21*SO22                                  
      REAL,save ::         SO212H          ! SO21*SO22*SO2H                             
      REAL,save ::         SO21H           ! SO21*SO2H                                  
      REAL,save ::         SO2L            ! SO2 conc in cloudwater (mol/liter)         
      REAL,save ::         SO3             ! SO3= conc in cloudwater (mol/liter)        
                                                                                
      REAL,save ::         PNO2            ! total NO2 partial pressure (atm)           
      REAL,save ::         NO2H            ! Henry's Law Constant for NO2 (M/atm)       
      REAL,save ::         NO2L            ! NO2 conc in cloudwater (mol/liter)         
                                                                        
      REAL,save ::         PNO             ! total NO partial pressure (atm)            
      REAL,save ::         NOH             ! Henry's Law Constant for NO (M/atm)        
      REAL,save ::         NOL             ! NO conc in cloudwater (mol/liter)          
                                                                    
      REAL,save ::         PCO20           ! total CO2 partial pressure (atm)           
      REAL,save ::         PCO2F           ! gas only CO2 partial pressure (atm)        
      REAL,save ::         CO2H            ! Henry's Law constant for CO2               
      REAL,save ::         CO21            ! First dissociation constant for CO2        
      REAL,save ::         CO22            ! Second dissociation constant for CO2       
      REAL,save ::         CO212           ! CO21*CO22                                  
      REAL,save ::         CO212H          ! CO2H*CO21*CO22                             
      REAL,save ::         CO21H           ! CO2H*CO21                                  
      REAL,save ::         CO2L            ! CO2 conc in cloudwater (mol/liter)         
                                                                           
      REAL,save ::         XL              ! conversion factor (liter-atm/mol)          
      REAL,save ::         ONE_OVER_XL     ! 1.0 / XL                                   
      REAL,save ::         PRES_ATM_OVER_XL     ! PRES_ATM / XL                         
      REAL,save ::         XLCO2           !                                            
      REAL,save ::         XLH2O2          !                                            
      REAL,save ::         XLHNO3          !                                            
      REAL,save ::         XLNH3           !                                            
      REAL,save ::         XLO3            !                                            
      REAL,save ::         XLSO2           !                                            
                                                                      
      REAL,save ::         HCO3            ! HCO3 conc in cloudwater (mol/liter)        
      REAL,save ::         HSO3            ! HSO3 conc in cloudwater (mol/liter)        
      REAL,save ::         HSO4            ! HSO4 concn in cloudwater (mol/liter)       
                                                                      
      REAL,save ::         MN              ! Mn++ conc in cloudwater (mol/liter)        
      REAL,save ::         MNA             ! initial Mn in cloudwater (mol/liter)       
      REAL,save ::         NH4             ! NH4+ conc in cloudwater (mol/liter)        
      REAL,save ::         NO3M            ! NO3- conc in cloudwater (mol/liter)        
      REAL,save ::         OHM             ! OH- conc in cloudwater (mol/liter)         
      REAL,save ::         SO4             ! SO4= conc in cloudwater (mol/liter)        
      REAL,save ::         CO3             ! CO3= conc in cloudwater (mol/liter)        
                                                                    
      REAL,save ::         A               ! iron's anion concentration                 
      REAL,save ::         B               ! manganese's anion concentration            
      REAL,save ::         FE              ! Fe+++ conc in cloudwater (mol/liter)       
      REAL,save ::         FEA             ! initial Fe in cloudwater (mol/liter)       
                                                                                
      INTEGER,save ::      I20C            ! loop counter for do loop 20                
      INTEGER,save ::      I30C            ! loop counter for do loop 30                
      INTEGER,save ::      ITERAT          ! # iterations of aqueous chemistry solver   
      INTEGER,save ::      I7777C          ! aqueous chem iteration counter             
      INTEGER,save ::      ICNTAQ          ! aqueous chem iteration counter             
      INTEGER,save ::      IOX             ! index over oxidation reactions             
                                                                                
      REAL,save ::         ACT1            ! activity correction factor!single ions     
      REAL,save ::         ACT2            ! activity factor correction!double ions     
      REAL,save ::         ACTB            !                                            
                                                                       
      REAL,save ::         FTST            !                                            
      REAL,save ::         GM1             !                                            
      REAL,save ::         GM1LOG          !                                            
      REAL,save ::         GM2             ! activity correction factor                 
      REAL,save ::         GM2LOG          !                                            
                                                                          
      REAL,save ::         FA              ! functional value ??                        
      REAL,save ::         FB              ! functional value ??                        
      REAL,save ::         AC              ! H+ concentration in cloudwater (mol/liter) 
      REAL,save ::         AE              ! guess for H+ conc in cloudwater (mol/l)    
      REAL,save ::         BB              ! lower limit guess of cloudwater pH         
      REAL,save ::         HA              !                                            
      REAL,save ::         HB              !                                            
      REAL,save ::         STION           ! ionic strength                             
                                                                           
      REAL,save ::         H2OW            !                                            
      REAL,save ::         HTST            !                                            
      REAL,save ::         RATE            !                                            
      REAL,save ::         RECIPA1         !                                            
      REAL,save ::         RECIPA2         !                                            
                                                                         
      REAL,save ::         RH2O2           !                                            
      REAL,save ::         RMHP            !                                            
      REAL,save ::         RPAA            !                                            
      REAL,save ::         SIV             ! dissolved so2 in cloudwater (mol/liter)    
      REAL,save ::         SK6             !                                            
      REAL,save ::         SK6TS6          !                                            
      REAL,save ::         DTS6            !                                            
      REAL,save ::         TAC             !                                            
      REAL,save ::         TEMP1           !                                            
      REAL,save ::         TIMEW           ! cloud chemistry clock (sec)                
      REAL,save ::         TOTOX           !                                            
      REAL,save ::         TS6             ! SO4 conc in cloudwater (mol/liter)         
      REAL,save ::         TSIV            !                                            
      REAL,save ::         TST             !                                            
                                                                        
      REAL,save ::         DSIVDT( 0:NUMOX ) ! rate of so2 oxid incloud (mol/liter/sec) 
      REAL,save ::         DS4   ( 0:NUMOX ) ! S(IV) oxidized over timestep DTW(0)      
      real,save ::         ds40old           ! total s(iv) oxidized                     
      REAL,save ::         DTW   ( 0:NUMOX ) ! cloud chemistry timestep (sec)           
                                                                     
      real,save ::         depfactor                                                    
                                                                                
!                                                                               
! --- Define background concs of iron and manganese (ug/m3)
      real cfe, cmn                             
      data cfe/0.01/, cmn/0.005/                                                
                                                                                
! ... Specify fraction of activation for particles.                            
      real fracma                                                                        
      data fracma/1.0/                                                         
                                                                                
!...........EXTERNAL FUNCTIONS and their descriptions:                          
                                                                                
      REAL HLCONST                                                              
      EXTERNAL HLCONST                                                          
                                                                                
!*********************************************************************          
!     begin body of subroutine AQRADM                                  


! --- wangzhm omp: critical
         !---$---OMP CRITICAL (call_aqradm) 
! --- wangzhm omp: set threadsprivate                                           
      !$OMP THREADPRIVATE(RT,RECIPAP1,ONE_OVER_TEMP,PH2O20,PH2O2F,H2O2H,H2O2L,&
      !$OMP& PHNO30,PHNO3F,HNO3H,HNO31,HNO31H,HNO3L,&
      !$OMP& PNH30,PNH3F,NH3H,NH31,NH3DH2O,NH31HDH,NH3L,&
      !$OMP& PO30,PO3F,O3H,O3L,PSO20,PSO2F,SO2H,SO21,SO22,SO212,SO212H,SO21H,SO2L,SO3,&
      !$OMP& PNO2,NO2H,NO2L,PNO,NOH,NOL,PCO20,PCO2F,CO2H,CO21,CO22,CO212,CO212H,CO21H,CO2L,&
      !$OMP& XL,ONE_OVER_XL,PRES_ATM_OVER_XL,XLCO2,XLH2O2,XLHNO3,XLNH3,XLO3,XLSO2,&
      !$OMP& HCO3,HSO3,HSO4,MN,MNA,NH4,NO3M,OHM,SO4,CO3,A,B,FE,FEA,&
      !$OMP& I20C,I30C,ITERAT,I7777C,ICNTAQ,IOX,ACT1,ACT2,ACTB,FTST,GM1,GM1LOG,GM2,GM2LOG,&
      !$OMP& FA,FB,AC,AE,BB,HA,HB,STION,H2OW,HTST,RATE,RECIPA1,RECIPA2,&
      !$OMP& RH2O2,RMHP,RPAA,SIV,SK6,SK6TS6,DTS6,TAC,TEMP1,TIMEW,TOTOX,TS6,TSIV,TST,&
      !$OMP& DSIVDT,DS4,ds40old,DTW,depfactor) 
                               
!...Check for bad temperature or pressure                                       
      if ( temp .LE. 0.0 .or. pres_atm .LE. 0.0 ) then                          
         write(io6,*)'Invalid temp and/or pressure; T,P: ',temp,pres_atm        
         stop 'Halted in AQRADM -- see list file.'                              
      end if                                                                    
                                                                                
      one_over_temp = 1.0 / temp                                                
                                                                                
!...Compute several conversion factors                                          
                                                                                
      icntaq = 0                                                                
      iterat = 0                                                                
      rt = ( MOLVOL / STDTEMP ) * temp    ! r * t (liter atm / mol)             
      xl   = wcavg * rt / H2ODENS         ! conversion factor (l-atm/mol)       
      one_over_xl = 1.0 / xl                                                    
      pres_atm_over_xl = pres_atm / xl                                          
      tst  = 0.999                                                              
      act1 = 1.0                                                                
      act2 = 1.0                                                                
      gm2  = 1.0                                                                
      timew = 0.0                                                               
      recipap1 = 1.0 / pres_atm                                                 
                                                                                
!...set equilibrium constants as a function of temperature                      
!...Henry's law constants                                                       
                                                                                
      so2h  = HLCONST( 'SO2', temp, .false., 0.0 )                              
      co2h  = HLCONST( 'CO2', temp, .false., 0.0 )                              
      nh3h  = HLCONST( 'NH3', temp, .false., 0.0 )                              
      h2o2h = HLCONST( 'H2O2', temp, .false., 0.0 )                             
      o3h   = HLCONST( 'O3', temp, .false., 0.0 )                               
      hno3h = HLCONST( 'HNO3', temp, .false., 0.0 )                             
      noh   = HLCONST( 'NO', temp, .false., 0.0 )                               
      no2h  = HLCONST( 'NO2', temp, .false., 0.0 )                              
                                                                                
!...Dissociation constants                                                      
                                                                                
      temp1 = one_over_temp - 1.0 / 298.0                                       
                                                                                
      sk6   = 1.02e-02 * EXP(  2.72e+03 * temp1 )    ! Smith and Martell (1976) 
      so21  = 1.30e-02 * EXP(  1.96e+03 * temp1 )    ! Smith and Martell (1976) 
      so22  = 6.60e-08 * EXP(  1.50e+03 * temp1 )    ! Smith and Martell (1976) 
      co21  = 4.30e-07 * EXP( -1.00e+03 * temp1 )    ! Smith and Martell (1976) 
      co22  = 4.68e-11 * EXP( -1.76e+03 * temp1 )    ! Smith and Martell (1976) 
      h2ow  = 1.00e-14 * EXP( -6.71e+03 * temp1 )    ! Smith and Martell (1976) 
      nh31  = 1.70e-05 * EXP( -4.50e+02 * temp1 )    ! Smith and Martell (1976) 
      hno31 = 1.54e+01 * EXP(  8.70e+03 * temp1 )    ! Schwartz (1984)          
                                                                                
!...Kinetic oxidation rates                                                     
!...   From Chamedies (1982)                                                    
                                                                                
      rh2o2 = 8.0e+04 * EXP( -3650.0 * temp1 )                                  
                                                                                
!...Make initializations                                                        
                                                                                
      do iox = 0, NUMOX                                                         
         dsivdt(iox) = 0.0                                                      
         dtw(iox)    = 0.0                                                      
         ds4(iox)    = 0.0                                                      
      end do                                                                    
                                                                                
      fea = cfe * 1.E-6 * pres_atm_over_xl / (rhoair * MW_FE)                   
      mna = cmn * 1.E-6 * pres_atm_over_xl / (rhoair * MW_MN)                   
                                                                                
!...Set constant factors that will be used in later multiplications (moles/atm) 
                                                                                
      xlh2o2  = h2o2h * xl                                                      
      xlo3    = o3h   * xl                                                      
      xlso2   = so2h  * xl                                                      
      xlnh3   = nh3h  * xl                                                      
      xlhno3  = hno3h * xl                                                      
      xlco2   = co2h  * xl                                                      
                                                                                
      so212   = so21  * so22                                                    
      so21h   = so21  * so2h                                                    
      so212h  = so212 * so2h                                                    
      co212   = co21  * co22                                                    
      co21h   = co21  * co2h                                                    
      co212h  = co22  * co21h                                                   
      nh3dh2o = nh31  / h2ow                                                    
      nh31hdh = nh3h  * nh3dh2o                                                 
      hno31h  = hno31 * hno3h                                                   
                                                                                
!...If kinetic calculations are made, return to this point                      
                                                                                
      i20c = 0                                                                  
20    continue                                                                  
                                                                                
      i20c = i20c + 1                                                           
      if ( i20c >= 1000 ) then                                                  
         write(io6,*) 'Excessive looping at I20C'                               
         stop 'Halted in AQRADM -- see list file.'                              
      end if                                                                    
                                                                                
!...Initial gas phase partial pressures (atm)                                   
      pnh30  = ctnh3 * pres_atm                                                 
      phno30 = ( conc(5)  + conc(6)*fracma ) * pres_atm                         
      pco20  = cco2 * pres_atm                                                  
                                                                                
! --- Reactive species                                                          
                                                                                
! --- H2O2                                                                      
      ph2o20 = ch2o2 * pres_atm + xl * ds4( 1 )                                 
! --- check if too much h2o2 has reacted (possible in plume with                
! --- high SO2 concs)                                                           
      if ( ph2o20 .LT. 0. ) then                                                
         ds4( 0 ) = ds4( 0 ) - ds4( 1 )                                         
         ds4( 1 ) = -ch2o2 * pres_atm * one_over_xl                             
         ds4( 0 ) = ds4( 0 ) + ds4( 1 )                                         
         ph2o20 = 0.                                                            
      end if                                                                    
                                                                                
! --- O3                                                                        
      po30   = coz * pres_atm + xl * ds4( 2 )                                   
! --- check if too much o3 has reacted                                          
      if ( po30 .LT. 0. ) then                                                  
         ds4( 0 ) = ds4( 0 ) - ds4( 2 )                                         
         ds4( 2 ) = -coz * pres_atm * one_over_xl                               
         ds4( 0 ) = ds4( 0 ) + ds4( 2 )                                         
         po30 = 0.                                                              
      end if                                                                    
                                                                                
! --- SO2                                                                       
      pso20  = conc( 1  ) * pres_atm + ds4( 0 ) * xl                            
! --- check if too much SO2 has reacted                                         
      if ( pso20 .LT. 0. ) then                                                 
         ds40old = ds4( 0 )                                                     
         ds4( 0 ) = -conc( 1 ) * pres_atm * one_over_xl                         
         do iox = 1, NUMOX                                                      
            ds4( iox ) = ds4( iox ) * ds4( 0 ) / ds40old                        
         end do                                                                 
         pso20 = 0.                                                             
         ph2o20 = ch2o2 * pres_atm + xl * ds4( 1 )                              
         po30   = coz * pres_atm + xl * ds4( 2 )                                
      end if                                                                    
                                                                                
!...Don't allow gas concentrations to go below zero                             
                                                                                
!     pso20  = MAX( pso20,  0.0 )                                               
!     ph2o20 = MAX( ph2o20, 0.0 )                                               
!     po30   = MAX( po30,   0.0 )                                               
      pnh30  = MAX( pnh30,  0.0 )                                               
      pco20  = MAX( pco20,  0.0 )                                               
      phno30 = MAX( phno30, 0.0 )                                               
                                                                                
!...Molar concentrations of soluble aerosols                                    
      ts6     = conc( 2 ) * fracma * pres_atm_over_xl                          &
     &        - ds4( 0 )   ! Sulfate                                            
                                                                                
      fe      = fea                                                             
      mn      = mna                                                             
      a       = 3.0 * fe                                                        
      b       = 2.0 * mn                                                        
                                                                                
!...Don't allow aerosol concentrations to go below zero                         
                                                                                
      ts6     = MAX( ts6,     0.0 )                                             
                                                                                
      sk6ts6 = sk6 * ts6                                                        
                                                                                
!...Find solution of the equation using a method of reiterative                 
!...bisections. Make initial guesses for pH between PHMIN  to  PHMAX.           
                                                                                
      ha = PHMIN                                                                
      hb = PHMAX                                                                
                                                                                
      i7777c = 0                                                                
7777  continue                                                                  
                                                                                
      i7777c = i7777c + 1                                                       
      if ( i7777c .GT. 1000 ) then                                              
         write(io6,*)'Excessive looping at I7777C'                              
         stop 'Halted in AQRADM -- see list file.'                              
      end if                                                                    
                                                                                
!     ha = MAX( ha - 0.8, 0.1 )                                                 
!     hb = MIN( hb + 0.8, 9.9 )                                                 
      ha = MAX( ha - 0.8, PHMIN )                                               
      hb = MIN( hb + 0.8, PHMAX )                                               
      ae = 10.0**( -ha )                                                        
                                                                                
      recipa1 = 1.0 / ( ae * act1 )                                             
      recipa2 = 1.0 / ( ae * ae * act2 )                                        
                                                                                
!...Calculate final gas phase partial pressure of SO2, NH3, HNO3 and            
!...CO2 (atm)                                                                   
                                                                                
      pso2f = pso20 / ( 1.0 + xlso2 * ( 1.0 + so21 * recipa1                   &
     &      + so212 * recipa2 ) )                                               
                                                                                
      pnh3f = pnh30 / ( 1.0 + xlnh3 * ( 1.0 + nh3dh2o * ae ) )                  
                                                                                
      phno3f = phno30 / ( 1.0 + xlhno3 * ( 1.0 + hno31 * recipa1 ) )            
                                                                                
      pco2f = pco20 / ( 1.0 + xlco2 * ( 1.0 + co21 * recipa1                   &
     &      + co212 * recipa2 ) )                                               
                                                                                
!...Calculate liquid phase concentrations (moles/liter)                         
                                                                                
      so4  = sk6ts6 / ( ae * gm2 + sk6 )                                        
      hso4 = ts6 - so4                                                          
      so3  = so212h  * pso2f  * recipa2                                         
      hso3 = so21h   * pso2f  * recipa1                                         
      co3  = co212h  * pco2f  * recipa2                                         
      hco3 = co21h   * pco2f  * recipa1                                         
      ohm  = h2ow    * recipa1                                                  
      nh4  = nh31hdh * pnh3f  * ae                                              
      no3m = hno31h  * phno3f * recipa1                                         
                                                                                
!...Compute functional value                                                    
                                                                                
      fa = ae + nh4 - 2.0 *  (co3 + so3 + so4 ) - ohm - hco3                   &
     &   - hso3 - no3m - hso4                                                   
                                                                                
!...Start iteration and bisection ****************<<<<<<<                       
                                                                                
      i30c = 0                                                                  
30    continue                                                                  
                                                                                
      i30c = i30c + 1                                                           
      if ( i30c .GT. 1000 ) then                                                
         write(io6,*)'Excessive looping at I30C'                                
         stop 'Halted in AQRADM -- see list file.'                              
      end if                                                                    
                                                                                
      bb = 0.5 * ( ha + hb )                                                    
      ae = 10.0**( -bb )                                                        
                                                                                
! --- don't solve for H+ if fa < 0 at first try                                 
      if ( i7777c .EQ. 1 .and. fa .LT. 0. ) then                                
                                                                                
         bb = ha                                                                
         hb = ha                                                                
         ae = 10.0**( -bb )                                                     
                                                                                
      end if                                                                    
                                                                                
      recipa1 = 1.0 / ( ae * act1 )                                             
      recipa2 = 1.0 / ( ae * ae * act2 )                                        
                                                                                
!...Calculate final gas phase partial pressure of SO2, NH3, HNO3 and            
!...CO2 (atm)                                                                   
                                                                                
      pso2f = pso20 / ( 1.0 + xlso2                                            &
     &	    * ( 1.0 + so21 * recipa1 + so212 * recipa2 ) )                       
                                                                                
      pnh3f = pnh30 / ( 1.0 + xlnh3 * ( 1.0 + nh3dh2o * ae ) )                  
                                                                                
      phno3f = phno30 / ( 1.0 + xlhno3 * ( 1.0 + hno31 * recipa1 ) )            
                                                                                
      pco2f = pco20 / ( 1.0 + xlco2 * ( 1.0 + co21 * recipa1                   &
     &      + co212 * recipa2 ) )                                               
                                                                                
!...Calculate liquid phase concentrations (moles/liter)                         
                                                                                
      so4  = sk6ts6 / ( ae * gm2 + sk6 )                                        
      hso4 = ts6 - so4                                                          
      so3  = so212h  * pso2f  * recipa2                                         
      hso3 = so21h   * pso2f  * recipa1                                         
      co3  = co212h  * pco2f  * recipa2                                         
      hco3 = co21h   * pco2f  * recipa1                                         
      ohm  = h2ow    * recipa1                                                  
      nh4  = nh31hdh * pnh3f  * ae                                              
      no3m = hno31h  * phno3f * recipa1                                         
                                                                                
!...compute functional value                                                    
                                                                                
      fb = ae + nh4 - 2.0 * ( co3 + so3 + so4 ) - ohm - hco3                   &
     &   - hso3 - no3m - hso4                                                   
                                                                                
!...Calculate and check the sign of the product of the two functional values    
                                                                                
      ftst = fa * fb                                                            
      if ( ftst .LE. 0.0 ) then                                                 
        hb = bb                                                                 
      else                                                                      
        ha = bb                                                                 
        fa = fb                                                                 
      end if                                                                    
                                                                                
!...Check convergence of solutions                                              
                                                                                
      htst = ha / hb                                                            
      if ( htst .LE. tst ) go to 30                                             
                                                                                
!...end of zero-finding routine ****************<<<<<<<<<<<<                    
                                                                                
!...compute Ionic strength and activity coefficient by the Davies equation      
                                                                                
      stion = 0.5 * (ae + nh4 + ohm + hco3 + hso3                              &
     &      + 4.0 * (so4 + co3 + so3)                                          &
     &      + no3m + hso4 + 9.0 * fe + a + b)                                   
      gm1log = -0.509 * ( SQRT( stion )                                        &
     &       / ( 1.0 + SQRT( stion ) ) - 0.2 * stion )                          
      gm2log = gm1log * 4.0                                                     
      gm1  = 10.0**gm1log                                                       
      gm2  = MAX( 10.0**gm2log, 1.0e-30 )                                       
      actb = act1                                                               
      act1 = MAX( gm1 * gm1, 1.0e-30 )                                          
      act2 = MAX( gm1 * gm1 * gm2, 1.0e-30 )                                    
                                                                                
!...check for convergence and possibly go to 7777, to recompute                 
!...  Gas and liquid phase concentrations                                       
                                                                                
! --- don't solve for H+ if fa < 0 at first try                                 
      if ( i7777c .EQ. 1 .and. fa .LT. 0. ) then                                
         actb = act1                                                            
      end if                                                                    
                                                                                
      tac = ABS( actb - act1 ) / actb                                           
      if ( tac >= 1.0e-2 ) then                                                 
                                                                                
         icntaq = icntaq + 1                                                    
         if ( icntaq .GT. 100 ) then                                            
            write(io6,*)'Maximum iterations for pH calculation exceeded'        
            write(io6,*)'Using last pH value'                                   
         else                                                                   
           go to 7777                                                           
         end if                                                                 
                                                                                
      end if                                                                    
                                                                                
!...return an error if the pH is not in range                                   
                                                                                
!cc      if ( ( ha .lt. 0.02 ) .or. ( ha .gt. 9.49 ) ) then                     
      if ( ( ha .LT. PHMIN ) .or. ( ha .GT. PHMAX ) ) then                      
         write(io6,*)'pH value out of range: ',ha                               
         stop 'Halted in AQRADM -- see list file.'                              
      end if                                                                    
                                                                                
!...Make those concentration calculations which can be made outside             
!...  of the function.                                                          
                                                                                
      so2l = so2h * pso2f                                                       
      ac = 10.0**( -bb )                                                        
      siv = so3 + hso3 + so2l                                                   
                                                                                
!...Calculate final gas phase concentrations of oxidants (atm)                  
                                                                                
      ph2o2f = ph2o20 / ( 1.0 + xlh2o2 )                                        
      po3f   = po30   / ( 1.0 + xlo3   )                                        
                                                                                
      ph2o2f = MAX( ph2o2f, 0.0 )                                               
      po3f   = MAX( po3f,   0.0 )                                               
                                                                                
!...Calculate liquid phase concentrations (moles/liter)                         
                                                                                
      h2o2l = ph2o2f * h2o2h                                                    
      o3l   = po3f   * o3h                                                      
      nh3l  = pnh3f  * nh3h                                                     
      co2l  = pco2f  * co2h                                                     
      hno3l = phno3f * hno3h                                                    
                                                                                
!...if the maximum cloud lifetime has not been reached, then compute            
!...the next timestep.                                                          
                                                                                
      if ( timew .LT. taucld ) then                                             
                                                                                
!...make kinetics calculations                                                  
!...  note: DS4(i) and DSIV(I) are negative numbers!                            
                                                                                
        iterat = iterat + 1                                                     
                                                                                
!...Define the total S(iv) available for oxidation                              
                                                                                
        tsiv = pso20 * one_over_xl                                              
                                                                                
!...Calculate sulfur iv oxidation rate due to H2O2                              
                                                                                
        dsivdt( 1 ) = -rh2o2 * h2o2l * so2l / ( 0.1 + ac )                      
        totox = ph2o20 * one_over_xl                                            
        if ( ( dsivdt( 1 ) .EQ. 0.0 ) .or.                                     &
     &       ( tsiv  .LE. CONCMIN ) .or.                                       &
     &       ( totox .LE. CONCMIN ) ) then                                      
          dtw(1) = taucld                                                       
        else                                                                    
          dtw( 1 ) = -0.05 * MIN( totox, tsiv ) / dsivdt( 1 )                   
        end if                                                                  
                                                                                
!...Calculate sulfur iv oxidation rate due to O3                                
                                                                                
        if ( bb .GE. 2.7 ) then                                                 
          dsivdt( 2 ) = -4.19e5 * ( 1.0 + 2.39e-4 / ac ) * o3l * siv            
        else                                                                    
          dsivdt( 2 ) = -1.9e4 * siv * o3l / SQRT( ac )                         
        end if                                                                  
        totox = po30 * one_over_xl                                              
        if ( ( dsivdt( 2 ) .EQ. 0.0 ) .or.                                     &
     &       ( tsiv  .LE. CONCMIN ) .or.                                       &
     &       ( totox .LE. CONCMIN ) ) then                                      
          dtw( 2 ) = taucld                                                     
        else                                                                    
          dtw( 2 ) = -0.01 * MIN( totox, tsiv ) / dsivdt( 2 )                   
        end if                                                                  
                                                                                
!...Calculate sulfur iv oxidation rate due to O2 catalyzed by Mn++              
!...  and Fe+++  See Table IV Walcek & Taylor ( 1986)                           
                                                                                
        if ( bb .GE. 4.0 )  then  ! 4.0  < ph                                   
	                                                                               
          if ( siv .LE. 1.0e-5 ) then                                           
            dsivdt( 3 ) = -5000.0 * mn * hso3                                   
          else                                                                  
            dsivdt( 3 ) = -( 4.7 * mn * mn / ac                                &
     &                  + 1.0e7 * fe * siv * siv )                              
          end if  ! end of first pass through siv conc.                         
                                                                                
        else          ! ph < 4.0                                                
                                                                                
          if ( siv .LE. 1.0e-5 ) then                                           
            dsivdt( 3 ) = -3.0 * ( 5000.0 * mn * hso3                          &
     &                  + 0.82 * fe * siv / ac )                                
          else                                                                  
            dsivdt( 3 ) = -( 4.7 * mn * mn / ac                                &
     &                  + ( 0.82 * fe * siv / ac )                             &
     &                  * ( 1.0 + 1.7e3 * mn**1.5 / ( 6.3e-6 + fe ) ) )         
          end if ! end of second pass through siv conc.                         
        end if  ! end of pass through ph                                        
                                                                                
        if ( ( dsivdt( 3 ) .EQ. 0.0 ) .or. ( tsiv .LE. CONCMIN ) ) then         
          dtw( 3 ) = taucld                                                     
        else                                                                    
          dtw( 3 ) = -0.1 * tsiv / dsivdt( 3 )                                  
        end if                                                                  
                                                                                
!...Calculate total sulfur iv oxidation rate                                    
                                                                                
        dsivdt( 0 ) = 0.0                                                       
        do iox = 1, NUMOX                                                       
          dsivdt( 0 ) = dsivdt( 0 ) + dsivdt( iox )                             
        end do                                                                  
                                                                                
!...Calculate a minimum time step required                                      
                                                                                
        dtw( 0 ) = MIN( dtw( 1 ), dtw( 2 ), dtw( 3 ) )                          
                                                                                
!...check for large time step                                                   
                                                                                
        if ( dtw( 0 ) .GT. 8.0e+37 ) then                                       
          write(io6,1001) dsivdt(0), ts6, dtw(0)                                
        else                                                                    
                                                                                
!...calculate the change in sulfur iv for this time step                        
                                                                                
60        continue                                                              
          dts6 = ABS( dtw( 0 ) * ( -dsivdt( 0 ) ) )                             
                                                                                
!...If DSIV(0), sulfur iv oxidized during this time step would be               
!...less than 5% of sulfur oxidized since time 0, then double DT                
                                                                                
          if ( dtw( 0 ) .LE. taucld ) then                                      
            if ( dts6 .LT. 0.05 * ts6 ) then                                    
              dtw( 0 ) = dtw( 0 ) * 2.0                                         
	      go to 60                                                                 
            end if                                                              
          end if                                                                
        end if                                                                  
        dtw( 0 ) = MIN( dtw( 0 ), taucld )                                      
                                                                                
!...If the total time after this time increment will be greater than            
!...  TAUCLD sec., then set DTW(0) so that total time will be TAUCLD            
                                                                                
        if ( timew + dtw( 0 ) .GT. taucld ) dtw( 0 ) = taucld - timew           
        if ( ts6 .LT. 1.0e-11 ) dtw( 0 ) = taucld - timew                       
        if ( iterat .GT. 100 ) dtw( 0 ) = taucld - timew                        
                                                                                
!...Set DSIV(I), I = 0,NUMOX, the amount of S(IV) oxidized by each              
!... individual oxidizing agent, as well as the total.                          
                                                                                
        do iox = 0, NUMOX                                                       
           ds4( iox ) = ds4( iox ) + dtw( 0 ) * dsivdt( iox )                   
        end do                                                                  
                                                                                
        timew = timew + dtw( 0 )                                                
                                                                                
!...Return to make additional calculations                                      
                                                                                
        go to 20                                                                
      end if                                                                    
                                                                                
! --- Calculate liquid-phase concentrations of other species                    
      pno    = conc(3) * pres_atm                                               
      pno2   = conc(4) * pres_atm                                               
                                                                                
      nol    = pno   * noh   / ( 1.0 + noh    * xl )                            
      no2l   = pno2  * no2h  / ( 1.0 + no2h   * xl )                            
                                                                                
!...Compute the output concentrations                                           
                                                                                
!...gas concentrations (mol/molV) (only for reactive species)                   
                                                                                
      conc(1) = (pso2f  + xl * siv)   * recipap1                                
      ch2o2   = (ph2o2f + xl * h2o2l) * recipap1                                
      coz     = (po3f   + xl * o3l)   * recipap1                                
                                                                                
! --- calculate scavenging coefficients for gases                               
      depfactor = prcrate * SEC2HR / ( rhoair * len )                           
      if ( conc( 1 ) > CONCMIN ) then                                           
         scav( 1 ) = siv * depfactor / conc( 1 )                                
      end if                                                                    
      if ( conc( 3 ) > CONCMIN ) then                                           
         scav( 3 ) = nol * depfactor / conc( 3 )                                
      end if                                                                    
      if ( conc( 4 ) > CONCMIN ) then                                           
         scav( 4 ) = no2l * depfactor / conc( 4 )                               
      end if                                                                    
      if ( conc( 5 ) > CONCMIN ) then                                           
         scav( 5 ) = hno3l * depfactor / conc( 5 )                              
      end if                                                                    
                                                                                
!...aerosol concentrations (mol/molV) (only reactive species-so4)               
                                                                                
      conc( 2 ) = conc( 2 ) * ( 1.0 - fracma ) + ts6 * xl * recipap1            
                                                                                
! --- calculate scavenging coefficients for aerosols                            
                                                                                
      if ( conc( 2 ) > CONCMIN ) then                                           
         scav( 2 ) = ts6 * depfactor / conc( 2 )                                
      end if                                                                    
      if ( conc( 6 ) > CONCMIN ) then                                           
         scav( 6 ) = no3m * depfactor / conc( 6 )                               
      end if                                                                    
! --- wangzhm omp: critical
         !---$---OMP END CRITICAL (call_aqradm)                                                                                 
      return                                                                    
                                                                                
!...formats                                                                     
                                                                                
1001  format (1X,'DSIVDT(0) =', F10.5,                                         &
     &       'TS6=', F10.5, 'DTW(0)=', F10.5)                                   


 
                                                                                
      end                                                                       
                                                                                
!----------------------------------------------------------------------         
      REAL FUNCTION HLCONST ( NAME, TEMP, EFFECTIVE, HPLUS )                    
!----------------------------------------------------------------------         
                                                                                
!                                                                               
!  FUNCTION: return the Henry's law constant for the specified substance        
!            at the given temperature                                           
!  Adapted for CALPUFF from CMAQ version, PK, AER, Feb 2007                     
                                                                                
      IMPLICIT NONE                                                             
                                                                                
!...........PARAMETERS and their descriptions:                                  
!      include 'params.puf'                                                   
! tem&orarily set io6=2 here since params.puf variables are not declared,     
! res&lting in compile errors.                                                
      INTEGER IO6                                                            
      PARAMETER (IO6 = 2) ! from params.puf                                  
                                                                              
! Num&er of species                                                           
      INTEGER MXSPCS                                                         
      PARAMETER (MXSPCS = 8)                                                 
! Num&er of dissociating species                                              
      INTEGER MXDSPCS                                                           
      PARAMETER (MXDSPCS = 8)                                                   
                                                                                
!...........ARGUMENTS and their descriptions                                    
                                                                                
      CHARACTER*(*) NAME                ! name of substance                     
      REAL          TEMP                ! temperature (K)                       
      LOGICAL       EFFECTIVE           ! true=compute the effective henry's law
      REAL          HPLUS               ! hydrogen ion concentration (mol/l)    
                                                                                
!...........SCRATCH LOCAL VARIABLES and their descriptions:                     
      CHARACTER*16 SUBNAME(MXSPCS)                                              
      SAVE SUBNAME                                                              
! --- wangzhm save                                                                                
      INTEGER,save ::       SPC                 ! species index                         
      INTEGER,save ::       LSO2                ! SO2 pointer                           
      INTEGER,save ::       LHSO3               ! HSO3 pointer                          
      INTEGER,save ::       LHNO3               ! HNO3 pointer                          
      INTEGER,save ::       LCO2                ! CO2 pointer                           
      INTEGER,save ::       LHCO3               ! HCO3 pointer                          
      INTEGER,save ::       LH2O2               ! H2O2 pointer                          
      INTEGER,save ::       LHO2                ! HO2 pointer                           
      INTEGER,save ::       LNH4OH              ! NH4OH pointer                         
      INTEGER,save ::       LH2O                ! H2O pointer                           
                                                                                
      REAL,save ::          HPLUSI              ! 1 / HPLUS                             
      REAL,save ::          HPLUS2I             ! 1 / HPLUS**2                          
      REAL,save ::          TFAC                ! (298-T)/(T*298)                       
      REAL,save ::          AKEQ1               ! temp var for dissociation constant    
      REAL,save ::          AKEQ2               ! temp var for dissociation constant    
      REAL,save ::          OHION               ! OH ion concentration                  
      REAL,save ::          KH                  ! temp var for henry's law constant     
                                                                               &
! Hen&y's law constants at 298.15K (M/atm) (taken from Rolf Sanders'           &
! Com&ilation of Henry's Law Constants for Inorganic and Organic Species       &
! of &otential Importance in Environment Chemistry, 1999)                       
      REAL A(MXSPCS)                                                            
      SAVE A                                                                    
                                                                               &
! Ent&alpy (like activation energy) (K) (taken from Rolf Sanders'              &
! Com&ilation of Henry's Law Constants for Inorganic and Organic Species       &
! of &otential Importance in Environment Chemistry, 1999)                       
      REAL E(MXSPCS)                                                            
      SAVE E                                                                    
                                                                               &
! Dis&ociation constants at 298.15K (M or M2) (taken from Table 6.A.1,         &
! Sei&feld and Pandis, Atmospheric Chemistry and Physics, 1997)                 
      REAL B(MXDSPCS)                                                           
      SAVE B                                                                    
                                                                               &
! -dH&R (K) (taken from Table 6.A.1,                                           &
! Sei&feld and Pandis, Atmospheric Chemistry and Physics, 1997)                 
      REAL D(MXDSPCS)                                                           
      SAVE D                                                                    
                                                                                
      DATA SUBNAME(1), A(1), E(1) / 'O3', 1.2E-02, 2.7E+03 /  ! Chameides 1984  
      DATA SUBNAME(2), A(2), E(2) / 'H2O2', 8.3E+04, 7.4E+03 /  ! O'Sullivan et 
      DATA SUBNAME(3), A(3), E(3) / 'NH3', 6.1E+01, 4.2E+03 /  ! Clegg and Brimb
      DATA SUBNAME(4), A(4), E(4) / 'NO', 1.9E-03, 1.4E+03 /  ! Lide and Frederi
      DATA SUBNAME(5), A(5), E(5) / 'NO2', 1.2E-02, 2.5E+03 /  ! Chameides 1984 
      DATA SUBNAME(6), A(6), E(6) / 'HNO3', 2.1E+05, 8.7E+03 /  ! Leieveld and C
      DATA SUBNAME(7), A(7), E(7) / 'SO2', 1.4E+00, 2.9E+03 /  ! Linde and Frede
      DATA SUBNAME(8), A(8), E(8) / 'CO2', 3.6E-02, 2.2E+03 /  ! Zheng et al. 19
                                                                                
      DATA LSO2,  B(1), D(1) /  1, 1.30E-02,  1.96E+03 /  ! SO2*H2O<=>HSO3+H    
      DATA LHSO3, B(2), D(2) /  2, 6.60E-08,  1.50E+03 /  ! HSO3<=>SO3+H        
      DATA LHNO3, B(3), D(3) /  3, 1.54E+01,  8.70E+03 /  ! HNO3(aq)<=>NO3+H    
      DATA LCO2,  B(4), D(4) /  4, 4.30E-07, -1.00E+03 /  ! CO2*H2O<=>HCO3+H    
      DATA LHCO3, B(5), D(5) /  5, 4.68E-11, -1.76E+03 /  ! HCO3<=>CO3+H        
      DATA LH2O2, B( 6), D(6) /  6, 2.20E-12, -3.73E+03 /  ! H2O2(aq)<=>HO2+H   
      DATA LNH4OH, B(7), D(7) /  7, 1.70E-05, -4.50E+02 /  ! NH4*OH<=>NH4+OH    
      DATA LH2O,   B(8), D(8) /  8, 1.00E-14, -6.71E+03 /  ! H2O<=>H+OH         
                                                                                
!...........EXTERNAL FUNCTIONS and their descriptions:                          
                                                                               &
! Fun&tion to look up name in table                                             
      INTEGER INDEX1                                                            
      EXTERNAL INDEX1                                                           

! --- wangzhm omp: set threadsprivate                                           
      !$OMP THREADPRIVATE(SUBNAME,SPC,LSO2,LHSO3,LHNO3,LCO2,LHCO3,LH2O2,LHO2,LNH4OH,LH2O,&
      !$OMP& HPLUSI,HPLUS2I,TFAC,AKEQ1,AKEQ2,OHION,KH,A,E,B,D) 
                                                                                
!-----------------------------------------------------------------------        
!  begin body of subroutine HLCONST                                             
                                                                                
      SPC = INDEX1( NAME, MXSPCS, SUBNAME )                                     
                                                                                
!...error if species not found in table                                         
                                                                                
      IF ( SPC <= 0 ) THEN                                                      
         write(io6,*)TRIM(NAME) //                                             &
     &         ' not found in Henrys Law Constant table'                        
         stop 'Halted in HLCONST -- see list file.'                             
      END IF                                                                    
                                                                                
!...compute the Henry's Law Constant                                            
      TFAC = (298.0 - TEMP) / (298.0 * TEMP)                                    
      KH = A(SPC) * EXP(E(SPC) * TFAC)                                          
      HLCONST = KH                                                              
                                                                                
!...compute the effective Henry's law constants                                 
                                                                                
      IF (EFFECTIVE) THEN                                                       
                                                                                
         IF ( HPLUS <= 0.0 ) THEN                                               
            write(io6,*)'Negative or Zero [H+] concentration specified '        
            stop 'Halted in HLCONST -- see list file.'                          
         END IF                                                                 
                                                                                
         HPLUSI = 1.0 / HPLUS                                                   
         HPLUS2I = HPLUSI * HPLUSI                                              
                                                                                
         IF (TRIM(NAME) .EQ. 'SO2') THEN                                        
                                                                                
            AKEQ1 = B(LSO2) * EXP(D(LSO2) * TFAC)     !SO2H2O <=> HSO3- + H+    
            AKEQ2 = B(LHSO3) * EXP(D(LHSO3) * TFAC)   !HSO3- <=> SO3= + H+      
            HLCONST = KH * (1.0 + AKEQ1*HPLUSI + AKEQ1*AKEQ2*HPLUS2I)           
                                                                                
         ELSE IF (TRIM(NAME) .EQ. 'HNO3') THEN                                  
                                                                                
            AKEQ1 = B(LHNO3) * EXP(D(LHNO3) * TFAC)   !HNO3(aq) <=> NO3- + H+   
            HLCONST = KH * (1.0 + AKEQ1*HPLUSI)                                 
                                                                                
         ELSE IF (TRIM(NAME) .EQ. 'CO2') THEN                                   
                                                                                
            AKEQ1 = B(LCO2) * EXP(D(LCO2)  * TFAC)    !CO2H2O <=> HCO3- + H+    
            AKEQ2 = B(LHCO3) * EXP(D(LHCO3) * TFAC)   !HCO3- <=> CO3= + H+      
            HLCONST = KH * (1.0 + AKEQ1*HPLUSI + AKEQ1*AKEQ2*HPLUS2I)           
                                                                                
         ELSE IF (TRIM(NAME) .EQ. 'H2O2') THEN                                  
                                                                                
            AKEQ1 = B(LH2O2) * EXP(D(LH2O2) * TFAC)   !H2O2(aq) <=> HO2- + H+   
            HLCONST = KH * (1.0 + AKEQ1*HPLUSI)                                 
                                                                                
         ELSE IF (TRIM(NAME) .EQ. 'NH3') THEN                                   
                                                                                
            AKEQ1 = B(LNH4OH) * EXP(D(LNH4OH) * TFAC) !NH4OH <=> NH4+ + OH-     
            AKEQ2 = B(LH2O) * EXP(D(LH2O) * TFAC)                               
            OHION = AKEQ2 * HPLUSI                                              
            HLCONST = KH * (1.0 + AKEQ1/OHION)                                  
                                                                                
         END IF                                                                 
                                                                                
      END IF                                                                    
                                                                                
      RETURN                                                                    
      END                                                                       
                                                                                
!----------------------------------------------------------------------         
      INTEGER FUNCTION INDEX1 (NAME, N, NLIST)                                  
!----------------------------------------------------------------------         
                                                                                
!***********************************************************************        
!  subroutine body starts at line 39                                            
!                                                                               
!  FUNCTION:                                                                    
!                                                                               
!    Searches for NAME in list NLIST and returns the subscript                  
!    (1...N) at which it is found, or returns 0 when NAME not                   
!    found in NLIST                                                             
!                                                                               
!  PRECONDITIONS REQUIRED:  none                                                
!                                                                               
!  SUBROUTINES AND FUNCTIONS CALLED:  none                                      
!                                                                               
!  Based on index1 routine from Models-3 I/O Library                            
!                                                                               
!***********************************************************************        
                                                                                
      IMPLICIT NONE                                                             
                                                                                
!.......   Arguments and their descriptions:                                    
                                                                                
      CHARACTER*(*) NAME        !  Character string being searched for          
      INTEGER       N           !  Length of array to be searched               
      CHARACTER*(*) NLIST(*)    !  array to be searched                         
                                                                                
!.......   Local variable:                                                      
                                                                                
      INTEGER       I   !  loop counter                                         
                                                                                
!.....................................................................          
!.......   begin body of INDEX1()                                               
                                                                                
      DO 100 I = 1, N                                                           
                                                                                
          IF (TRIM(NAME) .EQ. TRIM(NLIST(I))) THEN    ! Found NAME in NLIST     
              INDEX1 = I                                                        
              RETURN                                                            
          END IF                                                                
                                                                                
100   CONTINUE                                                                  
                                                                                
      INDEX1 = 0        !  not found                                            
      RETURN                                                                    
                                                                                
      END                                                                       
