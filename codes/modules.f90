! --- This group of modules takes on the role of data declarations and          
! --- included common files;  initializations are also done here using          
! --- data statements much like a BLOCK DATA structure.                         
! --- Individual modules are:                                                   
!      module mroad1                                                            
!      module mroad2                                                            
!      module mqscale                                                           
                                                                                
!-----------------------------------------------------------------------        
      module mroad1                                                             
!-----------------------------------------------------------------------        
! --- CALPUFF    Version: 7.2.1       Level: 141201           |MROAD1|          
! ---            Constant/Scaled Road-source data                               
!-----------------------------------------------------------------------        
      integer              :: nrd1,nrdseg1,nsfrds                               
      integer, allocatable :: nptrd1(:)                                 ! (nrd1)
      integer, allocatable :: iroad1(:),newrd1(:)                       ! (nrdse
      integer, allocatable :: idsfrds(:,:)                              ! (mxspe
      integer, allocatable :: ixrefrds(:)                               ! (nsfrd
                                                                                
      character(len=16), allocatable :: srcnamrd1(:)                    ! (nrd1)
      character(len=40), allocatable :: csfrds(:)                       ! (nsfrd
                                                                                
      real, allocatable    :: htrd1(:),sz0rd1(:),sy0rd1(:)              ! (nrd1)
      real, allocatable    :: qrd1(:,:)                                 ! (mxspe
      real, allocatable    :: rdlen1(:)                                 ! (nrdse
      real, allocatable    :: xrd1grd(:,:),yrd1grd(:,:),elrd1(:,:)      ! (2,nrd
                                                                                
! --- Variables:                                                                
! ---------------                                                               
!                                                                               
! --- Variables for named roads                                                 
!                     NRD1 - integer - Number of roads                          
!          SRCNAMRD1(nrd1) - char*16 - Road names                               
!              HTRD1(nrd1) - real    - Effective release height (m)             
!             SZ0RD1(nrd1) - real    - Initial sigma z (m)                      
!             SY0RD1(nrd1) - real    - Initial sigma y (m)                      
!        QRD1(mxspec,nrd1) - real    - Emission rate (g/s/m) for each           
!                                      pollutant                                
!             NPTRD1(nrd1) - real    - Number of points defining road           
!                                                                               
! --- Variables for road-species pairs with scaled emissions                    
!                   NSFRDS - integer - Number of road-species pairs             
!                                      with emissions scaling factors           
!     IDSFRDS(mxspec,nrd1) - integer - Pointer to road-species pair             
!                                      index, 0 to NSFRDS                       
!                                      (0 if no scaling)                        
!           CSFRDS(nsfrds) - char*40 - List of scale-factor table names         
!                                      for road-species pairs                   
!         IXREFRDS(nsfrds) - integer - Cross-reference pointer from             
!                                      road-species pairs to                    
!                                      scale-factor tables                      
!                                                                               
! --- Variables for road segments that emit puffs/slugs                         
!                  NRDSEG1 - integer - Number of emitting road segments         
!                                      (Total over all roads)                   
!          IROAD1(nrdseg1) - integer - Road number for this segment             
!          RDLEN1(nrdseg1) - real    - Road length (m) for this segment         
!       XRD1GRD(2,nrdseg1) - real    - X coordinate of the ends of road         
!                                      segments in grid units                   
!                                      (i.e., origin at (0.0,0.0))              
!       YRD1GRD(2,nrdseg1) - real    - Y coordinate of the ends of road         
!                                      segments in grid units                   
!                                      (i.e., origin at (0.0,0.0))              
!         ELRD1(2,nrdseg1) - real    - Ground elevation of the ends of          
!                                      road segments (m MSL)                    
!          NEWRD1(nrdseg1) - integer - Number of puffs/slugs released           
!                                      by each road during the current          
!                                      step                                     
!-----------------------------------------------------------------------        
      end module mroad1                                                         
!-----------------------------------------------------------------------        
                                                                                
                                                                                
!-----------------------------------------------------------------------        
      module mroad2                                                             
!-----------------------------------------------------------------------        
! --- CALPUFF    Version: TNG-7.1.0    Level: 141201           |MROAD2|         
! ---            Time-varying Road-source data                                  
!-----------------------------------------------------------------------        
      integer              :: nrd2,nrdseg2,nse7,nrddat                          
                                                                                
      integer, allocatable :: ixrem7(:)                                 ! (mxspe
      integer, allocatable :: nptrd2(:)                                 ! (nrd2)
      integer, allocatable :: iroad2(:),newrd2(:)                       ! (nrdse
                                                                                
      character(len=12), allocatable :: cslst7(:)                       ! (mxspe
      character(len=16), allocatable :: cid7(:)                         ! (nrd2)
                                                                                
      real,    allocatable :: xmwem7(:)                                 ! (mxspe
      real,    allocatable :: rdlen2(:)                                 ! (nrdse
      real,    allocatable :: xrd2grd(:,:),yrd2grd(:,:),elrd2(:,:)      ! (2,nrd
      real,    allocatable :: htrd2(:,:),sz0rd2(:,:),sy0rd2(:,:)        ! (mxqst
      real,    allocatable :: qrd2(:,:,:)                               ! (mxspe
                                                                                
! --- Arrays for data stored for each RDEMARB.DAT file (nrddat files)           
                                                                                
      integer, allocatable :: ibsrc7(:),iesrc7(:),ibdathr7(:),ibsec7(:) ! (nrdda
      integer, allocatable :: iedathr7(:),iesec7(:)                     ! (nrdda
      integer, allocatable :: iutmznrd2(:)                              ! (nrdda
      integer, allocatable :: nstep7(:),mfrd2(:)                        ! (nrdda
      integer, allocatable :: ndhrqb7(:,:),nsecqb7(:,:)                 ! (mxqst
      integer, allocatable :: ndhrqe7(:,:),nsecqe7(:,:)                 ! (mxqst
                                                                                
      real,    allocatable :: xtz7(:),t2btz7(:)                         ! (nrdda
      real,    allocatable :: feastrd2(:),fnorthrd2(:)                  ! (nrdda
      real,    allocatable :: rnlat0rd2(:),relon0rd2(:)                 ! (nrdda
      real,    allocatable :: rnlat1rd2(:),rnlat2rd2(:)                 ! (nrdda
                                                                                
      character(len=8),  allocatable :: pmaprd2(:),datumrd2(:)          ! (nrdda
      character(len=4),  allocatable :: utmhemrd2(:),xyunitrd2(:)       ! (nrdda
      character(len=12), allocatable :: datenrd2(:)                     ! (nrdda
      character(len=16), allocatable :: verrdarb(:)                     ! (nrdda
      character(len=132),allocatable :: rddat(:)                        ! (nrdda
                                                                                
! --- Variables:                                                                
! ---------------                                                               
!                                                                               
! --- Variables for named roads                                                 
!             NSE7 - integer  - Number of emitted species                       
!   CSLST7(mxspec) - char*12  - Species identifiers                             
!   XMWEM7(mxspec) - real     - Molecular weight for each species               
!   IXREM7(mxspec) - integer  - Cross referencing array of NSE7                 
!                               values relating species ordering                
!                               in the emissions file to the                    
!                               ordering in the main conc. array                
!             NRD2 - integer  - Total number of roads                           
!       CID7(nrd2) - char*16  - Road names                                      
!     NPTRD2(nrd2) - real     - Number of points defining road                  
!                                                                               
! --- Variables for each file                                                   
!           NRDDAT - integer  - Total number of RDEMARB.DAT files               
!    RDDAT(nrddat) - char*132 - Path & filename for the input CALPUFF           
!                               file(s) containing ROAD sources with            
!                               arbitrarily-varying location and                
!                               emissions                                       
!                               (default: RDEMARB.DAT, for 1 file)              
!    MFRD2(nrddat) - integer  - Flag for file type                              
!                                 0: UNFORMATTED (not supported!)               
!                                 1: FORMATTED                                  
! VERRDARB(nrddat) - char*16  - Version of the input CALPUFF                    
!                               file(s) containing road sources                 
!                               with arbitrarily-varying location and           
!                               emissions                                       
!                               (RDEMARB.DAT)                                   
!   IBSRC7(nrddat) - integer  - Index for first source in a RDEMARB.DAT         
!                               file                                            
!   IESRC7(nrddat) - integer  - Index for last source in a RDEMARB.DAT          
!                               file                                            
!  IBDATHR7(nrddat)- integer  - Date/hour at beginning of period for            
!                               the first data record in the file               
!                               (YYYYJJJHH, where YYYY=year,                    
!                               JJJ=Julian day, HH=hour [00-23 LST])            
!   IBSEC7(nrddat) - integer  - Seconds of the first data record in the         
!                               file  (0000-3599)                               
!  IEDATHR7(nrddat)- integer  - Date/hour at end of period for                  
!                               the last data record in the file                
!                               (YYYYJJJHH, where YYYY=year,                    
!                               JJJ=Julian day, HH=hour [00-23 LST])            
!   IESEC7(nrddat) - integer  - Seconds of the last data record in the          
!                               file  (0000-3599)                               
!     XTZ7(nrddat) - real     - Time zone (UTC=LST+XTZ7)                        
!   T2BTZ7(nrddat) - real     - Hours to ADD to Local Time to obtain            
!                               Base Time (xtz7-xbtz)                           
!                                                                               
! --- MAP Projection                                                            
!   PMAPRD2(nrddat) -char*8    - Character code for map projection              
!                                 UTM :  Universal Transverse Mercator          
!                                 LCC :  Lambert Conformal Conic                
!                                 PS  :  Polar Stereographic                    
!                                 EM  :  Equatorial Mercator                    
!                                 LAZA:  Lambert Azimuthal Equal Area           
!                                 TTM :  Tangential Transverse Mercator         
! UTMHEMRD2(nrddat) -char*4    - Base hemisphere for UTM projection             
!                                (S=southern, N=northern)                       
!  DATUMRD2(nrddat) -char*8    - Datum-Region for grid coordinates              
!  DATENRD2(nrddat) -char*12   - NIMA date for datum parameters                 
!                                 (MM-DD-YYYY  )                                
! XYUNITRD2(nrddat) -char*4    - Units for coordinates (e.g., KM)               
!                                                                               
!  IUTMZNRD2(nrddat) -integer  - UTM zone for UTM projection                    
!  FEASTRD2(nrddat)  -real     - False Easting (km) at projection origin        
!  FNORTHRD2(nrddat) -real     - False Northing (km) at projection origin       
!  RNLAT0RD2(nrddat) -real     - N. latitude & E. longitude of x=0 and y=0      
!  RELON0RD2(nrddat)  (deg)      of map projection (Used only if PMAP =         
!                                LCC, PS, EM, TTM or LAZA)                      
!                                NOTE: longitude neg in western hemisphere      
!  RNLAT1RD2(nrddat) - real    - Matching N. latitude(s) for projection         
!  RNLAT2RD2(nrddat)  (deg)      (Used only if PMAP3= LCC, PS, or EM)           
!                            LCC :  Projection cone slices through              
!                                   Earth's surface at XLAT1 and XLAT2          
!                            PS  :  Projection plane slices through             
!                                   Earth at XLAT1                              
!                            EM  :  Projection cylinder slices through          
!                                   Earth at [+/-] XLAT1                        
!                                                                               
! --- Variables for road-segments that emit puffs/slugs                         
! --- (other properties are taken from the (nrd2) arrays)                       
!                  NRDSEG2 - integer - Number of emitting road segments         
!                                      (Total over all roads)                   
!          IROAD2(nrdseg2) - integer - Road number for this segment             
!          RDLEN2(nrdseg2) - real    - Road length (m) for this segment         
!       XRD2GRD(2,nrdseg2) - real    - X coordinate of the ends of road         
!                                      segments in grid units                   
!                                      (i.e., origin at (0.0,0.0))              
!       YRD2GRD(2,nrdseg2) - real    - Y coordinate of the ends of road         
!                                      segments in grid units                   
!                                      (i.e., origin at (0.0,0.0))              
!         ELRD2(2,nrdseg2) - real    - Ground elevation of the ends of          
!                                      road segments (m MSL)                    
!                                                                               
! ---  Variable data  ---                                                       
!                                                                               
!    NSTEP7(nrddat)  - integer  - Number of emission steps in                   
!                                 current timestep for each file                
! NDHRQB7(mxqstep,nrddat) & NSECQB7(mxqstep,nrddat)                             
!                    - integer  - Starting time for which                       
!                                 emissions data in current set of              
!                                 records is valid                              
!                                 (YYYYJJJHH & SSSS)                            
! NDHRQE7(mxqstep,nrddat) & NSECQE7(mxqstep,nrddat)                             
!                    - integer  - Ending time for which                         
!                                 emissions data in current set of              
!                                 records is valid                              
!                                 (YYYYJJJHH & SSSS)                            
!    HTRD2(mxqstep,nrd2) - real - Effective height (mAGL)                       
!   SY0RD2(mxqstep,nrd2) - real - Initial sigma y (m)                           
!   SZ0RD2(mxqstep,nrd2) - real - Initial sigma z (m)                           
!  QRD2(mxspec,mxqstep,nrd2)                                                    
!                        - real - Emission rate (g/m/s) for each                
!                                 pollutant                                     
!    NEWRD2(nrdseg2) - integer  - Number of puffs/slugs released                
!                                 by each road during the current               
!                                 step                                          
!                                                                               
!-----------------------------------------------------------------------        
      end module mroad2                                                         
!-----------------------------------------------------------------------        
                                                                                
                                                                                
!-----------------------------------------------------------------------        
      module mqscale                                                            
!-----------------------------------------------------------------------        
! --- CALPUFF    Version: TNG-7.1.0    Level: 141201          |MQSCALE|         
! ---            Emission-Rate Scaling Factors (control file sources)           
!-----------------------------------------------------------------------        
      integer, parameter   :: nqsftype = 9                                      
      integer              :: nqsfval(nqsftype)                                 
      integer              :: nqsfcol(nqsftype),nqsfrow(nqsftype)               
      integer              :: mapivary(6)                                       
      character(len=24)    :: cqsftype(nqsftype)                                
      real                 :: wqsf(5,13),tqsf(11,13)                            
                                                                                
      integer              :: nqsftab                                           
      integer, allocatable :: iqsftype(:)                                ! (nqsf
      real, allocatable    :: qsftab(:,:)                                ! (mxqs
      character(len=40), allocatable :: cqsfname(:)                      ! (nqsf
                                                                                
! --- Assignments:                                                              
! -----------------                                                             
      data nqsfval/1,   12,   7,                                               &
     &             24, 168, 288,                                               &
     &             6,   36,  12/                                                
      data nqsfcol/1,   12,   7,                                               &
     &             24,  24,  24,                                               &
     &             6,    6,  12/                                                
      data nqsfrow/1,    1,   1,                                               &
     &             1,    7,  12,                                               &
     &             1,    6,   1/                                                
      data mapivary/1, 4, 2, 6, 8, 9/                                           
      data cqsftype/                'CONSTANT1               ',                &
     &   'MONTH12                 ','DAY7                    ',                &
     &   'HOUR24                  ','HOUR24_DAY7             ',                &
     &   'HOUR24_MONTH12          ','WSP6                    ',                &
     &   'WSP6_PGCLASS6           ','TEMPERATURE12           '/                 
! --- NOTE ---------                                                            
!           CONSTANT1        1   scaling factor                                 
!           MONTH12          12  scaling factors: months 1-12                   
!           DAY7             7   scaling factors: days 1-7                      
!                              [SUNDAY,MONDAY, ... FRIDAY,SATURDAY]             
!           HOUR24           24  scaling factors: hours 1-24                    
!           HOUR24_DAY7      168 scaling factors: hours 1-24,                   
!                              repeated  7 times:                               
!                              [SUNDAY,MONDAY, ... FRIDAY,SATURDAY]             
!           HOUR24_MONTH12   288 scaling factors: hours 1-24,                   
!                              repeated 12 times: months 1-12                   
!           WSP6             6   scaling factors: wind speed classes 1-6        
!                              [speed classes (WSCAT)]                          
!           WSP6_PGCLASS6    36  scaling factors: wind speed classes 1-6        
!                              repeated  6 times: PG classes A,B,C,D,E,F        
!                              [speed classes (WSCAT)]                          
!           TEMPERATURE12    12  scaling factors: temp(K) classes 1-12          
!                              [temperature classes (TKCAT)]                    
! -----------------                                                             
!                                                                               
! --- Variables:                                                                
! ---------------                                                               
!                                                                               
! --- Variables for defining emission-rate scaling factors                      
!                 NQSFTAB - integer - Number of tables of                       
!                                     emissions scaling factors                 
!       IQSFTYPE(nqsftab) - integer - Index of scale-factor type of             
!                                     each table                                
!       CQSFNAME(nqsftab) - char*40 - Name of each scale-factor table           
!   QSFTAB(mxqsf,nqsftab) - real    - Emission scale-factors                    
!                NQSFTYPE - integer - Number of types of                        
!                                     emissions scaling factors                 
!      CQSFTYPE(nqsftype) - char*24 - Name of each scale-factor type            
!             MAPIVARY(6) - integer - Map pointer from the 6 IVARY              
!                                     choices to the corresponding              
!                                     CQSFTYPE() index                          
!       NQSFVAL(nqsftype) - integer - Number of scaling factors for             
!                                     each type                                 
!                                     (Max must = MXQSF in /params/)            
!       NQSFCOL(nqsftype) - integer - Number of print columns for each          
!       NQSFROW(nqsftype) - integer - Number of print rows for each             
!                                                                               
! --- Temperature and wind speed classes by source type (13)                    
!              WQSF(5,13) - real    - Wind speed class boundaries (m/s)         
!                                     (boundary is upper limit of class)        
!             TQSF(11,13) - real    - Temperature class boundaries (K)          
!                                     (boundary is upper limit of class)        
!     Source Types are:                                                         
!            1 = Point         Constant Emissions                               
!            2 = Point         Variable Emissions (no WS/T class used)          
!            3 = Poly. Area    Constant Emissions                               
!            4 = Poly. Area    Variable Emissions (no WS/T class used)          
!            5 = Line          Constant Emissions                               
!            6 = Line          Variable Emissions (no WS/T class used)          
!            7 = Volume        Constant Emissions                               
!            8 = Grid Volume   Variable Emissions (no WS/T class used)          
!            9 = Boundary Condition                                             
!          (10)= Flare         Constant Emissions                               
!           11 = Flare         Variable Emissions (no WS/T class used)          
!           12 = Road          Constant Emissions                               
!           13 = Road          Variable Emissions (no WS/T class used)          
!                                                                               
!-----------------------------------------------------------------------        
      end module mqscale                                                        
!-----------------------------------------------------------------------        
