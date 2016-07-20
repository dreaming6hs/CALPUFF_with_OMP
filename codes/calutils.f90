!------------------------------------------------------------------------------ 
! --- CALUTILS -- CALPUFF SYSTEM UTILITIES                                      
!------------------------------------------------------------------------------ 
!                                                                               
! --- CALUTILS   Version: 7.0.0     Level: 141010                               
!                                                                               
!     Copyright (c) 2014 by Exponent, Inc.                                      
!                                                                               
! -----------------------------                                                 
! --- CONTENT:                                                                  
! -----------------------------                                                 
! --- Coordinates                                                               
!      subroutine xtractll                                                      
! --- Year 2000                                                                 
!      subroutine yr4                                                           
!      subroutine yr4c                                                          
!      subroutine qayr4                                                         
! --- Date/Time                                                                 
!      subroutine julday                                                        
!      subroutine grday                                                         
!      subroutine dedat                                                         
!      subroutine deltt                                                         
!      subroutine incr                                                          
!      subroutine indecr                                                        
!      subroutine incrs                                                         
!      subroutine deltsec                                                       
!      subroutine midnite                                                       
!      subroutine basrutc                                                       
!      subroutine utcbasr                                                       
! --- Control file                                                              
!      subroutine filcase                                                       
!      subroutine readin                                                        
!      subroutine altonu                                                        
!      subroutine deblnk                                                        
!      subroutine deplus                                                        
!      subroutine tright                                                        
!      subroutine tleft                                                         
!      subroutine setvar                                                        
!      subroutine allcap                                                        
! --- System                                                                    
!      subroutine datetm                                                        
!      subroutine fmt_date                                                      
!      subroutine etime                                                         
!      subroutine undrflw                                                       
!      subroutine comline                                                       
! --- Error                                                                     
!      subroutine open_err                                                      
! -----------------------------                                                 
!                                                                               
! --- UPDATE                                                                    
! --- V2.6.0-V7.0.0 141010    :Add error-report for file-open                   
!                              New     : OPEN_ERR                               
! --- V2.58-V2.6.0 140318(MBN):Use F95 intrinsic procedures for date and time.  
!                              Modified: DATETM                                 
!                              Removed obsolete Compaq, Microsoft, and HP       
!                              compiler codes, and removed getcl                
!                              Modified: COMLINE                                
! --- V2.571-V2.58 110225(DGS):Add variable type 5 to control file processor    
!                              to allow character array variables               
!                              Modified: READIN, ALTONU, SETVAR                 
! --- V2.57-V2.571 090511(DGS):Add routine to reformat a date string            
!                              New     : FMT_DATE                               
! --- V2.56-V2.57 090202(DGS): Increase control file line length to 200         
!                              characters                                       
!                              Modified: PARAMS.CAL, READIN                     
!                              Activate CPU clock using F95 system routine      
!                              Modified: DATETM                                 
! --- V2.55-V2.56 080407(DGS): Exponential notation processing in ALTONU did    
!                              not properly interpret an entry without a        
!                              decimal point.                                   
! --- V2.54-V2.55 070327(DGS): Format for output time zone stringin BASRUTC     
!                              wrote zone zero as 'UTC+0  0' instead of         
!                              'UTC+0000'                                       
!                              Add RETURN statement to BASRUTC and UTCBASR      
! --- V2.53-V2.54 061020(DGS): Allow negative increments in INCRS               
! --- V2.52-V2.53 060626(DGS): Remove routine GLOBE1 (move to COORDLIB)         
! --- V2.51-V2.52 060519(DGS): Modify search for '=' in READIN to allow         
!                              for blanks between c*12 variable name and        
!                              the '=' sign (internal blanks are not removed    
!                              after V2.2)                                      
! --- V2.5-V2.51 051019 (KAM): Add Albers Conical Equal Area projection         
!                              in GLOBE1                                        
! --- V2.4-V2.5  041123 (FRR): add subroutine BASRUTC to convert real           
!                              base time zone to character UTC time zone        
!                              and UTCBASR for the backward conversion          
! --- V2.3-V2.4  041029 (DGS): Add routine INCRS to change time by a            
!                              number of seconds                                
!                              Add routine MIDNITE - converts timestamp         
!                              from day N, time 0000                            
!                              to day N-1, time 2400                            
! --- V2.2-V2.3  040330 (DGS): Replace filename strings c*70 with c*132         
!                              (FILCASE, COMLINE)                               
!                              Allow for spaces within pathnames by adding      
!                              new TLEFT and TRIGHT trim subroutines            
! --- V2.1-V2.2  030528 (DGS): Screen for valid UTM zone using                  
!                              absolute value (S. Hem. zones are                
!                              negative) in GLOBE1                              
! --- V2.0-V2.1  030402 (DGS): Remove routine GLOBE                             
!                              Split DEBLNK action (removes ' ', '+')           
!                              into DEBLNK and DEPLUS                           
!                              Add routine UNDRFLW                              
!                              Add false Easting and Northing (GLOBE1)          
!                              Add TYPE argument to XTRACTLL                    
!                              Change format XTRACTLL (f16) to (f16.0)          
! --- V1.1-V2.0  021018 (DGS): Add routines for new COORDS                      
! --- V1.0-V1.1  020828 (DGS): Add check for YYYY on input   (YR4C)             
!                                                                               
!                                                                               
!----------------------------------------------------------------------         
      subroutine xtractll(io,type,clatlon,rlatlon)                              
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 030402               XTRACTLL         
!                D. Strimaitis   EarthTech                                      
!                                                                               
! --- PURPOSE:  Extract the real latitude or longitude from a character         
!               string that contains the N/S or E/W convention                  
!               character, and express result as either North Latitude          
!               or East Longitude                                               
!                                                                               
! --- UPDATE                                                                    
! --- V2.1 (030402) from V2.0 (010713) (DGS)                                    
!               - Add TYPE argument for QA                                      
!               - Change format (f16) to (f16.0) to satisfy different           
!                 compilers                                                     
!                                                                               
!                                                                               
! --- INPUTS:                                                                   
!               IO - integer    - Unit number for list file output              
!             TYPE - char*4     - LAT or LON                                    
!          CLATLON - char*16    - Latitude or longitude (degrees), with         
!                                 1 character that denotes convention           
!                                 (e.g. 'N  45.222' or  '-35.999s')             
!                                                                               
! --- OUTPUT:                                                                   
!          RLATLON - real       - North Latitude or East Longitude              
!                                 (degrees)                                     
!                                                                               
! --- XTRACTLL called by: (utility)                                             
! --- XTRACTLL calls:     DEBLNK, ALLCAP                                        
!----------------------------------------------------------------------         
!                                                                               
! --- Include parameter statements                                              
      include 'params.cal'                                                      
                                                                                
      character*1 cstor1(mxcol),cstor2(mxcol)                                   
      character*16 clatlon, clatlon2                                            
      character*4 type                                                          
      logical ltype                                                             
                                                                                
      ltype=.FALSE.                                                             
                                                                                
! --- Initialize character variables for output                                 
      clatlon2='                '                                               
      do i=1,20                                                                 
         cstor2(i)=' '                                                          
      enddo                                                                     
                                                                                
! --- Was valid type provided?                                                  
      if(type.NE.'LAT ' .AND. type.NE.'LON ') then                              
         write(io,*) 'XTRACTLL:  FATAL ERROR reported when ',                  &
     &               'extracting Latitude/Longitude'                            
         write(io,*) 'Invalid type:  ',type                                     
         write(io,*) 'Expected LAT or LON'                                      
         write(*,*)                                                             
         stop 'Halted in XTRACTLL -- see list file'                             
      endif                                                                     
                                                                                
! --- Pass c*16 string into storage array 1                                     
      do i=1,16                                                                 
         cstor1(i)=clatlon(i:i)                                                 
      enddo                                                                     
! --- Pad out to position 20                                                    
      do i=17,20                                                                
         cstor1(i)=' '                                                          
      enddo                                                                     
                                                                                
! --- Remove blank characters from string, place in storage array 2             
! --- (Use a 20-character field here for a margin at end of string)             
      call DEBLNK(cstor1,1,20,cstor2,nlim)                                      
!                                                                               
! --- Convert lower case letters to upper case                                  
      call ALLCAP(cstor2,nlim)                                                  
                                                                                
! --- Interpret valid convention character (N,S,E,W)                            
      nchar=0                                                                   
      ichar=0                                                                   
      ilat=0                                                                    
      ilon=0                                                                    
                                                                                
      do i=1,nlim                                                               
         if(cstor2(i).EQ.'N') then                                              
            ilat=1                                                              
            ichar=i                                                             
            nchar=nchar+1                                                       
         elseif(cstor2(i).EQ.'S') then                                          
            ilat=2                                                              
            ichar=i                                                             
            nchar=nchar+1                                                       
         elseif(cstor2(i).EQ.'W') then                                          
            ilon=1                                                              
            ichar=i                                                             
            nchar=nchar+1                                                       
         elseif(cstor2(i).EQ.'E') then                                          
            ilon=2                                                              
            ichar=i                                                             
            nchar=nchar+1                                                       
         endif                                                                  
      enddo                                                                     
                                                                                
! --- Was 1 valid character found?                                              
      if(nchar.NE.1) then                                                       
         write(io,*) 'XTRACTLL:  FATAL ERROR reported when ',                  &
     &               'extracting Latitude/Longitude'                            
         write(io,*) 'N,S,E,W character is missing or repeated'                 
         write(io,*) 'Lat/Lon = ',clatlon                                       
         write(*,*)                                                             
         stop 'Halted in XTRACTLL -- see list file'                             
      endif                                                                     
                                                                                
! --- Was valid character the right type?                                       
      if(type.EQ.'LAT ' .AND. ilat.EQ.0) ltype=.TRUE.                           
      if(type.EQ.'LON ' .AND. ilon.EQ.0) ltype=.TRUE.                           
      if(LTYPE) then                                                            
         write(io,*) 'XTRACTLL:  FATAL ERROR reported when ',                  &
     &               'extracting Latitude/Longitude'                            
         write(io,*) 'N,S,E,W character does not match type'                    
         write(io,*) 'Lat/Lon = ',clatlon                                       
         write(io,*) 'type    = ',type                                          
         write(*,*)                                                             
         stop 'Halted in XTRACTLL -- see list file'                             
      endif                                                                     
                                                                                
! --- Remove character from string                                              
      do i=ichar,nlim                                                           
         cstor2(i)=cstor2(i+1)                                                  
      enddo                                                                     
                                                                                
! --- Search for position of decimal point                                      
      ipt=0                                                                     
      do i=1,nlim                                                               
         if(cstor2(i).EQ.'.') ipt=i                                             
      enddo                                                                     
                                                                                
! --- Add a decimal point if needed                                             
      if(ipt.EQ.0) then                                                         
         cstor2(nlim)='.'                                                       
      endif                                                                     
                                                                                
! --- Pass resulting "number" back into c*16 variable                           
      do i=1,nlim                                                               
         clatlon2(i:i)=cstor2(i)                                                
      enddo                                                                     
                                                                                
! --- Get real part                                                             
      read(clatlon2,'(f16.0)') rlatlon                                          
                                                                                
! --- Convert to either N. Lat. or E. Lon., if needed                           
      if(ilat.EQ.2) then                                                        
         rlatlon=-rlatlon                                                       
      elseif(ilon.EQ.1) then                                                    
         rlatlon=-rlatlon                                                       
      endif                                                                     
                                                                                
! --- Condition longitude to be -180 to +180                                    
      if(ilon.GT.0) then                                                        
         if(rlatlon.GT.180.) then                                               
            rlatlon=rlatlon-360.                                                
         elseif(rlatlon.LT.-180.) then                                          
            rlatlon=rlatlon+360.                                                
         endif                                                                  
      endif                                                                     
                                                                                
      return                                                                    
      end                                                                       
!----------------------------------------------------------------------         
      subroutine yr4(io,iyr,ierr)                                               
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 991104                    YR4         
! ---            D. Strimaitis,   Earth Tech                                    
!                                                                               
! --- PURPOSE: Checks/converts 2-digit year to 4-digit year                     
!                                                                               
! --- INPUTS:                                                                   
!                IO - integer    - Unit number for list file output             
!               IYR - integer    - Year (YYYY or YY)                            
!                                                                               
!     Common block /Y2K/:                                                       
!             IYYLO - integer    - Smallest 2-digit year for which              
!                                  'old' century marker is used                 
!             ICCLO - integer    - 2-digit ('old') century                      
!                                                                               
! --- OUTPUT:                                                                   
!               IYR - integer    - Year (YYYY)                                  
!              IERR - integer    - Error code: 0=OK, 1=FATAL                    
!                                                                               
! --- YR4 called by:  Input routines reading 'year' data                        
! --- YR4 calls:      none                                                      
!----------------------------------------------------------------------         
!                                                                               
      common/y2k/iyylo,icclo                                                    
                                                                                
      ierr=0                                                                    
                                                                                
! --- Test for 4-digit year (must exceed 1000)                                  
      if(iyr.GT.1000) then                                                      
! ---    Passes 11th Century test (large year not trapped)                      
         return                                                                 
      elseif(iyr.LT.100 .AND. iyr.GE.0) then                                    
! ---    2-digit year                                                           
! ---    Construct 4-digit year                                                 
         if(iyr.LT.iyylo) then                                                  
            iyr=(icclo+1)*100+iyr                                               
         else                                                                   
            iyr=icclo*100+iyr                                                   
         endif                                                                  
      else                                                                      
! ---    Year not recognized                                                    
         ierr=1                                                                 
         write(io,*)'ERROR in YR4 --- Year not recognized: ',iyr                
         write(*,*)'ERROR in YR4 --- Year not recognized: ',iyr                 
      endif                                                                     
                                                                                
      return                                                                    
      end                                                                       
!----------------------------------------------------------------------         
      subroutine yr4c(iyr)                                                      
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 020828                   YR4C         
! ---            D. Strimaitis,  Earth Tech                                     
!                                                                               
! --- PURPOSE: Checks/converts 2-digit year to 4-digit year (CURRENT)           
!                                                                               
! --- UPDATE                                                                    
! --- V1.0-V1.1     020828  (DGS): Add check for YYYY on input                  
!                                                                               
! --- INPUTS:                                                                   
!               IYR - integer    - Year (YYYY or YY)                            
!                                                                               
! --- OUTPUT:                                                                   
!               IYR - integer    - Year (YYYY)                                  
!                                                                               
! --- YR4C called by:  host subroutines                                         
! --- YR4C calls:      none                                                     
!----------------------------------------------------------------------         
! --- Set parameters for converting a current year (1999 - 2098)                
! --- Use KCCLO as century digits for years GE KYYLO                            
      data kyylo/99/, kcclo/19/                                                 
                                                                                
! --- Test for 4-digit year (must exceed 1000)                                  
      if(iyr.GT.1000) then                                                      
! ---    Passes 11th Century test (large year not trapped)                      
         return                                                                 
      elseif(iyr.LT.100 .AND. iyr.GE.0) then                                    
! ---    2-digit year                                                           
! ---    Construct 4-digit year                                                 
         if(iyr.LT.kyylo) then                                                  
            iyr=(kcclo+1)*100+iyr                                               
         else                                                                   
            iyr=kcclo*100+iyr                                                   
         endif                                                                  
      else                                                                      
! ---    Year not recognized                                                    
         write(*,*)'ERROR in YR4C --- Year not recognized: ',iyr                
      endif                                                                     
                                                                                
      return                                                                    
      end                                                                       
!----------------------------------------------------------------------         
      subroutine qayr4(io,iyr,metrun,ierr)                                      
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 991104                  QAYR4         
! ---            D. Strimaitis,   Earth Tech                                    
!                                                                               
! --- PURPOSE: Defines century and year markers to use in converting            
! ---          2-digit year to 4-digit year                                     
! ---          The IBYR (YYYY) must be provided in the control file             
!                                                                               
! --- INPUTS:                                                                   
!                IO - integer    - Unit number for list file output             
!               IYR - integer    - Year provided for start of run               
!            METRUN - integer    - Flag to run period in met file               
!                                  0 = do not run period                        
!                                  1 = run period                               
!                                                                               
! --- OUTPUT:                                                                   
!              IERR - integer    - Error code: 0=OK, 1=FATAL                    
!                                                                               
!     Common block /Y2K/:                                                       
!             IYYLO - integer    - Smallest 2-digit year for which              
!                                  'old' century marker is used                 
!             ICCLO - integer    - 2-digit ('old') century                      
!                                                                               
! --- QAYR4 called by:  host subroutines                                        
! --- QAYR4 calls:      none                                                    
!----------------------------------------------------------------------         
!                                                                               
      common/y2k/iyylo,icclo                                                    
                                                                                
! --- Sets parameters for the starting century marker (CC) and the              
! --- 2-digit year (YY) used as the marker between the starting century         
! --- and the next century.  For example, if CC=19 and YY=30, then a            
! --- year less than 30 (say 15) is assumed to be 2015.  Any year               
! --- greater than or equal to 30 (say 56) is assumed to be 1956.               
                                                                                
! --- Set number of years prior to start of simulation that must not            
! --- be placed in the next century                                             
      data ibackyr/50/                                                          
                                                                                
      ierr=0                                                                    
                                                                                
! --- Expect explicit starting year (YYYY)                                      
! --- Test for 4-digit year (must exceed 1000)                                  
      if(iyr.GT.1000) then                                                      
! ---    Passes 11th Century test (large year not trapped)                      
! ---    Back up IBACKYR years to set IYYLO                                     
         kyr=iyr-ibackyr                                                        
! ---    Extract starting 2-digit century and 2-digit year                      
         icclo=kyr/100                                                          
         iyylo=kyr-icclo*100                                                    
                                                                                
! ---    Warn user that control file input is used to convert to YYYY           
         iyr1=icclo*100+iyylo                                                   
         iyr2=(icclo+1)*100+iyylo-1                                             
         write(io,*)                                                            
         write(io,*)'-------------------------------------------------'         
         write(io,*)'NOTICE: Starting year in control file sets the'            
         write(io,*)'        expected century for the simulation.  All'         
         write(io,*)'        YY years are converted to YYYY years in'           
         write(io,*)'        the range: ',iyr1,iyr2                             
         write(io,*)'-------------------------------------------------'         
         write(io,*)                                                            
      else                                                                      
         ierr=1                                                                 
         write(*,*)                                                             
         write(*,*)'--------------------------------------------'               
         write(*,*)'QAYR4 -- Start year must be 4-digits!: ',iyr                
         if(metrun.EQ.1) then                                                   
            write(*,*)'         and must always be provided'                    
         endif                                                                  
         write(*,*)'--------------------------------------------'               
         write(*,*)                                                             
         write(io,*)                                                            
         write(io,*)'-------------------------------------------'               
         write(io,*)'QAYR4 -- Start year must be 4-digits!: ',iyr               
         if(metrun.EQ.1) then                                                   
            write(io,*)'         and must always be provided'                   
         endif                                                                  
         write(io,*)'-------------------------------------------'               
         write(io,*)                                                            
      endif                                                                     
                                                                                
      return                                                                    
      end                                                                       
!----------------------------------------------------------------------         
      subroutine julday(io,iyr,imo,iday,ijuldy)                                 
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 000602                 JULDAY         
! ---            J. Scire, SRC                                                  
!                                                                               
! --- PURPOSE:  Compute the Julian day number from the Gregorian                
!               date (month, day)                                               
!                                                                               
! --- UPDATE                                                                    
! ---               000602  (DGS): YYYY format for year                         
!                                                                               
! --- INPUTS:                                                                   
!            IO - integer      - Unit number for list file output               
!           IYR - integer      - Year                                           
!           IMO - integer      - Month                                          
!          IDAY - integer      - Day                                            
!                                                                               
! --- OUTPUT:                                                                   
!          IJUL - integer      - Julian day                                     
!                                                                               
! --- JULDAY called by:  host subroutines                                       
! --- JULDAY calls:      none                                                   
!----------------------------------------------------------------------         
!                                                                               
      integer kday(12)                                                          
      data kday/0,31,59,90,120,151,181,212,243,273,304,334/                     
!                                                                               
! --- Check for valid input data                                                
      ierr=0                                                                    
! --- Check for valid month                                                     
      if(imo.lt.1.or.imo.gt.12)ierr=1                                           
! --- Check for valid day in 30-day months                                      
      if(imo.eq.4.or.imo.eq.6.or.imo.eq.9.or.imo.eq.11)then                     
         if(iday.gt.30)ierr=1                                                   
      else if(imo.eq.2)then                                                     
         if(mod(iyr,4).eq.0)then                                                
! ---       February in a leap year                                             
            if(iday.gt.29)ierr=1                                                
         else                                                                   
! ---       February in a non-leap year                                         
            if(iday.gt.28)ierr=1                                                
         endif                                                                  
      else                                                                      
! ---    Check for valid day in 31-day months                                   
         if(iday.gt.31)ierr=1                                                   
      endif                                                                     
!                                                                               
      if(ierr.eq.1)then                                                         
         write(io,*)                                                            
         write(io,*)'ERROR in SUBR. JULDAY'                                     
         write(io,*)'Invalid date - IYR = ',iyr,' IMO = ',                     &
     &    imo,' IDAY = ',iday                                                   
         write(*,*)                                                             
         stop 'Halted in JULDAY -- see list file.'                              
      endif                                                                     
!                                                                               
! --- Compute the Julian day                                                    
      ijuldy=kday(imo)+iday                                                     
      if(imo.le.2)return                                                        
      if(mod(iyr,4).EQ.0)ijuldy=ijuldy+1                                        
!                                                                               
      return                                                                    
      end                                                                       
!----------------------------------------------------------------------         
      subroutine grday(io,iyr,ijul,imo,iday)                                    
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 000602                  GRDAY         
!                J. Scire, SRC                                                  
!                                                                               
! --- PURPOSE:  Compute the Gregorian date (month, day) from the                
!               Julian day                                                      
!                                                                               
! --- UPDATE                                                                    
! ---               000602  (DGS): YYYY format for year                         
!                                                                               
! --- INPUTS:                                                                   
!            IO - integer      - Unit number for list file output               
!           IYR - integer      - Year                                           
!          IJUL - integer      - Julian day                                     
!                                                                               
! --- OUTPUT:                                                                   
!           IMO - integer      - Month                                          
!          IDAY - integer      - Day                                            
!                                                                               
! --- GRDAY called by:  host subroutines                                        
! --- GRDAY calls:      none                                                    
!----------------------------------------------------------------------         
!                                                                               
      integer kday(12,2)                                                        
      data kday/31,59,90,120,151,181,212,243,273,304,334,365,                  &
     &          31,60,91,121,152,182,213,244,274,305,335,366/                   
!                                                                               
!                                                                               
      ileap=1                                                                   
      if(mod(iyr,4).eq.0)ileap=2                                                
      if(ijul.lt.1.or.ijul.gt.kday(12,ileap))go to 11                           
!                                                                               
      do 10 i=1,12                                                              
      if(ijul.gt.kday(i,ileap))go to 10                                         
      imo=i                                                                     
      iday=ijul                                                                 
      if(imo.ne.1)iday=ijul-kday(imo-1,ileap)                                   
      return                                                                    
10    continue                                                                  
!                                                                               
11    continue                                                                  
      write(io,12)iyr,ijul                                                      
12    format(//2x,'ERROR in SUBR. GRDAY -- invalid Julian day '//2x,           &
     & 'iyr = ',i5,3x,'ijul = ',i5)                                             
      write(*,*)                                                                
      stop 'Halted in GRDAY -- see list file.'                                  
      end                                                                       
!------------------------------------------------------------------------------ 
      subroutine dedat(idathr,iyr,ijul,ihr)                                     
!------------------------------------------------------------------------------ 
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 941215                  DEDAT         
! ---            J. Scire, SRC                                                  
!                                                                               
! --- Decode a date-time variable                                               
!                                                                               
! --- INPUTS:                                                                   
!            IDATHR - integer    - Date-time variable (YYYYJJJHH)               
!                                                                               
! --- OUTPUT:                                                                   
!               IYR - integer    - Year of precip. data (4 digits)              
!              IJUL - integer    - Julian day number of precip. data            
!               IHR - integer    - Ending hour (1-24) of precip. data           
!                                                                               
! --- DEDAT called by:  host subroutines                                        
! --- DEDAT calls:      none                                                    
!------------------------------------------------------------------------------ 
!
! --- wangzhm declared
      integer,intent(in) :: idathr
      integer,intent(out) :: iyr,ijul,ihr

! --- decode date and time                                                      
      iyr=idathr/100000                                                         
      ijul=idathr/100-iyr*1000                                                  
      ihr=idathr-iyr*100000-ijul*100                                            
!                                                                               
      return                                                                    
      end
! --- wangzhm declared
!     integer,intent(in) :: idathr
!     integer,intent(out) :: iyr,ijul,ihr
!------------------------------------------------------------------------------ 
      subroutine deltt(j1yr,j1jul,j1hr,j2yr,j2jul,j2hr,jleng)                   
!------------------------------------------------------------------------------ 
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 941215                  DELTT         
! ---            J. Scire, SRC                                                  
!                                                                               
! --- Compute the difference (in hours) between two dates & times               
! ---    (time #2 - time #1)                                                    
!                                                                               
! --- INPUTS:                                                                   
!              J1YR - integer    - Year of date/time #1                         
!             J1JUL - integer    - Julian day of date/time #1                   
!              J1HR - integer    - Hour of date/time #1                         
!              J2YR - integer    - Year of date/time #2                         
!             J2JUL - integer    - Julian day of date/time #2                   
!              J2HR - integer    - Hour of date/time #2                         
!                                                                               
! --- OUTPUT:                                                                   
!             JLENG - integer    - Difference (#2 - #1) in hours                
!                                                                               
! --- DELTT called by:  host subroutines                                        
! --- DELTT calls:      none                                                    
!------------------------------------------------------------------------------ 
!                                                                               
      jmin=min0(j1yr,j2yr)                                                      
!                                                                               
! --- find the number of hours between Jan. 1 of the "base" year and            
! --- the first date/hour                                                       
      if(j1yr.eq.jmin)then                                                      
         j1=0                                                                   
      else                                                                      
         j1=0                                                                   
         j1yrm1=j1yr-1                                                          
         do 10 i=jmin,j1yrm1                                                    
         if(mod(i,4).eq.0)then                                                  
            j1=j1+8784                                                          
         else                                                                   
            j1=j1+8760                                                          
         endif                                                                  
10       continue                                                               
      endif                                                                     
      j1=j1+(j1jul-1)*24+j1hr                                                   
!                                                                               
! --- find the number of hours between Jan. 1 of the "base" year and            
! --- the second date/hour                                                      
      if(j2yr.eq.jmin)then                                                      
         j2=0                                                                   
      else                                                                      
         j2=0                                                                   
         j2yrm1=j2yr-1                                                          
         do 20 i=jmin,j2yrm1                                                    
         if(mod(i,4).eq.0)then                                                  
            j2=j2+8784                                                          
         else                                                                   
            j2=j2+8760                                                          
         endif                                                                  
20       continue                                                               
      endif                                                                     
      j2=j2+(j2jul-1)*24+j2hr                                                   
!                                                                               
! --- compute the time difference (in hours)                                    
      jleng=j2-j1                                                               
!                                                                               
      return                                                                    
      end                                                                       
!----------------------------------------------------------------------         
      subroutine incr(io,iyr,ijul,ihr,nhrinc)                                   
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 000602                   INCR         
!                J. Scire, SRC                                                  
!                                                                               
! --- PURPOSE:  Increment the time and date by "NHRINC" hours                   
!                                                                               
! --- UPDATE                                                                    
! ---               000602  (DGS): add message to "stop"                        
! ---               980304  (DGS): Allow for a negative "increment" of          
!                                  up to 24 hours                               
! ---               980304  (DGS): Allow for arbitrarily large nhrinc           
!                                                                               
! --- INPUTS:                                                                   
!       IO     - integer - Unit number for list file output                     
!       IYR    - integer - Current year                                         
!       IJUL   - integer - Current Julian day                                   
!       IHR    - integer - Current hour (00-23)                                 
!       NHRINC - integer - Time increment (hours)                               
!                                                                               
!               NOTE: "NHRINC" must >= -24                                      
!                      Hour is between 00-23                                    
!                                                                               
! --- OUTPUT:                                                                   
!       IYR    - integer - Updated year                                         
!       IJUL   - integer - Updated Julian day                                   
!       IHR    - integer - Updated hour (00-23)                                 
!                                                                               
! --- INCR called by: host subroutines                                          
! --- INCR calls:     none                                                      
!----------------------------------------------------------------------         
!                                                                               
! --- Check nhrinc                                                              
      if(nhrinc.lt.-24) then                                                    
         write(io,*)'ERROR IN SUBR. INCR -- Invalid value of NHRINC ',         &
     &   '-- NHRINC = ',nhrinc                                                  
         write(*,*)                                                             
         stop 'Halted in INCR -- see list file.'                                
      endif                                                                     
                                                                                
! --- Save increment remaining (needed if nhrinc > 8760)                        
      nleft=nhrinc                                                              
!                                                                               
! --- Process change in hour                                                    
      if(nhrinc.gt.0)then                                                       
!                                                                               
10       ninc=MIN0(nleft,8760)                                                  
         nleft=nleft-ninc                                                       
!                                                                               
! ---    Increment time                                                         
         ihr=ihr+ninc                                                           
         if(ihr.le.23)return                                                    
!                                                                               
! ---    Increment day                                                          
         ijul=ijul+ihr/24                                                       
         ihr=mod(ihr,24)                                                        
!                                                                               
! ---    ILEAP = 0 (non-leap year) or 1 (leap year)                             
         if(mod(iyr,4).eq.0)then                                                
            ileap=1                                                             
         else                                                                   
            ileap=0                                                             
         endif                                                                  
!                                                                               
         if(ijul.gt.365+ileap) then                                             
! ---       Update year                                                         
            iyr=iyr+1                                                           
            ijul=ijul-(365+ileap)                                               
         endif                                                                  
!                                                                               
! ---    Repeat if more hours need to be added                                  
         if(nleft.GT.0) goto 10                                                 
!                                                                               
      elseif(nhrinc.lt.0)then                                                   
! ---    Decrement time                                                         
         ihr=ihr+nhrinc                                                         
         if(ihr.lt.0)then                                                       
            ihr=ihr+24                                                          
            ijul=ijul-1                                                         
            if(ijul.lt.1)then                                                   
               iyr=iyr-1                                                        
               if(mod(iyr,4).eq.0)then                                          
                  ijul=366                                                      
               else                                                             
                  ijul=365                                                      
               endif                                                            
            endif                                                               
         endif                                                                  
      endif                                                                     
!                                                                               
      return                                                                    
      end                                                                       
!------------------------------------------------------------------------------ 
      subroutine indecr(io,iyr,ijul,ihr,idelt,ihrmin,ihrmax)                    
!------------------------------------------------------------------------------ 
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 961014                 INDECR         
! ---            J. Scire, SRC                                                  
!                                                                               
! --- Increment or decrement a date/time by "IDELT" hours                       
! --- (-24 <= IDELT <= 24)                                                      
! --- Allows specification of 0-23 or 1-24 hour clock                           
!                                                                               
! --- INPUTS:                                                                   
!                IO - integer    - Unit number for list file output             
!               IYR - integer    - Input Year                                   
!              IJUL - integer    - Input Julian day                             
!               IHR - integer    - Input hour (ihrmin <= IHR <= ihrmax)         
!             IDELT - integer    - Change in time (hours) -- must be            
!                                  between -24 to +24, inclusive                
!            IHRMIN - integer    - Minimum hour (i.e., either  0 or  1)         
!            IHRMAX - integer    - Maximum hour (i.e., either 23 or 24)         
!                                                                               
! --- OUTPUT:                                                                   
!               IYR - integer    - Year after change of "IDELT" hours           
!              IJUL - integer    - Julian day after change of "IDELT" hours     
!               IHR - integer    - Hour after change of "IDELT" hours           
!                                                                               
! --- INDECR called by:  host subroutines                                       
! --- INDECR calls:      none                                                   
!------------------------------------------------------------------------------ 
!                                                                               
      if(iabs(idelt).gt.24)then                                                 
         write(io,10)'IDELT',iyr,ijul,ihr,idelt,ihrmin,ihrmax                   
10       format(/1x,'ERROR in subr. INDECR -- invalid "',a,'" -- ',            &
     &   ' iyr,ijul,ihr,idelt,ihrmin,ihrmax = ',6i10)                           
         write(*,987)                                                           
987      format(1x,'ERROR in run - see the .LST file')                          
         stop                                                                   
      endif                                                                     
      if(ihr.lt.ihrmin.or.ihr.gt.ihrmax)then                                    
         write(io,10)'IHR',iyr,ijul,ihr,idelt,ihrmin,ihrmax                     
         write(*,987)                                                           
         stop                                                                   
      endif                                                                     
!                                                                               
      if(idelt.lt.0)then                                                        
! ---    idelt is negative                                                      
         ihr=ihr+idelt                                                          
         if(ihr.lt.ihrmin)then                                                  
            ihr=ihr+24                                                          
            ijul=ijul-1                                                         
            if(ijul.lt.1)then                                                   
               iyr=iyr-1                                                        
               if(mod(iyr,4).eq.0)then                                          
                  ijul=366                                                      
               else                                                             
                  ijul=365                                                      
               endif                                                            
            endif                                                               
         endif                                                                  
      else                                                                      
! ---    idelt is positive or zero                                              
         ihr=ihr+idelt                                                          
         if(ihr.gt.ihrmax)then                                                  
            ihr=ihr-24                                                          
            ijul=ijul+1                                                         
            if(mod(iyr,4).eq.0)then                                             
               ndays=366                                                        
            else                                                                
               ndays=365                                                        
            endif                                                               
            if(ijul.gt.ndays)then                                               
               ijul=1                                                           
               iyr=iyr+1                                                        
            endif                                                               
         endif                                                                  
      endif                                                                     
!                                                                               
      return                                                                    
      end                                                                       
!----------------------------------------------------------------------         
      subroutine incrs(io,iyr,ijul,ihr,isec,nsec)                               
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 061020                  INCRS         
!                D. Strimaitis, EARTH TECH                                      
!                                                                               
! --- PURPOSE:  Increment the time and date by "NSEC" seconds                   
!                                                                               
! --- UPDATE                                                                    
! --- V2.54 (061020) from V2.4 (041029) (DGS)                                   
!               - Allow negative increment                                      
!                                                                               
! --- INPUTS:                                                                   
!       IO     - integer - Unit number for list file output                     
!       IYR    - integer - Current year (YYYY)                                  
!       IJUL   - integer - Current Julian day (JJJ)                             
!       IHR    - integer - Current hour (00-23)                                 
!       ISEC   - integer - Current second (0000-3599)                           
!       NSEC   - integer - Time increment (seconds)                             
!       Parameters: IO6                                                         
!                                                                               
! --- OUTPUT:                                                                   
!       IYR    - integer - Updated year                                         
!       IJUL   - integer - Updated Julian day                                   
!       IHR    - integer - Updated hour (00-23)                                 
!       ISEC   - integer - Updated seconds (0000-3599)                          
!                                                                               
! --- INCRS called by: host subroutines                                         
! --- INCRS calls:     INCR                                                     
!----------------------------------------------------------------------         
! --- wangzhm declared
     integer,intent(in) :: io,nsec
     integer,intent(inout) :: iyr,ijul,ihr,isec
! --- wangzhm omp: critical
         !---$---OMP CRITICAL (call_incrs)                                                                                 
      if(nsec.GE.0) then                                                        
! ---    Increment seconds                                                      
         isec=isec+nsec                                                         
         if(isec.GE.3600) then                                                  
            nhrinc=isec/3600                                                    
            isec=MOD(isec,3600)                                                 
            call INCR(io,iyr,ijul,ihr,nhrinc)                                   
         endif                                                                  
                                                                                
      else                                                                      
! ---   Decrement seconds                                                       
         isec=isec+nsec                                                         
         if(isec.LT.0) then                                                     
! ---       Earlier hour                                                        
            ksec=-isec                                                          
            if(ksec.GE.3600) then                                               
! ---          Back up at least 1 hour                                          
               nhrinc=ksec/3600                                                 
               ksec=MOD(ksec,3600)                                              
               nhrinc=-nhrinc                                                   
               call INCR(io,iyr,ijul,ihr,nhrinc)                                
            endif                                                               
            isec=-ksec                                                          
            if(isec.LT.0) then                                                  
! ---          Back up 1 more hour                                              
               nhrinc=-1                                                        
               isec=3600+isec                                                   
               call INCR(io,iyr,ijul,ihr,nhrinc)                                
            endif                                                               
         endif                                                                  
                                                                                
      endif                                                                     
! --- wangzhm omp: critical
         !---$---OMP END CRITICAL (call_incrs)                                                                                 
      return                                                                    
      end
! --- wangzhm declared
!     integer,intent(in) :: io,nsec
!     integer,intent(inout) :: iyr,ijul,ihr,isec
!----------------------------------------------------------------------         
      subroutine deltsec(ndhrb,nsecb,ndhre,nsece,ndelsec)                       
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 041029                DELTSEC         
! ---            D. Strimaitis, Earth Tech                                      
!                                                                               
! --- PURPOSE: Compute the difference (in seconds) between two dates &          
!              times (timeE - timeB)                                            
!                                                                               
! --- INPUTS:                                                                   
!             NDHRB - integer    - Beginning year & hour (YYYYJJJHH)            
!             NSECB - integer    - Beginning second (SSSS)                      
!             NDHRE - integer    - Ending year & hour (YYYYJJJHH)               
!             NSECE - integer    - Ending second (SSSS)                         
!                                                                               
! --- OUTPUT:                                                                   
!           NDELSEC - integer    - Length of interval (seconds)                 
!                                                                               
! --- DELTSEC called by: host subroutines                                       
! --- DELTSEC calls:     DELTT                                                  
!----------------------------------------------------------------------         
!                                                                               
! --- Extract year, Julian day, and hour from date-time variables               
! --- Beginning                                                                 
      j1yr=ndhrb/100000                                                         
      iyyjjj=ndhrb/100                                                          
      j1jul=iyyjjj-j1yr*1000                                                    
      j1hr=ndhrb-iyyjjj*100                                                     
! --- Ending                                                                    
      j2yr=ndhre/100000                                                         
      iyyjjj=ndhre/100                                                          
      j2jul=iyyjjj-j2yr*1000                                                    
      j2hr=ndhre-iyyjjj*100                                                     
                                                                                
! --- Find difference between hours (in seconds)                                
      call DELTT(j1yr,j1jul,j1hr,j2yr,j2jul,j2hr,jdelhr)                        
      ndelsec=jdelhr*3600                                                       
                                                                                
! --- Add difference between seconds                                            
      ndelsec=ndelsec+(nsece-nsecb)                                             
                                                                                
      return                                                                    
      end                                                                       
!----------------------------------------------------------------------         
      subroutine midnite(io,ctrans,iyr,imo,iday,ijul,                          &
     &                             kyr,kmo,kday,kjul)                           
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 041029                MIDNITE         
! ---            D. Strimaitis, Earth Tech                                      
!                                                                               
! --- PURPOSE:  Converts date/time at midnight between day N, 0000              
!               and day N-1, 2400.  Direction is determined by the              
!               CTRANS instruction.                                             
!                                                                               
! --- INPUTS:                                                                   
!            IO - integer      - Unit number for list file output               
!        CTRANS - character    - Instruction 'TO 24h' or 'TO 00h'               
!           IYR - integer      - Year                                           
!           IMO - integer      - Month                                          
!          IDAY - integer      - Day                                            
!          IJUL - integer      - Julian day                                     
!                                                                               
! --- OUTPUT:                                                                   
!           KYR - integer      - Year                                           
!           KMO - integer      - Month                                          
!          KDAY - integer      - Day                                            
!          KJUL - integer      - Julian day                                     
!                                                                               
! --- MIDNITE called by:  host subroutines                                      
! --- MIDNITE calls:      JULDAY, INCR, GRDAY                                   
!----------------------------------------------------------------------         
      character*6 ctrans                                                        
                                                                                
      ierr =0                                                                   
                                                                                
! --- Get Julian day from month/day if needed                                   
      if(ijul.LE.0) call JULDAY(io,iyr,imo,iday,ijul)                           
                                                                                
      kyr=iyr                                                                   
      kmo=imo                                                                   
      kday=iday                                                                 
      kjul=ijul                                                                 
                                                                                
      if(ctrans.EQ.'TO 24h') then                                               
! ---    Convert from 0000 on ijul to 2400 on kjul                              
         ihr=0                                                                  
         nhr=-1                                                                 
         call INCR(io,kyr,kjul,ihr,nhr)                                         
         call GRDAY(io,kyr,kjul,kmo,kday)                                       
      elseif(ctrans.EQ.'TO 00h') then                                           
! ---    Convert from 2400 on ijul to 0000 on kjul                              
         ihr=23                                                                 
         nhr=1                                                                  
         call INCR(io,kyr,kjul,ihr,nhr)                                         
         call GRDAY(io,kyr,kjul,kmo,kday)                                       
      else                                                                      
         ierr=1                                                                 
      endif                                                                     
                                                                                
      if(ierr.eq.1)then                                                         
         write(io,*)                                                            
         write(io,*)'ERROR in SUBR. MIDNITE'                                    
         write(io,*)'Invalid instruction: ',ctrans                              
         write(io,*)'           Expected: TO 24h'                               
         write(io,*)'              OR   : TO 00h'                               
         write(*,*)                                                             
         stop 'Halted in MIDNITE -- see list file.'                             
      endif                                                                     
                                                                                
      return                                                                    
      end                                                                       
!----------------------------------------------------------------------         
      subroutine utcbasr(axtz,xbtz)                                             
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 070327                UTCBASR         
! ---            F.Robe, Earth Tech                                             
!                                                                               
! --- PURPOSE:  Converts character string UTC time zone                         
!               to real base time zone                                          
!                                                                               
! --- V2.55 (070327) from V2.5 (041123) (DGS)                                   
!               - Add RETURN statement                                          
!                                                                               
! --- INPUT:                                                                    
!          AXTZ - char*8    - time zone (international convention:              
!                             relative to UTC/GMT)UTC-HHMM                      
! --- OUTPUT:                                                                   
!          XBTZ - real      - base time zone (old convention: positive          
!                             in North America i.e. opposite to UTC)            
!                                                                               
! --- UTCBASR called by:  host subroutines                                      
! --- UTCBASR calls:      none                                                  
!----------------------------------------------------------------------         
      character*8 axtz                                                          
                                                                                
      read(axtz(4:6),'(i3)')ihr                                                 
      read(axtz(7:8),'(i2)')imin                                                
      if(ihr.lt.0)imin=-imin                                                    
                                                                                
      xbtz=ihr+imin/60.                                                         
                                                                                
! --- Flip sign as base time convention is opposite UTC/GMT                     
      xbtz=-xbtz                                                                
                                                                                
      return                                                                    
      end                                                                       
!----------------------------------------------------------------------         
      subroutine basrutc(xbtz,axtz)                                             
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 070327                BASRUTC         
! ---            F.Robe, Earth Tech                                             
!                                                                               
! --- PURPOSE:  Converts real base time zone  to character string               
!               UTC time zone                                                   
!                                                                               
! --- UPDATE                                                                    
! --- V2.55 (070327) from V2.5 (041123) (DGS)                                   
!               - Fix output format of time zone string for zone=0              
!               - Add RETURN statement                                          
!                                                                               
! --- INPUT:                                                                    
!          XBTZ - real      - base time zone (old convention: positive          
!                             in North America i.e. opposite to UTC)            
                                                                                
! --- OUTPUT:                                                                   
!          AXTZ - real      - time zone (international convention:              
!                             relative to UTC/GMT)UTC-HHMM                      
!                                                                               
! --- BASRUTC called by:  host subroutines                                      
! --- BASRUTC calls:      none                                                  
!----------------------------------------------------------------------         
      character*8 axtz                                                          
                                                                                
      ixbtz=int(xbtz)                                                           
!     convert fractional real to minutes                                        
      imin=(xbtz-ixbtz)*60                                                      
      ixbtz=ixbtz*100+imin                                                      
                                                                                
! --- Define time as "UTC-HHMM" (hours/minutes)                                 
      axtz(1:3)="UTC"                                                           
                                                                                
! --- Flip sign as base time zone is minus UTC zone                             
      if (xbtz.gt.0.) then                                                      
         axtz(4:4)="-"                                                          
      else                                                                      
         axtz(4:4)="+"                                                          
      endif                                                                     
! --- Make sure time zone is written as 4 digits                                
      write(axtz(5:8),'(i4.4)')abs(ixbtz)                                       
                                                                                
      return                                                                    
      end                                                                       
!----------------------------------------------------------------------         
      subroutine filcase(lcfiles,cfile)                                         
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 040330                FILCASE         
! ---            J. Scire, SRC                                                  
!                                                                               
! --- PURPOSE:  Convert all characters within a file name to lower              
!               case (if LCFILES=T) or UPPER CASE (if LCFILES=F).               
!                                                                               
! --- UPDATE                                                                    
! --- V2.2 (950610) to V2.3 (040330)  DGS                                       
!               - Replace filename strings c*70 with c*132                      
!                                                                               
! --- INPUTS:                                                                   
!                                                                               
!         LCFILES - logical - Switch indicating if all characters in the        
!                             filenames are to be converted to lower case       
!                             letters (LCFILES=T) or converted to UPPER         
!                             CASE letters (LCFILES=F).                         
!           CFILE - char*132- Input character string                            
!                                                                               
! --- OUTPUT:                                                                   
!                                                                               
!           CFILE - char*132- Output character string with                      
!                             letters converted                                 
!                                                                               
! --- FILCASE called by:  READFN                                                
! --- FILCASE calls:      none                                                  
!----------------------------------------------------------------------         
!                                                                               
      character*132 cfile                                                       
      character*1 cchar,clc(29),cuc(29)                                         
      logical lcfiles                                                           
!                                                                               
      data clc/'i','n','x','a','e','o','u','b','c','d','f','g','h',            &
     & 'j','k','l','m','p','q','r','s','t','v','w','y','z','-','.',            &
     & '*'/                                                                     
      data cuc/'I','N','X','A','E','O','U','B','C','D','F','G','H',            &
     & 'J','K','L','M','P','Q','R','S','T','V','W','Y','Z','-','.',            &
     & '*'/                                                                     
!                                                                               
      if(lcfiles)then                                                           
!                                                                               
! ---    Convert file name to lower case letters                                
         do i=1,132                                                             
            cchar=cfile(i:i)                                                    
!                                                                               
            do j=1,29                                                           
               if(cchar.eq.cuc(j))then                                          
                  cfile(i:i)=clc(j)                                             
                  go to 52                                                      
               endif                                                            
            enddo                                                               
52          continue                                                            
         enddo                                                                  
      else                                                                      
!                                                                               
! ---    Convert file name to UPPER CASE letters                                
         do i=1,132                                                             
            cchar=cfile(i:i)                                                    
!                                                                               
            do j=1,29                                                           
               if(cchar.eq.clc(j))then                                          
                  cfile(i:i)=cuc(j)                                             
                  go to 62                                                      
               endif                                                            
            enddo                                                               
62          continue                                                            
         enddo                                                                  
      endif                                                                     
!                                                                               
      return                                                                    
      end                                                                       
!----------------------------------------------------------------------         
      subroutine readin(cvdic,ivleng,ivtype,ioin,ioout,lecho,                  &
     & i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,         &
     & i19,i20,i21,i22,i23,i24,i25,i26,i27,i28,i29,i30,i31,i32,i33,i34,        &
     & i35,i36,i37,i38,i39,i40,i41,i42,i43,i44,i45,i46,i47,i48,i49,i50,        &
     & i51,i52,i53,i54,i55,i56,i57,i58,i59,i60)                                 
!----------------------------------------------------------------------         
! *** Change number of characters in line from 150 to 200 ***                   
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 110225                 READIN         
!                J. Scire                                                       
!                                                                               
! --- PURPOSE:  Read one input group of the free formatted control              
!               file -- allows comments within the input file --                
!               ignores all text except that within delimiters                  
!                                                                               
! ---           NOTE:  All variables (real, integer, logical,                   
!                      or character) must be 4 bytes                            
! ---           NOTE:  Character*4 array uses only one character                
!                      per word -- it must be dimensioned large                 
!                      enough to accommodate the number of characters           
!                      in the variable field                                    
!                                                                               
! --- UPDATE                                                                    
! --- V2.58 (110225) from V2.57 (090202) (DGS)                                  
!               - Add IVTYPE=5 (char*4 array with commas retained               
!                 as delimiters for parsing)                                    
! --- V2.57 (090202) from V2.52 (060519) (DGS)                                  
!               - Increase max line length from 150 to 200                      
!                 (requires MXCOL=200)                                          
! --- V2.52 (060519) from V2.3 (040330) (DGS)                                   
!               - Search for '=' beyond position 14 because blanks are          
!                 not automatically removed within string                       
! --- V2.3 (040330) from V2.1 (030402) (DGS)                                    
!               - Preserve spaces within character variables                    
! --- V2.1 (030402) from V2.0 (000602) (DGS)                                    
!               - Split DEBLNK action (removes ' ', '+') into                   
!                 DEBLNK and DEPLUS(new)                                        
!                                                                               
!                                                                               
! --- INPUTS:                                                                   
!                                                                               
!     CVDIC(mxvar) - character*12 array - Variable dictionary                   
!                                         containing up to "MXVAR"              
!                                         variable names                        
!    IVLENG(mxvar) - integer array      - Dimension of each variable            
!                                         (dim. of scalars = 1)                 
!    IVTYPE(mxvar) - integer array      - Type of each variable                 
!                                           1 = real,                           
!                                           2 = integer,                        
!                                           3 = logical,                        
!                                           4 = character*4                     
!                                           5 = character*4 with commas         
!             IOIN - integer            - Fortran unit of control file          
!                                         input                                 
!            IOOUT - integer            - Fortran unit of list file             
!                                         output                                
!            LECHO - logical            - Control variable determining          
!                                         if input data are echoed to           
!                                         list file (IOOUT)                     
!        Parameters: MXVAR, MXCOL                                               
!                                                                               
! --- OUTPUT:                                                                   
!                                                                               
!      I1, I2, ... - integer arrays     - Variables being read                  
!                    (integer array locally, but can be a real,                 
!                     integer, logical, or character*4 array in                 
!                     the calling routine)                                      
!                                                                               
! --- READIN called by:  host subroutines                                       
! --- READIN calls:      DEBLNK, ALTONU, SETVAR, ALLCAP, DEPLUS,                
!                        TRIGHT, TLEFT                                          
!----------------------------------------------------------------------         
!                                                                               
! --- Include parameter statements                                              
      include 'params.cal'                                                      
!                                                                               
      integer*4 i1(*),i2(*),i3(*),i4(*),i5(*),i6(*),i7(*),i8(*),i9(*),         &
     & i10(*),i11(*),i12(*),i13(*),i14(*),i15(*),i16(*),i17(*),i18(*),         &
     & i19(*),i20(*),i21(*),i22(*),i23(*),i24(*),i25(*),i26(*),i27(*),         &
     & i28(*),i29(*),i30(*),i31(*),i32(*),i33(*),i34(*),i35(*),i36(*),         &
     & i37(*),i38(*),i39(*),i40(*),i41(*),i42(*),i43(*),i44(*),i45(*),         &
     & i46(*),i47(*),i48(*),i49(*),i50(*),i51(*),i52(*),i53(*),i54(*),         &
     & i55(*),i56(*),i57(*),i58(*),i59(*),i60(*)                                
      integer*4 ivleng(mxvar),jdex(mxvar),ivtype(mxvar)                         
!                                                                               
      logical*4 lv                                                              
      logical lecho                                                             
!                                                                               
      character*12 cvdic(mxvar),cvar,cblank                                     
      character*4 cv(mxcol)                                                     
      character*1 cstor1(mxcol),cstor2(mxcol)                                   
! --- Intermediate scratch arrays                                               
      character*1 cstor3(mxcol),cstor4(mxcol)                                   
      character*1 cdelim,ceqls,ce,cn,cd,comma,cblnk                             
!                                                                               
      data cblank/'            '/                                               
      data cdelim/'!'/,ceqls/'='/,ce/'E'/,cn/'N'/,cd/'D'/,comma/','/            
      data cblnk/' '/                                                           
!                                                                               
      ilim2=99                                                                  
      do 2 i=1,mxvar                                                            
      jdex(i)=1                                                                 
2     continue                                                                  
!                                                                               
! --- begin loop over lines                                                     
!                                                                               
! --- read a line of input                                                      
5     continue                                                                  
      read(ioin,10)cstor1                                                       
10    format(200a1)                                                             
      if(lecho)write(ioout,7)cstor1                                             
7     format(1x,200a1)                                                          
!                                                                               
! --- check if this is a continuation line                                      
      if(ilim2.gt.0)go to 16                                                    
!                                                                               
! --- continuation line -- find the second delimiter                            
      do 12 i=1,mxcol                                                           
      if(cstor1(i).eq.cdelim)then                                               
         ilim2=i                                                                
         go to 14                                                               
      endif                                                                     
12    continue                                                                  
14    continue                                                                  
      il2=ilim2                                                                 
      if(il2.eq.0)il2=mxcol                                                     
!                                                                               
! --- Trim blanks from left and right sides of string within delimiters         
! -----------------------                                                       
!c --- remove blank characters from string within delimiters                    
!      call deblnk(cstor1,1,il2,cstor2,nlim)                                    
!c --- Remove '+' characters as well (is this needed?)                          
!      if(nlim.gt.0) then                                                       
!         do k=1,mxcol                                                          
!            cstor3(k)=cstor2(k)                                                
!         enddo                                                                 
!         il3=nlim                                                              
!         call deplus(cstor3,1,il3,cstor2,nlim)                                 
!      endif                                                                    
! -----------------------                                                       
! --- Remove blank characters on right side                                     
      call TRIGHT(cstor1,1,il2,cstor2,nlim)                                     
! --- Remove blank characters on left side                                      
      if(nlim.gt.0) then                                                        
         do k=1,mxcol                                                           
            cstor3(k)=cstor2(k)                                                 
         enddo                                                                  
         il3=nlim                                                               
         call TLEFT(cstor3,1,il3,cstor2,nlim)                                   
      endif                                                                     
! -----------------------                                                       
      icom=0                                                                    
!                                                                               
! --- convert lower case letters to upper case                                  
      call allcap(cstor2,nlim)                                                  
      go to 55                                                                  
!                                                                               
16    continue                                                                  
      ibs=1                                                                     
!                                                                               
! --- begin loop over delimiter pairs                                           
17    continue                                                                  
      if(ibs.ge.mxcol)go to 5                                                   
!                                                                               
! --- find location of delimiters                                               
      do 20 i=ibs,mxcol                                                         
      if(cstor1(i).eq.cdelim)then                                               
         ilim1=i                                                                
         if(ilim1.eq.mxcol)go to 22                                             
         ip1=ilim1+1                                                            
         do 18 j=ip1,mxcol                                                      
         if(cstor1(j).eq.cdelim)then                                            
            ilim2=j                                                             
            go to 22                                                            
         endif                                                                  
18       continue                                                               
!                                                                               
! ---    second delimiter not on this line                                      
         ilim2=0                                                                
         go to 22                                                               
      endif                                                                     
20    continue                                                                  
!                                                                               
! --- no delimiters found -- skip line and read next line of text               
      go to 5                                                                   
22    continue                                                                  
      ibs=ilim2+1                                                               
      if(ilim2.eq.0)ibs=mxcol+1                                                 
!                                                                               
! --- Trim blanks from left and right sides of string within delimiters         
! -----------------------                                                       
!c --- remove blanks from string within delimiters                              
!      il2=ilim2                                                                
!      if(il2.eq.0)il2=mxcol                                                    
!      call deblnk(cstor1,ilim1,il2,cstor2,nlim)                                
!c --- Remove '+' characters as well (is this needed?)                          
!      if(nlim.gt.0) then                                                       
!         do k=1,mxcol                                                          
!            cstor3(k)=cstor2(k)                                                
!         enddo                                                                 
!         il3=nlim                                                              
!         call deplus(cstor3,1,il3,cstor2,nlim)                                 
!      endif                                                                    
! -----------------------                                                       
      il2=ilim2                                                                 
      if(il2.eq.0)il2=mxcol                                                     
! --- Remove blank characters on right side                                     
      call TRIGHT(cstor1,ilim1,il2,cstor2,nlim)                                 
! --- Remove blank characters on left side                                      
      if(nlim.gt.0) then                                                        
         do k=1,mxcol                                                           
            cstor3(k)=cstor2(k)                                                 
         enddo                                                                  
         il3=nlim                                                               
         call TLEFT(cstor3,1,il3,cstor2,nlim)                                   
      endif                                                                     
! -----------------------                                                       
!                                                                               
! --- convert lower case letters to upper case                                  
      call allcap(cstor2,nlim)                                                  
!                                                                               
! --- search for equals sign (cstor2(1) is delimiter; cstor2(2) is              
! --- first letter of variable; cstor2(3) is earliest '=' can occur)            
! --- (060519)  Search entire string as now there may be blanks before '='      
!      do 30 i=3,14                                                             
      do 30 i=3,nlim                                                            
      if(cstor2(i).eq.ceqls)then                                                
         ieq=i                                                                  
         go to 32                                                               
      endif                                                                     
30    continue                                                                  
!                                                                               
! --- "END" within delimiters signifies the end of the read for                 
! --- this input group                                                          
      if(cstor2(2).eq.ce.and.cstor2(3).eq.cn.and.cstor2(4).eq.cd)return         
      write(ioout,31)(cstor2(n),n=1,nlim)                                       
31    format(/1x,'ERROR IN SUBR. READIN -- Error in input data -- '/           &
     & 1x,'Variable too long (Equals sign not found in string) -- ',           &
     & 'CSTOR2 = ',200a1)                                                       
      write(*,*)                                                                
      stop 'Halted in READIN -- see list file.'                                 
!                                                                               
! --- CVAR is character*12 variable name                                        
32    continue                                                                  
      cvar=cblank                                                               
      ieqm1=ieq-1                                                               
! --- Grab string to left of '=', and remove blanks                             
      call deblnk(cstor2,1,ieqm1,cstor3,keqm1)                                  
! --- Pass string to variable name                                              
      do 40 i=2,keqm1                                                           
      il=i-1                                                                    
      cvar(il:il)=cstor3(i)                                                     
40    continue                                                                  
!                                                                               
! --- find the variable name in the variable dictionary                         
      do 50 i=1,mxvar                                                           
      if(cvar.eq.cvdic(i))then                                                  
         nvar=i                                                                 
         go to 52                                                               
      endif                                                                     
50    continue                                                                  
      write(ioout,51)cvar,(cvdic(n),n=1,mxvar)                                  
51    format(/1x,'ERROR IN SUBR. READIN -- Error in input data -- '/           &
     & 1x,'Variable not found in variable dictionary'/                         &
     & 1x,'Variable: ',a12/                                                    &
     & 1x,'Variable Dictionary: ',9(a12,1x)/                                   &
     & 10(22x,9(a12,1x)/))                                                      
      write(*,*)                                                                
      stop 'Halted in READIN -- see list file.'                                 
!                                                                               
52    continue                                                                  
! --- Assign current variable type                                              
      itype=ivtype(nvar)                                                        
!                                                                               
! --- Check for invalid value of variable type                                  
      if(itype.le.0.or.itype.ge.6)then                                          
         write(ioout,53)itype,nvar,ivtype(nvar),cvdic(nvar)                     
53       format(/1x,'ERROR IN SUBR. READIN -- Error in input data -- '/        &
     &   1x,'Invalid value of variable type -- ITYPE must be 1, 2, 3, ',       &
     &   '4, or 5'/1x,'ITYPE = ',i10/1x,'NVAR = ',i10/1x,                      &
     &   'IVTYPE(nvar) = ',i10/1x,'CVDIC(nvar) = ',a12)                         
      write(*,*)                                                                
      stop 'Halted in READIN -- see list file.'                                 
      endif                                                                     
!                                                                               
! --- search for comma                                                          
      icom=ieq                                                                  
!                                                                               
! --- beginning of loop over values within delimiters                           
55    continue                                                                  
      ivb=icom+1                                                                
!                                                                               
! --- if reaches end of line, read next line                                    
      if(ivb.gt.nlim)go to 5                                                    
      do 60 i=ivb,nlim                                                          
      if(cstor2(i).eq.comma)then                                                
         icom=i                                                                 
         go to 64                                                               
      endif                                                                     
60    continue                                                                  
!                                                                               
! --- no comma found                                                            
      icom=0                                                                    
      ive=nlim-1                                                                
!                                                                               
! --- comma between last value and delimiter is allowed                         
      if(cstor2(ivb).eq.cdelim.and.cstor2(ive).eq.comma)go to 17                
!                                                                               
! --- if no comma & last non-blank character is not a delimiter,                
! --- then the input is in error                                                
      if(cstor2(nlim).eq.cdelim)go to 66                                        
      write(ioout,63)cstor1                                                     
63    format(/1x,'ERROR IN SUBR. READIN -- Error in input data -- '/           &
     & 1x,'If a string within delimiters covers more than one line, ',         &
     & 'the last character in the line must be a comma'/                       &
     & 1x,'Input line: ',200a1)                                                 
      write(*,*)                                                                
      stop 'Halted in READIN -- see list file.'                                 
64    continue                                                                  
!                                                                               
! --- value of variable is contained in elements IVB to IVE of                  
! --- CSTOR2 array                                                              
! --- Include comma for variable type 5 (character array) so that it            
! --- can be used outside of READIN to parse the array values from the          
! --- single string that is returned                                            
      if(itype.EQ.5) then                                                       
         ive=icom                                                               
      else                                                                      
         ive=icom-1                                                             
      endif                                                                     
66    continue                                                                  
!      ncar=ive-ivb+1                                                           
      index=jdex(nvar)                                                          
!                                                                               
! --- Convert character string to numeric or logical value                      
!     (if ITYPE = 1,2, or 3) -- If 4 or 5 transfer characters to the            
!     work array CV)                                                            
                                                                                
! --- Remove all blanks from variable string if type is numeric or              
! --- logical;  otherwise, trim left and right side of string                   
      if(itype.LT.4) then                                                       
         call deblnk(cstor2,ivb,ive,cstor4,nv)                                  
! ---    Remove '+' characters as well (is this needed?)                        
         if(nv.gt.0) then                                                       
            do k=1,mxcol                                                        
               cstor3(k)=cstor4(k)                                              
            enddo                                                               
            il3=nv                                                              
            call deplus(cstor3,1,il3,cstor4,nv)                                 
         endif                                                                  
         call altonu(ioout,cstor4(1),nv,itype,irep,rlno,ino,lv,cv)              
      else                                                                      
! ---    Pass variable string into cstor4                                       
         nv=ive-ivb+1                                                           
         do k=1,nv                                                              
            cstor4(k)=cstor2(ivb+k-1)                                           
         enddo                                                                  
         do k=nv+1,mxcol                                                        
            cstor4(k)=cblnk                                                     
         enddo                                                                  
! ---    Remove blank characters on right side of character variable            
! ---    if last character is either a blank or comma                           
         if(cstor4(nv).EQ.cblnk .OR.                                           &
     &      cstor4(nv).EQ.comma) call TRIGHT(cstor2,ivb,ive,cstor4,nv)          
! ---    Remove blank characters on left side of character variable             
         if(nv.GT.0 .AND. cstor4(1).EQ.cblnk) then                              
            do k=1,mxcol                                                        
               cstor3(k)=cstor4(k)                                              
            enddo                                                               
            il3=nv                                                              
            call TLEFT(cstor3,1,il3,cstor4,nv)                                  
         endif                                                                  
         call altonu(ioout,cstor4(1),nv,itype,irep,rlno,ino,lv,cv)              
      endif                                                                     
!                                                                               
! --- check that array bounds are not exceeded                                  
      if(index+irep-1.gt.ivleng(nvar))go to 201                                 
!                                                                               
      go to (101,102,103,104,105,106,107,108,109,110,111,112,113,114,          &
     & 115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,        &
     & 131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,        &
     & 147,148,149,150,151,152,153,154,155,156,157,158,159,160),nvar            
!                                                                               
! --- code currently set up to handle up to 60 variables/source group           
      write(ioout,71)nvar,(cstor2(n),n=1,nlim)                                  
71    format(/1x,'ERROR IN SUBR. READIN -- Current code ',                     &
     & 'configuration allows up to 60 variables per source group'/             &
     & 1x,'No. variables (NVAR) = ',i10/                                       &
     & 1x,'Input data (CSTOR2)  = ',200a1)                                      
      write(*,*)                                                                
      stop 'Halted in READIN -- see list file.'                                 
!                                                                               
! --- transfer value into output variable                                       
101   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i1(index),i1(index),               &
     & i1(index),i1(index))                                                     
      go to 161                                                                 
102   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i2(index),i2(index),               &
     & i2(index),i2(index))                                                     
      go to 161                                                                 
103   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i3(index),i3(index),               &
     & i3(index),i3(index))                                                     
      go to 161                                                                 
104   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i4(index),i4(index),               &
     & i4(index),i4(index))                                                     
      go to 161                                                                 
105   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i5(index),i5(index),               &
     & i5(index),i5(index))                                                     
      go to 161                                                                 
106   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i6(index),i6(index),               &
     & i6(index),i6(index))                                                     
      go to 161                                                                 
107   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i7(index),i7(index),               &
     & i7(index),i7(index))                                                     
      go to 161                                                                 
108   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i8(index),i8(index),               &
     & i8(index),i8(index))                                                     
      go to 161                                                                 
109   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i9(index),i9(index),               &
     & i9(index),i9(index))                                                     
      go to 161                                                                 
110   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i10(index),i10(index),             &
     & i10(index),i10(index))                                                   
      go to 161                                                                 
111   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i11(index),i11(index),             &
     & i11(index),i11(index))                                                   
      go to 161                                                                 
112   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i12(index),i12(index),             &
     & i12(index),i12(index))                                                   
      go to 161                                                                 
113   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i13(index),i13(index),             &
     & i13(index),i13(index))                                                   
      go to 161                                                                 
114   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i14(index),i14(index),             &
     & i14(index),i14(index))                                                   
      go to 161                                                                 
115   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i15(index),i15(index),             &
     & i15(index),i15(index))                                                   
      go to 161                                                                 
116   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i16(index),i16(index),             &
     & i16(index),i16(index))                                                   
      go to 161                                                                 
117   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i17(index),i17(index),             &
     & i17(index),i17(index))                                                   
      go to 161                                                                 
118   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i18(index),i18(index),             &
     & i18(index),i18(index))                                                   
      go to 161                                                                 
119   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i19(index),i19(index),             &
     & i19(index),i19(index))                                                   
      go to 161                                                                 
120   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i20(index),i20(index),             &
     & i20(index),i20(index))                                                   
      go to 161                                                                 
121   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i21(index),i21(index),             &
     & i21(index),i21(index))                                                   
      go to 161                                                                 
122   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i22(index),i22(index),             &
     & i22(index),i22(index))                                                   
      go to 161                                                                 
123   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i23(index),i23(index),             &
     & i23(index),i23(index))                                                   
      go to 161                                                                 
124   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i24(index),i24(index),             &
     & i24(index),i24(index))                                                   
      go to 161                                                                 
125   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i25(index),i25(index),             &
     & i25(index),i25(index))                                                   
      go to 161                                                                 
126   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i26(index),i26(index),             &
     & i26(index),i26(index))                                                   
      go to 161                                                                 
127   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i27(index),i27(index),             &
     & i27(index),i27(index))                                                   
      go to 161                                                                 
128   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i28(index),i28(index),             &
     & i28(index),i28(index))                                                   
      go to 161                                                                 
129   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i29(index),i29(index),             &
     & i29(index),i29(index))                                                   
      go to 161                                                                 
130   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i30(index),i30(index),             &
     & i30(index),i30(index))                                                   
      go to 161                                                                 
131   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i31(index),i31(index),             &
     & i31(index),i31(index))                                                   
      go to 161                                                                 
132   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i32(index),i32(index),             &
     & i32(index),i32(index))                                                   
      go to 161                                                                 
133   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i33(index),i33(index),             &
     & i33(index),i33(index))                                                   
      go to 161                                                                 
134   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i34(index),i34(index),             &
     & i34(index),i34(index))                                                   
      go to 161                                                                 
135   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i35(index),i35(index),             &
     & i35(index),i35(index))                                                   
      go to 161                                                                 
136   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i36(index),i36(index),             &
     & i36(index),i36(index))                                                   
      go to 161                                                                 
137   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i37(index),i37(index),             &
     & i37(index),i37(index))                                                   
      go to 161                                                                 
138   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i38(index),i38(index),             &
     & i38(index),i38(index))                                                   
      go to 161                                                                 
139   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i39(index),i39(index),             &
     & i39(index),i39(index))                                                   
      go to 161                                                                 
140   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i40(index),i40(index),             &
     & i40(index),i40(index))                                                   
      go to 161                                                                 
141   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i41(index),i41(index),             &
     & i41(index),i41(index))                                                   
      go to 161                                                                 
142   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i42(index),i42(index),             &
     & i42(index),i42(index))                                                   
      go to 161                                                                 
143   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i43(index),i43(index),             &
     & i43(index),i43(index))                                                   
      go to 161                                                                 
144   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i44(index),i44(index),             &
     & i44(index),i44(index))                                                   
      go to 161                                                                 
145   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i45(index),i45(index),             &
     & i45(index),i45(index))                                                   
      go to 161                                                                 
146   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i46(index),i46(index),             &
     & i46(index),i46(index))                                                   
      go to 161                                                                 
147   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i47(index),i47(index),             &
     & i47(index),i47(index))                                                   
      go to 161                                                                 
148   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i48(index),i48(index),             &
     & i48(index),i48(index))                                                   
      go to 161                                                                 
149   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i49(index),i49(index),             &
     & i49(index),i49(index))                                                   
      go to 161                                                                 
150   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i50(index),i50(index),             &
     & i50(index),i50(index))                                                   
      go to 161                                                                 
151   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i51(index),i51(index),             &
     & i51(index),i51(index))                                                   
      go to 161                                                                 
152   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i52(index),i52(index),             &
     & i52(index),i52(index))                                                   
      go to 161                                                                 
153   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i53(index),i53(index),             &
     & i53(index),i53(index))                                                   
      go to 161                                                                 
154   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i54(index),i54(index),             &
     & i54(index),i54(index))                                                   
      go to 161                                                                 
155   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i55(index),i55(index),             &
     & i55(index),i55(index))                                                   
      go to 161                                                                 
156   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i56(index),i56(index),             &
     & i56(index),i56(index))                                                   
      go to 161                                                                 
157   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i57(index),i57(index),             &
     & i57(index),i57(index))                                                   
      go to 161                                                                 
158   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i58(index),i58(index),             &
     & i58(index),i58(index))                                                   
      go to 161                                                                 
159   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i59(index),i59(index),             &
     & i59(index),i59(index))                                                   
      go to 161                                                                 
160   continue                                                                  
      call setvar(itype,irep,rlno,ino,lv,cv,i60(index),i60(index),             &
     & i60(index),i60(index))                                                   
!                                                                               
161   continue                                                                  
      jdex(nvar)=jdex(nvar)+irep                                                
!                                                                               
! --- continue reading values for this array until array is filled              
! --- or delimiter is reached                                                   
      if(icom.ne.0.and.jdex(nvar).le.ivleng(nvar))go to 55                      
      go to 17                                                                  
201   continue                                                                  
      iatt=index+irep-1                                                         
      write(ioout,202)cvdic(nvar),ivleng(nvar),iatt,cstor1                      
202   format(/1x,'ERROR IN SUBR. READIN -- Error in input data',               &
     & 1x,'Array bounds exceeded -- Variable: ',a12,3x,' Declared ',           &
     & 'dimension = ',i8/1x,'Input attempted to element ',i8/1x,               &
     & 'Input line: ',200a1)                                                    
      write(*,*)                                                                
      stop 'Halted in READIN -- see list file.'                                 
      end                                                                       
!----------------------------------------------------------------------         
      subroutine altonu(ioout,alp,ncar,itype,irep,rlno,ino,lv,cv)               
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 110225                 ALTONU         
! ---            J. Scire                                                       
!                                                                               
! --- PURPOSE:  Convert a character string into a real, integer or              
!               logical variable -- also compute the repetition factor          
!               for the variable                                                
!                                                                               
! --- UPDATES                                                                   
! --- V2.58 (110225) from V2.56 (080407) (DGS)                                  
!               - Add ITYPE=5 (char*4 array with commas retained                
!                 as delimiters for parsing)                                    
! --- V2.56 (080407) from V1.0 (000602) (DGS)                                   
!               - Treat case in which exponential notation is used              
!                 without a decimal point.  Pointer had been left at            
!                 'zero' which placed the decimal location in front of          
!                 a number so that 2e02 became 0.2e02 instead of 2.0e02         
!               - Trap case where no number appears in front the E or D         
!                 in exponential notation                                       
!                                                                               
! ---             000602  (DGS): add message to "stop"                          
!                                                                               
! --- INPUTS:                                                                   
!            IOOUT - integer           - Fortran unit of list file              
!                                        output                                 
!        ALP(ncar) - character*1 array - Characters to be converted             
!             NCAR - integer           - Number of characters                   
!            ITYPE - integer           - Type of each variable                  
!                                           1 = real,                           
!                                           2 = integer,                        
!                                           3 = logical,                        
!                                           4 = character*4                     
!                                           5 = character*4 with commas         
!                                                                               
!       Parameter:   MXCOL                                                      
!                                                                               
! --- OUTPUT:                                                                   
!             IREP - integer           - Repetition factor for value            
!             RLNO - real              - Real variable produced from            
!                                        character string                       
!              INO - integer           - Integer variable produced from         
!                                        character string                       
!               LV - logical*4         - Logical variable produced from         
!                                        character string                       
!        CV(mxcol) - character*4       - Character*4 variable produced          
!                                        from character string                  
!                                        (NOTE: Only 1 (NOT 4)                  
!                                        character(s) per word)                 
!                                                                               
! --- ALTONU called by:  READIN                                                 
! --- ALTONU calls:      none                                                   
!----------------------------------------------------------------------         
!                                                                               
! --- Include parameter statements                                              
      include 'params.cal'                                                      
!                                                                               
      real*8 rno,xmult,ten                                                      
      integer num2(mxcol)                                                       
      logical*4 lv                                                              
      character*4 cv(mxcol)                                                     
      character*1 alp(ncar),alpsv,ad(17),astar,adec                             
!                                                                               
      data ad/'0','1','2','3','4','5','6','7','8','9','-',                   &   
! ---   num2 = 0   1   2   3   4   5   6   7   8   9  11                       
     &        '*','.','E','D','T','F'/                                          
! ---   num2 = 12  13  14  15  16  17                                           
      data astar/'*'/,adec/'.'/,ten/10.0d0/                                     
!                                                                               
! --- If dealing with a character*4 variable, transfer characters               
!     into the work array CV (ONE character per 4-byte word)                    
      if(itype.eq.4 .OR. itype.eq.5)then                                        
         do 5 i=1,ncar                                                          
         cv(i)(1:1)=alp(i)                                                      
5        continue                                                               
!                                                                               
! ---    NOTE: Repetition factor refers to the number of                        
!              characters in the field, if ITYPE = 4, 5                         
         irep=ncar                                                              
         return                                                                 
      endif                                                                     
!                                                                               
! --- Convert character array elements into numeric codes                       
      do 30 i=1,ncar                                                            
      alpsv=alp(i)                                                              
      do 20 j=1,17                                                              
      if(alpsv.eq.ad(j))then                                                    
         num2(i)=j                                                              
         if(j.lt.11)num2(i)=j-1                                                 
         go to 30                                                               
      endif                                                                     
20    continue                                                                  
      write(ioout,21)(alp(n),n=1,ncar)                                          
21    format(/1x,'ERROR IN SUBR. ALTONU -- Unrecognizable character ',         &
     & 'in input -- Character string (ALP) = ',15a1)                            
      write(*,*)                                                                
      stop 'Halted in ALTONU -- see list file.'                                 
30    continue                                                                  
!                                                                               
! --- Locally classify variable type (1=real, 2=integer, 3=logical)             
      do 40 i=1,ncar                                                            
      if(num2(i).le.12)go to 40                                                 
      if(num2(i).ge.16)then                                                     
!                                                                               
! ---    logical variable ("T", "F")                                            
         jtype=3                                                                
         go to 41                                                               
      else                                                                      
!                                                                               
! ---    real variable (".", "E", "D")                                          
         jtype=1                                                                
         go to 41                                                               
      endif                                                                     
40    continue                                                                  
!                                                                               
! --- integer variable                                                          
      jtype=2                                                                   
41    continue                                                                  
!                                                                               
! --- determine if repetition factor "*" is used                                
      do 50 i=1,ncar                                                            
      if(alp(i).eq.astar)then                                                   
         istar=i                                                                
         go to 51                                                               
      endif                                                                     
50    continue                                                                  
      istar=0                                                                   
51    continue                                                                  
      if(istar.ne.0)go to 400                                                   
      irep=1                                                                    
      go to (101,201,301),jtype                                                 
      write(ioout,55)jtype,(alp(n),n=1,ncar)                                    
55    format(/1x,'ERROR IN SUBR. ALTONU -- JTYPE must be 1, 2, or 3 ',         &
     & '-- JTYPE = ',i3/3x,'Text string (ALP) = ',15a1)                         
      write(*,*)                                                                
      stop 'Halted in ALTONU -- see list file.'                                 
!                                                                               
! --------------------------------------------------------------------          
! --- REAL number w/o "*"                                                       
! --------------------------------------------------------------------          
! --- Determine sign -- ISTAR is position of array containing "*"               
!                       (ISTAR = 0 if no repetition factor)                     
101   continue                                                                  
      if(num2(1+istar).eq.11)then                                               
         isgn=-1                                                                
         istart=istar+2                                                         
      else                                                                      
         isgn=1                                                                 
         istart=istar+1                                                         
      endif                                                                     
!                                                                               
! --- Locate decimal point                                                      
      idec=0                                                                    
      do 109 i=istart,ncar                                                      
      if(alp(i).eq.adec)then                                                    
         if(idec.eq.0)then                                                      
            idec=i                                                              
            go to 109                                                           
         endif                                                                  
!                                                                               
! ---    More than one decimal point found                                      
         write(ioout,120)(alp(n),n=1,ncar)                                      
120      format(/1x,'ERROR IN SUBR. ALTONU -- Invalid real variable ',         &
     &   'entry'/5x,'Input text (ALP) = ',15a1)                                 
         write(*,*)                                                             
         stop 'Halted in ALTONU -- see list file.'                              
      endif                                                                     
109   continue                                                                  
!                                                                               
! --- Search for E or D                                                         
      do 110 i=istart,ncar                                                      
      if(num2(i).eq.14.or.num2(i).eq.15)then                                    
         istop=i-1                                                              
         go to 111                                                              
      endif                                                                     
110   continue                                                                  
      istop=ncar                                                                
111   continue                                                                  
                                                                                
! --- 080407 Update:                                                            
! --- Correct for missing decimal point before decoding                         
      if(idec.EQ.0) idec=istop+1                                                
! --- Trap missing number in front of E,D                                       
      if(istop.LT.1 .OR. istart.GT.istop) then                                  
         write(ioout,120)(alp(n),n=1,ncar)                                      
         write(*,*)                                                             
         write(*,*)'Missing number!'                                            
         stop 'Halted in ALTONU -- see list file.'                              
      endif                                                                     
!                                                                               
! --- Convert integer numerics to real number                                   
      rno=0.0                                                                   
      do 130 i=istart,istop                                                     
      if(i.eq.idec)go to 130                                                    
      if(num2(i).ge.10)then                                                     
         write(ioout,120)(alp(n),n=1,ncar)                                      
         write(*,*)                                                             
         stop 'Halted in ALTONU -- see list file.'                              
      endif                                                                     
      iexp=idec-i                                                               
      if(iexp.gt.0)iexp=iexp-1                                                  
      xmult=1.0                                                                 
      if(iexp.ne.0)xmult=ten**iexp                                              
      rno=rno+xmult*num2(i)                                                     
                                                                                
130   continue                                                                  
!                                                                               
! --- Account for minus sign (if present)                                       
      rno=isgn*rno                                                              
      rlno=rno                                                                  
! --- Also set integer variable in case of improper input                       
      if(rlno.lt.0.0)then                                                       
         ino=rlno-0.0001                                                        
      else                                                                      
         ino=rlno+0.0001                                                        
      endif                                                                     
      if(istop.eq.ncar)return                                                   
!                                                                               
! --- Find exponent (istop+1 is position in array containing E or D)            
      isgn=1                                                                    
      istart=istop+2                                                            
      if(num2(istart).ne.11)go to 135                                           
      isgn=-1                                                                   
      istart=istart+1                                                           
135   continue                                                                  
      if(istart.gt.ncar)then                                                    
         write(ioout,120)(alp(n),n=1,ncar)                                      
         write(*,*)                                                             
         stop 'Halted in ALTONU -- see list file.'                              
      endif                                                                     
      rexp=0.0                                                                  
      do 140 i=istart,ncar                                                      
      if(num2(i).ge.10)then                                                     
         write(ioout,120)(alp(n),n=1,ncar)                                      
         write(*,*)                                                             
         stop 'Halted in ALTONU -- see list file.'                              
      endif                                                                     
      iexp=ncar-i                                                               
      xmult=1.0                                                                 
      if(iexp.ne.0)xmult=ten**iexp                                              
      rexp=rexp+xmult*num2(i)                                                   
140   continue                                                                  
      xmult=1.0                                                                 
      if(rexp.ne.0.0)xmult=ten**(isgn*rexp)                                     
      rno=rno*xmult                                                             
      rlno=rno                                                                  
!                                                                               
! --- Also set integer variable in case of improper input                       
      if(rlno.lt.0.0)then                                                       
         ino=rlno-0.0001                                                        
      else                                                                      
         ino=rlno+0.0001                                                        
      endif                                                                     
      return                                                                    
!                                                                               
! --------------------------------------------------------------------          
! --- INTEGER variables                                                         
! --------------------------------------------------------------------          
201   continue                                                                  
      if(num2(1+istar).ne.11)go to 228                                          
      isgn=-1                                                                   
      istart=istar+2                                                            
      go to 229                                                                 
228   continue                                                                  
      isgn=1                                                                    
      istart=istar+1                                                            
229   continue                                                                  
      ino=0                                                                     
      do 230 i=istart,ncar                                                      
      if(num2(i).ge.10)go to 208                                                
      iexp=ncar-i                                                               
      xmult=1.0                                                                 
      if(iexp.ne.10)xmult=ten**iexp                                             
      ino=ino+xmult*num2(i)+0.5                                                 
230   continue                                                                  
      ino=isgn*ino                                                              
!                                                                               
! --- Also set real variable in case of improper input                          
      rlno=ino                                                                  
      return                                                                    
208   continue                                                                  
      write(ioout,220)(alp(n),n=1,ncar)                                         
220   format(/1x,'ERROR IN SUBR. ALTONU -- Invalid integer variable ',         &
     & 'entry'/5x,'Input text (ALP) = ',15a1)                                   
      write(*,*)                                                                
      stop 'Halted in ALTONU -- see list file.'                                 
!                                                                               
! --------------------------------------------------------------------          
! --- LOGICAL variables                                                         
! --------------------------------------------------------------------          
301   continue                                                                  
      if(ncar-istar.ne.1)go to 308                                              
      if(num2(istar+1).eq.16)then                                               
!                                                                               
! ---    Variable = T                                                           
         lv=.true.                                                              
         return                                                                 
      else if(num2(istar+1).eq.17)then                                          
!                                                                               
! ---    Variable = F                                                           
         lv=.false.                                                             
         return                                                                 
      endif                                                                     
308   continue                                                                  
      write(ioout,320)(alp(n),n=1,ncar)                                         
320   format(/1x,'ERROR IN SUBR. ALTONU -- Invalid logical variable ',         &
     & 'entry'/5x,'Input text (ALP) = ',15a1)                                   
      write(*,*)                                                                
      stop 'Halted in ALTONU -- see list file.'                                 
!                                                                               
! --- Determine repetition factor                                               
400   continue                                                                  
      irep=0                                                                    
!                                                                               
! --- ISTAR is the position of array containing "*"                             
      istrm1=istar-1                                                            
      do 430 i=1,istrm1                                                         
      if(num2(i).ge.10)go to 408                                                
      iexp=istrm1-i                                                             
      xmult=1.0                                                                 
      if(iexp.ne.0)xmult=ten**iexp                                              
      irep=irep+xmult*num2(i)+0.5                                               
430   continue                                                                  
      go to(101,201,301),jtype                                                  
      write(ioout,55)jtype,(alp(n),n=1,ncar)                                    
      write(*,*)                                                                
      stop 'Halted in ALTONU -- see list file.'                                 
408   continue                                                                  
      write(ioout,420)(alp(n),n=1,ncar)                                         
420   format(/1x,'ERROR IN SUBR. ALTONU -- Invalid repetition factor ',        &
     & 'entry'/5x,'Input text (ALP) = ',15a1)                                   
      write(*,*)                                                                
      stop 'Halted in ALTONU -- see list file.'                                 
      end                                                                       
!----------------------------------------------------------------------         
      subroutine deblnk(cstor1,ilim1,il2,cstor2,nlim)                           
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 030402                 DEBLNK         
! ---            J. Scire, Earth Tech, Inc.                                     
!                                                                               
! --- PURPOSE:  Remove all blank or "+" characters from the character           
!               string within delimiters                                        
!               Only characters in the range ilim1 to il2 may be                
!               written to output array                                         
!                                                                               
! --- UPDATE                                                                    
! --- V2.1 (030402) from V2.0 (980918) (DGS)                                    
!               - Split DEBLNK action (removes ' ', '+') into                   
!                 DEBLNK and DEPLUS(new)                                        
!                                                                               
! --- INPUTS:                                                                   
!                                                                               
!    CSTOR1(mxcol) - character*1 array - Input character string                 
!            ILIM1 - integer           - Array element at which search          
!                                        for blanks begins                      
!              IL2 - integer           - Array element at which search          
!                                        for blanks ends                        
!        Parameters: MXCOL                                                      
!                                                                               
! --- OUTPUT:                                                                   
!                                                                               
!    CSTOR2(mxcol) - character*1 array - Output character string                
!                                        (without blanks within text)           
!             NLIM - integer           - Length of output string                
!                                        (characters)                           
!                                                                               
! --- DEBLNK called by:  (utility)                                              
! --- DEBLNK calls:      none                                                   
!----------------------------------------------------------------------         
!                                                                               
! --- Include parameter statements                                              
      include 'params.cal'                                                      
!                                                                               
      character*1 cstor1(mxcol),cstor2(mxcol),cblnk                             
      data cblnk/' '/                                                           
!                                                                               
      ind=0                                                                     
      do 10 i=ilim1,il2                                                         
      if(cstor1(i).eq.cblnk)go to 10                                            
!                                                                               
! --- transfer non-blank character into output array                            
      ind=ind+1                                                                 
      cstor2(ind)=cstor1(i)                                                     
10    continue                                                                  
      nlim=ind                                                                  
      if(ind.eq.mxcol)return                                                    
!                                                                               
! --- pad rest of output array                                                  
      indp1=ind+1                                                               
      do 20 i=indp1,mxcol                                                       
      cstor2(i)=cblnk                                                           
20    continue                                                                  
      return                                                                    
      end                                                                       
!----------------------------------------------------------------------         
      subroutine deplus(cstor1,ilim1,il2,cstor2,nlim)                           
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 030402                 DEPLUS         
! ---            J. Scire, Earth Tech, Inc.                                     
!                                                                               
! --- PURPOSE:  Remove all "+" characters from the character                    
!               string within delimiters                                        
!               Only characters in the range ilim1 to il2 may be                
!               written to output array                                         
!                                                                               
! --- INPUTS:                                                                   
!                                                                               
!    CSTOR1(mxcol) - character*1 array - Input character string                 
!            ILIM1 - integer           - Array element at which search          
!                                        for plus begins                        
!              IL2 - integer           - Array element at which search          
!                                        for plus ends                          
!        Parameters: MXCOL                                                      
!                                                                               
! --- OUTPUT:                                                                   
!                                                                               
!    CSTOR2(mxcol) - character*1 array - Output character string                
!                                        (without plus within text)             
!             NLIM - integer           - Length of output string                
!                                        (characters)                           
!                                                                               
! --- DEPLUS called by:  (utility)                                              
! --- DEPLUS calls:      none                                                   
!----------------------------------------------------------------------         
!                                                                               
! --- Include parameter statements                                              
      include 'params.cal'                                                      
!                                                                               
      character*1 cstor1(mxcol),cstor2(mxcol),cblnk,cplus                       
      data cblnk/' '/,cplus/'+'/                                                
!                                                                               
      ind=0                                                                     
      do 10 i=ilim1,il2                                                         
      if(cstor1(i).eq.cplus)go to 10                                            
!                                                                               
! --- transfer non-plus character into output array                             
      ind=ind+1                                                                 
      cstor2(ind)=cstor1(i)                                                     
10    continue                                                                  
      nlim=ind                                                                  
      if(ind.eq.mxcol)return                                                    
!                                                                               
! --- pad rest of output array                                                  
      indp1=ind+1                                                               
      do 20 i=indp1,mxcol                                                       
      cstor2(i)=cblnk                                                           
20    continue                                                                  
      return                                                                    
      end                                                                       
!----------------------------------------------------------------------         
      subroutine tright(cstor1,ilim1,il2,cstor2,nlim)                           
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 040330                 TRIGHT         
! ---            D. Strimaitis, Earth Tech, Inc.                                
!                                                                               
! --- PURPOSE:  Remove all blank characters in the range ilim1 to il2           
!               that lie to the RIGHT of the last non-blank character           
!               in the string before il2. Also remove the character             
!               at il2 if it is blank.                                          
!               Only characters in the range ilim1 to il2 may be                
!               written to the output array.                                    
!                                                                               
!               Example --                                                      
!               Range    : ilim1=3, il2=21                                      
!               CSTOR1   :  2   for this run   !                                
!               Position : 000000000111111111122                                
!                          123456789012345678901                                
!               CSTOR2   :    for this run!                                     
!                                                                               
! --- INPUTS:                                                                   
!                                                                               
!    CSTOR1(mxcol) - character*1 array - Input character string                 
!            ILIM1 - integer           - Array element at which search          
!                                        for blanks begins                      
!              IL2 - integer           - Array element at which search          
!                                        for blanks ends                        
!        Parameters: MXCOL                                                      
!                                                                               
! --- OUTPUT:                                                                   
!                                                                               
!    CSTOR2(mxcol) - character*1 array - Output character string                
!                                        (with right-blanks removed)            
!             NLIM - integer           - Length of output string                
!                                        (characters)                           
!                                                                               
! --- TRIGHT called by:  (utility)                                              
! --- TRIGHT calls:      none                                                   
!----------------------------------------------------------------------         
!                                                                               
! --- Include parameter statements                                              
      include 'params.cal'                                                      
!                                                                               
      character*1 cstor1(mxcol),cstor2(mxcol),cblnk                             
      data cblnk/' '/                                                           
                                                                                
! --- Position of last non-blank character                                      
      klast=0                                                                   
      il2m1=il2-1                                                               
      do k=ilim1,il2m1                                                          
         if(cstor1(k).NE.cblnk) klast=k                                         
      enddo                                                                     
                                                                                
! --- Transfer all characters in range up to klast                              
      ind=0                                                                     
      if(klast.GT.0) then                                                       
         do k=ilim1,klast                                                       
            ind=ind+1                                                           
            cstor2(ind)=cstor1(k)                                               
         enddo                                                                  
      endif                                                                     
! --- Add last character in range if non-blank                                  
      if(cstor1(il2).NE.cblnk) then                                             
         ind=ind+1                                                              
         cstor2(ind)=cstor1(il2)                                                
      endif                                                                     
      nlim=ind                                                                  
      if(ind.EQ.mxcol) return                                                   
                                                                                
! --- Pad rest of output array                                                  
      indp1=ind+1                                                               
      do i=indp1,mxcol                                                          
         cstor2(i)=cblnk                                                        
      enddo                                                                     
                                                                                
      return                                                                    
      end                                                                       
!----------------------------------------------------------------------         
      subroutine tleft(cstor1,ilim1,il2,cstor2,nlim)                            
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 040330                  TLEFT         
! ---            D. Strimaitis, Earth Tech, Inc.                                
!                                                                               
! --- PURPOSE:  Remove all blank characters in the range ilim1 to il2           
!               that lie to the LEFT of the first non-blank character           
!               in the string after ilim1. Also remove the character            
!               at ilim1 if it is blank.                                        
!               Only characters in the range ilim1 to il2 may be                
!               written to the output array.                                    
!                                                                               
!               Example --                                                      
!               Range    : ilim1=2, il2=19                                      
!               CSTOR1   :  2   for this run   !                                
!               Position : 123456789111111111122                                
!                                   012345678901                                
!               CSTOR2   : 2for this run                                        
!                                                                               
! --- INPUTS:                                                                   
!                                                                               
!    CSTOR1(mxcol) - character*1 array - Input character string                 
!            ILIM1 - integer           - Array element at which search          
!                                        for blanks begins                      
!              IL2 - integer           - Array element at which search          
!                                        for blanks ends                        
!        Parameters: MXCOL                                                      
!                                                                               
! --- OUTPUT:                                                                   
!                                                                               
!    CSTOR2(mxcol) - character*1 array - Output character string                
!                                        (with left-blanks removed)             
!             NLIM - integer           - Length of output string                
!                                        (characters)                           
!                                                                               
! --- TLEFT called by:  (utility)                                               
! --- TLEFT calls:      none                                                    
!----------------------------------------------------------------------         
!                                                                               
! --- Include parameter statements                                              
      include 'params.cal'                                                      
!                                                                               
      character*1 cstor1(mxcol),cstor2(mxcol),cblnk                             
      data cblnk/' '/                                                           
                                                                                
! --- Position of first non-blank character                                     
      kfrst=0                                                                   
      ilim1p1=ilim1+1                                                           
      do k=il2,ilim1p1,-1                                                       
         if(cstor1(k).NE.cblnk) kfrst=k                                         
      enddo                                                                     
                                                                                
      ind=0                                                                     
! --- Pass first character in range if non-blank                                
      if(cstor1(ilim1).NE.cblnk) then                                           
         ind=ind+1                                                              
         cstor2(ind)=cstor1(ilim1)                                              
      endif                                                                     
                                                                                
! --- Transfer all characters in range from kfrst                               
      if(kfrst.GT.0) then                                                       
         do k=kfrst,il2                                                         
            ind=ind+1                                                           
            cstor2(ind)=cstor1(k)                                               
         enddo                                                                  
      endif                                                                     
      nlim=ind                                                                  
      if(ind.EQ.mxcol) return                                                   
                                                                                
! --- Pad rest of output array                                                  
      indp1=ind+1                                                               
      do i=indp1,mxcol                                                          
         cstor2(i)=cblnk                                                        
      enddo                                                                     
                                                                                
      return                                                                    
      end                                                                       
!----------------------------------------------------------------------         
      subroutine setvar(itype,irep,xx,jj,ll,cv,xarr,jarr,larr,carr)             
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 110225                 SETVAR         
! ---            J. Scire                                                       
!                                                                               
! --- PURPOSE:  Fill the output variable or array with the value read           
!               from the input file                                             
!                                                                               
! --- UPDATE                                                                    
! --- V2.58 (110225) from V1.0 (950122) (DGS)                                   
!               - Add IVTYPE=5 (char*4 array with commas retained               
!                 as delimiters for parsing)                                    
!                                                                               
! --- INPUTS:                                                                   
!                                                                               
!            ITYPE - integer        - Variable type (1=real, 2=integer,         
!                                     3=logical, 4=character*4,                 
!                                     5=character*4 includes commas)            
!             IREP - integer        - Repetition factor                         
!                                     If ITYPE = 4, IREP refers to the          
!                                     number of characters in the field)        
!               XX - real           - Real value read from input                
!                                     file (Used only if ITYPE=1)               
!               JJ - integer        - Integer value read from input             
!                                     file (Used only if ITYPE=2)               
!               LL - logical*4      - Logical value read from input             
!                                     file (Used only if ITYPE=3)               
!        CV(mxcol) - character*4    - Character*4 values read from input        
!                                     file (Used only if ITYPE=4)               
!                                                                               
!         PARAMETER:  MXCOL                                                     
!                                                                               
! --- OUTPUT:                                                                   
!                                                                               
!          XARR(*) - real array     - Output real array (or scalar if           
!                                     IREP=1) -- Used only if ITYPE=1           
!          JARR(*) - integer array  - Output integer array (or scalar if        
!                                     IREP=1) -- Used only if ITYPE=2           
!          LARR(*) - logical array  - Output logical array (or scalar if        
!                                     IREP=1) -- Used only if ITYPE=3           
!          CARR(*) - character*4    - Output character*4 array (or              
!                                     scalar if IREP=1) -- Used only if         
!                                     ITYPE=4                                   
!                                                                               
! --- SETVAR called by:  READIN                                                 
! --- SETVAR calls:      none                                                   
!----------------------------------------------------------------------         
!                                                                               
! --- Include parameter statements                                              
      include 'params.cal'                                                      
!                                                                               
      real xarr(*)                                                              
      integer jarr(*)                                                           
      logical*4 larr(*),ll                                                      
      character*4 carr(*),cv(mxcol)                                             
!                                                                               
      go to(10,20,30,40,50),itype                                               
!                                                                               
! --- real variable                                                             
10    continue                                                                  
      do 15 i=1,irep                                                            
      xarr(i)=xx                                                                
15    continue                                                                  
      return                                                                    
!                                                                               
! --- integer variable                                                          
20    continue                                                                  
      do 25 i=1,irep                                                            
      jarr(i)=jj                                                                
25    continue                                                                  
      return                                                                    
!                                                                               
! --- logical variable                                                          
30    continue                                                                  
      do 35 i=1,irep                                                            
      larr(i)=ll                                                                
35    continue                                                                  
      return                                                                    
!                                                                               
! --- character*4 variable string                                               
40    continue                                                                  
      do 45 i=1,irep                                                            
      carr(i)=cv(i)                                                             
45    continue                                                                  
      return                                                                    
!                                                                               
! --- character*4 variable string                                               
50    continue                                                                  
      do 55 i=1,irep                                                            
      carr(i)=cv(i)                                                             
55    continue                                                                  
      return                                                                    
                                                                                
      end                                                                       
!----------------------------------------------------------------------         
      subroutine allcap(cstor2,nlim)                                            
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 950122                 ALLCAP         
! ---            J. Scire, SRC                                                  
!                                                                               
! --- PURPOSE:  Convert all lower case letters within a character               
!               string to upper case                                            
!                                                                               
! --- INPUTS:                                                                   
!                                                                               
!    CSTOR2(mxcol) - character*1 array - Input character string                 
!             NLIM - integer           - Length of string (characters)          
!        Parameters: MXCOL                                                      
!                                                                               
! --- OUTPUT:                                                                   
!                                                                               
!    CSTOR2(mxcol) - character*1 array - Output character string with           
!                                        lower case letters converted           
!                                        to upper case                          
!                                                                               
! --- ALLCAP called by:  READIN                                                 
! --- ALLCAP calls:      none                                                   
!----------------------------------------------------------------------         
!                                                                               
! --- Include parameter statements                                              
      include 'params.cal'                                                      
!                                                                               
      character*1 cstor2(mxcol),cchar,clc(29),cuc(29)                           
!                                                                               
      data clc/'i','n','x','a','e','o','u','b','c','d','f','g','h',            &
     & 'j','k','l','m','p','q','r','s','t','v','w','y','z','-','.',            &
     & '*'/                                                                     
      data cuc/'I','N','X','A','E','O','U','B','C','D','F','G','H',            &
     & 'J','K','L','M','P','Q','R','S','T','V','W','Y','Z','-','.',            &
     & '*'/                                                                     
!                                                                               
      do 100 i=1,nlim                                                           
      cchar=cstor2(i)                                                           
!                                                                               
      do 50 j=1,29                                                              
      if(cchar.eq.clc(j))then                                                   
         cstor2(i)=cuc(j)                                                       
         go to 52                                                               
      endif                                                                     
50    continue                                                                  
52    continue                                                                  
100   continue                                                                  
!                                                                               
      return                                                                    
      end                                                                       
!----------------------------------------------------------------------         
      subroutine datetm(rdate,rtime,rcpu)                                       
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 140318                 DATETM         
! ---            J. Scire                                                       
!                                                                               
! --- PURPOSE:  Get system date and time from system clock, and                 
!               elapsed CPU time                                                
! --- UPDATES                                                                   
! --- V2.57-V2.6.0 140318(MBN):Remove obsolete Lahey F77L code,                 
!                              and etime calls.                                 
! --- V1.0-V2.57  090202 (DGS): Activate CPU time (F95 call)                    
!                                                                               
! --- INPUTS:  none                                                             
!                                                                               
! --- OUTPUT:  rdate  - C*10 - Current system date (MM-DD-YYYY)                 
!              rtime  - C*8  - Current system time (HH:MM:SS)                   
!               rcpu  - real - CPU time (sec) from system utility               
!                                                                               
! --- DATETM called by:  SETUP, FIN                                             
! --- DATETM calls:      DATE_AND_TIME (F95)                                    
!                        CPU_TIME (F95)                                         
!                        YR4C                                                   
!----------------------------------------------------------------------         
      character*8  rtime                                                        
      character*10 rdate                                                        
                                                                                
! --- Local store                                                               
      character*11 stime                                                        
      character*8 sdate                                                         
                                                                                
! --- Set initial base CPU time to -1.                                          
      data rcpu0/-1./                                                           
      SAVE rcpu0                                                                
                                                                                
! --- System date in CCYYMMDD                                                   
! --- System clock in HHMMSS.sss, where sss = thousandths of seconds            
      call DATE_AND_TIME(sdate,stime)                                           
! --- Pass to output formats (MM-DD-YYYY) and (HH:MM:SS)                        
      rdate='  -  -    '                                                        
      rdate(1:2)=sdate(5:6)                                                     
      rdate(4:5)=sdate(7:8)                                                     
      rdate(7:10)=sdate(1:4)                                                    
      rtime='  :  :  '                                                          
      rtime(1:2)=stime(1:2)                                                     
      rtime(4:5)=stime(3:4)                                                     
      rtime(7:8)=stime(5:6)                                                     
! --- Get CPU time from F95 intrinsic procedure                                 
      call CPU_TIME(rcpu1)                                                      
                                                                                
! --- Construct 4-digit year from current 2-digit year (if found)               
      read(rdate(7:10),'(i4)') iyr                                              
      call YR4C(iyr)                                                            
      write(rdate(7:10),'(i4)') iyr                                             
                                                                                
! --- Update base CPU time on first call                                        
      if(rcpu0.LT.0.0) rcpu0=rcpu1                                              
                                                                                
! --- Return CPU time difference from base                                      
      rcpu=rcpu1-rcpu0                                                          
                                                                                
!c --- DEBUG                                                                    
!      write(*,*)'DATETM: stime,rcpu0,rcpu1,rcpu = ',                           
!     &                   stime,rcpu0,rcpu1,rcpu                                
                                                                                
      return                                                                    
      end                                                                       
!----------------------------------------------------------------------         
      subroutine fmt_date(io,fmt1,fmt2,sdate)                                   
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 090511               FMT_DATE         
!                D. Strimaitis                                                  
!                                                                               
! --- PURPOSE:  Change the format of a date string                              
!                                                                               
! --- INPUTS:                                                                   
!            io - integer       - Listfile output unit number                   
!          fmt1 - character*12  - Input date format                             
!                 MM-DD-YYYY                                                    
!                 DD-MM-YYYY                                                    
!                 YYYY-MM-DD                                                    
!                 YYYY-DD-MM                                                    
!                 DD-MMM-YYYY                                                   
!                 MMM-DD-YYYY                                                   
!         sdate - character*12  - Date string to convert                        
!          fmt2 - character*12  - Output date format                            
!                 MM-DD-YYYY                                                    
!                 DD-MM-YYYY                                                    
!                 YYYY-MM-DD                                                    
!                 YYYY-DD-MM                                                    
!                 DD-MMM-YYYY                                                   
!                 MMM-DD-YYYY                                                   
!                                                                               
! --- OUTPUT:                                                                   
!         sdate - character*12  - Converted date string                         
!                                                                               
! --- FMT_DATE called by:  (any)                                                
! --- FMT_DATE calls:      ALLCAP                                               
!----------------------------------------------------------------------         
      character*12 fmt1,fmt2,sdate                                              
      character*3 month3(12),month3uc(12),amon3                                 
      character*1 amon(3)                                                       
      integer io                                                                
                                                                                
! --- Set abbreviation names for months                                         
      data month3/'Jan','Feb','Mar','Apr','May','Jun',                         &
     &            'Jul','Aug','Sep','Oct','Nov','Dec'/                          
      data month3uc/'JAN','FEB','MAR','APR','MAY','JUN',                       &
     &            'JUL','AUG','SEP','OCT','NOV','DEC'/                          
                                                                                
! --- Extract input month, day and year                                         
      if(fmt1(1:10).EQ.'MM-DD-YYYY') then                                       
         read(sdate(1:2),'(i2)') imon                                           
         read(sdate(4:5),'(i2)') iday                                           
         read(sdate(7:10),'(i4)') iyear                                         
      elseif(fmt1(1:10).EQ.'DD-MM-YYYY') then                                   
         read(sdate(1:2),'(i2)') iday                                           
         read(sdate(4:5),'(i2)') imon                                           
         read(sdate(7:10),'(i4)') iyear                                         
      elseif(fmt1(1:10).EQ.'YYYY-MM-DD') then                                   
         read(sdate(1:4),'(i4)') iyear                                          
         read(sdate(6:7),'(i2)') imon                                           
         read(sdate(9:10),'(i4)') iday                                          
      elseif(fmt1(1:10).EQ.'YYYY-DD-MM') then                                   
         read(sdate(1:4),'(i4)') iyear                                          
         read(sdate(6:7),'(i2)') iday                                           
         read(sdate(9:10),'(i4)') imon                                          
      elseif(fmt1(1:11).EQ.'DD-MMM-YYYY') then                                  
         read(sdate(1:2),'(i2)') iday                                           
         read(sdate(4:6),'(3a1)') amon                                          
         read(sdate(8:11),'(i4)') iyear                                         
         call ALLCAP(amon,3)                                                    
         amon3=amon(1)//amon(2)//amon(3)                                        
         imon=0                                                                 
         do k=1,12                                                              
            if(amon3.EQ.month3uc(k)) imon=k                                     
         enddo                                                                  
      elseif(fmt1(1:11).EQ.'MMM-DD-YYYY') then                                  
         read(sdate(1:3),'(3a1)') amon                                          
         read(sdate(5:6),'(i2)') iday                                           
         read(sdate(8:11),'(i4)') iyear                                         
         call ALLCAP(amon,3)                                                    
         amon3=amon(1)//amon(2)//amon(3)                                        
         imon=0                                                                 
         do k=1,12                                                              
            if(amon3.EQ.month3uc(k)) imon=k                                     
         enddo                                                                  
      else                                                                      
         write(io,*)'FMT_DATE:  Invalid input format = ',fmt1                   
         write(io,*)'Expected: MM-DD-YYYY, DD-MM-YYYY,  YYYY-MM-DD'             
         write(io,*)'          YYYY-DD-MM, DD-MMM-YYYY, MMM-DD-YYYY'            
         stop 'Halted in FMT_DATE --- see list file'                            
      endif                                                                     
                                                                                
! --- Check for valid month index                                               
      if(imon.LT.1 .OR. imon.GT.12) then                                        
         write(io,*)'FMT_DATE:  Invalid month in date = ',sdate                 
         write(io,*)'           for input format      = ',fmt1                  
         stop 'Halted in FMT_DATE --- see list file'                            
      endif                                                                     
                                                                                
! --- Create output date string                                                 
      if(fmt2(1:10).EQ.'MM-DD-YYYY') then                                       
         sdate='MM-DD-YYYY  '                                                   
         write(sdate(1:2),'(i2.2)') imon                                        
         write(sdate(4:5),'(i2.2)') iday                                        
         write(sdate(7:10),'(i4.4)') iyear                                      
      elseif(fmt2(1:10).EQ.'DD-MM-YYYY') then                                   
         sdate='DD-MM-YYYY  '                                                   
         write(sdate(1:2),'(i2.2)') iday                                        
         write(sdate(4:5),'(i2.2)') imon                                        
         write(sdate(7:10),'(i4.4)') iyear                                      
      elseif(fmt2(1:10).EQ.'YYYY-MM-DD') then                                   
         sdate='YYYY-MM-DD  '                                                   
         write(sdate(1:4),'(i4.4)') iyear                                       
         write(sdate(6:7),'(i2.2)') imon                                        
         write(sdate(9:10),'(i2.2)') iday                                       
      elseif(fmt2(1:10).EQ.'YYYY-DD-MM') then                                   
         sdate='YYYY-DD-MM  '                                                   
         write(sdate(1:4),'(i4.4)') iyear                                       
         write(sdate(6:7),'(i2.2)') iday                                        
         write(sdate(9:10),'(i2.2)') imon                                       
      elseif(fmt2(1:11).EQ.'DD-MMM-YYYY') then                                  
         sdate='DD-MMM-YYYY '                                                   
         write(sdate(1:2),'(i2.2)') iday                                        
         sdate(4:6)=month3(imon)                                                
         write(sdate(8:11),'(i4.4)') iyear                                      
      elseif(fmt2(1:11).EQ.'MMM-DD-YYYY') then                                  
         sdate='MMM-DD-YYYY '                                                   
         sdate(1:3)=month3(imon)                                                
         write(sdate(5:6),'(i2.2)') iday                                        
         write(sdate(8:11),'(i4.4)') iyear                                      
      else                                                                      
         write(io,*)'FMT_DATE:  Invalid output format = ',fmt2                  
         write(io,*)'Expected: MM-DD-YYYY, DD-MM-YYYY,  YYYY-MM-DD'             
         write(io,*)'          YYYY-DD-MM, DD-MMM-YYYY, MMM-DD-YYYY'            
         stop 'Halted in FMT_DATE --- see list file'                            
      endif                                                                     
                                                                                
      return                                                                    
      end                                                                       
!----------------------------------------------------------------------         
      subroutine etime(rcpu)                                                    
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 941215                  ETIME         
! ---            J. Scire, SRC                                                  
!                                                                               
! --- PURPOSE:  Dummy system CPU time routine for PC                            
!               DO NOT USE THIS ROUTINE ON SUNs                                 
!                                                                               
! --- INPUTS:  none                                                             
!                                                                               
! --- OUTPUT:  RCPU  - real - CPU time (sec) -- set to zero for PC              
!                                                                               
! --- ETIME called by:  DATETM                                                  
! --- ETIME calls:      none                                                    
!----------------------------------------------------------------------         
      rcpu=0.0                                                                  
!                                                                               
      return                                                                    
      end                                                                       
!----------------------------------------------------------------------         
      subroutine undrflw(lflag)                                                 
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 030402                UNDRFLW         
!                D. Strimaitis,  Earth Tech Inc.                                
!                                                                               
! --- PURPOSE:  This routine takes advantage of the Lahey F77L routine          
!               UNDER0 to set underflows to zero.  When other compilers         
!               are used, there may be a similar routine.  If none              
!               exists, place a dummy statement here and use compiler           
!               switches to configure the NDP response to an underflow.         
!                                                                               
!               This routine contains calls for several different               
!               compilers, but only one should be active at any one             
!               time.                                                           
!                                                                               
!----------------------------------------------------------------------         
      logical lflag                                                             
                                                                                
!c --- Lahey F77L Compiler (begin)                                              
!c -------------------------------                                              
!c --- Lahey F77 compiler -- set underflows ( < 10**-38 ) to zero               
!      call UNDER0(lflag)                                                       
!c --- Lahey F77L Compiler (end)                                                
                                                                                
! --- Dummy (no action on underflows)                                           
! -----------------------------------                                           
      lflag=.TRUE.                                                              
! --- Dummy (end)                                                               
                                                                                
      return                                                                    
      end                                                                       
!----------------------------------------------------------------------         
      subroutine comline(ctext)                                                 
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 040330                COMLINE         
!                J. Scire, SRC                                                  
!                                                                               
! --- PURPOSE:  Call the compiler-specific system routine that will             
!               pass back the command line argument after the text              
!               that executed the program                                       
!                                                                               
!               This routine contains calls for several different               
!               compilers, but only one should be active at any one             
!               time.                                                           
!                                                                               
! --- UPDATE                                                                    
! --- V2.3 (040330) to V2.6.0 (040330)  MBN                                     
!               - Removed obsolete Compaq, Microsoft, and HP compiler codes     
!               - Removed getcl (Lahey-only function not needed)                
! --- V2.2 (960521) to V2.3 (040330)  DGS                                       
!               - Replace strings c*70 with c*132                               
!                                                                               
! --- INPUTS:                                                                   
!                                                                               
!          CTEXT - character*132 - Default command line argument #1             
!                                                                               
! --- OUTPUT:                                                                   
!                                                                               
!          CTEXT - character*132 - Command line argument #1                     
!                                  If command line argument is                  
!                                  missing, CTEXT is not changed                
!                                                                               
! --- COMLINE called by:  SETUP                                                 
! --- COMLINE calls:      IARGC, GETARG - compiler routines                     
!                                                                               
!----------------------------------------------------------------------         
!                                                                               
      character*132 ctext,cdeflt                                                
!                                                                               
! --- The following is for any system without a command line routine            
! --- and is also used as a default                                             
      cdeflt=ctext                                                              
!                                                                               
! ----------------                                                              
! --- Intel ifort, Lahey lf95, and GNU gfortran compilers:                      
! ----------------                                                              
      numargs=IARGC()                                                           
      if(numargs.ge.1)then                                                      
         call GETARG(1,ctext)                                                   
      endif                                                                     
!                                                                               
! --- If no command line arguments, use default                                 
      if(ctext(1:1).eq.' ')ctext=cdeflt                                         
                                                                                
      return                                                                    
      end                                                                       
                                                                                
!----------------------------------------------------------------------         
      subroutine open_err(iolst,cfrom,cftype,cfname,iunit)                      
!----------------------------------------------------------------------         
!                                                                               
! --- CALUTILS   Version: 7.0.0    Level: 141010              OPEN_ERR          
!               D. Strimaitis, Exponent Inc.                                    
!                                                                               
! --- PURPOSE:  Report error in opening a file                                  
!                                                                               
! --- INPUTS:                                                                   
!               IOLST - integer   - Unit number of output list file             
!                                   (<0 if not available)                       
!               CFROM - char*     - Called-From string to report error          
!              CFTYPE - char*     - File-type string                            
!              CFNAME - char*     - File-name string                            
!               IUNIT - integer   - File unit number                            
!                                                                               
! --- OUTPUT:                                                                   
!                                                                               
! --- OPEN_ERR called by:  ()                                                   
! --- OPEN_ERR calls:                                                           
!----------------------------------------------------------------------         
      implicit none                                                             
                                                                                
! --- Declare arguments                                                         
      character(len=*) :: cfrom,cftype,cfname                                   
      integer          :: iolst, iunit                                          
                                                                                
      if(iolst.GT.0) then                                                       
         write(iolst,*)                                                         
         write(iolst,*)'ERROR opening '//TRIM(cftype)                           
         write(iolst,*)'   File Name: '//TRIM(cfname)                           
         write(iolst,*)'   File Unit: ',iunit                                   
         write(iolst,*)'Problem reported from '//TRIM(cfrom)                    
         write(iolst,*)                                                         
         write(iolst,*)'The file may not exist in this location'                
         write(iolst,*)'Check the spelling of the name and the location'        
         write(*,*)                                                             
         stop 'ERROR: File not found -- see list file'                          
      else                                                                      
         write(*,*)                                                             
         write(*,*)'ERROR opening '//TRIM(cftype)                               
         write(*,*)'   File Name: '//TRIM(cfname)                               
         write(*,*)'   File Unit: ',iunit                                       
         write(*,*)'Problem reported from '//TRIM(cfrom)                        
         write(*,*)                                                             
         write(*,*)'The file may not exist in this location'                    
         write(*,*)'Check the spelling of the name, and the location'           
         stop                                                                   
      endif                                                                     
                                                                                
      end                                                                       
                                                                                
