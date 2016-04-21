
     
      subroutine wooddec (ttswatlyr)

      use parm
      
      implicit none
      
      integer j, idf
      integer ttswatlyr
      

!! ...... Argument declarations
!!      real             dtm
!!      double precision newminrl

!! ...... Wood decomposition       written by vek 04/91

!! ...... defac  = decomposition factor based on water and temperature
!! ......          (computed in prelim and in cycle)
!! ...... pligst = fixed parameter that represents the effect of
!! ......          of lignin-to-structural-ratio on structural
!! ......          decomposition

!! ...... Function declarations


!! ...... Local variables
      real      tcflow, pheff, ph
      real      tfunc   !! tempature function
      real      agdefac !! temperature effect on surface litter decomposition
      real      normalizer
      real      stem !! soil temperature
      real      agwfunc   !! water effect on surface decomposition
      real      rel_wc    !! relative ater content
      real      swclimit
      real      bgdefac
      real      bgwfunc
      real      krainwfunc
      real      avg_rel_wc
      integer   lyr
      real      drain
      real      rprpet
      real      aminrl(3)
      real      ps1co2(2)
      real      ratnew1(3,2)
      real      ratnew2(3,2)
      real      varat11(3,3)
      real      varat12(3,3)
      real      varat21(3,3)
      real      varat22(3,3)
      real      rsplig
      real      anerb
     
  !! function declaration
      real      agdrat, bgdrat, anerob    
      integer   iel
      real :: sol_thick(sol_nly(ihru))
  !! initialize variables    
      j = 0
      idf = 0
    
 
      j = ihru
      idf = idplt(j)
      
      aminrl = 0. 
      
      ps1co2(1) = 0.45000           
      ps1co2(2) = 0.55000   
      
      
      varat11(1,1) = 16.00000          !!'VARAT11(1,1)'
      varat11(2,1) = 8.00000           !!'VARAT11(2,1)'
      varat11(3,1) = 2.00000           !!'VARAT11(3,1)'
      varat11(1,2) = 150.0000          !!'VARAT11(1,2)'
      varat11(2,2) = 30.00000          !!'VARAT11(2,2)'
      varat11(3,2) = 2.00000           !!'VARAT11(3,2)'
      varat11(1,3) = 200.0000          !! 'VARAT11(1,3)'
      varat11(2,3) = 50.00000          !!'VARAT11(2,3)'
      varat11(3,3) = 2.00000           !!'VARAT11(3,3)'
      varat12(1,1) = 16.00000          !!'VARAT12(1,1)'
      varat12(2,1) = 8.00000           !!'VARAT12(2,1)'
      varat12(3,1) = 2.00000           !!'VARAT12(3,1)'
      varat12(1,2) = 150.0000          !!'VARAT12(1,2)'
      varat12(2,2) = 30.00000          !!'VARAT12(2,2)'
      varat12(3,2) = 2.00000           !!'VARAT12(3,2)'
      varat12(1,3) = 200.0000          !!'VARAT12(1,3)'
      varat12(2,3) = 50.00000          !!'VARAT12(2,3)'
      varat12(3,3) = 2.00000           !!'VARAT12(3,3)'
      varat21(1,1) = 20.0000           !!'VARAT21(1,1)'
      varat21(2,1) = 10.0000           !!'VARAT21(2,1)'
      varat21(3,1) = 2.00000           !!'VARAT21(3,1)'
      varat21(1,2) = 400.000           !!'VARAT21(1,2)'
      varat21(2,2) = 100.000           !!'VARAT21(2,2)'
      varat21(3,2) = 2.00000           !!'VARAT21(3,2)'
      varat21(1,3) = 400.000           !!'VARAT21(1,3)'
      varat21(2,3) = 100.000           !!'VARAT21(2,3)'
      varat21(3,3) = 2.00000           !!'VARAT21(3,3)'
      varat22(1,1) = 20.00000          !!'VARAT22(1,1)'
      varat22(2,1) = 10.00000          !!'VARAT22(2,1)'
      varat22(3,1) = 2.00000           !!'VARAT22(3,1)'
      varat22(1,2) = 400.00000         !!'VARAT22(1,2)'
      varat22(2,2) = 100.0000          !!'VARAT22(2,2)'
      varat22(3,2) = 2.00000           !!'VARAT22(3,2)'
      varat22(1,3) = 400.00000         !!'VARAT22(1,3)'
      varat22(2,3) = 100.0000          !!'VARAT22(2,3)'
      varat22(3,3) = 2.00000           !!'VARAT22(3,3)'        
      
      rsplig = 0.3
      

      
      do lyr = 1, sol_nly(j)
        
           if (lyr == 1) then
	              sol_thick(lyr) = sol_z(lyr,j)
	        else	
	              sol_thick(lyr) = sol_z(lyr,j) - sol_z(lyr-1,j)
	        end if
       end do 

      
!! ...... FINE BRANCHES
!! ......   wood1c       = C in dead fine branch component of forest system (g/m2)
!! ......   decw1        = intrinsic rate of decomposition of dead fine branches
!! ......   wdlig(FBRCH) = lignin fraction for fine branches
      stem = sol_tmp(1,j) !! use surface temperature
      normalizer = 11.75+(29.7/DPI)*atan(DPI*0.031*(30.-15.4))

      tfunc = max(0.01, (11.75+(29.7/3.14)*atan(3.14*0.031*(stem-15.4)))
     &  /normalizer)
           

        !! Compute water effect for surface decomposition using the */
        !! top soil layer (first layer)
       swclimit = sol_wpmm(1,j)*0.1-0.01
       rel_wc = ((sol_st(1,j) + sol_wpmm(1,j))*0.1)/(sol_thick(1)*0.1) -
     &                    swclimit/
     &                  ((sol_fc(1, j)+sol_wpmm(1,j))*0.1 + 0.01)
                       
      if (rel_wc > 1.0) then
            agwfunc = 1.0
         else 
         if (rel_wc< 0.0)  rel_wc = 0.0
          
            agwfunc = 1.0/(1.0 + 30.0 * exp(-9.0 * rel_wc))
      endif
       avg_rel_wc = 0 
!!       
      do lyr = 2, 3  !! check, what if there are less than 3 layers
 !!    average relative water content
        
       rel_wc = ((sol_st(lyr,j) + sol_wpmm(lyr,j))*0.1)/(sol_thick(lyr)
     &                    *0.1) -(sol_wpmm(lyr,j)*0.1-0.01)/
     &                  ((sol_fc(lyr, j)+sol_wpmm(lyr,j))*0.1 + 0.01)
 
       avg_rel_wc = avg_rel_wc + rel_wc*sol_thick(lyr)*0.1 
       
                       
       end do

       
      !! calculate et related variables

       rprpet = (precipday + sol_st(3,j))/ pet_day
    
      
        
        avg_rel_wc = avg_rel_wc / 
     &                     ((sol_thick(2) + sol_thick(3))*0.1)
      !!   calcualte below ground water scalor                 
                                  
          if (avg_rel_wc > 1.0) then
            bgwfunc = 1.0
           else 
            bgwfunc = 1.0/(1.0 + 30.0 * exp(-9.0 * avg_rel_wc))
          endif
     
 !! calculate aboveground scaler
      if (precipday> 5.) then
        agdefac = max(0.000001, tfunc * 2.0)
        else 
        agdefac = max(0.000001, tfunc * agwfunc)
      endif
      
 
       
      !! calculate belowground scaler
      !! bgdefac = max(0.000001, tfunc **bgwfunc*krainwfunc) !! this is the original function, but do not know what is krainwfunc
       bgdefac = max(0.000001, tfunc **bgwfunc)
      
!! ... Calculate the effect impact of anerobic conditions on decomposition
!! ... Last parameter = 1 if microcosm
      drain = 1.
      anerb = anerob(drain,rprpet,pet_day,0)
      
      ph = sol_ph(1,j)
      
      
      
      if (woodc(j,1) .gt. 1.e-07) then

!! ......Compute pH effect on decomposition
        pheff = 0.5 + (1.1/DPI)*atan(DPI*0.7*(ph-4))
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)
        
        
      !! calculate mineralize nutrient 
       
       
       !!availabel nutrient for surface litter
      
       
       aminrl (N) = (sol_no3(1,j) + sol_nh3(1,j))/10. !! P and S are not accounted for yet here we assume that mineral N are stored in the first three layers.
       
       aminrl (P) = sol_solp(1,j) / 10.
       
       
       
!! ... Determine C/E ratios for flows from structural material to
!! ... surface som1 and surface som2
!! ... The second index of ratnew(iel,*) should be 1 or 2
!! ... (not SRFC or SOIL) for som1c or som2c. -MDH 1/23/2012
!! ... Create ratnew1 and ratnew2 for surface and soil. -MDH 9/24/2012
      do 30 iel=1,nelem
!! ..... ratnew1: SRFC som1 and som2
        ratnew1(iel,1) = agdrat(aminrl,varat11,iel)
        ratnew1(iel,2) = agdrat(aminrl,varat21,iel)
!! ..... ratnew2: SOIL som1 and som2
        ratnew2(iel,1) = bgdrat(aminrl,varat12,iel)
        ratnew2(iel,2) = bgdrat(aminrl,varat22,iel)
30    continue
        

!! ......Compute total C flow out of fine branches
!! ......Add pH effect on decomposition to calculation, cak - 08/02/02
!!        tcflow = wood1c * defac * decw1 * exp(-pligst(SRFC) * 
!!     &           wdlig(FBRCH)) * dtm
        tcflow = woodc(j,1) * agdefac * DECWf(idf,1)/365.0 * 
     &      exp(-pligst(SRFC) * WDLIGf(idf,FBRCH)) * pheff

!! ......Decompose fine branches into som1 and som2 with CO2 loss.
!! ......Use surface decomposition ratios, MDH Nov 2012



       call declig(1,aminrl,WDLIGf(idf,FBRCH),SRFC,nelem,1,ps1co2,
     &             ratnew1, rsplig,tcflow)
      endif

!! ... LARGE WOOD
!! ...   wood2c       = C in dead large wood component of forest system (g/m2)
!1 ...   decw2        = intrinsic rate of decomposition of dead large wood
!! ...   wdlig(LWOOD) = lignin fraction for large wood

      if (woodc(j,2) .gt. 1.e-07) then

!! ......Compute pH effect on decomposition
        !!pheff = catanf(ph, 4.0, 0.5, 1.1, 0.7)
        !!pheff = min(pheff, 1.0)
        !!pheff = max(pheff, 0.0)

!! ......Compute total C flow out of large wood
!! ......Add pH effect on decomposition to calculation, cak - 08/02/02
!!        tcflow = wood2c * defac * decw2 * exp(-pligst(SRFC) * 
!!     &           wdlig(LWOOD)) * dtm
        tcflow = woodc(j,2) * agdefac * DECWf(idf,2) *
     &    exp(-pligst(SRFC) * WDLIGf(idf,LWOOD)) *  pheff

!! ......Decompose large wood into som1 and som2 with CO2 loss.
!! ......Use surface decomposition ratios, MDH Nov 2012
        call declig(2,aminrl,WDLIGf(idf,LWOOD),SRFC,nelem,1,ps1co2,
     &              ratnew1,rsplig,tcflow)
            
      endif


 
!! ... COARSE ROOTS
!! ...   wood3c       = C in dead coarse root component of forest system (g/m2)
!! ...   decw3        = intrinsic rate of decomposition of dead coarse roots
!! ...   wdlig(CROOT) = lignin fraction for coarse roots

      if (woodc(j,3) .gt. 1.e-07) then

!! ......Compute pH effect on decomposition
        !!pheff = catanf(ph, 4.0, 0.5, 1.1, 0.7)
        !!pheff = min(pheff, 1.0)
        !!pheff = max(pheff, 0.0)

!! ......Compute total C flow out of coarse roots.
!! ......Add pH effect on decomposition to calculation, cak - 08/02/02
!!        tcflow = wood3c * defac * decw3 * exp(-pligst(SOIL) *
!!     &           wdlig(CROOT)) *  anerb * dtm
        tcflow = woodc(j,3) * bgdefac * DECWf(idf,3)/365. * 
     &   exp(-pligst(SOIL) * WDLIGf(idf,CROOT)) *  anerb *  pheff

!! ......Decompose coarse roots into som1 and som2 with CO2 loss.
!! ......Use soil decomposition ratios, MDH Nov 2012

      !!availabel nutrient for litter in soils
      aminrl = 0.
      
      do lyr = 2, sol_nly(j)
      
       aminrl (N) = aminrl (N) + (sol_no3(lyr,j) + sol_nh3(lyr,j))/10. 
       
       aminrl (P) = aminrl (P) + sol_solp(lyr,j) / 10.
     
      end do  


        call declig(3,aminrl,WDLIGf(idf,CROOT),SOIL,nelem,1,ps1co2,
     &              ratnew2,rsplig,tcflow)


      endif

      return
      end
