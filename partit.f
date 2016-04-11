
      subroutine partit(cpart,recres,lyr,cdonor,edonor,frlign)


      use parm
      
      implicit none
      
      integer j, idf
      

!! ...... Argument declarations
      integer lyr
      real    cpart, recres(MAXIEL), cdonor(ISOS), edonor(MAXIEL), 
     &        frlign, friso

!! ...... Partition residue from compartments cdonor and edonor
!! ...... into layer lyr of structural and metabolic.
!! ...... cpart is the amount of carbon in the residue.
!! ...... recres contains the n/c, p/c, and s/c ratios in the residue.
!! ...... frlign is the fraction of the incoming material which is lignin.
!! ...... friso is the fraction of cpart which is labeled.

!! ...... Local variables
      integer   iel, clyr
      real      accum, caddm, cadds, dirabs(MAXIEL),
     &          eaddm, eadds, epart(MAXIEL),
     &          fligst, frmet, frn, rcetot, rlnres, delin,
     &          dellig, c13c12r, c13frac, c13lig, c13nlig,
     &          c13struc, c12struc, c13met, c12met
      real      namt
      real      no3nh4suf, no3nh4sol
      double precision frac_nh4, frac_no3
      character subname*10
      real      damr(2,3)
      real      pabres
      real      damrmn (3)
      real      frac_nh4suf, frac_no3suf
      real      frac_nh4sol, frac_no3sol
      real      spl(2)
      
   !!  fake variables defined here to explain the patition process to xuesong. Real varialbes need to be find from xuesong's soc process
   
         
      real metcis  !! metbolic carbon pool
      real strcis  !! structural carbon pool
      real rcestr(3)  !! carbon to elementy ratio in structural pool
      real struce(3,3)  !! structural element pool
      real metabe(3,3)  !! metabolic element pool
      
      j = 0
      idf = 0
      

      
      j = ihru
      idf = idplt(j)
      
      rcestr = 1. !! give fake values to this array
      
      no3nh4suf = sol_no3(1,j) + sol_nh3(1,j)  !! need to account for P later
      no3nh4sol = sol_no3(2,j) + sol_nh3(2,j)
      
      
      frac_nh4suf = sol_nh3(1,j) / no3nh4suf
      frac_no3suf = 1 - frac_nh4suf
      frac_nh4sol = sol_nh3(2,j) / no3nh4sol 
      frac_no3sol = 1 - frac_nh4sol 
     
      accum = 0.0
      
      damr(1,1) = 0.00000           
      damr(1,3) = 0.00000           
      damr(1,3) = 0.01000          
      damr(2,1) = 0.02000           
      damr(2,2) = 0.02000           
      damr(2,3) = 0.04000    
      
      
      damrmn (1) = 15.00000          
      damrmn (2) = 150.00000         
      damrmn (3) = 150.00000
      
      spl(1) = 0.85000          !! 'SPL(1)'    
      spl(2) = 0.01300           !!'SPL(2)'     
          
      pabres = 100.        
     

      if (cpart .lt. 1.e-07) then
        goto 999
      endif


!! ...... For each mineral element...
      do 10 iel = 1, nelem

!! ........ Compute amount of element in residue.
        epart(iel) = cpart * recres(iel)

!! ........ Direct absorption of mineral element by residue
!! ........ (mineral will be transferred to donor compartment
!! ........ and then partitioned into structural and metabolic
!! ........ using flow routines.)

!! ........ If minerl(SRFC,iel) is negative then dirabs = zero.

        if (no3nh4suf .lt. 0.) then
          dirabs(iel) = 0.0
        else
          dirabs(iel)=damr(lyr,iel)*no3nh4suf*
     &                amax1(cpart/pabres,1.)
        endif

!! ........ If C/E ratio is too low, transfer just enough to make
!! ........ C/E of residue = damrmn
        if (epart(iel)+dirabs(iel) .le. 0.0) then
          rcetot = 0.0
        else
          rcetot = cpart/(epart(iel)+dirabs(iel))
        endif

        if (rcetot .lt. damrmn(iel)) then
          dirabs(iel) = cpart/damrmn(iel) - epart(iel)
        endif
        if (dirabs(iel) .lt. 0.) then
          dirabs(iel) = 0.
        endif
        
        if (iel .eq. N) then
          namt = -1.0*dirabs(iel)
          clyr = 1
         !! call cmpnfrac(clyr,ammonium,nitrate,minerl,frac_nh4,frac_no3)
         !! call update_npool(clyr, namt, frac_nh4, frac_no3, ammonium,
    !! &                      nitrate, subname)
    
    
       !! for this part, daycent have very complicated processes. we assume nutrient exchange between soil and detritus only occure between litter and the first soil layer    
    
    
        sol_no3 (1,j) = sol_no3 (1,j) - namt * frac_no3suf 
        sol_nh3 (1,j) = sol_nh3 (1,j) - namt * frac_nh4suf 
    
        endif
        !! call flow(minerl(1,iel),edonor(iel),time,dirabs(iel))
10    continue

!! ...... Partition carbon into structural and metabolic fraction of
!! ...... residue (including direct absorption) which is nitrogen
      frn = (epart(1)+dirabs(1)) / (cpart*2.5)

!! ...... Lignin/nitrogen ratio of residue
      rlnres = frlign/frn

!! ...... Carbon added to metabolic
!! ...... Compute the fraction of carbon that goes to metabolic.
      frmet = spl(INTCPT)-spl(SLOPE)*rlnres

!! ...... Make sure the fraction of residue which is lignin isn't
!! ...... greater than the fraction which goes to structural.  -rm 12/91
      if (frlign .gt. (1.0 - frmet)) then
        frmet = (1.0 - frlign)
      endif

!! ...... Make sure at least 1% goes to metabolic
      if (frmet .lt. 0.20) then
        frmet = .20
      endif

!! ...... Compute amounts to flow
      caddm = cpart * frmet !! carbon added to metabolic pool
      if (caddm .lt. 0) then
        caddm = 0.0
      endif
      cadds = cpart-caddm   !! carbon added to structural pool
!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX STOP here because xuesong may have done the particition
!! ...... Adjust lignin content of structural.
!! ...... fligst is the fraction of incoming structural residue
!! ...... which is lignin; restricting it to a maximum of .8
      fligst = frlign/(cadds/cpart)

!! ...... Changed allowable maximum from .8 to .6 -rm 5/92
!! ...... Changed maximum fraction from .6 to 1.0  -lh 1/93
      if (fligst .gt. 1.0) then
         fligst = 1.0
      endif

!! ...... Determine what type of labeling is to be done


!! ........ Carbon added to metabolic
           cdonor = cdonor -caddm     !! carbon transfer from donor to metobalic pool
           metcis = metcis + caddm

           cdonor = cdonor -cadds     !! carbon transfer from donor to metobalic pool
           strcis = strcis + cadds

   
       

!! ...... Adjust lignin
    !1  call adjlig(strucc(lyr),fligst,cadds,strlig(lyr))

!! ...... Partition mineral elements into structural and metabolic
      do 20 iel = 1, nelem

!! ........ Flow into structural
        eadds = cadds/rcestr(iel)
       !! call flow(edonor(iel),struce(lyr,iel),time,eadds)
         edonor(iel) = edonor(iel) - eadds
         struce(lyr,iel) = struce(lyr,iel) + eadds
!! ........ Flow into metabolic
        eaddm = epart(iel)+dirabs(iel)-eadds
       !! call flow(edonor(iel),metabe(lyr,iel),time,eaddm)
        edonor(iel) = edonor(iel) - eaddm
        metabe(lyr,iel) = metabe(lyr,iel) + eaddm
        
        
 20    continue

 999   continue

      return
      end
