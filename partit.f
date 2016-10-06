
      subroutine partit(cpart,recres,lyr,frlign)


      use parm
      
      implicit none
      
      integer j, idf
      

!! ...... Argument declarations
      integer lyr
      real    cpart, recres(MAXIEL), !! cdonor(ISOS), edonor(MAXIEL), 
     &        frlign, friso

!! ...... Partition residue from compartments cdonor and edonor
!! ...... into layer lyr of structural and metabolic.
!! ...... cpart is the amount of carbon in the residue.
!! ...... recres contains the n/c, p/c, and s/c ratios in the residue.
!! ...... frlign is the fraction of the incoming material which is lignin.
!! ...... friso is the fraction of cpart which is labeled.

!! ...... Local variables
      integer   iel, 
      real      accum, caddm, cadds, dirabs(MAXIEL),
     &          eaddm, eadds, epart(MAXIEL),
     &          fligst, frmet, frn, rcetot, rlnres, delin,
     &          dellig, c13c12r, c13frac, c13lig, c13nlig,
     &          c13struc, c12struc, c13met, c12met
      real      namt
      real      no3nh4
      double precision frac_nh4, frac_no3
      character subname*10
      real      damr(2,3)
      real      pabres
      real      damrmn (3)
      !!real      frac_nh4suf, frac_no3suf
      !!real      frac_nh4sol, frac_no3sol
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
      
      rcestr = 1. !! give fake values to this array, SWAT should have have this variable. Need to check with Xuesong 
      
      no3nh4 = sol_no3(lyr,j) + sol_nh3(lyr,j)  !! daycent only consider impacts of nitrogen
      
      
      
      frac_nh4 = sol_nh3(lyr,j) / no3nh4
      frac_no3 = 1 - frac_nh4

      accum = 0.0
      
      damr(1,1) = 0.00000           
      damr(1,2) = 0.00000           
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
      do 10 iel = 1, 1 !! only consider N here

!! ........ Compute amount of element in residue.
        epart(iel) = cpart * recres(iel)

!! ........ Direct absorption of mineral element by residue
!! ........ (mineral will be transferred to donor compartment
!! ........ and then partitioned into structural and metabolic
!! ........ using flow routines.)

!! ........ If minerl(SRFC,iel) is negative then dirabs = zero.

        if (no3nh4.lt. 0.) then
          dirabs(iel) = 0.0
        else
          dirabs(iel)=damr(lyr,iel)*0.1*no3nh4*
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
          
   
        sol_no3 (lyr,j) = sol_no3 (lyr,j) - 10.*namt * frac_no3 
        sol_nh3 (lyr,j) = sol_nh3 (lyr,j) - 10.*namt * frac_nh4 
    
        endif
      
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

!! ...... Make sure at least 20% goes to metabolic
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

c ... Adjust lignin
!!XXXXX      call adjlig(strcis,fligst,cadds,strlig(lyr))   !! need to check with Xuesong if this process is needed, the key is if SWAT has a lignin fraction for the structural pool


!! ........ Carbon added to metabolic

           metcis = metcis + caddm   !! need to be changed according to SWAT setting


           strcis = strcis + cadds   !! need to be changed according to SWAT setting

   
       

!! ...... Adjust lignin
    !1  call adjlig(strucc(lyr),fligst,cadds,strlig(lyr))

!! ...... Partition mineral elements into structural and metabolic
      do 20 iel = 1, nelem

!! ........ Flow into structural
        eadds = cadds/rcestr(iel)
       !! call flow(edonor(iel),struce(lyr,iel),time,eadds)
       !!  edonor(iel) = edonor(iel) - eadds
         struce(lyr,iel) = struce(lyr,iel) + eadds
!! ........ Flow into metabolic
        eaddm = epart(iel)+dirabs(iel)-eadds
       !! call flow(edonor(iel),metabe(lyr,iel),time,eaddm)
        !! edonor(iel) = edonor(iel) - eaddm
        metabe(lyr,iel) = metabe(lyr,iel) + eaddm
        
        
 20    continue

 999   continue

      return
      end
