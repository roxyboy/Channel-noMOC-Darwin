C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/darwin/DARWIN_PARAMS.h,v 1.1 2011/04/13 18:56:23 jahn Exp $
C $Name:  $

c DARWIN_PARAMS.h 
c description: runtime configurable parameters for darwin package
c

c Iron chemistry parameters
c    scav           - chemical scavenging rate for iron
c    ligand_tot     - total organic ligand (assumed constant)
c    ligand_stab    - ligand stability
c    freefemax      - maximum free iron
      COMMON/darwin_ironparam/ alpfe, scav, ligand_tot, ligand_stab,
     &       freefemax
      _RL alpfe, scav, ligand_tot, ligand_stab
      _RL freefemax

#ifdef DAR_DIAG_CHL
      COMMON/DAR_DIAG_CHL_PARAMS/
     &       Geider_Bigalphachl, Geider_smallalphachl,
     &       Geider_Bigchl2cmax, Geider_smallchl2cmax,
     &       Geider_Bigchl2cmin, Geider_smallchl2cmin,
     &       Doney_Bmin, Doney_Bmax, Doney_PARstar,
     &       Cloern_A, Cloern_B, Cloern_C, Cloern_chl2cmin
      _RL Geider_Bigalphachl
      _RL Geider_Smallalphachl
      _RL Geider_Bigchl2cmax
      _RL Geider_Smallchl2cmax
      _RL Geider_Bigchl2cmin
      _RL Geider_smallchl2cmin
      _RL Doney_Bmin
      _RL Doney_Bmax
      _RL Doney_PARstar
      _RL Cloern_chl2cmin
      _RL Cloern_A
      _RL Cloern_B
      _RL Cloern_C
#endif /* DAR_DIAG_CHL */

