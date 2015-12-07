/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// Several special cases are defined:
//   NBTYPE: exclusion method (NBPAIR, NBSELF -- mutually exclusive)
//   FULLELECT full electrostatics calculation?

#ifdef ARCH_POWERPC
#include <builtins.h>
#endif

#if defined(__SSE2__) && ! defined(NAMD_DISABLE_SSE)
#include <emmintrin.h>  // We're using SSE2 intrinsics
#if defined(__INTEL_COMPILER)
#define __align(X) __declspec(align(X) )
#elif defined(__GNUC__) || defined(__PGI)
#define __align(X)  __attribute__((aligned(X) ))
#else
#define __align(X) __declspec(align(X) )
#endif
#endif

#ifdef DEFINITION // (
  #include "LJTable.h"
  #include "Molecule.h"
  #include "ComputeNonbondedUtil.h"
#endif // )
  #include "Parameters.h"
#if NAMD_ComputeNonbonded_SortAtoms != 0
  #include "PatchMap.h"
#endif

// determining class name
#undef NAME
#undef CLASS
#undef CLASSNAME
#define NAME CLASSNAME(calc)

#undef PAIR
#if NBTYPE == NBPAIR
  #define PAIR(X) X
  #define CLASS ComputeNonbondedPair
  #define CLASSNAME(X) ENERGYNAME( X ## _pair )
#else
  #define PAIR(X)
#endif

#undef SELF
#if NBTYPE == NBSELF
  #define SELF(X) X
  #define CLASS ComputeNonbondedSelf
  #define CLASSNAME(X) ENERGYNAME( X ## _self )
#else
  #define SELF(X)
#endif

#undef ENERGYNAME
#undef ENERGY
#undef NOENERGY
#ifdef CALCENERGY
  #define ENERGY(X) X
  #define NOENERGY(X)
  #define ENERGYNAME(X) SLOWONLYNAME( X ## _energy )
#else
  #define ENERGY(X)
  #define NOENERGY(X) X
  #define ENERGYNAME(X) SLOWONLYNAME( X )
#endif

#undef SLOWONLYNAME
#undef FAST
#ifdef SLOWONLY
  #define FAST(X)
  #define SLOWONLYNAME(X) MERGEELECTNAME( X ## _slow )
#else
  #define FAST(X) X
  #define SLOWONLYNAME(X) MERGEELECTNAME( X )
#endif

#undef MERGEELECTNAME
#undef SHORT
#undef NOSHORT
#ifdef MERGEELECT
  #define SHORT(X)
  #define NOSHORT(X) X
  #define MERGEELECTNAME(X) FULLELECTNAME( X ## _merge )
#else
  #define SHORT(X) X
  #define NOSHORT(X)
  #define MERGEELECTNAME(X) FULLELECTNAME( X )
#endif

#undef FULLELECTNAME
#undef FULL
#undef NOFULL
#ifdef FULLELECT
  #define FULLELECTNAME(X) TABENERGYNAME( X ## _fullelect )
  #define FULL(X) X
  #define NOFULL(X)
#else
  #define FULLELECTNAME(X) TABENERGYNAME( X )
  #define FULL(X)
  #define NOFULL(X) X
#endif

#undef TABENERGYNAME
#undef TABENERGY
#undef NOTABENERGY
#ifdef TABENERGYFLAG
  #define TABENERGYNAME(X) FEPNAME( X ## _tabener )
  #define TABENERGY(X) X
  #define NOTABENERGY(X)
#else
  #define TABENERGYNAME(X) FEPNAME( X )
  #define TABENERGY(X)
  #define NOTABENERGY(X) X
#endif

// Here are what these macros stand for:
// FEP/NOT_FEP: FEP free energy perturbation is active/inactive 
//      (does NOT use LAM)
// LES: locally-enhanced sampling is active
// LAM: scale the potential by a factor lambda (except FEP)
// INT: measure interaction energies
// PPROF: pressure profiling

#undef FEPNAME
#undef FEP
#undef LES
#undef SPT
#undef INT
#undef PPROF
#undef LAM
#undef CUDA
#undef ALCH
#undef TI
#undef GO
#define FEPNAME(X) LAST( X )
#define FEP(X)
#define ALCHPAIR(X)
#define NOT_ALCHPAIR(X) X
#define LES(X)
#define SPT(X)
#define INT(X)
#define PPROF(X)
#define LAM(X)
#define CUDA(X)
#define ALCH(X)
#define TI(X)
#define GO(X)
#ifdef FEPFLAG
  #undef FEPNAME
  #undef FEP
  #undef ALCH
  #define FEPNAME(X) LAST( X ## _fep )
  #define FEP(X) X
  #define ALCH(X) X
#endif
#ifdef TIFLAG
  #undef FEPNAME
  #undef TI
  #undef ALCH
  #define FEPNAME(X) LAST( X ## _ti )
  #define TI(X) X
  #define ALCH(X) X
#endif
#ifdef LESFLAG
  #undef FEPNAME
  #undef LES
  #undef LAM
  #define FEPNAME(X) LAST( X ## _les )
  #define LES(X) X
  #define LAM(X) X
#endif
#ifdef SPTFLAG
  #undef FEPNAME
  #undef SPT
  #undef LAM
  #define FEPNAME(X) LAST( X ## _spt )
  #define SPT(X) X
  #define LAM(X) X
#endif
#ifdef INTFLAG
  #undef FEPNAME
  #undef INT
  #define FEPNAME(X) LAST( X ## _int )
  #define INT(X) X
#endif
#ifdef PPROFFLAG
  #undef FEPNAME
  #undef INT
  #undef PPROF
  #define FEPNAME(X) LAST( X ## _pprof )
  #define INT(X) X
  #define PPROF(X) X
#endif
#ifdef GOFORCES
  #undef FEPNAME
  #undef GO
  #define FEPNAME(X) LAST( X ## _go )
  #define GO(X) X
#endif
#ifdef NAMD_CUDA
  #undef CUDA
  #define CUDA(X) X
#endif

#define LAST(X) X

// see if things are really messed up
SELF( PAIR( foo bar ) )
LES( FEP( foo bar ) )
LES( INT( foo bar ) )
SPT( FEP( foo bar ) )
SPT( INT( foo bar ) )
FEP( INT( foo bar ) )
LAM( INT( foo bar ) )
FEP( NOENERGY( foo bar ) )
ENERGY( NOENERGY( foo bar ) )
TABENERGY(NOTABENERGY( foo bar ) )

#if NAMD_ComputeNonbonded_SortAtoms != 0 && ( 0 PAIR( + 1 ) )
  #define COMPONENT_DOTPRODUCT(A,B)  ((A##_x * B##_x) + (A##_y * B##_y) + (A##_z * B##_z))
#endif


// ************************************************************
// function header
void ComputeNonbondedUtil :: NAME
  ( nonbonded *params )

// function body
{
  // int NAME;  // to track errors via unused variable warnings

  if ( ComputeNonbondedUtil::commOnly ) return;

#ifdef NAMD_CUDA
  NOENERGY(return;)
#endif

  // speedup variables
  BigReal *reduction = params->reduction;
  SimParameters *simParams = Node::Object()->simParameters;
  Node *node = Node::Object();

  PPROF(
  BigReal *pressureProfileReduction = params->pressureProfileReduction;
  const BigReal invThickness = 1.0 / pressureProfileThickness;
  )

  Pairlists &pairlists = *(params->pairlists);
#ifdef NAMD_CUDA
  int savePairlists = 0;
  int usePairlists = 0;
#else
  int savePairlists = params->savePairlists;
  int usePairlists = params->usePairlists;
#endif
  pairlists.reset();
  // PAIR(iout << "--------\n" << endi;)
  
  // BEGIN LA
  const int doLoweAndersen = params->doLoweAndersen;
  // END LA

  // local variables
  int exclChecksum = 0;
  FAST
  (
   // JLai
  ENERGY( BigReal vdwEnergy = 0;
	  GO( BigReal groLJEnergy = 0;
	      BigReal groGaussEnergy = 0;
	      BigReal goEnergyNative = 0;
	      BigReal goEnergyNonnative = 0; ) )
  SHORT
  (
  ENERGY( BigReal electEnergy = 0; )
  )

  FEP
  (
  ENERGY( BigReal vdwEnergy_s = 0; BigReal vdwEnergy_s_Left = 0; )
  SHORT
  (
  ENERGY( BigReal electEnergy_s = 0; )
  )
  )
  )
    
#ifndef A2_QPX
  FAST 
  (
  SHORT
  (
  BigReal virial_xx = 0;
  BigReal virial_xy = 0;
  BigReal virial_xz = 0;
  BigReal virial_yy = 0;
  BigReal virial_yz = 0;
  BigReal virial_zz = 0;
  )
  )
  FULL
  (
  ENERGY( BigReal fullElectEnergy = 0; )
  FEP
  (
  ENERGY( BigReal fullElectEnergy_s = 0; )
  )
  BigReal fullElectVirial_xx = 0;
  BigReal fullElectVirial_xy = 0;
  BigReal fullElectVirial_xz = 0;
  BigReal fullElectVirial_yy = 0;
  BigReal fullElectVirial_yz = 0;
  BigReal fullElectVirial_zz = 0;
  )
#else 
  vector4double virial_v0   = (vector4double)(0.0);
  vector4double virial_v1   = (vector4double)(0.0);
  vector4double virial_v2   = (vector4double)(0.0);
  vector4double fullvirial_v0   = (vector4double) (0.0);
  vector4double fullvirial_v1   = (vector4double) (0.0);
  vector4double fullvirial_v2   = (vector4double) (0.0);
  BigReal fullElectEnergy = 0;
  BigReal fullElectEnergy_s = 0;
#endif

  // Bringing stuff into local namespace for speed.
  const BigReal offset_x = params->offset.x;
  const BigReal offset_y = params->offset.y;
  const BigReal offset_z = params->offset.z;

  register const BigReal plcutoff2 = \
 			params->plcutoff * params->plcutoff;
  register const BigReal groupplcutoff2 = \
	 		params->groupplcutoff * params->groupplcutoff;
  const BigReal dielectric_1 = ComputeNonbondedUtil:: dielectric_1;
  const LJTable* const ljTable = ComputeNonbondedUtil:: ljTable;
  LJTable::TableEntry ljNull;  ljNull.A = 0; ljNull.B = 0;
  const LJTable::TableEntry* const lj_null_pars = &ljNull;
  const Molecule* const mol = ComputeNonbondedUtil:: mol;
  SHORT
  (
  const BigReal* const table_four = ComputeNonbondedUtil:: table_short;
  )
  FULL
  (
  SHORT
  (
  const BigReal* const slow_table = ComputeNonbondedUtil:: slow_table;
  )
  NOSHORT
  (
//#if 1 ALCH(-1)
  const BigReal* const table_four = ComputeNonbondedUtil:: table_noshort;
//#else  // have to switch this for ALCH
//  BigReal* table_four = ComputeNonbondedUtil:: table_noshort;
//#endif
  )
  )
  BigReal scaling = ComputeNonbondedUtil:: scaling;
//#ifdef  A2_QPX
//  vector4double scalingv = vec_splats(scaling);
//#endif
  const BigReal modf_mod = 1.0 - scale14;
  FAST
  (
  const BigReal switchOn2 = ComputeNonbondedUtil:: switchOn2;
  const BigReal c1 = ComputeNonbondedUtil:: c1;
  const BigReal c3 = ComputeNonbondedUtil:: c3;
  )
  const BigReal r2_delta = ComputeNonbondedUtil:: r2_delta;
  const int r2_delta_exp = ComputeNonbondedUtil:: r2_delta_exp;
  // const int r2_delta_expc = 64 * (r2_delta_exp - 127);
  const int r2_delta_expc = 64 * (r2_delta_exp - 1023);

  ALCH(
    const BigReal switchdist2 = ComputeNonbondedUtil::switchOn2;
    const BigReal cutoff2 = ComputeNonbondedUtil::cutoff2;
    const BigReal switchfactor = 1./((cutoff2 - switchdist2)*(cutoff2 - switchdist2)*(cutoff2 - switchdist2));
    const BigReal alchElecLambdaStart = ComputeNonbondedUtil::alchElecLambdaStart;
    const BigReal alchVdwLambdaEnd = ComputeNonbondedUtil::alchVdwLambdaEnd;
    const BigReal alchVdwShiftCoeff = ComputeNonbondedUtil::alchVdwShiftCoeff;
    const Bool Fep_WCA_repuOn = ComputeNonbondedUtil::Fep_WCA_repuOn;
    const Bool Fep_WCA_dispOn = ComputeNonbondedUtil::Fep_WCA_dispOn;
    const Bool Fep_ElecOn = ComputeNonbondedUtil::Fep_ElecOn;
    const Bool Fep_Wham = ComputeNonbondedUtil::Fep_Wham;
    const BigReal WCA_rcut1 = ComputeNonbondedUtil::WCA_rcut1;
    const BigReal WCA_rcut2 = ComputeNonbondedUtil::WCA_rcut2;
    const BigReal WCA_rcut3 = ComputeNonbondedUtil::WCA_rcut3;

    /*lambda values 'up' are for atoms scaled up with lambda (partition 1)*/
    BigReal lambdaUp = ComputeNonbondedUtil::alchLambda;
    BigReal elecLambdaUp =  (lambdaUp <= alchElecLambdaStart)? 0. : \
              (lambdaUp - alchElecLambdaStart) / (1. - alchElecLambdaStart);
    BigReal vdwLambdaUp = 
        (lambdaUp >= alchVdwLambdaEnd)? 1. : lambdaUp / alchVdwLambdaEnd; 
    BigReal vdwShiftUp = ComputeNonbondedUtil::alchVdwShiftCoeff * (1-vdwLambdaUp);
    FEP(BigReal lambda2Up = ComputeNonbondedUtil::alchLambda2;)
    FEP(BigReal elecLambda2Up = (lambda2Up <= alchElecLambdaStart)? 0. : \
              (lambda2Up - alchElecLambdaStart) / (1. - alchElecLambdaStart);)
    FEP(BigReal vdwLambda2Up = 
        (lambda2Up >= alchVdwLambdaEnd)? 1. : lambda2Up / alchVdwLambdaEnd;) 
    FEP(BigReal vdwShift2Up = ComputeNonbondedUtil::alchVdwShiftCoeff * (1-vdwLambda2Up);)


    FEP( if( (Fep_Wham) && (Fep_WCA_repuOn) ) {
    	elecLambdaUp=0.0; 
    	vdwLambdaUp=ComputeNonbondedUtil::alchRepLambda;
    })
    FEP( if( (Fep_Wham) && (Fep_WCA_dispOn) ) {
    	elecLambdaUp=0.0; 
    	vdwLambdaUp=ComputeNonbondedUtil::alchDispLambda;
    })
    FEP( if( (Fep_Wham) && (Fep_ElecOn) ) {
    	elecLambdaUp=ComputeNonbondedUtil::alchElecLambda; 
    	vdwLambdaUp=1.0;
    	vdwLambda2Up=1.0;
	    vdwShiftUp = 0.0;
	    vdwShift2Up = 0.0;
    })
        
    /*lambda values 'down' are for atoms scaled down with lambda (partition 2)*/
    BigReal lambdaDown = 1 - ComputeNonbondedUtil::alchLambda;
    BigReal elecLambdaDown =  (lambdaDown <= alchElecLambdaStart)? 0. : \
              (lambdaDown - alchElecLambdaStart) / (1. - alchElecLambdaStart);
    BigReal vdwLambdaDown = 
        (lambdaDown >= alchVdwLambdaEnd)? 1. : lambdaDown / alchVdwLambdaEnd; 
    BigReal vdwShiftDown = ComputeNonbondedUtil::alchVdwShiftCoeff * (1-vdwLambdaDown);
    FEP(BigReal lambda2Down = 1 - ComputeNonbondedUtil::alchLambda2;)
    FEP(BigReal elecLambda2Down = (lambda2Down <= alchElecLambdaStart)? 0. : \
        (lambda2Down - alchElecLambdaStart) / (1. - alchElecLambdaStart); )
    FEP(BigReal vdwLambda2Down = 
        (lambda2Down >= alchVdwLambdaEnd)? 1. : lambda2Down / alchVdwLambdaEnd; )
    FEP(BigReal vdwShift2Down = ComputeNonbondedUtil::alchVdwShiftCoeff * (1-vdwLambda2Down);)


    FEP( if( (Fep_Wham) && (Fep_WCA_repuOn) ) {
    	elecLambdaDown=0.0; 
    	vdwLambdaDown = 1.0 - ComputeNonbondedUtil::alchRepLambda;
    })
    FEP( if( (Fep_Wham) && (Fep_WCA_dispOn) ) {
    	elecLambdaDown=0.0; 
    	vdwLambdaDown = 1.0 - ComputeNonbondedUtil::alchDispLambda;
    })
    FEP( if( (Fep_Wham) && (Fep_ElecOn) ) {
    	elecLambdaDown = 1.0 - ComputeNonbondedUtil::alchElecLambda; 
    	vdwLambdaDown = 1.0;
    	vdwLambda2Down = 1.0;
	    vdwShiftDown = 0.0;
	    vdwShift2Down = 0.0;
    })
  
  // Thermodynamic Integration Notes: 
  // Separation of atom pairs into different pairlist according to lambda
  // coupling, for alchemical free energy calculations. Normal pairlists
  // (pairlist[nxm]_save) are re-used for non-lambda-coupled pairs; new ones
  // (pairlist[nxm][12]_save are created for atoms switched up or down with
  // lambda respectively.
  // This makes TI coding far easier and more readable, since it's necessary 
  // to store dU/dlambda in different variables depending on whether atoms are
  // being switched up or down. Further, it allows the prior explicit coding of 
  // the separation-shifted vdW potential to stay in place without a big 
  // performance hit, since it's no longer necessary to calculate vdW potentials
  // explicity for the bulk (non-alchemical) interactions - so that part of the 
  // free energy code stays readable and easy to modify. Finally there should
  // be some performance gain over the old FEP implementation as we no
  // longer have to check the partitions of each atom pair and switch
  // parameters accordingly for every single NonbondedBase2.h loop - this is 
  // done at the pairlist level
  
  int pswitchTable[3*3];
  // determines which pairlist alchemical pairings get sent to
  // 0: non-alchemical pairs, partition 0 <-> partition 0
  // 1: atoms scaled up as lambda increases, p0<->p1
  // 2: atoms scaled down as lambda increases, p0<->p2
  // all p1<->p2 interactions to be dropped (99)
  // in general, 'up' refers to 1, 'down' refers to 2
  for (int ip=0; ip<3; ++ip) {
    for (int jp=0; jp<3; ++jp ) {
      pswitchTable[ip+3*jp] = 0;
      if ((ip==1 && jp==0) || (ip==0 && jp==1)) pswitchTable[ip+3*jp] = 1;
      if ((ip==2 && jp==0) || (ip==0 && jp==2)) pswitchTable[ip+3*jp] = 2;
      if (ip + jp == 3) pswitchTable[ip+3*jp] = 99; // no interaction

      if (! ComputeNonbondedUtil::alchDecouple) {
        // no decoupling: interactions within a partition are treated like
        // normal alchemical pairs
        if (ip == 1 && jp == 1) pswitchTable[ip+3*jp] = 1;
        if (ip == 2 && jp == 2) pswitchTable[ip+3*jp] = 2;
      }
      if (ComputeNonbondedUtil::alchDecouple) {
        // decoupling: PME calculates extra grids so that while PME 
        // interaction with the full system is switched off, a new PME grid
        // containing only alchemical atoms is switched on. Full interactions 
        // between alchemical atoms are maintained; potentials within one 
        // partition need not be scaled here.
        if (ip == 1 && jp == 1) pswitchTable[ip+3*jp] = 0;
        if (ip == 2 && jp == 2) pswitchTable[ip+3*jp] = 0;
      }
    }
  }

  // dU/dlambda variables for thermodynamic integration
  TI(
      BigReal vdwEnergy_ti_1 = 0;
      BigReal vdwEnergy_ti_2 = 0;
      SHORT(BigReal electEnergy_ti_1 = 0;
      BigReal electEnergy_ti_2 = 0;)
      FULL(BigReal fullElectEnergy_ti_1 = 0; 
      BigReal fullElectEnergy_ti_2 = 0;) 
   )
  )
        
        
  const int i_upper = params->numAtoms[0];
  register const int j_upper = params->numAtoms[1];
  register int j;
  register int i;
  const CompAtom *p_0 = params->p[0];
  const CompAtom *p_1 = params->p[1];
  const CompAtomExt *pExt_0 = params->pExt[0];
  const CompAtomExt *pExt_1 = params->pExt[1];

  char * excl_flags_buff = 0;
  const int32 * full_excl = 0;
  const int32 * mod_excl = 0;

  plint *pairlistn_save;  int npairn;
  plint *pairlistx_save;  int npairx;
  plint *pairlistm_save;  int npairm;

  ALCH(
  plint *pairlistnA1_save;  int npairnA1;
  plint *pairlistxA1_save;  int npairxA1;
  plint *pairlistmA1_save;  int npairmA1;
  plint *pairlistnA2_save;  int npairnA2;
  plint *pairlistxA2_save;  int npairxA2;
  plint *pairlistmA2_save;  int npairmA2;
  )

  NBWORKARRAYSINIT(params->workArrays);

  int arraysize = j_upper+5;

  NBWORKARRAY(plint,pairlisti,arraysize)
  NBWORKARRAY(BigReal,r2list,arraysize)

  union { double f; int32 i[2]; } byte_order_test;
  byte_order_test.f = 1.0;  // should occupy high-order bits only
  int32 *r2iilist = (int32*)r2list + ( byte_order_test.i[0] ? 0 : 1 );

  if ( ! ( savePairlists || ! usePairlists ) ) arraysize = 0;


  // DMK - Atom Sort
  // 
  // Basic Idea: Determine the line between the center of masses of
  //   the two patches.  Project and then sort the lists of atoms
  //   along this line.  Then, as the pairlists are being generated
  //   for the atoms in the first atom list, use the sorted
  //   list to only add atoms from the second list that are between
  //   +/- ~cutoff from the atoms position on the line.
  #if NAMD_ComputeNonbonded_SortAtoms != 0 && ( 0 PAIR( + 1 ) )

  // NOTE: For the first atom list, just the sort values themselves need to be
  // calculated (a BigReal vs. SortEntry data structure).  However, a second
  // SortEntry array is needed to perform the merge sort on the second list of
  // atoms.  Use the atomSort_[0/1]_sortValues__ arrays to sort the second
  // list of atoms and then use the left over array to make a version of the
  // list that only includes non-fixed groups (fixg).

  NBWORKARRAY(SortEntry, atomSort_0_sortValues__, arraysize);
  NBWORKARRAY(SortEntry, atomSort_1_sortValues__, arraysize);
  NBWORKARRAY(BigReal, p_0_sortValues, arraysize);

  register SortEntry* p_1_sortValues = atomSort_0_sortValues__;
  register SortEntry* p_1_sortValues_fixg = atomSort_1_sortValues__;

  int p_0_sortValues_len = 0;
  int p_1_sortValues_len = 0;
  int p_1_sortValues_fixg_len = 0;

  // Calculate the distance between to projected points on the line that
  //   represents the cutoff distance.
  BigReal atomSort_windowRadius = sqrt(groupplcutoff2);

  if (savePairlists || !usePairlists) {

    register const BigReal projLineVec_x = params->projLineVec.x;
    register const BigReal projLineVec_y = params->projLineVec.y;
    register const BigReal projLineVec_z = params->projLineVec.z;

    // Calculate the sort values for the atoms in patch 1
    {
      register int nbgs = p_1->nonbondedGroupSize;
      register BigReal p_x = p_1->position.x;
      register BigReal p_y = p_1->position.y;
      register BigReal p_z = p_1->position.z;
      register int index = 0;

      for (register int j = nbgs; j < j_upper; j += nbgs) {

        // Set p_j_next to point to the atom for the next iteration and begin
        //   loading the 'nbgs' value for that atom.
        register const CompAtom* p_j_next = p_1 + j;
        nbgs = p_j_next->nonbondedGroupSize;

        // Calculate the distance along the projection vector
        // NOTE: If the vector from the origin to the point is 'A' and the vector
        //   between the patches is 'B' then to project 'A' onto 'B' we take the dot
        //   product of the two vectors 'A dot B' divided by the length of 'B'.
        // So... projection of A onto B = (A dot B) / length(B), however, note
        //   that length(B) == 1, therefore the sort value is simply (A dot B)
        register BigReal sortVal = COMPONENT_DOTPRODUCT(p,projLineVec);

        // Start loading the next iteration's atom's position
        p_x = p_j_next->position.x;
        p_y = p_j_next->position.y;
        p_z = p_j_next->position.z;

        // Store the caclulated sort value into the array of sort values
        register SortEntry* p_1_sortValStorePtr = p_1_sortValues + p_1_sortValues_len;
        p_1_sortValStorePtr->index = index;
        p_1_sortValStorePtr->sortValue = sortVal;
        p_1_sortValues_len++;

        // Update index for the next iteration
        index = j;       

      } // end while (j < j_upper)

      register BigReal sortVal = COMPONENT_DOTPRODUCT(p,projLineVec);

      register SortEntry* p_1_sortValStorePtr = p_1_sortValues + p_1_sortValues_len;
      p_1_sortValStorePtr->index = index;
      p_1_sortValStorePtr->sortValue = sortVal;
      p_1_sortValues_len++;
    }

    // NOTE: This list and another version of it with only non-fixed
    //   atoms will be used in place of grouplist and fixglist.
    #if 0   // Selection Sort
      sortEntries_selectionSort(p_1_sortValues, p_1_sortValues_len);
    #elif 0   // Bubble Sort
      sortEntries_bubbleSort(p_1_sortValues, p_1_sortValues_len);
    #else   // Merge Sort
      #if NAMD_ComputeNonbonded_SortAtoms_LessBranches == 0
        sortEntries_mergeSort_v1(p_1_sortValues, p_1_sortValues_fixg, p_1_sortValues_len);
      #else
        sortEntries_mergeSort_v2(p_1_sortValues, p_1_sortValues_fixg, p_1_sortValues_len);
      #endif
    #endif

    // Calculate the sort values for the atoms in patch 0
    {
      register int nbgs = p_0->nonbondedGroupSize;
      register BigReal p_x = p_0->position.x + offset_x;
      register BigReal p_y = p_0->position.y + offset_y;
      register BigReal p_z = p_0->position.z + offset_z;
      register int index = 0;

      for (register int i = nbgs; i < i_upper; i += nbgs) {

        // Set p_i_next to point to the atom for the next iteration and begin
        //   loading the 'nbgs' value for that atom.
        register const CompAtom* p_i_next = p_0 + i;
        nbgs = p_i_next->nonbondedGroupSize;

        // Calculate the distance along the projection vector
        register BigReal sortVal = COMPONENT_DOTPRODUCT(p,projLineVec);

        // Start loading the next iteration's atom's position
        p_x = p_i_next->position.x + offset_x;
        p_y = p_i_next->position.y + offset_y;
        p_z = p_i_next->position.z + offset_z;

        // Store the calculated sort value into the array of sort values
	register BigReal* p_0_sortValStorePtr = p_0_sortValues + index;
        *p_0_sortValStorePtr = sortVal;

        // Update index for the next iteration
        index = i;
      }

      register BigReal sortVal = COMPONENT_DOTPRODUCT(p,projLineVec);

      register BigReal* p_0_sortValStorePtr = p_0_sortValues + index;
      *p_0_sortValStorePtr = sortVal;

      p_0_sortValues_len = i_upper;
    }

  }  // end if (savePairlists || !usePairlists)

  #endif  // NAMD_ComputeNonbonded_SortAtoms != 0 && ( 0 PAIR( + 1 ) )

  // Atom Sort : The grouplist and fixglist arrays are not needed when the
  //   the atom sorting code is in use.
  #if ! (NAMD_ComputeNonbonded_SortAtoms != 0 && ( 0 PAIR( + 1 ) ) )
    NBWORKARRAY(plint,grouplist,arraysize);
    NBWORKARRAY(plint,fixglist,arraysize);
  #endif

  NBWORKARRAY(plint,goodglist,arraysize);
  NBWORKARRAY(plint,pairlistx,arraysize);
  NBWORKARRAY(plint,pairlistm,arraysize);
  NBWORKARRAY(plint,pairlist,arraysize);
  NBWORKARRAY(plint,pairlist2,arraysize);
  ALCH(
  NBWORKARRAY(plint,pairlistnA1,arraysize);
  NBWORKARRAY(plint,pairlistxA1,arraysize);
  NBWORKARRAY(plint,pairlistmA1,arraysize);
  NBWORKARRAY(plint,pairlistnA2,arraysize);
  NBWORKARRAY(plint,pairlistxA2,arraysize);
  NBWORKARRAY(plint,pairlistmA2,arraysize);
  )

  int fixg_upper = 0;
  int g_upper = 0;

  if ( savePairlists || ! usePairlists ) {

    #if NAMD_ComputeNonbonded_SortAtoms != 0 && ( 0 PAIR ( + 1 ) )

      // Create a sorted list of non-fixed groups
      register int fixg = 0;
      for (int tmpI = 0; tmpI < p_1_sortValues_len; tmpI++) {
        register SortEntry* p_1_sortEntry = p_1_sortValues + tmpI;
        register int p_1_index = p_1_sortEntry->index;
        if (!pExt_1[p_1_index].groupFixed) {
          register SortEntry* p_1_sortEntry_fixg = p_1_sortValues_fixg + p_1_sortValues_fixg_len;
          p_1_sortEntry_fixg->index = p_1_sortEntry->index;
          p_1_sortEntry_fixg->sortValue = p_1_sortEntry->sortValue;
          p_1_sortValues_fixg_len++;
	}
      }

    #else

      register int g = 0;
      for ( j = 0; j < j_upper; ++j ) {
        if ( p_1[j].nonbondedGroupSize ) {
          grouplist[g++] = j;
        }
      }
      g_upper = g;
      if ( g_upper ) grouplist[g_upper] = grouplist[g_upper-1];
      int fixg = 0;

      if ( fixedAtomsOn ) {
        for ( g = 0; g < g_upper; ++g ) {
          j = grouplist[g];
          if ( ! pExt_1[j].groupFixed ) {
            fixglist[fixg++] = j;
          }
        }
      }

      fixg_upper = fixg;
      if ( fixg_upper ) fixglist[fixg_upper] = fixglist[fixg_upper-1];

    #endif // NAMD_ComputeNonbonded_SortAtoms != 0

    pairlists.addIndex();
    pairlists.setIndexValue(i_upper);

  } else { // if ( savePairlists || ! usePairlists )

    if ( pairlists.getIndexValue() != i_upper )
      NAMD_bug("pairlist i_upper mismatch!");

  } // if ( savePairlists || ! usePairlists )

  SELF(
  int j_hgroup = 0;
  int g_lower = 0;
  int fixg_lower = 0;
  )
  int pairlistindex=0;
  int pairlistoffset=0;

  

#if ( SHORT( FAST( 1+ ) ) 0 )
#if ( PAIR( 1+ ) 0 )
    Force *f_0 = params->ff[0];
    Force *f_1 = params->ff[1];
#else
#define f_1 f_0
    NBWORKARRAY(Force,f_0,i_upper)
    memset( (void*) f_0, 0, i_upper * sizeof(Force) );
#endif
#endif
#if ( FULL( 1+ ) 0 )
#if ( PAIR( 1+ ) 0 )
    Force *fullf_0 = params->fullf[0];
    Force *fullf_1 = params->fullf[1];
#else
#define fullf_1 fullf_0
    NBWORKARRAY(Force,fullf_0,i_upper);
    memset( (void*) fullf_0, 0, i_upper * sizeof(Force) );
#endif
#endif
    

  int maxPart = params->numParts - 1;
  int groupCount = params->minPart;
  PAIR( 
    if ( savePairlists || ! usePairlists ) { 
      i = 0;
      pairlists.addIndex(); 
    } else {
      i = pairlists.getIndexValue();
    }
  )
   PAIR(for ( ; i < (i_upper);)) SELF(for ( i=0; i < (i_upper- 1);i++))
    {
    const CompAtom &p_i = p_0[i];
    const CompAtomExt &pExt_i = pExt_0[i];

    PAIR(if (savePairlists || ! usePairlists){)
    if ( p_i.hydrogenGroupSize ) {
      if ( groupCount ) {  // skip this group
        --groupCount;
        i += p_i.hydrogenGroupSize;
#ifdef ARCH_POWERPC
        __dcbt((void *) &(p_0[i]));
#endif
        SELF(--i;)
        continue;
      } else {  // compute this group
        groupCount = maxPart;
      }
    }
    PAIR(})

    register const BigReal p_i_x = p_i.position.x + offset_x;
    register const BigReal p_i_y = p_i.position.y + offset_y;
    register const BigReal p_i_z = p_i.position.z + offset_z;
#ifdef A2_QPX
    vector4double p_i_v = {p_i_x, p_i_y, p_i_z, 0.0};
#endif

    ALCH(const int p_i_partition = p_i.partition;)

    PPROF(
        const int p_i_partition = p_i.partition;
        int n1 = (int)floor((p_i.position.z-pressureProfileMin)*invThickness);
        pp_clamp(n1, pressureProfileSlabs);
        )

  SELF ( if ( p_i.hydrogenGroupSize ) j_hgroup = i + p_i.hydrogenGroupSize; )

  if ( savePairlists || ! usePairlists ) {

    #ifdef MEM_OPT_VERSION
    const ExclusionCheck *exclcheck = mol->get_excl_check_for_idx(pExt_i.exclId);        
    const int excl_min = pExt_i.id + exclcheck->min;
    const int excl_max = pExt_i.id + exclcheck->max;
    #else
    const ExclusionCheck *exclcheck = mol->get_excl_check_for_atom(pExt_i.id);
    const int excl_min = exclcheck->min;
    const int excl_max = exclcheck->max;
    #endif
    const char * excl_flags_var;
    if ( exclcheck->flags ) excl_flags_var = exclcheck->flags - excl_min;
    else {  // need to build list on the fly

    //TODO: Should change later!!!!!!!!!! --Chao Mei
    //Now just for the sake of passing compilation
    #ifndef MEM_OPT_VERSION 
      if ( excl_flags_buff ) {
        int nl,l;
        nl = full_excl[0] + 1;
        for ( l=1; l<nl; ++l ) excl_flags_buff[full_excl[l]] = 0;
        nl = mod_excl[0] + 1;
        for ( l=1; l<nl; ++l ) excl_flags_buff[mod_excl[l]] = 0;
      } else {
        excl_flags_buff = new char[mol->numAtoms];
        memset( (void*) excl_flags_buff, 0, mol->numAtoms);
      }
      int nl,l;
      full_excl = mol->get_full_exclusions_for_atom(pExt_i.id);
      nl = full_excl[0] + 1;
      for ( l=1; l<nl; ++l ) excl_flags_buff[full_excl[l]] = EXCHCK_FULL;
      mod_excl = mol->get_mod_exclusions_for_atom(pExt_i.id);
      nl = mod_excl[0] + 1;
      for ( l=1; l<nl; ++l ) excl_flags_buff[mod_excl[l]] = EXCHCK_MOD;
      excl_flags_var = excl_flags_buff;
    #endif

    }
    const char * const excl_flags = excl_flags_var;

  if ( p_i.nonbondedGroupSize ) {

    pairlistindex = 0;	// initialize with 0 elements
    pairlistoffset=0;
    const int groupfixed = ( fixedAtomsOn && pExt_i.groupFixed );

    // If patch divisions are not made by hydrogen groups, then
    // hydrogenGroupSize is set to 1 for all atoms.  Thus we can
    // carry on as if we did have groups - only less efficiently.
    // An optimization in this case is to not rebuild the temporary
    // pairlist but to include every atom in it.  This should be a
    // a very minor expense.

    SELF
    (
      if ( p_i.hydrogenGroupSize ) {
        // exclude child hydrogens of i
        // j_hgroup = i + p_i.hydrogenGroupSize;  (moved above)
        while ( g_lower < g_upper &&
                grouplist[g_lower] < j_hgroup ) ++g_lower;
        while ( fixg_lower < fixg_upper &&
                fixglist[fixg_lower] < j_hgroup ) ++fixg_lower;
      }
      // add all child or sister hydrogens of i
      for ( j = i + 1; j < j_hgroup; ++j ) {
	pairlist[pairlistindex++] = j;
      }
    )

    // add remaining atoms to pairlist via hydrogen groups
    register plint *pli = pairlist + pairlistindex;

    {
      // Atom Sort : Modify the values of g and gu based on the added information
      //   of the linear projections (sort values) information.
      #if NAMD_ComputeNonbonded_SortAtoms != 0 && ( 0 PAIR ( + 1 ) )

        register int g = 0;
        register BigReal p_i_sortValue = p_0_sortValues[i];
        const BigReal p_i_sortValue_plus_windowRadius = p_i_sortValue + atomSort_windowRadius;
        register SortEntry* sortValues = ( groupfixed ? p_1_sortValues_fixg : p_1_sortValues );

        // Find the actual gu (upper bound in sorted list for this outer-loop atom) based on the sort values
        register int lower = 0;
        register int upper = (groupfixed ? p_1_sortValues_fixg_len : p_1_sortValues_len);
        while ((upper - lower) > 1) {
          register int j = ((lower + upper) >> 1);
          register BigReal jSortVal = sortValues[j].sortValue;
          if (jSortVal < p_i_sortValue_plus_windowRadius) {
            lower = j;
          } else {
            upper = j;
	  }
	}
        const int gu = (sortValues[lower].sortValue >= p_i_sortValue_plus_windowRadius) ? lower : upper;

      #else

        register plint *gli = goodglist;
        const plint *glist = ( groupfixed ? fixglist : grouplist );
        SELF( const int gl = ( groupfixed ? fixg_lower : g_lower ); )
        const int gu = ( groupfixed ? fixg_upper : g_upper );
        register int g = PAIR(0) SELF(gl);

      #endif


      if ( g < gu ) {
	int hu = 0;
#if defined(__SSE2__) && ! defined(NAMD_DISABLE_SSE)
	if ( gu - g  >  6 ) { 

          #if NAMD_ComputeNonbonded_SortAtoms != 0 && ( 0 PAIR ( + 1 ) )
	    register SortEntry* sortEntry0 = sortValues + g;
	    register SortEntry* sortEntry1 = sortValues + g + 1;
            register int jprev0 = sortEntry0->index;
	    register int jprev1 = sortEntry1->index;
          #else
            register int jprev0 = glist[g    ];
	    register int jprev1 = glist[g + 1];
	  #endif
	  
	  register int j0;
	  register int j1;

          __m128d PJ_X_01 = _mm_set_pd(p_1[jprev1].position.x, p_1[jprev0].position.x);
          __m128d PJ_Y_01 = _mm_set_pd(p_1[jprev1].position.y, p_1[jprev0].position.y);
          __m128d PJ_Z_01 = _mm_set_pd(p_1[jprev1].position.z, p_1[jprev0].position.z);

          // these don't change here, so we could move them into outer scope
          const __m128d P_I_X = _mm_set1_pd(p_i_x);
          const __m128d P_I_Y = _mm_set1_pd(p_i_y);
          const __m128d P_I_Z = _mm_set1_pd(p_i_z);
 
	  g += 2;
	  for ( ; g < gu - 2; g +=2 ) {
	    // compute 1d distance, 2-way parallel	 
	    j0     =  jprev0;
	    j1     =  jprev1;

            __m128d T_01 = _mm_sub_pd(P_I_X, PJ_X_01);
            __m128d R2_01 = _mm_mul_pd(T_01, T_01);
            T_01 = _mm_sub_pd(P_I_Y, PJ_Y_01);
            R2_01 = _mm_add_pd(R2_01, _mm_mul_pd(T_01, T_01));
            T_01 = _mm_sub_pd(P_I_Z, PJ_Z_01);
            R2_01 = _mm_add_pd(R2_01, _mm_mul_pd(T_01, T_01));
	    
            #if NAMD_ComputeNonbonded_SortAtoms != 0 && ( 0 PAIR ( + 1 ) )
	      sortEntry0 = sortValues + g;
	      sortEntry1 = sortValues + g + 1;
              jprev0 = sortEntry0->index;
	      jprev1 = sortEntry1->index;
            #else
	      jprev0     =  glist[g  ];
	      jprev1     =  glist[g+1];
	    #endif
	   
            PJ_X_01 = _mm_set_pd(p_1[jprev1].position.x, p_1[jprev0].position.x);
            PJ_Y_01 = _mm_set_pd(p_1[jprev1].position.y, p_1[jprev0].position.y);
            PJ_Z_01 = _mm_set_pd(p_1[jprev1].position.z, p_1[jprev0].position.z);

            __align(16) double r2_01[2];
            _mm_store_pd(r2_01, R2_01); // 16-byte-aligned store

            // XXX these could possibly benefit from SSE-based conditionals
	    bool test0 = ( r2_01[0] < groupplcutoff2 );
	    bool test1 = ( r2_01[1] < groupplcutoff2 ); 
	    
	    //removing ifs benefits on many architectures
	    //as the extra stores will only warm the cache up
	    goodglist [ hu         ] = j0;
	    goodglist [ hu + test0 ] = j1;
	    
	    hu += test0 + test1;
	  }
	  g-=2;
	}
#elif defined (A2_QPX)
	if ( gu - g  >  6 ) { 
#if NAMD_ComputeNonbonded_SortAtoms != 0 && ( 0 PAIR ( + 1 ) )
	  register SortEntry* sortEntry0 = sortValues + g;
	  register SortEntry* sortEntry1 = sortValues + g + 1;
	  register int jprev0 = sortEntry0->index;
	  register int jprev1 = sortEntry1->index;
          #else
	  register int jprev0 = glist[g    ];
	  register int jprev1 = glist[g + 1];
          #endif
          
	  __dcbt ((void*)(p_1 + jprev0));
          register  int j0; 
          register  int j1;           
          vector4double    pj_v_0, pj_v_1; 
          vector4double    v_0, v_1;
          register BigReal r2_0, r2_1;
          
          pj_v_0 = vec_ld(jprev0 * sizeof(CompAtom), (BigReal *)p_1);
          pj_v_1 = vec_ld(jprev1 * sizeof(CompAtom), (BigReal *)p_1);  
          
          g += 2;
          for ( ; g < gu - 2; g +=2 ) {
            // compute 1d distance, 2-way parallel       
            j0     =  jprev0;
            j1     =  jprev1;

#if NAMD_ComputeNonbonded_SortAtoms != 0 && ( 0 PAIR ( + 1 ) )
	    sortEntry0 = sortValues + g;
	    sortEntry1 = sortValues + g + 1;
	    jprev0 = sortEntry0->index;
	    jprev1 = sortEntry1->index;
            #else
	    jprev0     =  glist[g  ];
	    jprev1     =  glist[g+1];
            #endif

            v_0 = vec_sub (p_i_v, pj_v_0);
            v_1 = vec_sub (p_i_v, pj_v_1);
            v_0 = vec_mul (v_0, v_0);
            v_1 = vec_mul (v_1, v_1);

            r2_0 = vec_extract(v_0, 0) + vec_extract(v_0, 1) + vec_extract(v_0, 2);
            r2_1 = vec_extract(v_1, 0) + vec_extract(v_1, 1) + vec_extract(v_1, 2);
            
            pj_v_0 = vec_ld(jprev0 * sizeof(CompAtom), (BigReal *)p_1);
            pj_v_1 = vec_ld(jprev1 * sizeof(CompAtom), (BigReal *)p_1);  
            
            size_t test0 = ( groupplcutoff2 >  r2_0 );
            size_t test1 = ( groupplcutoff2 >  r2_1 ); 
            
            //removing ifs benefits on many architectures
            //as the extra stores will only warm the cache up
            goodglist [ hu         ] = j0;
            goodglist [ hu + test0 ] = j1;
            
            hu += test0 + test1;
          }
          g-=2;
        }
#else
	if ( gu - g  >  6 ) { 

          #if NAMD_ComputeNonbonded_SortAtoms != 0 && ( 0 PAIR ( + 1 ) )
	    register SortEntry* sortEntry0 = sortValues + g;
	    register SortEntry* sortEntry1 = sortValues + g + 1;
            register int jprev0 = sortEntry0->index;
	    register int jprev1 = sortEntry1->index;
          #else
            register int jprev0 = glist[g    ];
	    register int jprev1 = glist[g + 1];
	  #endif
	  
	  register  int j0; 
	  register  int j1; 
	  
	  register  BigReal pj_x_0, pj_x_1; 
	  register  BigReal pj_y_0, pj_y_1; 
	  register  BigReal pj_z_0, pj_z_1; 
	  register  BigReal t_0, t_1, r2_0, r2_1;
	  
	  pj_x_0 = p_1[jprev0].position.x;
	  pj_x_1 = p_1[jprev1].position.x;  
	  
	  pj_y_0 = p_1[jprev0].position.y; 
	  pj_y_1 = p_1[jprev1].position.y;  
	  
	  pj_z_0 = p_1[jprev0].position.z; 
	  pj_z_1 = p_1[jprev1].position.z;
	  
	  g += 2;
	  for ( ; g < gu - 2; g +=2 ) {
	    // compute 1d distance, 2-way parallel	 
	    j0     =  jprev0;
	    j1     =  jprev1;
	    
	    t_0    =  p_i_x - pj_x_0;
	    t_1    =  p_i_x - pj_x_1;
	    r2_0   =  t_0 * t_0;
	    r2_1   =  t_1 * t_1;
	    
	    t_0    =  p_i_y - pj_y_0;
	    t_1    =  p_i_y - pj_y_1;
	    r2_0  +=  t_0 * t_0;
	    r2_1  +=  t_1 * t_1;
	    
	    t_0    =  p_i_z - pj_z_0;
	    t_1    =  p_i_z - pj_z_1;
	    r2_0  +=  t_0 * t_0;
	    r2_1  +=  t_1 * t_1;
	    
            #if NAMD_ComputeNonbonded_SortAtoms != 0 && ( 0 PAIR ( + 1 ) )
	      sortEntry0 = sortValues + g;
	      sortEntry1 = sortValues + g + 1;
              jprev0 = sortEntry0->index;
	      jprev1 = sortEntry1->index;
            #else
	      jprev0     =  glist[g  ];
	      jprev1     =  glist[g+1];
	    #endif
	    
	    pj_x_0     =  p_1[jprev0].position.x;
	    pj_x_1     =  p_1[jprev1].position.x;
	    pj_y_0     =  p_1[jprev0].position.y; 
	    pj_y_1     =  p_1[jprev1].position.y;
	    pj_z_0     =  p_1[jprev0].position.z; 
	    pj_z_1     =  p_1[jprev1].position.z;
	    
	    bool test0 = ( r2_0 < groupplcutoff2 );
	    bool test1 = ( r2_1 < groupplcutoff2 ); 
	    
	    //removing ifs benefits on many architectures
	    //as the extra stores will only warm the cache up
	    goodglist [ hu         ] = j0;
	    goodglist [ hu + test0 ] = j1;
	    
	    hu += test0 + test1;
	  }
	  g-=2;
	}
#endif
	
	for (; g < gu; g++) {

          #if NAMD_ComputeNonbonded_SortAtoms != 0 && ( 0 PAIR ( + 1 ) )
	    register SortEntry* sortEntry = sortValues + g;
            register int j = sortEntry->index;
          #else
	    int j = glist[g];
	  #endif

	  BigReal p_j_x = p_1[j].position.x;
	  BigReal p_j_y = p_1[j].position.y;
	  BigReal p_j_z = p_1[j].position.z;
	  
	  BigReal r2 = p_i_x - p_j_x;
	  r2 *= r2;
	  BigReal t2 = p_i_y - p_j_y;
	  r2 += t2 * t2;
	  t2 = p_i_z - p_j_z;
	  r2 += t2 * t2;

	  if ( r2 <= groupplcutoff2 ) 
	    goodglist[hu ++] = j; 
	}

	for ( int h=0; h<hu; ++h ) {
          int j = goodglist[h];
          int nbgs = p_1[j].nonbondedGroupSize;
	  pli[0] = j;   // copy over the next three in any case
	  pli[1] = j+1;
	  pli[2] = j+2;
          if ( nbgs & 4 ) {  // if nbgs > 3, assume nbgs <= 5
	    pli[3] = j+3;
	    pli[4] = j+4;
          }
          pli += nbgs;
	}
      }
    }

    pairlistindex = pli - pairlist;
    // make sure padded element on pairlist points to real data
    if ( pairlistindex ) {
       pairlist[pairlistindex] = pairlist[pairlistindex-1];
    } PAIR( else {  // skip empty loops if no pairs were found
       i += p_i.nonbondedGroupSize;
       continue;
    } )
  } // if i is nonbonded group parent
  SELF
    (
    // self-comparisions require list to be incremented
    // pair-comparisions use entire list (pairlistoffset is 0)
    else pairlistoffset++;
    )

    const int atomfixed = ( fixedAtomsOn && pExt_i.atomFixed );

    register plint *pli = pairlist2;

    plint *pairlistn = pairlists.newlist(j_upper + 5 + 5 ALCH( + 20 ));
    register plint *plin = pairlistn;

    INT(
    if ( pairInteractionOn ) {
      const int ifep_type = p_i.partition;
      if (pairInteractionSelf) {
        if (ifep_type != 1) { PAIR(++i;) continue; }
        for (int k=pairlistoffset; k<pairlistindex; k++) {
          j = pairlist[k];
          const int jfep_type = p_1[j].partition;
          // for pair-self, both atoms must be in group 1.
          if (jfep_type == 1) {
            *(pli++) = j;
          }
        }
      } else {
        if (ifep_type != 1 && ifep_type != 2) { PAIR(++i;) continue; }
        for (int k=pairlistoffset; k<pairlistindex; k++) {
          j = pairlist[k];
          const int jfep_type = p_1[j].partition;
          // for pair, must have one from each group.
          if (ifep_type + jfep_type == 3) {
            *(pli++) = j;
          }
        }
      }
      int npair2_int = pli - pairlist2;
      pli = pairlist2;
      for (int k=0; k<npair2_int; k++) {
        j = pairlist2[k];
        BigReal p_j_x = p_1[j].position.x;
	BigReal r2 = p_i_x - p_j_x;
	r2 *= r2;
        BigReal p_j_y = p_1[j].position.y;
	BigReal t2 = p_i_y - p_j_y;
	r2 += t2 * t2;
        BigReal p_j_z = p_1[j].position.z;
	t2 = p_i_z - p_j_z;
	r2 += t2 * t2;
	if ( ( ! (atomfixed && pExt_1[j].atomFixed) ) && (r2 <= plcutoff2) ) {
          int atom2 = pExt_1[j].id;
          if ( atom2 >= excl_min && atom2 <= excl_max ) *(pli++) = j;
          else *(plin++) = j;
        }
      }
    } else
    )
    if ( atomfixed ) {
      for (int k=pairlistoffset; k<pairlistindex; k++) {
        j = pairlist[k];
        BigReal p_j_x = p_1[j].position.x;
	BigReal r2 = p_i_x - p_j_x;
	r2 *= r2;
        BigReal p_j_y = p_1[j].position.y;
	BigReal t2 = p_i_y - p_j_y;
	r2 += t2 * t2;
        BigReal p_j_z = p_1[j].position.z;
	t2 = p_i_z - p_j_z;
	r2 += t2 * t2;
	if ( (! pExt_1[j].atomFixed) && (r2 <= plcutoff2) ) {
          int atom2 = pExt_1[j].id;
          if ( atom2 >= excl_min && atom2 <= excl_max ) *(pli++) = j;
          else *(plin++) = j;
        }
      }
    } else {
      int k = pairlistoffset;
      int ku = pairlistindex;
      if ( k < ku ) {
#if defined(__SSE2__) && ! defined(NAMD_DISABLE_SSE)
	if ( ku - k  >  6 ) { 	   
	  register  int jprev0 = pairlist [k    ];
	  register  int jprev1 = pairlist [k + 1];
	  
	  register  int j0; 
	  register  int j1; 

          __m128d PJ_X_01 = _mm_set_pd(p_1[jprev1].position.x, p_1[jprev0].position.x);
          __m128d PJ_Y_01 = _mm_set_pd(p_1[jprev1].position.y, p_1[jprev0].position.y);
          __m128d PJ_Z_01 = _mm_set_pd(p_1[jprev1].position.z, p_1[jprev0].position.z);

          // these don't change here, so we could move them into outer scope
          const __m128d P_I_X = _mm_set1_pd(p_i_x);
          const __m128d P_I_Y = _mm_set1_pd(p_i_y);
          const __m128d P_I_Z = _mm_set1_pd(p_i_z);
	  
	  int atom2_0 = pExt_1[jprev0].id;
	  int atom2_1 = pExt_1[jprev1].id;
	  
	  k += 2;
	  for ( ; k < ku - 2; k +=2 ) {
	    // compute 1d distance, 2-way parallel	 
	    j0     =  jprev0;
	    j1     =  jprev1;
	    
            __m128d T_01 = _mm_sub_pd(P_I_X, PJ_X_01);
            __m128d R2_01 = _mm_mul_pd(T_01, T_01);
            T_01 = _mm_sub_pd(P_I_Y, PJ_Y_01);
            R2_01 = _mm_add_pd(R2_01, _mm_mul_pd(T_01, T_01));
            T_01 = _mm_sub_pd(P_I_Z, PJ_Z_01);
            R2_01 = _mm_add_pd(R2_01, _mm_mul_pd(T_01, T_01));
	    
	    jprev0     =  pairlist[k];
	    jprev1     =  pairlist[k+1];
	    
            PJ_X_01 = _mm_set_pd(p_1[jprev1].position.x, p_1[jprev0].position.x);
            PJ_Y_01 = _mm_set_pd(p_1[jprev1].position.y, p_1[jprev0].position.y);
            PJ_Z_01 = _mm_set_pd(p_1[jprev1].position.z, p_1[jprev0].position.z);

            __align(16) double r2_01[2];
            _mm_store_pd(r2_01, R2_01); // 16-byte-aligned store
	    
	    if (r2_01[0] <= plcutoff2) {
	      if ( atom2_0 >= excl_min && atom2_0 <= excl_max ) 
		*(pli++) = j0;
	      else 
		*(plin++) = j0;
	    }
	    atom2_0 = pExt_1[jprev0].id;
	    
	    if (r2_01[1] <= plcutoff2) {
	      if ( atom2_1 >= excl_min && atom2_1 <= excl_max ) 
		*(pli++) = j1;
	      else 
		*(plin++) = j1;
	    }
	    atom2_1 = pExt_1[jprev1].id;	    
	  }
	  k-=2;
	}       
#elif defined(A2_QPX)
        if ( ku - k  >  6 ) {      
          register  int jprev0 = pairlist [k];
          register  int jprev1 = pairlist [k + 1];
          
          register  int j0; 
          register  int j1; 
          vector4double    pj_v_0, pj_v_1; 
          vector4double    v_0, v_1;
          BigReal          r2_0, r2_1;

          pj_v_0 = vec_ld(jprev0 * sizeof(CompAtom), (BigReal *)p_1);
          pj_v_1 = vec_ld(jprev1 * sizeof(CompAtom), (BigReal *)p_1);  
          
          int atom2_0 = pExt_1[jprev0].id;
          int atom2_1 = pExt_1[jprev1].id;
          
          k += 2;
          for ( ; k < ku - 2; k +=2 ) {
            // compute 1d distance, 2-way parallel       
            j0     =  jprev0;
            j1     =  jprev1;
          
            v_0 = vec_sub (p_i_v, pj_v_0);
            v_1 = vec_sub (p_i_v, pj_v_1);
            v_0 = vec_mul (v_0, v_0);
            v_1 = vec_mul (v_1, v_1);

            r2_0 = vec_extract(v_0, 0) + vec_extract(v_0, 1) + vec_extract(v_0, 2);
            r2_1 = vec_extract(v_1, 0) + vec_extract(v_1, 1) + vec_extract(v_1, 2);

            jprev0     =  pairlist[k];
            jprev1     =  pairlist[k+1];
            
            pj_v_0 = vec_ld(jprev0 * sizeof(CompAtom), (BigReal *)p_1);
            pj_v_1 = vec_ld(jprev1 * sizeof(CompAtom), (BigReal *)p_1);  
                    
            if (r2_0 <= plcutoff2) {
              if ( atom2_0 >= excl_min && atom2_0 <= excl_max ) 
                *(pli++) = j0;
              else 
                *(plin++) = j0;
            }
            atom2_0 = pExt_1[jprev0].id;
            
            if (r2_1 <= plcutoff2) {
              if ( atom2_1 >= excl_min && atom2_1 <= excl_max ) 
                *(pli++) = j1;
              else 
                *(plin++) = j1;
	    }
            atom2_1 = pExt_1[jprev1].id;            
          }
          k-=2;
        }      
#else
	if ( ku - k  >  6 ) { 	   
	  register  int jprev0 = pairlist [k];
	  register  int jprev1 = pairlist [k + 1];
	  
	  register  int j0; 
	  register  int j1; 
	  
	  register  BigReal pj_x_0, pj_x_1; 
	  register  BigReal pj_y_0, pj_y_1; 
	  register  BigReal pj_z_0, pj_z_1; 
	  register  BigReal t_0, t_1, r2_0, r2_1;
	  
	  pj_x_0 = p_1[jprev0].position.x;
	  pj_x_1 = p_1[jprev1].position.x;  
	  
	  pj_y_0 = p_1[jprev0].position.y; 
	  pj_y_1 = p_1[jprev1].position.y;  
	  
	  pj_z_0 = p_1[jprev0].position.z; 
	  pj_z_1 = p_1[jprev1].position.z;
	  
	  int atom2_0 = pExt_1[jprev0].id;
	  int atom2_1 = pExt_1[jprev1].id;
	  
	  k += 2;
	  for ( ; k < ku - 2; k +=2 ) {
	    // compute 1d distance, 2-way parallel	 
	    j0     =  jprev0;
	    j1     =  jprev1;
	    
	    t_0    =  p_i_x - pj_x_0;
	    t_1    =  p_i_x - pj_x_1;
	    r2_0   =  t_0 * t_0;
	    r2_1   =  t_1 * t_1;
	    
	    t_0    =  p_i_y - pj_y_0;
	    t_1    =  p_i_y - pj_y_1;
	    r2_0  +=  t_0 * t_0;
	    r2_1  +=  t_1 * t_1;
	    
	    t_0    =  p_i_z - pj_z_0;
	    t_1    =  p_i_z - pj_z_1;
	    r2_0  +=  t_0 * t_0;
	    r2_1  +=  t_1 * t_1;
	    
	    jprev0     =  pairlist[k];
	    jprev1     =  pairlist[k+1];
	    
	    pj_x_0     =  p_1[jprev0].position.x;
	    pj_x_1     =  p_1[jprev1].position.x;
	    pj_y_0     =  p_1[jprev0].position.y; 
	    pj_y_1     =  p_1[jprev1].position.y;
	    pj_z_0     =  p_1[jprev0].position.z; 
	    pj_z_1     =  p_1[jprev1].position.z;
	    
	    if (r2_0 <= plcutoff2) {
	      if ( atom2_0 >= excl_min && atom2_0 <= excl_max ) 
		*(pli++) = j0;
	      else 
		*(plin++) = j0;
	    }
	    atom2_0 = pExt_1[jprev0].id;
	    
	    if (r2_1 <= plcutoff2) {
	      if ( atom2_1 >= excl_min && atom2_1 <= excl_max ) 
		*(pli++) = j1;
	      else 
		*(plin++) = j1;
	     }
	    atom2_1 = pExt_1[jprev1].id;	    
	  }
	  k-=2;
	}       
#endif

	for (; k < ku; k++) {
	  int j = pairlist[k];
	  int atom2 = pExt_1[j].id;
	  
	  BigReal p_j_x = p_1[j].position.x;
	  BigReal p_j_y = p_1[j].position.y;
	  BigReal p_j_z = p_1[j].position.z;
	  
	  BigReal r2 = p_i_x - p_j_x;
	  r2 *= r2;
	  BigReal t2 = p_i_y - p_j_y;
	  r2 += t2 * t2;
	  t2 = p_i_z - p_j_z;
	  r2 += t2 * t2;
	  
	  if (r2 <= plcutoff2) {
	    if ( atom2 >= excl_min && atom2 <= excl_max ) 
	      *(pli++) = j;
	    else 
	      *(plin++) = j;
	  }
	}
      }
    }

    PAIR(
    if ( plin == pairlistn && pli == pairlist2 ) {
      ++i;
      continue;
    }
    pairlists.setIndexValue(i); 
    )

    plint *plix = pairlistx;
    plint *plim = pairlistm;
    ALCH(
    plint *plinA1 = pairlistnA1;
    plint *plixA1 = pairlistxA1;
    plint *plimA1 = pairlistmA1;
    plint *plinA2 = pairlistnA2;
    plint *plixA2 = pairlistxA2;
    plint *plimA2 = pairlistmA2;
    )

    int k;

#if 0 ALCH(+1)
    int unsortedNpairn = plin - pairlistn;
    plin = pairlistn;
    for ( k=0; k<unsortedNpairn; ++k ) {
      int j = pairlistn[k];
      switch(pswitchTable[p_i_partition + 3*(p_1[j].partition)]) {
        case 0:  *(plin++) = j; break;
        case 1:  *(plinA1++) = j; break;
        case 2:  *(plinA2++) = j; break;
      }
    }
#endif

    int npair2 = pli - pairlist2;
    // if ( npair2 ) pairlist2[npair2] = pairlist2[npair2-1];
    // removed code for implicit exclusions within hydrogen groups -JCP
    for (k=0; k < npair2; ++k ) {
      int j = pairlist2[k];
      int atom2 = pExt_1[j].id;
      int excl_flag = excl_flags[atom2];
      ALCH(int pswitch = pswitchTable[p_i_partition + 3*(p_1[j].partition)];)
      switch ( excl_flag ALCH( + 3 * pswitch)) {
      case 0:  *(plin++) = j;  break;
      case 1:  *(plix++) = j;  break;
      case 2:  *(plim++) = j;  break;
      ALCH(
      case 3:  *(plinA1++) = j; break;
      case 6:  *(plinA2++) = j; break;
      case 5:  *(plimA1++) = j; break;
      case 8:  *(plimA2++) = j; break;
      case 4:  *(plixA1++) = j; break;
      case 7:  *(plixA2++) = j; break;
      )
      }
    }

    npairn = plin - pairlistn;
    pairlistn_save = pairlistn;
    pairlistn_save[npairn] = npairn ? pairlistn_save[npairn-1] : -1;
    pairlists.newsize(npairn + 1);

    npairx = plix - pairlistx;
    pairlistx_save = pairlists.newlist();
    for ( k=0; k<npairx; ++k ) {
      pairlistx_save[k] = pairlistx[k];
    }
    pairlistx_save[k] = k ? pairlistx_save[k-1] : -1;
    pairlists.newsize(npairx + 1);

    npairm = plim - pairlistm;
    pairlistm_save = pairlists.newlist();
    for ( k=0; k<npairm; ++k ) {
      pairlistm_save[k] = pairlistm[k];
    }
    pairlistm_save[k] = k ? pairlistm_save[k-1] : -1;
    pairlists.newsize(npairm + 1);


#if 0 ALCH(+1)
#define PAIRLISTFROMARRAY(NPAIRS,PL1,PL2,PLSAVE) \
{ \
  (NPAIRS) = (PL2) - (PL1); \
  (PLSAVE) = pairlists.newlist(); \
  for ( k=0; k<(NPAIRS); ++k ) { \
    (PLSAVE)[k] = (PL1)[k]; \
  } \
  (PLSAVE)[k] = k ? (PLSAVE)[k-1] : -1; \
  pairlists.newsize((NPAIRS) + 1); \
}
  
    PAIRLISTFROMARRAY(npairnA1,pairlistnA1,plinA1,pairlistnA1_save);
    PAIRLISTFROMARRAY(npairxA1,pairlistxA1,plixA1,pairlistxA1_save);
    PAIRLISTFROMARRAY(npairmA1,pairlistmA1,plimA1,pairlistmA1_save);
    PAIRLISTFROMARRAY(npairnA2,pairlistnA2,plinA2,pairlistnA2_save);
    PAIRLISTFROMARRAY(npairxA2,pairlistxA2,plixA2,pairlistxA2_save);
    PAIRLISTFROMARRAY(npairmA2,pairlistmA2,plimA2,pairlistmA2_save);
#undef PAIRLISTFROMARRAY

#endif
    
    if ( ! savePairlists ) pairlists.reset();  // limit space usage
    PAIR( pairlists.addIndex(); )
    
	// PAIR( iout << i << " " << i_upper << " save\n" << endi;)
  } else { // if ( savePairlists || ! usePairlists )
	// PAIR( iout << i << " " << i_upper << " use\n" << endi;)

    pairlists.nextlist(&pairlistn_save,&npairn);  --npairn;
    pairlists.nextlist(&pairlistx_save,&npairx);  --npairx;
    pairlists.nextlist(&pairlistm_save,&npairm);  --npairm;
    ALCH(
    pairlists.nextlist(&pairlistnA1_save,&npairnA1);  --npairnA1;
    pairlists.nextlist(&pairlistxA1_save,&npairxA1);  --npairxA1;
    pairlists.nextlist(&pairlistmA1_save,&npairmA1);  --npairmA1;
    pairlists.nextlist(&pairlistnA2_save,&npairnA2);  --npairnA2;
    pairlists.nextlist(&pairlistxA2_save,&npairxA2);  --npairxA2;
    pairlists.nextlist(&pairlistmA2_save,&npairmA2);  --npairmA2;
    )

  } // if ( savePairlists || ! usePairlists )

    LES( BigReal *lambda_table_i =
			lambda_table + (lesFactor+1) * p_i.partition; )
    SPT( BigReal *lambda_table_i =
                        lambda_table + p_i.partition;)

    INT(
    const BigReal force_sign = (
      ( pairInteractionOn && ! pairInteractionSelf ) ?
        ( ( p_i.partition == 1 ) ? 1. : -1. ) : 0.
    );
      
    )
    BigReal scaling_qi = scaling;
    if (sptOn && node->molecule->get_spt_type(pExt_i.id) && sptScaleFactor2 == -1) {
       scaling_qi = sqrt(sptScaleFactor) * scaling; }
    const BigReal kq_i = COULOMB * p_i.charge * scaling_qi * dielectric_1;
    const LJTable::TableEntry * const lj_row =
		ljTable->table_row(p_i.vdwType);

#ifdef A2_QPX
    vector4double kq_iv = (vector4double)(0.0);
    vector4double f_i_v = (vector4double)(0.0);
    vector4double fullf_i_v = (vector4double)(0.0);
    vector4double full_cnst = (vector4double)(0.);
#define f_i_x       vec_extract(f_i_v, 0)
#define f_i_y       vec_extract(f_i_v, 1)
#define f_i_z       vec_extract(f_i_v, 2)
#define fullf_i_x   vec_extract(fullf_i_v, 0)
#define fullf_i_y   vec_extract(fullf_i_v, 1)
#define fullf_i_z   vec_extract(fullf_i_v, 2)
#else
    SHORT( FAST( BigReal f_i_x = 0.; ) )
    SHORT( FAST( BigReal f_i_y = 0.; ) )
    SHORT( FAST( BigReal f_i_z = 0.; ) )
    FULL( BigReal fullf_i_x = 0.; )
    FULL( BigReal fullf_i_y = 0.; )
    FULL( BigReal fullf_i_z = 0.; )
#endif

    int npairi;
    int k;

    npairi = pairlist_from_pairlist(ComputeNonbondedUtil::cutoff2,
	p_i_x, p_i_y, p_i_z, p_1, pairlistn_save, npairn, pairlisti,
	r2_delta, r2list);

// BEGIN NBTHOLE OF DRUDE MODEL
#if (FAST(1+)0)
    if (drudeNbthole && p_i.hydrogenGroupSize > 1) {
  
      Parameters *parameters = Node::Object()->parameters;
      const NbtholePairValue * const nbthole_array = parameters->nbthole_array;
      const int NumNbtholePairParams = parameters->NumNbtholePairParams;
      BigReal drudeNbtholeCut = simParams -> drudeNbtholeCut;
      BigReal drudeNbtholeCut2 = (drudeNbtholeCut * drudeNbtholeCut) + r2_delta;
      BigReal CC = COULOMB * scaling * dielectric_1;
      int kk;

      for (k = 0; k < npairi; k++) {
            if (r2list[k] > drudeNbtholeCut2) { continue; }
  
            const int j = pairlisti[k];
            const CompAtom& p_j = p_1[j];

            if ( p_j.hydrogenGroupSize < 2 ) { continue; }
  
        for (kk = 0;kk < NumNbtholePairParams; kk++) {

               if (((nbthole_array[kk].ind1 == p_i.vdwType) &&
                  (nbthole_array[kk].ind2 == p_j.vdwType)) ||
                  ((nbthole_array[kk].ind2 == p_i.vdwType) &&
                   (nbthole_array[kk].ind1 == p_j.vdwType))) {
                 break;
               }
        }
        if ( kk < NumNbtholePairParams ) {
    
                BigReal aprod = mol->GetAtomAlpha(pExt_0[i].id) * mol->GetAtomAlpha(pExt_1[j].id);
                const BigReal aa = nbthole_array[kk].tholeij * powf(aprod, -1.f/6);
                const BigReal qqaa = CC * p_0[i].charge * p_1[j].charge;
                const BigReal qqad = CC * p_0[i].charge * p_1[j+1].charge;
                const BigReal qqda = CC * p_0[i+1].charge * p_1[j].charge;
                const BigReal qqdd = CC * p_0[i+1].charge * p_1[j+1].charge;

                Vector raa = (p_0[i].position + params->offset) - p_1[j].position;
                Vector rad = (p_0[i].position + params->offset) - p_1[j+1].position;
                Vector rda = (p_0[i+1].position + params->offset) - p_1[j].position;
                Vector rdd = (p_0[i+1].position + params->offset) - p_1[j+1].position;

                BigReal raa_invlen = raa.rlength();  // reciprocal of length
                BigReal rad_invlen = rad.rlength();
                BigReal rda_invlen = rda.rlength();
                BigReal rdd_invlen = rdd.rlength();

                BigReal auaa = aa / raa_invlen;
                BigReal auad = aa / rad_invlen;
                BigReal auda = aa / rda_invlen;
                BigReal audd = aa / rdd_invlen;

                BigReal expauaa = exp(-auaa);
                BigReal expauad = exp(-auad);
                BigReal expauda = exp(-auda);
                BigReal expaudd = exp(-audd);

                BigReal polyauaa = 1 + 0.5*auaa;
                BigReal polyauad = 1 + 0.5*auad;
                BigReal polyauda = 1 + 0.5*auda;
                BigReal polyaudd = 1 + 0.5*audd;

                BigReal ethole = 0;
                ethole += qqaa * raa_invlen * (- polyauaa * expauaa);
                ethole += qqad * rad_invlen * (- polyauad * expauad);
                ethole += qqda * rda_invlen * (- polyauda * expauda);
                ethole += qqdd * rdd_invlen * (- polyaudd * expaudd);

                polyauaa = 1 + auaa*polyauaa;
                polyauad = 1 + auad*polyauad;
                polyauda = 1 + auda*polyauda;
                polyaudd = 1 + audd*polyaudd;

                BigReal raa_invlen3 = raa_invlen * raa_invlen * raa_invlen;
                BigReal rad_invlen3 = rad_invlen * rad_invlen * rad_invlen;
                BigReal rda_invlen3 = rda_invlen * rda_invlen * rda_invlen;
                BigReal rdd_invlen3 = rdd_invlen * rdd_invlen * rdd_invlen;

                BigReal dfaa = qqaa * raa_invlen3 * (polyauaa*expauaa);
                BigReal dfad = qqad * rad_invlen3 * (polyauad*expauad);
                BigReal dfda = qqda * rda_invlen3 * (polyauda*expauda);
                BigReal dfdd = qqdd * rdd_invlen3 * (polyaudd*expaudd);

                Vector faa = -dfaa * raa;
                Vector fad = -dfad * rad;
                Vector fda = -dfda * rda;
                Vector fdd = -dfdd * rdd;

                params->ff[0][i] += faa + fad;
                params->ff[0][i+1] += fda + fdd;
                params->ff[1][j] -= faa + fda;
                params->ff[1][j+1] -= fad + fdd;

                reduction[virialIndex_XX] += faa.x * raa.x + fad.x * rad.x
                + fda.x * rda.x + fdd.x * rdd.x;
                reduction[virialIndex_XY] += faa.x * raa.y + fad.x * rad.y
                + fda.x * rda.y + fdd.x * rdd.y;
                reduction[virialIndex_XZ] += faa.x * raa.z + fad.x * rad.z
                + fda.x * rda.z + fdd.x * rdd.z;
                reduction[virialIndex_YX] += faa.y * raa.x + fad.y * rad.x
                + fda.y * rda.x + fdd.y * rdd.x;
                reduction[virialIndex_YY] += faa.y * raa.y + fad.y * rad.y
                + fda.y * rda.y + fdd.y * rdd.y;
                reduction[virialIndex_YZ] += faa.y * raa.z + fad.y * rad.z
                + fda.y * rda.z + fdd.y * rdd.z;
                reduction[virialIndex_ZX] += faa.z * raa.x + fad.z * rad.x
                + fda.z * rda.x + fdd.z * rdd.x;
                reduction[virialIndex_ZY] += faa.z * raa.y + fad.z * rad.y
                + fda.z * rda.y + fdd.z * rdd.y;
                reduction[virialIndex_ZZ] += faa.z * raa.z + fad.z * rad.z
                + fda.z * rda.z + fdd.z * rdd.z;
                reduction[electEnergyIndex] += ethole;
          }
       }
    }
#endif
// END NBTHOLE OF DRUDE MODEL
    
    // BEGIN LA
#if (FAST(1+)0)
    if (doLoweAndersen && p_i.hydrogenGroupSize) {
	BigReal loweAndersenCutoff = simParams->loweAndersenCutoff;
	BigReal loweAndersenCutoff2 = (loweAndersenCutoff * loweAndersenCutoff) + r2_delta;
	BigReal loweAndersenProb = simParams->loweAndersenRate * (simParams->dt * simParams->nonbondedFrequency) * 0.001; // loweAndersenRate is in 1/ps
	const bool loweAndersenUseCOMvelocity = (simParams->rigidBonds != RIGID_NONE);
	const BigReal kbT = BOLTZMANN * (simParams->loweAndersenTemp);
	const BigReal dt_inv = TIMEFACTOR / (simParams->dt * simParams->nonbondedFrequency);
	//const BigReal dt_inv = 1.0 / simParams->dt;
	//BigReal kbT = 8.3145e-7 * (simParams->loweAndersenTemp); // in A^2/fs^2 * K
	
	const CompAtom* v_0 = params->v[0];
	const CompAtom* v_1 = params->v[1];
	const CompAtom& v_i = v_0[i];
	Mass mass_i = v_i.charge;
	Velocity vel_i = v_i.position;
	Position pos_i = p_i.position;
	int atom_i = pExt_0[i].id;
	
	if (loweAndersenUseCOMvelocity) {
	    Mass mass = 0;
	    Velocity vel = 0;
	    Position pos = 0;
	    for (int l = 0; l < p_i.hydrogenGroupSize; l++) {
		vel += v_0[i+l].charge * v_0[i+l].position;
		pos += v_0[i+l].charge * p_0[i+l].position;
		mass += v_0[i+l].charge;
	    }
	    vel_i = vel / mass;
	    pos_i = pos / mass;
	    mass_i = mass;
	}
	
	// Note v[0].charge is actually mass
	//Random rand(CkMyPe()); // ?? OK ?? NO!!!!
	Random *rand = Node::Object()->rand;
	for (k = 0; k < npairi; k++) {
	    if (r2list[k] > loweAndersenCutoff2) { continue; }
		
	    const int j = pairlisti[k];
	    const CompAtom& v_j = v_1[j];
	    const CompAtom& p_j = p_1[j];
		
	    if (!p_j.hydrogenGroupSize) { continue; }
	    if (rand->uniform() > loweAndersenProb) { continue; }
		
	    Mass mass_j = v_j.charge;
	    Velocity vel_j = v_j.position;
	    Position pos_j = p_j.position;
	    int atom_j = pExt_1[j].id;
		
	    if (loweAndersenUseCOMvelocity) {
		Mass mass = 0;
		Velocity vel = 0;
		Position pos = 0;
		for (int l = 0; l < p_j.hydrogenGroupSize; l++) {
		    vel += v_1[j+l].charge * v_1[j+l].position;
		    pos += v_1[j+l].charge * p_1[j+l].position;
		    mass += v_1[j+l].charge;
		}
		vel_j = vel / mass;
		pos_j = pos / mass;
		mass_j = mass;
	    }
		
	    //Velocity deltaV = v_i.position - v_j.position;
	    Velocity deltaV = vel_i - vel_j;
	    Mass mu_ij = (mass_i * mass_j)/(mass_i + mass_j);
	    //Vector sep = (p_i.position + params->offset) - p_j.position;
	    Vector sep = (pos_i + params->offset) - pos_j;
	    Vector sigma_ij = sep.unit();
	    BigReal lambda = rand->gaussian() * sqrt(kbT / mu_ij);
	    Force force = mu_ij * dt_inv * (lambda - (deltaV * sigma_ij)) * sigma_ij;
		
	    //DebugM(1, "atom1 atom2 = " << atom1 << " " << atom2 << " lambda = " << lambda << " force = " << force << " deltaP = " << sep.length() << " sqrt(r2) = " << sqrt(r2list[k]) << "\n");
	    //DebugM(1, "atom1 atom2 = " << atom_i << " " << atom_j << " mass1 mass2 = " << mass_i << " " << mass_j << " hgrp1 hgrp2 " << p_i.hydrogenGroupSize << " " << p_j.hydrogenGroupSize << "\n");
		
	    if (loweAndersenUseCOMvelocity) {
		BigReal inv_mass_i = 1.0 / mass_i;
		BigReal inv_mass_j = 1.0 / mass_j;
		for (int l = 0; l < p_i.hydrogenGroupSize; l++) {
		    params->ff[0][i+l] += (v_0[i+l].charge * inv_mass_i) * force;
		}
		for (int l = 0; l < p_j.hydrogenGroupSize; l++) {
		    params->ff[1][j+l] -= (v_1[j+l].charge * inv_mass_j) * force;
		}
	    } else {
		params->ff[0][i] += force;
		params->ff[1][j] -= force;
	    }
		
	    reduction[virialIndex_XX] += force.x * sep.x;
	    reduction[virialIndex_XY] += force.x * sep.y;
	    reduction[virialIndex_XZ] += force.x * sep.z;
	    reduction[virialIndex_YX] += force.y * sep.x;
	    reduction[virialIndex_YY] += force.y * sep.y;
	    reduction[virialIndex_YZ] += force.y * sep.z;
	    reduction[virialIndex_ZX] += force.z * sep.x;
	    reduction[virialIndex_ZY] += force.z * sep.y;
	    reduction[virialIndex_ZZ] += force.z * sep.z;
	}
    }
#endif
    // END LA

#define NORMAL(X) X
#define EXCLUDED(X)
#define MODIFIED(X)
#include  "ComputeNonbondedBase2.h"
#undef NORMAL
#undef EXCLUDED
#undef MODIFIED

    npairi = pairlist_from_pairlist(ComputeNonbondedUtil::cutoff2,
	p_i_x, p_i_y, p_i_z, p_1, pairlistm_save, npairm, pairlisti,
	r2_delta, r2list);
    exclChecksum += npairi;

#define NORMAL(X)
#define EXCLUDED(X)
#define MODIFIED(X) X
#include  "ComputeNonbondedBase2.h"
#undef NORMAL
#undef EXCLUDED
#undef MODIFIED

#ifdef FULLELECT
    npairi = pairlist_from_pairlist(ComputeNonbondedUtil::cutoff2,
	p_i_x, p_i_y, p_i_z, p_1, pairlistx_save, npairx, pairlisti,
	r2_delta, r2list);
    exclChecksum += npairi;

#undef FAST
#define FAST(X)
#define NORMAL(X)
#define EXCLUDED(X) X
#define MODIFIED(X)
#include  "ComputeNonbondedBase2.h"
#undef FAST
#ifdef SLOWONLY
  #define FAST(X)
#else
  #define FAST(X) X
#endif
#undef NORMAL
#undef EXCLUDED
#undef MODIFIED
#else
    exclChecksum += npairx;
#endif

#if 0 ALCH(+1)   // nonbondedbase2 for alchemical-separeted pairlists

    #undef ALCHPAIR
    #define ALCHPAIR(X) X
    #undef NOT_ALCHPAIR
    #define NOT_ALCHPAIR(X)
    BigReal myLambda; FEP(BigReal myLambda2;)
    BigReal myElecLambda;  FEP(BigReal myElecLambda2;)
    BigReal myVdwLambda; FEP(BigReal myVdwLambda2;)
    BigReal myVdwShift; FEP(BigReal myVdwShift2;)
    BigReal alch_vdw_energy; BigReal alch_vdw_force;
    FEP(BigReal alch_vdw_energy_2; BigReal alch_vdw_energy_2_Left;) TI(BigReal alch_vdw_dUdl;)
    BigReal shiftedElec; BigReal shiftedElecForce;
    
    /********************************************************************/
    /*******NONBONDEDBASE2 FOR NORMAL INTERACTIONS SCALED BY LAMBDA******/
    /********************************************************************/
    #define NORMAL(X) X
    #define EXCLUDED(X)
    #define MODIFIED(X)

    #define ALCH1(X) X
    #define ALCH2(X)
    npairi = pairlist_from_pairlist(ComputeNonbondedUtil::cutoff2,
            p_i_x, p_i_y, p_i_z, p_1, pairlistnA1_save, npairnA1, pairlisti,
            r2_delta, r2list);
    #include  "ComputeNonbondedBase2.h" // normal, direction 'up'
    #undef ALCH1
    #undef ALCH2

    #define ALCH1(X)
    #define ALCH2(X) X
    npairi = pairlist_from_pairlist(ComputeNonbondedUtil::cutoff2,
            p_i_x, p_i_y, p_i_z, p_1, pairlistnA2_save, npairnA2, pairlisti,
            r2_delta, r2list);
    #include  "ComputeNonbondedBase2.h" // normal, direction 'down'
    #undef ALCH1
    #undef ALCH2

    #undef NORMAL
    #undef EXCLUDED
    #undef MODIFIED
    
    /********************************************************************/
    /******NONBONDEDBASE2 FOR MODIFIED INTERACTIONS SCALED BY LAMBDA*****/
    /********************************************************************/
    #define NORMAL(X)
    #define EXCLUDED(X)
    #define MODIFIED(X) X

    #define ALCH1(X) X
    #define ALCH2(X)
    npairi = pairlist_from_pairlist(ComputeNonbondedUtil::cutoff2,
            p_i_x, p_i_y, p_i_z, p_1, pairlistmA1_save, npairmA1, pairlisti,
            r2_delta, r2list);
        exclChecksum += npairi;
    #include  "ComputeNonbondedBase2.h" // modified, direction 'up'
    #undef ALCH1
    #undef ALCH2

    #define ALCH1(X)
    #define ALCH2(X) X
    npairi = pairlist_from_pairlist(ComputeNonbondedUtil::cutoff2,
            p_i_x, p_i_y, p_i_z, p_1, pairlistmA2_save, npairmA2, pairlisti,
            r2_delta, r2list);
        exclChecksum += npairi;
    #include  "ComputeNonbondedBase2.h" // modified, direction 'down'
    #undef ALCH1
    #undef ALCH2

    #undef NORMAL
    #undef EXCLUDED
    #undef MODIFIED
    
    /********************************************************************/
    /******NONBONDEDBASE2 FOR EXCLUDED INTERACTIONS SCALED BY LAMBDA*****/
    /********************************************************************/
    #ifdef FULLELECT
    #undef FAST
    #define FAST(X)
    #define NORMAL(X)
    #define EXCLUDED(X) X
    #define MODIFIED(X)

    #define ALCH1(X) X
    #define ALCH2(X)
    npairi = pairlist_from_pairlist(ComputeNonbondedUtil::cutoff2,
            p_i_x, p_i_y, p_i_z, p_1, pairlistxA1_save, npairxA1, pairlisti,
            r2_delta, r2list);
        exclChecksum += npairi;
    #include  "ComputeNonbondedBase2.h"  //excluded, direction 'up'
    #undef ALCH1
    #undef ALCH2

    #define ALCH1(X)
    #define ALCH2(X) X
    npairi = pairlist_from_pairlist(ComputeNonbondedUtil::cutoff2,
            p_i_x, p_i_y, p_i_z, p_1, pairlistxA2_save, npairxA2, pairlisti,
            r2_delta, r2list);
        exclChecksum += npairi;
    #include  "ComputeNonbondedBase2.h"  //excluded, direction 'down'
    #undef ALCH1
    #undef ALCH2

    #undef FAST
    #ifdef SLOWONLY
      #define FAST(X)
    #else
      #define FAST(X) X
    #endif
    #undef NORMAL
    #undef EXCLUDED
    #undef MODIFIED
    #else
        exclChecksum += npairxA1 + npairxA2;
    #endif

    #undef ALCHPAIR
    #define ALCHPAIR(X)
    #undef NOT_ALCHPAIR
    #define NOT_ALCHPAIR(X) X

#endif // end nonbondedbase2.h includes for alchemical pairlists

#ifdef NETWORK_PROGRESS
    CkNetworkProgress();
#endif

#ifdef ARCH_POWERPC
    //data cache block touch the position structure
    __dcbt ((void *) &(p_0[i+1]));
    //__prefetch_by_load ((void *)&(groupCount));
#endif

#ifndef NAMD_CUDA
    SHORT( FAST( f_0[i].x += f_i_x; ) )
    SHORT( FAST( f_0[i].y += f_i_y; ) )
    SHORT( FAST( f_0[i].z += f_i_z; ) )
    FULL( fullf_0[i].x += fullf_i_x; )
    FULL( fullf_0[i].y += fullf_i_y; )
    FULL( fullf_0[i].z += fullf_i_z; )
#endif
PAIR(
    if ( savePairlists || ! usePairlists ) {
      i++;
    } else {
      i = pairlists.getIndexValue();
    }
)

	// PAIR( iout << i << " " << i_upper << " end\n" << endi;)
  } // for i

  // PAIR(iout << "++++++++\n" << endi;)
  PAIR( if ( savePairlists ) { pairlists.setIndexValue(i); } )

#ifdef A2_QPX
    BigReal  virial_xx   =  vec_extract (virial_v0, 0);
    BigReal  virial_xy   =  vec_extract (virial_v0, 1);
    BigReal  virial_xz   =  vec_extract (virial_v0, 2);
    BigReal  virial_yy   =  vec_extract (virial_v1, 1);
    BigReal  virial_yz   =  vec_extract (virial_v1, 2);
    BigReal  virial_zz   =  vec_extract (virial_v2, 2);
#endif
#ifdef f_1
#undef f_1
#endif
#if ( SHORT( FAST( 1+ ) ) 0 )
#if ( SELF( 1+ ) 0 )
  {
    Force *patch_f_0 = params->ff[0];
#ifndef NAMD_CUDA
#ifndef ARCH_POWERPC 
#pragma ivdep
#endif
    for ( int i = 0; i < i_upper; ++i ) {
      patch_f_0[i].x += f_0[i].x;
      patch_f_0[i].y += f_0[i].y;
      patch_f_0[i].z += f_0[i].z;
      virial_xx += f_0[i].x * p_0[i].position.x;
      virial_xy += f_0[i].x * p_0[i].position.y;
      virial_xz += f_0[i].x * p_0[i].position.z;
      virial_yy += f_0[i].y * p_0[i].position.y;
      virial_yz += f_0[i].y * p_0[i].position.z;
      virial_zz += f_0[i].z * p_0[i].position.z;
    }
#endif
  }
#endif
#endif

#ifdef A2_QPX
    BigReal  fullElectVirial_xx  =  vec_extract(fullvirial_v0, 0);
    BigReal  fullElectVirial_xy  =  vec_extract(fullvirial_v0, 1);
    BigReal  fullElectVirial_xz  =  vec_extract(fullvirial_v0, 2);    
    BigReal  fullElectVirial_yy  =  vec_extract(fullvirial_v1, 1);
    BigReal  fullElectVirial_yz  =  vec_extract(fullvirial_v1, 2);
    BigReal  fullElectVirial_zz  =  vec_extract(fullvirial_v2, 2);
#endif
#ifdef fullf_1
#undef fullf_1
#endif
#if ( FULL( 1+ ) 0 )
#if ( SELF( 1+ ) 0 )
  {
    Force *patch_fullf_0 = params->fullf[0];
#ifndef NAMD_CUDA
#ifndef ARCH_POWERPC 
#pragma ivdep
#endif
    for ( int i = 0; i < i_upper; ++i ) {
      patch_fullf_0[i].x += fullf_0[i].x;
      patch_fullf_0[i].y += fullf_0[i].y;
      patch_fullf_0[i].z += fullf_0[i].z;
      fullElectVirial_xx += fullf_0[i].x * p_0[i].position.x;
      fullElectVirial_xy += fullf_0[i].x * p_0[i].position.y;
      fullElectVirial_xz += fullf_0[i].x * p_0[i].position.z;
      fullElectVirial_yy += fullf_0[i].y * p_0[i].position.y;
      fullElectVirial_yz += fullf_0[i].y * p_0[i].position.z;
      fullElectVirial_zz += fullf_0[i].z * p_0[i].position.z;
    }
#endif
  }
#endif
#endif

#ifndef NAMD_CUDA
  reduction[exclChecksumIndex] += exclChecksum;
#endif
  FAST
  (
   // JLai
  ENERGY( reduction[vdwEnergyIndex] += vdwEnergy;
	  GO( reduction[groLJEnergyIndex] += groLJEnergy;
	      reduction[groGaussEnergyIndex] += groGaussEnergy;
	      reduction[goNativeEnergyIndex] += goEnergyNative;
	      reduction[goNonnativeEnergyIndex] += goEnergyNonnative; ) )
  SHORT
  (
  ENERGY( reduction[electEnergyIndex] += electEnergy; )
  )
  ALCH
  (
    TI(reduction[vdwEnergyIndex_ti_1] += vdwEnergy_ti_1;) 
    TI(reduction[vdwEnergyIndex_ti_2] += vdwEnergy_ti_2;) 
    FEP( reduction[vdwEnergyIndex_s] += vdwEnergy_s; )
    FEP( reduction[vdwEnergyIndex_s_Left] += vdwEnergy_s_Left; )
  SHORT
  (
    FEP( reduction[electEnergyIndex_s] += electEnergy_s; )
    TI(reduction[electEnergyIndex_ti_1] += electEnergy_ti_1;) 
    TI(reduction[electEnergyIndex_ti_2] += electEnergy_ti_2;) 
  )
  )
  SHORT
  (
  reduction[virialIndex_XX] += virial_xx;
  reduction[virialIndex_XY] += virial_xy;
  reduction[virialIndex_XZ] += virial_xz;
  reduction[virialIndex_YX] += virial_xy;
  reduction[virialIndex_YY] += virial_yy;
  reduction[virialIndex_YZ] += virial_yz;
  reduction[virialIndex_ZX] += virial_xz;
  reduction[virialIndex_ZY] += virial_yz;
  reduction[virialIndex_ZZ] += virial_zz;
  )
  )
  FULL
  (
  ENERGY( reduction[fullElectEnergyIndex] += fullElectEnergy; )
  ALCH
  (
    FEP( reduction[fullElectEnergyIndex_s] += fullElectEnergy_s; )
    TI(reduction[fullElectEnergyIndex_ti_1] += fullElectEnergy_ti_1;) 
    TI(reduction[fullElectEnergyIndex_ti_2] += fullElectEnergy_ti_2;) 
  )
  reduction[fullElectVirialIndex_XX] += fullElectVirial_xx;
  reduction[fullElectVirialIndex_XY] += fullElectVirial_xy;
  reduction[fullElectVirialIndex_XZ] += fullElectVirial_xz;
  reduction[fullElectVirialIndex_YX] += fullElectVirial_xy;
  reduction[fullElectVirialIndex_YY] += fullElectVirial_yy;
  reduction[fullElectVirialIndex_YZ] += fullElectVirial_yz;
  reduction[fullElectVirialIndex_ZX] += fullElectVirial_xz;
  reduction[fullElectVirialIndex_ZY] += fullElectVirial_yz;
  reduction[fullElectVirialIndex_ZZ] += fullElectVirial_zz;
  )

  delete [] excl_flags_buff;

}

