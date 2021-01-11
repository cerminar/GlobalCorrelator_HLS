#ifndef FIRMWARE_PFTKEGALGO_H
#define FIRMWARE_PFTKEGALGO_H

// FIXME: move to data.h
// FIXME: tune on the basis of the region occupancy
// FIXME: this needs to go somewhere else
#define SELHWQUAL 4

#ifndef REG_HGCal

  #ifndef CMSSW_GIT_HASH
    #warning "REG_HGCal is not #defined, but this algorithm has only been tested there"
  #endif
#endif

// #include "pfalgo_common.h"
#include "data.h"

void pftkegalgo(const EmCaloObj emcalo[NEMCALO], const TkObj track[NTRACK],
  EGIsoParticle photons[NEM_EGOUT],
  EGIsoEleParticle eles[NEM_EGOUT]
) ;

void link_emCalo2tk(const EmCaloObj emcalo[NEMCALOSEL_EGIN], const TkObj track[NTRACK], ap_uint<NTRACK> emCalo2tk_bit[NEMCALOSEL_EGIN]);

void link_emCalo2emCalo(const EmCaloObj emcalo[NEMCALOSEL_EGIN], ap_uint<NEMCALOSEL_EGIN> emCalo2emcalo_bit[NEMCALOSEL_EGIN]);

void sel_emCalo(const EmCaloObj emcalo[NEMCALO], EmCaloObj emcalo_sel[NEMCALOSEL_EGIN]);


#endif
