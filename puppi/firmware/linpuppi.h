#ifndef FIRMWARE_LINPUPPI_H
#define FIRMWARE_LINPUPPI_H

#include <cmath>
#ifdef CMSSW_GIT_HASH
#include "../dataformmats/pf.h"
#else
#include "../../dataformats/pf.h"
#endif

#if defined(PACKING_DATA_SIZE) && defined(PACKING_NCHANN)
#include "../../dataformats/l1pf_encoding.h"
#endif

int dr2_int(eta_t eta1, phi_t phi1, eta_t eta2, eta_t phi2);

// charged
void linpuppi_chs(z0_t pvZ0, const PFChargedObj pfch[NTRACK], PuppiObj outallch[NTRACK]) ;

// neutrals, in the tracker
void linpuppiNoCrop(const TkObj track[NTRACK], z0_t pvZ0, const PFNeutralObj pfallne[NALLNEUTRALS], PuppiObj outallne[NALLNEUTRALS]) ;
void linpuppi(const TkObj track[NTRACK], z0_t pvZ0, const PFNeutralObj pfallne[NALLNEUTRALS], PuppiObj outselne[NNEUTRALS]) ;

// streaming versions, taking one object at a time 
struct linpuppi_refobj { ap_uint<17> pt2_shift; eta_t hwEta; phi_t hwPhi; };
linpuppi_refobj linpuppi_prepare_track(const TkObj & track, z0_t pvZ0);
PuppiObj linpuppi_one(const PFNeutralObj & in, const linpuppi_refobj sel_track[NTRACK]);
PuppiObj linpuppi_chs_one(const PFChargedObj pfch, z0_t pvZ0) ;
void linpuppiNoCrop_streamed(const TkObj track[NTRACK], z0_t pvZ0, const PFNeutralObj pfallne[NALLNEUTRALS], PuppiObj outallne[NALLNEUTRALS]) ;
void linpuppi_chs_streamed(z0_t pvZ0, const PFChargedObj pfch[NTRACK], PuppiObj outallch[NTRACK]) ;

// neutrals, forward
void fwdlinpuppi(const HadCaloObj caloin[NCALO], PuppiObj pfselne[NNEUTRALS]);
void fwdlinpuppiNoCrop(const HadCaloObj caloin[NCALO], PuppiObj pfallne[NCALO]);

#define LINPUPPI_DATA_SIZE_IN 72
#define LINPUPPI_DATA_SIZE_OUT 64
#define LINPUPPI_DATA_SIZE_FWD 64
#define LINPUPPI_NCHANN_IN (1+NTRACK+NALLNEUTRALS)
#define LINPUPPI_NCHANN_OUTNC (NALLNEUTRALS)
#define LINPUPPI_NCHANN_OUT (NNEUTRALS)
#define LINPUPPI_CHS_NCHANN_IN (1+NTRACK)
#define LINPUPPI_CHS_NCHANN_OUT (NTRACK)
#define LINPUPPI_NCHANN_FWDNC (NCALO)
#define LINPUPPI_NCHANN_FWD (NNEUTRALS)

void packed_fwdlinpuppi(const ap_uint<LINPUPPI_DATA_SIZE_FWD> input[LINPUPPI_NCHANN_FWDNC], ap_uint<LINPUPPI_DATA_SIZE_FWD> output[LINPUPPI_NCHANN_FWD]) ;
void packed_fwdlinpuppiNoCrop(const ap_uint<LINPUPPI_DATA_SIZE_FWD> input[LINPUPPI_NCHANN_FWDNC], ap_uint<LINPUPPI_DATA_SIZE_FWD> output[LINPUPPI_NCHANN_FWDNC]) ;

void packed_linpuppi_chs(const ap_uint<LINPUPPI_DATA_SIZE_IN> input[LINPUPPI_CHS_NCHANN_IN], ap_uint<LINPUPPI_DATA_SIZE_OUT> output[LINPUPPI_CHS_NCHANN_OUT]);
void packed_linpuppi(const ap_uint<LINPUPPI_DATA_SIZE_IN> input[LINPUPPI_NCHANN_IN], ap_uint<LINPUPPI_DATA_SIZE_OUT> output[LINPUPPI_NCHANN_OUT]);
void packed_linpuppiNoCrop(const ap_uint<LINPUPPI_DATA_SIZE_IN> input[LINPUPPI_NCHANN_IN], ap_uint<LINPUPPI_DATA_SIZE_OUT> output[LINPUPPI_NCHANN_OUTNC]);

void linpuppi_pack_in(const TkObj track[NTRACK], z0_t pvZ0, const PFNeutralObj pfallne[NALLNEUTRALS], ap_uint<LINPUPPI_DATA_SIZE_IN> input[LINPUPPI_CHS_NCHANN_IN]);
void linpuppi_unpack_in(const ap_uint<LINPUPPI_DATA_SIZE_IN> input[LINPUPPI_CHS_NCHANN_IN], TkObj track[NTRACK], z0_t & pvZ0, PFNeutralObj pfallne[NALLNEUTRALS]);
void linpuppi_chs_pack_in(z0_t pvZ0, const PFChargedObj pfch[NTRACK], ap_uint<LINPUPPI_DATA_SIZE_IN> input[LINPUPPI_CHS_NCHANN_IN]);
void linpuppi_chs_unpack_in(const ap_uint<LINPUPPI_DATA_SIZE_IN> input[LINPUPPI_CHS_NCHANN_IN], z0_t & pvZ0, PFChargedObj pfch[NTRACK]);

typedef ap_uint<17+eta_t::width+phi_t::width> packed_linpuppi_refobj;
inline linpuppi_refobj linpuppi_refobj_unpack(const packed_linpuppi_refobj & src) {
    linpuppi_refobj ret;
    ret.pt2_shift(16,0)        = src(16,0);
    ret.hwEta(eta_t::width-1, 0) = src(17+eta_t::width-1,17);
    ret.hwPhi(phi_t::width-1, 0) = src(17+eta_t::width+phi_t::width-1,17+eta_t::width);
    return ret;
}
inline packed_linpuppi_refobj linpuppi_refobj_pack(const linpuppi_refobj & src) {
    packed_linpuppi_refobj ret;
    ret(16,0) = src.pt2_shift;
    ret(17+eta_t::width-1,17) = src.hwEta(eta_t::width-1, 0);
    ret(17+eta_t::width+phi_t::width-1,17+eta_t::width) = src.hwPhi(phi_t::width-1, 0);
    return ret;
}

packed_linpuppi_refobj packed_linpuppi_prepare_track(const ap_uint<TkObj::BITWIDTH> & track, const z0_t & pvZ0);
ap_uint<PuppiObj::BITWIDTH> packed_linpuppi_one(const ap_uint<PFNeutralObj::BITWIDTH> & in, const packed_linpuppi_refobj sel_tracks[NTRACK]);
ap_uint<PuppiObj::BITWIDTH> packed_linpuppi_chs_one(const ap_uint<PFChargedObj::BITWIDTH> & pfch, const z0_t & pvZ0) ;

// these two call the packed versions internally, and are used for valiation (they are not synthethized directly)
void packed_linpuppiNoCrop_streamed(const TkObj track[NTRACK], z0_t pvZ0, const PFNeutralObj pfallne[NALLNEUTRALS], PuppiObj outallne[NALLNEUTRALS]) ;
void packed_linpuppi_chs_streamed(z0_t pvZ0, const PFChargedObj pfch[NTRACK], PuppiObj outallch[NTRACK]) ;

void linpuppi_set_debug(bool debug);

#define LINPUPPI_ptLSB 0.25
#define LINPUPPI_DR2LSB 1.9e-5
#define LINPUPPI_dzLSB  0.05
#define LINPUPPI_pt2LSB LINPUPPI_ptLSB*LINPUPPI_ptLSB
#define LINPUPPI_pt2DR2_scale LINPUPPI_ptLSB*LINPUPPI_ptLSB/LINPUPPI_DR2LSB

#define LINPUPPI_sum_bitShift  15
#define LINPUPPI_x2_bits  6    // decimal bits the discriminator values
#define LINPUPPI_alpha_bits  5 // decimal bits of the alpha values
#define LINPUPPI_alphaSlope_bits  5 // decimal bits of the alphaSlope values
#define LINPUPPI_ptSlope_bits  6    // decimal bits of the ptSlope values 
#define LINPUPPI_weight_bits  8


//=================================================
#if defined(REG_Barrel)

#define LINPUPPI_DR2MAX  4727 // 0.3 cone
#define LINPUPPI_DR2MIN   257 // 0.07 cone
#define LINPUPPI_dzCut     10
#define LINPUPPI_ptMax    200 // 50.0/LINPUPPI_ptLSB 

#define LINPUPPI_ptSlopeNe  0.3
#define LINPUPPI_ptSlopePh  0.3
#define LINPUPPI_ptZeroNe   4.0
#define LINPUPPI_ptZeroPh   2.5
#define LINPUPPI_alphaSlope 0.7
#define LINPUPPI_alphaZero  6.0
#define LINPUPPI_alphaCrop  4.0
#define LINPUPPI_priorNe    5.0
#define LINPUPPI_priorPh    1.0

#define LINPUPPI_ptCut        4 // 1.0/LINPUPPI_ptLSB

//=================================================
#elif defined(REG_HGCal) 

#define LINPUPPI_etaBins 2
#define LINPUPPI_etaCut  0 // assuming the region spans [1.5,2.5] (or 1.25,2.75 with the overlaps), 
                           // the cut is exactly at the region center (2.0), so it's at 0 in integer coordinates
#define LINPUPPI_invertEta 0 // 0 if we're building the FW for the positive eta endcap.

#define LINPUPPI_DR2MAX  4727 // 0.3 cone
#define LINPUPPI_DR2MIN    84 // 0.04 cone
#define LINPUPPI_dzCut     40
#define LINPUPPI_ptMax    200 // 50.0/LINPUPPI_ptLSB 

#define LINPUPPI_ptSlopeNe  0.3 
#define LINPUPPI_ptSlopePh  0.4 
#define LINPUPPI_ptZeroNe   5.0 
#define LINPUPPI_ptZeroPh   3.0 
#define LINPUPPI_alphaSlope 1.5 
#define LINPUPPI_alphaZero  6.0 
#define LINPUPPI_alphaCrop  3.0 
#define LINPUPPI_priorNe    5.0 
#define LINPUPPI_priorPh    1.5 

#define LINPUPPI_ptSlopeNe_1  0.3 
#define LINPUPPI_ptSlopePh_1  0.4 
#define LINPUPPI_ptZeroNe_1   7.0 
#define LINPUPPI_ptZeroPh_1   4.0 
#define LINPUPPI_alphaSlope_1 1.5 
#define LINPUPPI_alphaZero_1  6.0 
#define LINPUPPI_alphaCrop_1  3.0 
#define LINPUPPI_priorNe_1    5.0 
#define LINPUPPI_priorPh_1    1.5 


#define LINPUPPI_ptCut        4 // 1.0/LINPUPPI_ptLSB
#define LINPUPPI_ptCut_1      8 // 2.0/LINPUPPI_ptLSB

//=================================================
#elif defined(REG_HGCalNoTK)

#define LINPUPPI_DR2MAX  4727 // 0.3 cone
#define LINPUPPI_DR2MIN    84 // 0.04 cone
#define LINPUPPI_dzCut     40 // unused
#define LINPUPPI_ptMax    200 // 50.0/LINPUPPI_ptLSB 

#define LINPUPPI_ptSlopeNe  0.3
#define LINPUPPI_ptSlopePh  0.4
#define LINPUPPI_ptZeroNe   9.0
#define LINPUPPI_ptZeroPh   5.0
#define LINPUPPI_alphaSlope 2.2
#define LINPUPPI_alphaZero  9.0
#define LINPUPPI_alphaCrop  4.0
#define LINPUPPI_priorNe    7.0
#define LINPUPPI_priorPh    5.0

#define LINPUPPI_ptCut       16 // 4.0/LINPUPPI_ptLSB

//=================================================
#elif defined(REG_HF)

#define LINPUPPI_DR2MAX  4727 // 0.3 cone
#define LINPUPPI_DR2MIN   525 // 0.1 cone
#define LINPUPPI_dzCut     40 // unused

#define LINPUPPI_ptMax    400 // 100.0/LINPUPPI_ptLSB 

#define LINPUPPI_ptSlopeNe  0.25
#define LINPUPPI_ptSlopePh  0.25
#define LINPUPPI_ptZeroNe   14.
#define LINPUPPI_ptZeroPh   14.
#define LINPUPPI_alphaSlope 0.6
#define LINPUPPI_alphaZero  9.0
#define LINPUPPI_alphaCrop  4.0
#define LINPUPPI_priorNe    6.0
#define LINPUPPI_priorPh    6.0

#define LINPUPPI_ptCut      40  // 10.0/LINPUPPI_ptLSB

#endif

#endif
