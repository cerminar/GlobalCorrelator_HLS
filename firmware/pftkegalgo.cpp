#include "pftkegalgo.h"
#include <cassert>

// #include "pfalgo_common.icc"




//----- BEGIN ---------------
// FIXME: get all these from common once merged

template<int N> struct ct_log2_ceil { enum { value = ct_log2_ceil<(N/2)+(N%2)>::value + 1 }; };
template<> struct ct_log2_ceil<2> { enum { value = 1 }; };
template<> struct ct_log2_ceil<1> { enum { value = 0 }; };

template<typename DR_T, typename INDEX_T> struct match {
    DR_T val; INDEX_T idx;
    match(DR_T dr, int i) : val(dr), idx(i) {}
    match() {}

    static match best(const match & m1, const match & m2) {
        #pragma HLS inline
        return (m1.val <= m2.val) ? m1 : m2;
    }
};

template<typename MATCH_T, int N> struct picker_of_closest {
    static MATCH_T pick(const MATCH_T in[N]) {
        #pragma HLS inline
        static constexpr int halfWidth = N / 2;
        static constexpr int reducedSize = halfWidth + N % 2;
        MATCH_T reduced[reducedSize];
        #pragma HLS ARRAY_PARTITION variable=reduced complete
        for (int i = 0; i < halfWidth; ++i) {
            #pragma HLS unroll
            reduced[i] = MATCH_T::best(in[2*i], in[2*i+1]);
        }
        if (halfWidth != reducedSize) {
            reduced[reducedSize - 1] = in[N - 1];
        }
        return picker_of_closest<MATCH_T,reducedSize>::pick(reduced);
    }
};




template<typename MATCH_T>
struct picker_of_closest<MATCH_T,2> {
    static MATCH_T pick(const MATCH_T vals[2]) {
        #pragma HLS inline
        return MATCH_T::best(vals[0],vals[1]);
    }
};

template<typename MATCH_T>
struct picker_of_closest<MATCH_T,1> {
    static MATCH_T pick(const MATCH_T vals[1]) {
        #pragma HLS inline
        return vals[0];
    }
};


template<int DR2MAX, int NCA, typename DR_T>
ap_uint<NCA>  pick_closest(const DR_T calo_track_drval[NCA]) {
    const DR_T eDR2MAX = DR2MAX;
    typedef match<DR_T,ap_uint<ct_log2_ceil<NCA>::value>> match_t;
    match_t candidates[NCA];
    #pragma HLS ARRAY_PARTITION variable=candidates complete
    for (int i = 0; i < NCA; ++i) {
        candidates[i] = match_t(calo_track_drval[i], i);
    }

    match_t best = picker_of_closest<match_t, NCA>::pick(candidates);
    ap_uint<NCA> calo_track_link_bit = 0;
    calo_track_link_bit[best.idx] = (best.val >= eDR2MAX) ? 0 : 1;
    return calo_track_link_bit;
}


template<typename T, int NIn, int NOut>
void ptsort_hwopt(const T in[NIn], T out[NOut]) {
    T tmp[NOut];
    #pragma HLS ARRAY_PARTITION variable=tmp complete

    for (int iout = 0; iout < NOut; ++iout) {
        #pragma HLS unroll
        tmp[iout].hwPt = 0;
    }

    for (int it = 0; it < NIn; ++it) {
        for (int iout = NOut-1; iout >= 0; --iout) {
            if (tmp[iout].hwPt <= in[it].hwPt) {
                if (iout == 0 || tmp[iout-1].hwPt > in[it].hwPt) {
                    tmp[iout] = in[it];
                } else {
                    tmp[iout] = tmp[iout-1];
                }
            }
        }

    }
    for (int iout = 0; iout < NOut; ++iout) {
        out[iout] = tmp[iout];
    }

}


//-------- END


ap_int<pt_t::width+1> ell_dpt_int_cap(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, pt_t pt1, pt_t pt2, ap_int<pt_t::width+1> max) {
#pragma HLS INLINE
   // FIXME: this should be configurable
    const ap_uint<10> cdeta = 16;
    const ap_uint<10> cm = 256;

    ap_int<etaphi_t::width+1> d_eta = (eta1-eta2);
    ap_int<etaphi_t::width+1> d_phi = (phi1-phi2);

    int ell = d_phi*d_phi + d_eta*d_eta*cdeta;
    // FIXME: should be ap_uint<pt_t::width> ?
    ap_int<pt_t::width+1> d_pt = abs(pt1 - pt2);
    return (ell <= int(cm)) ? d_pt : max;
}



template<int DPTMAX>
void calo2tk_ellipticdptvals(const EmCaloObj &em, const TkObj track[NTRACK], ap_int<pt_t::width+1> calo_track_dptval[NTRACK]) {
#pragma HLS INLINE
    // ap_int<pt_t::width+1>  ?
    // FIXME: uint<pt_t::width>
    // FIXME: some of this should be configurable
    const ap_int<pt_t::width+1> eDPTMAX = DPTMAX;
    const pt_t trkQualityPtMin_ = 40; // 10 GeV

    track_loop: for (int itk = 0; itk < NTRACK; ++itk) {
      if (track[itk].hwPt < trkQualityPtMin_ || em.hwPt == 0) {
        calo_track_dptval[itk] = eDPTMAX;
      } else {
        calo_track_dptval[itk] = ell_dpt_int_cap(em.hwEta, em.hwPhi, track[itk].hwEta, track[itk].hwPhi, em.hwPt, track[itk].hwPt, eDPTMAX);
        // std::cout << "[" << itk << "] dpt: " << calo_track_dptval[itk] << std::endl;
      }
    }

}



void link_emCalo2emCalo(const EmCaloObj emcalo[NEMCALOSEL_EGIN], ap_uint<NEMCALOSEL_EGIN> emCalo2emcalo_bit[NEMCALOSEL_EGIN]) {
  #pragma HLS ARRAY_PARTITION variable=emcalo complete dim=1
  #pragma HLS ARRAY_PARTITION variable=emCalo2emcalo_bit complete dim=1
  #pragma HLS INLINE

  // FIXME: configuration
  const ap_int<etaphi_t::width+1>  dEtaMaxBrem_ = 5; // 0.02; -> round(0.02*4*180/3.14)
  const ap_int<etaphi_t::width+1>  dPhiMaxBrem_ = 23; // 0.1; -> round(0.1*4*180/3.14)

  // NOTE: we assume the input to be sorted!!!
  brem_reco_outer_loop: for (int ic = 0; ic < NEMCALOSEL_EGIN; ++ic) {
    auto &calo = emcalo[ic];
    brem_reco_inner_loop: for (int jc = ic + 1; jc < NEMCALOSEL_EGIN; ++jc) {
      auto &otherCalo = emcalo[jc];
      if (calo.hwPt != 0 && otherCalo.hwPt != 0 &&
        abs(otherCalo.hwEta - calo.hwEta) < dEtaMaxBrem_ &&
          abs(otherCalo.hwPhi - calo.hwPhi) < dPhiMaxBrem_) {
            emCalo2emcalo_bit[ic][jc] = 1;
            emCalo2emcalo_bit[jc][jc] = 1; // use diagonal bit to mark the cluster as already used
      }
    }
  }
}





void link_emCalo2tk(const EmCaloObj emcalo[NEMCALOSEL_EGIN],
                    const TkObj track[NTRACK],
                    ap_uint<NTRACK> emCalo2tk_bit[NEMCALOSEL_EGIN]) {
  #pragma HLS INLINE

  #pragma HLS ARRAY_PARTITION variable=emcalo complete dim=1
  #pragma HLS ARRAY_PARTITION variable=track complete dim=1
  #pragma HLS ARRAY_PARTITION variable=emCalo2tk_bit complete dim=1

  // FIXME: parametrize
  const int DPTMAX = 65535; // DPT = 16+1 bits -> max = 2^16-1 = ((1<<16)-1)
  calo_loop: for (int ic = 0; ic < NEMCALOSEL_EGIN; ++ic) {
    ap_int<pt_t::width+1> dptvals[NTRACK];
    calo2tk_ellipticdptvals<DPTMAX>(emcalo[ic], track, dptvals);
    emCalo2tk_bit[ic] = pick_closest<DPTMAX,NTRACK,ap_int<pt_t::width+1>>(dptvals);
  }

}


void sel_emCalo(const EmCaloObj emcalo[NEMCALO], EmCaloObj emcalo_sel[NEMCALOSEL_EGIN]) {
  #pragma HLS INLINE

  EmCaloObj emcalo_selt[NEMCALO];
  #pragma HLS ARRAY_PARTITION variable=emcalo_sel complete dim=1
  #pragma HLS ARRAY_PARTITION variable=emcalo_selt complete dim=1
  EmCaloObj emcalo_zero;
  clear(emcalo_zero);
  in_select_loop: for(int ic = 0; ic < NEMCALO; ++ic) {
    emcalo_selt[ic] = (emcalo[ic].hwFlags == SELHWQUAL) ? emcalo[ic] : emcalo_zero;
  }
  ptsort_hwopt<EmCaloObj,NEMCALO,NEMCALOSEL_EGIN>(emcalo_selt, emcalo_sel);
}


void pftkegalgo(const EmCaloObj emcalo[NCALO], const TkObj track[NTRACK],
  // ap_uint<NTRACK> emCalo2tk_bit[NEMCALOSEL_EGIN], ap_uint<NEMCALOSEL_EGIN> emCalo2emcalo_bit[NEMCALOSEL_EGIN],
  EGIsoParticle photons[NEM_EGOUT], EGIsoEleParticle eles[NEM_EGOUT]) {
  #pragma HLS PIPELINE II=1
  #pragma HLS ARRAY_PARTITION variable=emcalo complete dim=1
  #pragma HLS ARRAY_PARTITION variable=track complete dim=1

  #pragma HLS ARRAY_PARTITION variable=photons complete dim=1
  #pragma HLS ARRAY_PARTITION variable=eles complete dim=1

  // FIXME: filter out qualities which are not useful for later processing (eg quality !=4 in HGC)
  // EmCaloObj emcalo_sel[NEMCALOSEL_EGIN];
  EmCaloObj emcalo_sel[NEMCALOSEL_EGIN];
  EmCaloObj emcalo_selt[NEMCALO];
  #pragma HLS ARRAY_PARTITION variable=emcalo_sel complete dim=1
  #pragma HLS ARRAY_PARTITION variable=emcalo_selt complete dim=1
  EmCaloObj emcalo_zero;
  clear(emcalo_zero);
  in_select_loop: for(int ic = 0; ic < NEMCALO; ++ic) {
    emcalo_selt[ic] = (emcalo[ic].hwFlags == SELHWQUAL) ? emcalo[ic] : emcalo_zero;
  }
  ptsort_hwopt<EmCaloObj,NEMCALO,NEMCALOSEL_EGIN>(emcalo_selt, emcalo_sel);

  // FIXME: shall we forseen the same selection step for tracks?
  ap_uint<NTRACK> emCalo2tk_bit[NEMCALOSEL_EGIN];
  ap_uint<NEMCALOSEL_EGIN> emCalo2emcalo_bit[NEMCALOSEL_EGIN];
  #pragma HLS ARRAY_PARTITION variable=emCalo2tk_bit complete dim=1
  #pragma HLS ARRAY_PARTITION variable=emCalo2emcalo_bit complete dim=1

  // initialize
  init_loop: for (int ic = 0; ic < NEMCALOSEL_EGIN; ++ic) {
    emCalo2tk_bit[ic] = 0;
    emCalo2emcalo_bit[ic] = 0;
  }

  #if defined(DOBREMRECOVERY)
    link_emCalo2emCalo(emcalo_sel, emCalo2emcalo_bit);
  #endif
  //
  link_emCalo2tk(emcalo_sel, track, emCalo2tk_bit);


  EGIsoParticle photons_temp[NEMCALOSEL_EGIN];
  EGIsoEleParticle eles_temp[NEMCALOSEL_EGIN];
  #pragma HLS ARRAY_PARTITION variable=photons_temp complete dim=1
  #pragma HLS ARRAY_PARTITION variable=eles_temp complete dim=1

  loop_calo: for (int ic = 0; ic < NEMCALOSEL_EGIN; ++ic) {
    loop_track_matched: for(int it = 0; it < NTRACK; ++it) {
      if(emCalo2tk_bit[ic][it]) break;
    }

    pt_t ptcorr = emcalo_sel[ic].hwPt;
    #if defined(DOBREMRECOVERY)
      if(emCalo2emcalo_bit[ic][ic] != 1) { // FIXME: we still save the multiclusters used as brem???
        // FIXME: we should set the quality bit as "brem-recovery performed"
        loop_calo_brem_reco: for (int ioc = 0; ioc < NEMCALOSEL_EGIN; ++ioc) {
          if(emCalo2emcalo_bit[ic][ioc]) {
            ptcorr += emcalo_sel[ioc].hwPt;
          }
        }
      }
    #endif

    photons_temp[ic].hwPt = ptcorr;
    photons_temp[ic].hwEta = emcalo_sel[ic].hwEta;
    photons_temp[ic].hwPhi = emcalo_sel[ic].hwPhi;

    if(emCalo2tk_bit[ic]) {
      eles_temp[ic].hwPt = ptcorr;
      eles_temp[ic].hwEta = emcalo_sel[ic].hwEta;
      eles_temp[ic].hwPhi = emcalo_sel[ic].hwPhi;
      // FIXME: add track properties @ vertex using
      //track[it]
    } else {
      clear(eles_temp[ic]);
    }

  }
  ptsort_hwopt<EGIsoParticle,NEMCALOSEL_EGIN,NEM_EGOUT>(photons_temp, photons);
  ptsort_hwopt<EGIsoEleParticle,NEMCALOSEL_EGIN,NEM_EGOUT>(eles_temp, eles);

}
