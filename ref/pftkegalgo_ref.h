#ifndef PFTKEGALGO_REF_H
#define PFTKEGALGO_REF_H
#ifndef CMSSW_GIT_HASH
  #include "../DiscretePFInputs.h"
#else
  #include "../../interface/DiscretePFInputs.h"
#endif

#include "../firmware/pfalgo2hgc.h"
#include "pfalgo_common_ref.h"



namespace PFTkEGAlgo {

  struct pftkegalgo_config {
      unsigned int nTRACK;
      unsigned int nEMCALO;
      unsigned int nEMCALOSEL_EGIN;
      unsigned int nEM_EGOUT;

      pftkegalgo_config(unsigned int nTrack, unsigned int nEmCalo, unsigned int nEmCaloSel_in, unsigned int nEmOut) :
          nTRACK(nTrack), nEMCALO(nEmCalo), nEMCALOSEL_EGIN(nEmCaloSel_in), nEM_EGOUT(nEmOut) {}
  };

  struct Region {
    std::vector<l1tpf_impl::CaloCluster> emcalo;
    std::vector<l1tpf_impl::PropagatedTrack> track;
    std::vector<EGIsoParticle> egphotons;
    std::vector<EGIsoEleParticle> egeles;
  };

  void pftkegalgo_ref_set_debug(int debug) ;

  void pftkegalgo_ref(const pftkegalgo_config &cfg,
                      const EmCaloObj emcalo[/*cfg.nCALO*/],
                      const TkObj track[/*cfg.nTRACK*/],
                      int emCalo2tk_ref[],
                      int emCalo2emcalo_ref[],
                      EGIsoParticle egphs_ref[],
                      EGIsoEleParticle egele_ref[]) ;

  void link_emCalo2emCalo(const std::vector<l1tpf_impl::CaloCluster>& emcalo,
                          std::vector<int> &emCalo2emCalo);

  void link_emCalo2tk(const std::vector<l1tpf_impl::CaloCluster> &emcalo,
                      const std::vector<l1tpf_impl::PropagatedTrack> &track,
                      std::vector<int> &emCalo2tk);

  float deltaPhi(float phi1, float phi2);

  bool hwpt_sort(const l1tpf_impl::CaloCluster &i, const l1tpf_impl::CaloCluster &j);

  void sel_emCalo_ref(const EmCaloObj emcalo[/*cfg.nCALO*/],
                      EmCaloObj emcalo_sel_ref[/*cfg.nCALO*/]);

  void eg_algo(Region &r,
               const std::vector<int> &emCalo2emCalo,
               const std::vector<int> &emCalo2tk);

   EGIsoParticle &addEGIsoToPF(std::vector<EGIsoParticle> &egobjs,
                               const l1tpf_impl::CaloCluster &calo,
                               const int hwQual,
                               const int ptCorr);

   EGIsoEleParticle &addEGIsoEleToPF(std::vector<EGIsoEleParticle> &egobjs,
                                     const l1tpf_impl::CaloCluster &calo,
                                     const l1tpf_impl::PropagatedTrack &track,
                                     const int hwQual,
                                     const int ptCorr);

   void addEgObjsToPF(Region &r, const int calo_idx, const int hwQual, const int ptCorr, const int tk_idx);
}

#endif
