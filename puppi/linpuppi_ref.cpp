#include "linpuppi_ref.h"
#include "firmware/linpuppi_bits.h"
#include <cmath>
#include <algorithm>

namespace l1ct {
    unsigned int linpuppi_ieta_ref(const linpuppi_config &cfg, eta_t eta) ;
    std::pair<pt_t,puppiWgt_t> linpuppi_ref_sum2puppiPt(const linpuppi_config &cfg, uint64_t sum, pt_t pt, unsigned int ieta, bool isEM, int icand, bool debug) ;
    std::pair<float,float> linpuppi_flt_sum2puppiPt(const linpuppi_config &cfg, float sum, float pt, unsigned int ieta, bool isEM, int icand, bool debug) ;
};

using namespace l1ct;

l1ct::linpuppi_config::linpuppi_config(unsigned int nTrack_, unsigned int nIn_, unsigned int nOut_,
                    unsigned int dR2Min_, unsigned int dR2Max_, unsigned int ptMax_, unsigned int dzCut_,
                    int etaCut, bool invertEta,
                    float ptSlopeNe_0, float ptSlopeNe_1, float ptSlopePh_0, float ptSlopePh_1, float ptZeroNe_0, float ptZeroNe_1, float ptZeroPh_0, float ptZeroPh_1, 
                    float alphaSlope_0, float alphaSlope_1, float alphaZero_0, float alphaZero_1, float alphaCrop_0, float alphaCrop_1, 
                    float priorNe_0, float priorNe_1, float priorPh_0, float priorPh_1, 
                    unsigned int ptCut_0, unsigned int ptCut_1) :
                nTrack(nTrack_), nIn(nIn_), nOut(nOut_),
                dR2Min(dR2Min_), dR2Max(dR2Max_), ptMax(ptMax_), dzCut(dzCut_),
                absEtaBins(1, etaCut), invertEtaBins(invertEta),
                ptSlopeNe(2), ptSlopePh(2), ptZeroNe(2), ptZeroPh(2), alphaSlope(2), alphaZero(2), alphaCrop(2), priorNe(2), priorPh(2), 
                ptCut(2) 
{
    ptSlopeNe[0] = ptSlopeNe_0; 
    ptSlopeNe[1] = ptSlopeNe_1;
    ptSlopePh[0] = ptSlopePh_0; 
    ptSlopePh[1] = ptSlopePh_1;
    ptZeroNe[0] = ptZeroNe_0; 
    ptZeroNe[1] = ptZeroNe_1;
    ptZeroPh[0] = ptZeroPh_0; 
    ptZeroPh[1] = ptZeroPh_1;
    alphaSlope[0] = alphaSlope_0; 
    alphaSlope[1] = alphaSlope_1;
    alphaZero[0] = alphaZero_0;
    alphaZero[1] = alphaZero_1;
    alphaCrop[0] = alphaCrop_0; 
    alphaCrop[1] = alphaCrop_1;
    priorNe[0] = priorNe_0; 
    priorNe[1] = priorNe_1;
    priorPh[0] = priorPh_0; 
    priorPh[1] = priorPh_1;
    ptCut[0] = ptCut_0; 
    ptCut[1] = ptCut_1;
}


void l1ct::puppisort_and_crop_ref(unsigned int nIn, unsigned int nOut, const PuppiObj in[/*nIn*/], PuppiObj out[/*nOut*/]) {
    std::vector<PuppiObj> tmp(nOut);

    for (unsigned int iout = 0; iout < nOut; ++iout) {
        tmp[iout].clear();
    }

    for (unsigned int it = 0; it < nIn; ++it) {
        for (int iout = int(nOut)-1; iout >= 0; --iout) {
            if (tmp[iout].hwPt <= in[it].hwPt) {
                if (iout == 0 || tmp[iout-1].hwPt > in[it].hwPt) {
                    tmp[iout] = in[it];
                } else {
                    tmp[iout] = tmp[iout-1];
                }
            }
        }
    }

    for (unsigned int iout = 0; iout < nOut; ++iout) {
        out[iout] = tmp[iout];
    }

}


void l1ct::linpuppi_chs_ref(const linpuppi_config &cfg, z0_t pvZ0, const PFChargedObj pfch[/*cfg.nTrack*/], PuppiObj outallch[/*cfg.nTrack*/], bool debug) {
    for (unsigned int i = 0; i < cfg.nTrack; ++i) {
        int z0diff = pfch[i].hwZ0 - pvZ0;
        if (std::abs(z0diff) <= cfg.dzCut || pfch[i].hwId.isMuon()) {
            outallch[i].fill(pfch[i]);
            if (debug && pfch[i].hwPt > 0) printf("ref candidate %02u pt %7.2f pid %1d   vz %+6d  dz %+6d (cut %5d) -> pass\n", i, 
                                    pfch[i].floatPt(), pfch[i].intId(), int(pfch[i].hwZ0), z0diff, cfg.dzCut);
        } else {
            outallch[i].clear();
            if (debug && pfch[i].hwPt > 0) printf("ref candidate %02u pt %7.2f pid %1d   vz %+6d  dz %+6d (cut %5d) -> fail\n", i, 
                                    pfch[i].floatPt(), pfch[i].intId(), int(pfch[i].hwZ0), z0diff, cfg.dzCut);
        }
    }
}

unsigned int l1ct::linpuppi_ieta_ref(const linpuppi_config &cfg, eta_t eta) {
    int n = cfg.absEtaBins.size();
    for (int i = 0; i < n; ++i) {
        if (int(eta) <= cfg.absEtaBins[i]) return (cfg.invertEtaBins ? n-i : i);
    }
    return cfg.invertEtaBins ? 0 : n;
}

std::pair<pt_t,puppiWgt_t> l1ct::linpuppi_ref_sum2puppiPt(const linpuppi_config &cfg, uint64_t sum, pt_t pt, unsigned int ieta, bool isEM, int icand, bool debug) {
    const int sum_bitShift = LINPUPPI_sum_bitShift;
    const int x2_bits = LINPUPPI_x2_bits;    // decimal bits the discriminator values
    const int alpha_bits = LINPUPPI_alpha_bits; // decimal bits of the alpha values
    const int alphaSlope_bits = LINPUPPI_alphaSlope_bits; // decimal bits of the alphaSlope values
    const int ptSlope_bits = LINPUPPI_ptSlope_bits;    // decimal bits of the ptSlope values 
    const int weight_bits = LINPUPPI_weight_bits;

    const int ptSlopeNe = cfg.ptSlopeNe[ieta] * (1 << ptSlope_bits);
    const int ptSlopePh = cfg.ptSlopePh[ieta] * (1 << ptSlope_bits);
    const int ptZeroNe = cfg.ptZeroNe[ieta] / LINPUPPI_ptLSB; // in pt scale
    const int ptZeroPh = cfg.ptZeroPh[ieta] / LINPUPPI_ptLSB; // in pt scale
    const int alphaCrop = cfg.alphaCrop[ieta] * (1 << x2_bits);
    const int alphaSlopeNe = cfg.alphaSlope[ieta] * std::log(2.) * (1 << alphaSlope_bits); // we put a log(2) here since we compute alpha as log2(sum) instead of ln(sum)
    const int alphaSlopePh = cfg.alphaSlope[ieta] * std::log(2.) * (1 << alphaSlope_bits);
    const int alphaZeroNe = cfg.alphaZero[ieta] / std::log(2.) * (1 << alpha_bits);
    const int alphaZeroPh = cfg.alphaZero[ieta] / std::log(2.) * (1 << alpha_bits);
    const int priorNe = cfg.priorNe[ieta] * (1 << x2_bits);
    const int priorPh = cfg.priorPh[ieta] * (1 << x2_bits);

    // -- simplest version
    //int alpha = sum > 0 ? int(std::log2(float(sum) * LINPUPPI_pt2DR2_scale / (1<<sum_bitShift)) * (1 << alpha_bits) + 0.5) :  0;
    // -- re-written bringing terms out of the log
    //int alpha = sum > 0 ? int(std::log2(float(sum))*(1 << alpha_bits) + (std::log2(LINPUPPI_pt2DR2_scale) - sum_bitShift)*(1 << alpha_bits) + 0.5 ) :  0;
    // -- re-written for a LUT implementation of the log2
    const int log2lut_bits = 10;
    int alpha = 0; uint64_t logarg = sum;
    if (logarg > 0) {
        alpha = int((std::log2(LINPUPPI_pt2DR2_scale) - sum_bitShift)*(1 << alpha_bits) + 0.5);
        while (logarg >= (1 << log2lut_bits)) { logarg = logarg >> 1; alpha += (1 << alpha_bits); }
        alpha += int(std::log2(float(logarg))*(1 << alpha_bits)); // the maximum value of this term is log2lut_bits * (1 << alpha_bits) ~ 10*16 = 160 => fits in ap_uint<4+alpha_bits>
    }
    int alphaZero  = (isEM ? alphaZeroPh : alphaZeroNe);
    int alphaSlope = (isEM ? alphaSlopePh : alphaSlopeNe);
    int x2a = std::min(std::max( alphaSlope * (alpha - alphaZero)  >> (alphaSlope_bits + alpha_bits - x2_bits), -alphaCrop), alphaCrop);

    // -- re-written to fit in a single LUT
    int x2a_lut = - alphaSlope * alphaZero; logarg = sum;
    if (logarg > 0) {
        x2a_lut += alphaSlope * int((std::log2(LINPUPPI_pt2DR2_scale) - sum_bitShift)*(1 << alpha_bits) + 0.5);
        while (logarg >= (1 << log2lut_bits)) { 
            logarg = logarg >> 1; x2a_lut += alphaSlope * (1 << alpha_bits); 
        }
        x2a_lut += alphaSlope * int(std::log2(float(logarg))*(1 << alpha_bits)); 
        /*if (in <= 3) printf("ref [%d]:  x2a(sum = %9lu): logarg = %9lu, sumterm = %9d, table[logarg] = %9d, ret pre-crop = %9d\n", 
          in, sum, logarg, 
          alphaSlope * int((std::log2(LINPUPPI_pt2DR2_scale) - sum_bitShift)*(1 << alpha_bits) + 0.5) - alphaSlope * alphaZero,
          alphaSlope * int(std::log2(float(logarg))*(1 << alpha_bits)), 
          x2a_lut); */
    } else {
        //if (in <= 3) printf("ref [%d]:  x2a(sum = %9lu): logarg = %9lu, ret pre-crop = %9d\n", 
        //        in, sum, logarg, x2a_lut); 
    }
    x2a_lut = std::min(std::max( x2a_lut >> (alphaSlope_bits + alpha_bits - x2_bits), -alphaCrop), alphaCrop);
    assert( x2a_lut == x2a );

    int ptZero  = (isEM ? ptZeroPh : ptZeroNe);
    int ptSlope = (isEM ? ptSlopePh : ptSlopeNe);
    int x2pt    = ptSlope * (Scales::ptToInt(pt) - ptZero) >> (ptSlope_bits + 2 - x2_bits);

    int prior  = (isEM ? priorPh : priorNe);

    int x2 = x2a + x2pt - prior;

    int weight = std::min<int>( 1.0/(1.0 + std::exp(- float(x2)/(1<<x2_bits))) * ( 1 << weight_bits ) + 0.5, (1 << weight_bits) );

    pt_t ptPuppi = Scales::makePt(( Scales::ptToInt(pt) * weight ) >> weight_bits);

    if (debug) printf("ref candidate %02d pt %7.2f  em %1d: alpha %+7.2f   x2a %+5d = %+7.3f  x2pt %+5d = %+7.3f   x2 %+5d = %+7.3f  --> weight %4d = %.4f  puppi pt %7.2f\n",
               icand, Scales::floatPt(pt), int(isEM), 
               std::max<float>(alpha/float(1<<alpha_bits)*std::log(2.),-99.99f), 
               x2a, x2a/float(1<<x2_bits), x2pt, x2pt/float(1<<x2_bits), x2, x2/float(1<<x2_bits), 
               weight, weight/float( 1 << weight_bits ), 
               Scales::floatPt(ptPuppi));

    return std::make_pair(ptPuppi, puppiWgt_t(weight));
}


void l1ct::fwdlinpuppi_ref(const linpuppi_config &cfg, const HadCaloObj caloin[/*cfg.nIn*/], PuppiObj outallne_nocut[/*cfg.nIn*/], PuppiObj outallne[/*cfg.nIn*/], PuppiObj outselne[/*cfg.nOut*/], bool debug) {
    const int DR2MAX = cfg.dR2Max; 
    const int DR2MIN = cfg.dR2Min; 
    const int PTMAX2 = (cfg.ptMax*cfg.ptMax);

    const int sum_bitShift = LINPUPPI_sum_bitShift;

    for (unsigned int in = 0; in < cfg.nIn; ++in) {
        outallne_nocut[in].clear(); outallne[in].clear();
        if (caloin[in].hwPt == 0) continue;
        uint64_t sum = 0; // 2 ^ sum_bitShift times (int pt^2)/(int dr2)
        for (unsigned int it = 0; it < cfg.nIn; ++it) {
            if (it == in || caloin[it].hwPt == 0) continue;
            int dr2 = dr2_int(caloin[it].hwEta, caloin[it].hwPhi, caloin[in].hwEta, caloin[in].hwPhi); // if dr is inside puppi cone
            if (dr2 <= DR2MAX) {
                ap_uint<9> dr2short = (dr2 >= DR2MIN ? dr2 : DR2MIN) >> 5; // reduce precision to make divide LUT cheaper
                uint64_t pt2 = Scales::ptToInt(caloin[it].hwPt)*Scales::ptToInt(caloin[it].hwPt);
                uint64_t term = std::min<uint64_t>(pt2 >> 5, PTMAX2 >> 5)*((1 << sum_bitShift)/int(dr2short));
                //      dr2short >= (DR2MIN >> 5) = 2
                //      num <= (PTMAX2 >> 5) << sum_bitShift = (2^11) << 15 = 2^26
                //      ==> term <= 2^25
                //printf("ref term [%2d,%2d]: dr = %8d  pt2_shift = %8lu  term = %12lu\n", in, it, dr2, std::min<uint64_t>(pt2 >> 5, PTMAX2 >> 5), term);
                assert(uint64_t(PTMAX2 << (sum_bitShift-5))/(DR2MIN >> 5) <= (1 << 25));
                assert(term < (1 << 25));
                sum += term;
                //printf("    pT cand %5.1f    pT item %5.1f    dR = %.3f   term = %.1f [dbl] = %lu [int]\n",
                //            caloin[in].floatPt(), caloin[it].floatPt(), std::sqrt(dr2*LINPUPPI_DR2LSB),
                //            double(std::min<uint64_t>(pt2 >> 5, 131071)<<15)/double(std::max<int>(dr2,DR2MIN) >> 5),
                //            term);
            }
        }
        unsigned int ieta = linpuppi_ieta_ref(cfg, caloin[in].hwEta);
        std::pair<pt_t,puppiWgt_t> ptAndW = linpuppi_ref_sum2puppiPt(cfg, sum, caloin[in].hwPt, ieta, caloin[in].hwIsEM, in, debug);

        outallne_nocut[in].fill(caloin[in], ptAndW.first, ptAndW.second);
        if (outallne_nocut[in].hwPt >= cfg.ptCut[ieta]) {
            outallne[in] = outallne_nocut[in];
        }
    }
    puppisort_and_crop_ref(cfg.nIn, cfg.nOut, outallne, outselne);
}

void l1ct::linpuppi_ref(const linpuppi_config &cfg, const TkObj track[/*cfg.nTrack*/], z0_t pvZ0, const PFNeutralObj pfallne[/*cfg.nIn*/], PuppiObj outallne_nocut[/*cfg.nIn*/], PuppiObj outallne[/*cfg.nIn*/], PuppiObj outselne[/*cfg.nOut*/], bool debug) {
    const int DR2MAX = cfg.dR2Max; 
    const int DR2MIN = cfg.dR2Min; 
    const int PTMAX2 = (cfg.ptMax*cfg.ptMax);

    const int sum_bitShift = LINPUPPI_sum_bitShift;

    for (unsigned int in = 0; in < cfg.nIn; ++in) {
        outallne_nocut[in].clear(); outallne[in].clear();
        if (pfallne[in].hwPt == 0) continue;
        uint64_t sum = 0; // 2 ^ sum_bitShift times (int pt^2)/(int dr2)
        for (unsigned int it = 0; it < cfg.nTrack; ++it) {
            if (track[it].hwPt == 0) continue;
            if (std::abs(int(track[it].hwZ0 - pvZ0)) > cfg.dzCut) continue;
            int dr2 = dr2_int(pfallne[in].hwEta, pfallne[in].hwPhi, track[it].hwEta, track[it].hwPhi); // if dr is inside puppi cone
            if (dr2 <= DR2MAX) {
                ap_uint<9> dr2short = (dr2 >= DR2MIN ? dr2 : DR2MIN) >> 5; // reduce precision to make divide LUT cheaper
                uint64_t pt2 = Scales::ptToInt(track[it].hwPt)*Scales::ptToInt(track[it].hwPt);
                uint64_t term = std::min<uint64_t>(pt2 >> 5, PTMAX2 >> 5)*((1 << sum_bitShift)/int(dr2short));
                //      dr2short >= (DR2MIN >> 5) = 2
                //      num <= (PTMAX2 >> 5) << sum_bitShift = (2^11) << 15 = 2^26
                //      ==> term <= 2^25
                //printf("ref term [%2d,%2d]: dr = %8d  pt2_shift = %8lu  term = %12lu\n", in, it, dr2, std::min<uint64_t>(pt2 >> 5, PTMAX2 >> 5), term);
                assert(uint64_t(PTMAX2 << (sum_bitShift-5))/(DR2MIN >> 5) <= (1 << 25));
                assert(term < (1 << 25));
                sum += term;
                //printf("    pT cand %5.1f    pT item %5.1f    dR = %.3f   term = %.1f [dbl] = %lu [int]\n",
                //            pfallne[in].floatPt(), track[it].floatPt(), std::sqrt(dr2*LINPUPPI_DR2LSB),
                //            double(std::min<uint64_t>(pt2 >> 5, 131071)<<15)/double(std::max<int>(dr2,DR2MIN) >> 5),
                //            term);
            }
        }

        unsigned int ieta = linpuppi_ieta_ref(cfg, pfallne[in].hwEta);
        bool isEM = (pfallne[in].hwId.isPhoton());
        std::pair<pt_t,puppiWgt_t> ptAndW = linpuppi_ref_sum2puppiPt(cfg, sum, pfallne[in].hwPt, ieta, isEM, in, debug);
        outallne_nocut[in].fill(pfallne[in], ptAndW.first, ptAndW.second);
        if (outallne_nocut[in].hwPt >= cfg.ptCut[ieta]) {
            outallne[in] = outallne_nocut[in];
        }
     }
    puppisort_and_crop_ref(cfg.nIn, cfg.nOut, outallne, outselne);

}

std::pair<float,float> l1ct::linpuppi_flt_sum2puppiPt(const linpuppi_config &cfg, float sum, float pt, unsigned int ieta, bool isEM, int icand, bool debug) {
    float alphaZero  = cfg.alphaZero[ieta], alphaSlope = cfg.alphaSlope[ieta], alphaCrop = cfg.alphaCrop[ieta];
    float alpha = sum > 0 ? std::log(sum) : -9e9;
    float x2a = std::min(std::max( alphaSlope * (alpha - alphaZero), -alphaCrop), alphaCrop);

    float ptZero  = (isEM ? cfg.ptZeroPh[ieta]  : cfg.ptZeroNe[ieta]);
    float ptSlope = (isEM ? cfg.ptSlopePh[ieta] : cfg.ptSlopeNe[ieta]);
    float x2pt    = ptSlope * (pt - ptZero);

    float prior  = (isEM ? cfg.priorPh[ieta] : cfg.priorNe[ieta]);

    float x2 = x2a + x2pt - prior;

    float weight = 1.0/(1.0 + std::exp(-x2));

    float puppiPt = pt *weight;
    if (debug) printf("flt candidate %02d pt %7.2f  em %1d: alpha %+7.2f   x2a         %+7.3f  x2pt         %+7.3f   x2         %+7.3f  --> weight        %.4f  puppi pt %7.2f\n",
                   icand, pt, int(isEM), std::max(alpha,-99.99f), x2a, x2pt, x2, weight, puppiPt);

    return std::make_pair(puppiPt,weight);
}


void l1ct::fwdlinpuppi_flt(const linpuppi_config &cfg, const HadCaloObj caloin[/*cfg.nIn*/], PuppiObj outallne_nocut[/*cfg.nIn*/], PuppiObj outallne[/*cfg.nIn*/], PuppiObj outselne[/*cfg.nOut*/], bool debug) {
    const int DR2MAX = cfg.dR2Max; 
    const int DR2MIN = cfg.dR2Min; 
    const float f_ptMax = Scales::floatPt(Scales::makePt(cfg.ptMax));

    for (unsigned int in = 0; in < cfg.nIn; ++in) {
        outallne_nocut[in].clear(); outallne[in].clear();
        if (caloin[in].hwPt == 0) continue;
        float sum = 0;
        for (unsigned int it = 0; it < cfg.nIn; ++it) {
            if (it == in || caloin[it].hwPt == 0) continue;
            int dr2 = dr2_int(caloin[it].hwEta, caloin[it].hwPhi, caloin[in].hwEta, caloin[in].hwPhi); // if dr is inside puppi cone
            if (dr2 <= DR2MAX) {
                sum += std::pow(std::min<float>(caloin[it].floatPt(),f_ptMax),2) / (std::max<int>(dr2,DR2MIN) * LINPUPPI_DR2LSB);
            }
        }

        unsigned int ieta = linpuppi_ieta_ref(cfg, caloin[in].hwEta);
        std::pair<float,float> ptAndW = linpuppi_flt_sum2puppiPt(cfg, sum, caloin[in].floatPt(), ieta, caloin[in].hwIsEM, in, debug);
        outallne_nocut[in].fill(caloin[in], Scales::makePtFromFloat(ptAndW.first), int(ptAndW.second*256) );
        if (outallne_nocut[in].hwPt >= cfg.ptCut[ieta]) {
            outallne[in] = outallne_nocut[in];
        }
    }

    puppisort_and_crop_ref(cfg.nIn, cfg.nOut, outallne, outselne);
}

void l1ct::linpuppi_flt(const linpuppi_config &cfg, const TkObj track[/*cfg.nTrack*/], z0_t pvZ0, const PFNeutralObj pfallne[/*cfg.nIn*/], PuppiObj outallne_nocut[/*cfg.nIn*/], PuppiObj outallne[/*cfg.nIn*/], PuppiObj outselne[/*cfg.nOut*/], bool debug) {
    const int DR2MAX = cfg.dR2Max; 
    const int DR2MIN = cfg.dR2Min; 
    const float f_ptMax = Scales::floatPt(Scales::makePt(cfg.ptMax));

    for (unsigned int in = 0; in < cfg.nIn; ++in) {
        outallne_nocut[in].clear(); outallne[in].clear();
        if (pfallne[in].hwPt == 0) continue;
        float sum = 0;
        for (unsigned int it = 0; it < cfg.nTrack; ++it) {
            if (track[it].hwPt == 0) continue;
            if (std::abs(int(track[it].hwZ0 - pvZ0)) > cfg.dzCut) continue;
            int dr2 = dr2_int(pfallne[in].hwEta, pfallne[in].hwPhi, track[it].hwEta, track[it].hwPhi); // if dr is inside puppi cone
            if (dr2 <= DR2MAX) {
                sum += std::pow(std::min<float>(track[it].floatPt(),f_ptMax),2) / (std::max<int>(dr2,DR2MIN) * LINPUPPI_DR2LSB);
            }
        }
        unsigned int ieta = linpuppi_ieta_ref(cfg, pfallne[in].hwEta);
        bool isEM = pfallne[in].hwId.isPhoton();
        std::pair<float,float> ptAndW = linpuppi_flt_sum2puppiPt(cfg, sum, pfallne[in].floatPt(), ieta, isEM, in, debug);
        outallne_nocut[in].fill(pfallne[in], Scales::makePtFromFloat(ptAndW.first), int(ptAndW.second*256) );
        if (outallne_nocut[in].hwPt >= cfg.ptCut[ieta]) {
            outallne[in] = outallne_nocut[in];
        }
    }
    puppisort_and_crop_ref(cfg.nIn, cfg.nOut, outallne, outselne);
}



