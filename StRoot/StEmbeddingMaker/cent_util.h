#ifndef YGHUANG_CENT_UTIL_HEADER
#define YGHUANG_CENT_UTIL_HEADER

/*
    Version: 1.1
*/

#include <iostream>
#include "TF1.h"

using std::cout;

class Corr{
    private:

        // pile up rejection parameters
        Double_t k11, k12, b11, b12;
        Double_t k21, k22, b21, b22;
        Double_t k31, k32, b31, b32;
        bool do_pile_up;

        static const int nTrg = 3; // this depends on your dataset
        // static const int nTrg = 1;
        Double_t k_lumi[nTrg];
        Double_t b_lumi[nTrg];
        bool do_lumi;

        // vz correction parameter
        Double_t par[nTrg][7];
        TF1* hfunc;
        bool do_vz;

        // centrality definition with centrality bin edge
        Double_t edge[9];
        bool do_split;

    public:
        Corr(){
            k11 = 0;
            k12 = 0;
            b11 = 0;
            b12 = 0;
            k21 = 0;
            k22 = 0;
            b21 = 0;
            b22 = 0;
            k31 = 0;
            k32 = 0;
            b31 = 0;
            b32 = 0;
            for (int i=0; i<nTrg; i++){
                k_lumi[i] = 0;
                b_lumi[i] = 0;
            }
            hfunc = new TF1("hfunc", "pol6", -70, 70);
            for (int i=0; i<9; i++){
                edge[i] = 0;
            }
            do_pile_up = true;
            do_lumi = true;
            do_vz = true;
            do_split = false;
        }
        ~Corr(){}

        Int_t convert_trg(Int_t trg){ // this depends on your dataset
            if (trg == 810010){
                return 0;
            } else if (trg == 810020){
                return 1;
            } else if (trg == 810030){
                return 2;
            } else if (trg == 810040){
                return 3;
            } else {
                return -1;
            }
        }

        void set_do_pile_up(bool do_){
            do_pile_up = do_;
            if (do_){
                cout << "[LOG] Pile up rejection: ON.\n";
            } else {
                cout << "[LOG] Pile up rejection: OFF.\n";
            }
        }

        void set_do_lumi(bool do_){
            do_lumi = do_;
            if (do_){
                cout << "[LOG] Luminosity correction: ON.\n";
            } else {
                cout << "[LOG] Luminosity correction: OFF.\n";
            }
        }

        void set_do_vz(bool do_){
            do_vz = do_;
            if (do_){
                cout << "[LOG] Vz correction: ON.\n";
            } else {
                cout << "[LOG] Vz correction: OFF.\n";
            }
        }

        void set_pile_up_BetaEta1_lower_par(Double_t k11, Double_t b11){
            this->k11 = k11;
            this->b11 = b11;
        }

        void set_pile_up_BetaEta1_upper_par(Double_t k12, Double_t b12){
            this->k12 = k12;
            this->b12 = b12;
        }

        void set_pile_up_nTofMatch_lower_par(Double_t k21, Double_t b21){
            this->k21 = k21;
            this->b21 = b21;
        }

        void set_pile_up_nTofMatch_upper_par(Double_t k22, Double_t b22){
            this->k22 = k22;
            this->b22 = b22;
        }

        void set_pile_up_TofMult_lower_par(Double_t k31, Double_t b31){
            this->k31 = k31;
            this->b31 = b31;
        }

        void set_pile_up_TofMult_upper_par(Double_t k32, Double_t b32){
            this->k32 = k32;
            this->b32 = b32;
        }

        void set_pile_up_par(
            Double_t k11, Double_t b11, Double_t k12, Double_t b12,
            Double_t k21, Double_t b21, Double_t k22, Double_t b22,
            Double_t k31, Double_t b31, Double_t k32, Double_t b32
        ){
            set_pile_up_BetaEta1_lower_par(k11, b11);
            set_pile_up_BetaEta1_upper_par(k12, b12);
            set_pile_up_nTofMatch_lower_par(k21, b21);
            set_pile_up_nTofMatch_upper_par(k22, b22);
            set_pile_up_TofMult_lower_par(k31, b31);
            set_pile_up_TofMult_upper_par(k32, b32);
        }

        bool is_BetaEta1_bad(Int_t ref3, Int_t beta_eta1){
            return (k11 != -999 && ref3 * k11 + b11 > beta_eta1) || (k12 != -999 && ref3 * k12 + b12 < beta_eta1);
        }

        bool is_nTofMatch_bad(Int_t ref3, Int_t nTofMatch){
            return (k21 != -999 && ref3 * k21 + b21 > nTofMatch) || (k22 != -999 && ref3 * k22 + b22 < nTofMatch);
        }

        bool is_TofMult_bad(Int_t ref3, Int_t tofMult){
            return (k31 != -999 && ref3 * k31 + b31 > tofMult) || (k32 != -999 && ref3 * k32 + b32 < tofMult);
        }

        bool is_pile_up(Int_t ref3, Int_t beta_eta1, Int_t nTofMatch, Int_t tofMult){
            return is_BetaEta1_bad(ref3, beta_eta1) || is_nTofMatch_bad(ref3, nTofMatch) || is_TofMult_bad(ref3, tofMult);
        }

        bool is_BetaEta1_good(Int_t ref3, Int_t beta_eta1){
            return !is_BetaEta1_bad(ref3, beta_eta1);
        }

        bool is_nTofMatch_good(Int_t ref3, Int_t nTofMatch){
            return !is_nTofMatch_bad(ref3, nTofMatch);
        }

        bool is_TofMult_good(Int_t ref3, Int_t tofMult){
            return !is_TofMult_bad(ref3, tofMult);
        }

        void set_lumi_par(Int_t trgid, Double_t b, Double_t k){
            Int_t trg = convert_trg(trgid);
            if (trg != -1){
                k_lumi[trg] = k;
                b_lumi[trg] = b;
            }
        }

        Double_t lumi_correction(Int_t trg, Int_t ref3, Double_t bbc){
            if (trg >= nTrg || trg < 0){
                return -1;
            }
            Double_t f = bbc * k_lumi[trg] + b_lumi[trg];
            Double_t factor = 0.0;
            if (f != 0){
                factor = b_lumi[trg] / f;
            } 
            return (Int_t)(factor * ref3);
        }

        void set_vz_par(Int_t trgid, Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4, Double_t p5, Double_t p6){
            Int_t trg = convert_trg(trgid);
            par[trg][0] = p0;
            par[trg][1] = p1;
            par[trg][2] = p2;
            par[trg][3] = p3;
            par[trg][4] = p4;
            par[trg][5] = p5;
            par[trg][6] = p6;
            hfunc->SetParameters(&par[trg][0]);
        }

        Double_t vz_correction(Int_t trg, Int_t refmult3, Double_t vz){
            Double_t factor = par[trg][0] / hfunc->Eval(vz);
            return (Int_t)(refmult3*factor);
        }

        Int_t get_corr_refmult3(Int_t ref3, Int_t beta_eta1, Int_t nTofMatch, Int_t tofMult, Double_t bbc, Double_t vz, Int_t trgid){
            Int_t trg = convert_trg(trgid);
            if (trg < 0){
                return -1;
            }
            if (fabs(vz) > 70.0){
                return -1;
            }
            if (do_pile_up && is_pile_up(ref3, beta_eta1, nTofMatch, tofMult)){
                return -1;
            }
            if (do_lumi){
                ref3 = lumi_correction(trg, ref3, bbc);
            }
            if (do_vz){
                ref3 = vz_correction(trg, ref3, vz);
            }
            return (Int_t)ref3;
        }

        void set_cent_edge(Int_t e0, Int_t e1, Int_t e2, Int_t e3, Int_t e4, Int_t e5, Int_t e6, Int_t e7, Int_t e8){
            edge[0] = e0;
            edge[1] = e1;
            edge[2] = e2;
            edge[3] = e3;
            edge[4] = e4;
            edge[5] = e5;
            edge[6] = e6;
            edge[7] = e7;
            edge[8] = e8;
            do_split = true;
            cout << "[LOG] Centrality bin edge specified.\n";
        }

        void set_cent_edge(Int_t* arr){
            for (int i=0; i<9; i++){
                edge[i] = *(arr+i);
            }
            do_split = true;
            cout << "[LOG] Centrality bin edge specified.\n";
        }

        Int_t get_centrality9(Int_t ref3){
            if (!do_split){
                cout << "[WARNING] Centrality bin edge not specified, can not get correct centrality bin.\n";
                return -1;
            }
            for (int i=0; i<9; i++){
                if (ref3 > edge[i]){
                    return i;
                }
            }
            return -1;
        }
};

#endif