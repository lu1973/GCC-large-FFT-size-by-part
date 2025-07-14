#pragma once

#include <string>
#include <Eigen/Dense>
#include <vector>

class ClongGCC
{
private:
    class CParams {
    public:
        int m_iFS;
        int m_iN;				        ///< length of the large FFT
        int m_iN1;				        ///< buffer length
        double m_dLambda;               ///< smoothing factor
        int m_iMaxlag;                  ///< maximum correlation lag for peak search (used as positive and as negative value)
        double rho;                     ///< whitening parameter for rho-CSP GCC
        double m;                       ///< root power for CSP-m GCC
        std::string m_sGccWeight;       ///< GCC method : CSPm / rhoCSP / PHAT

        int m_iN2;
        Eigen::MatrixXcd m_N2_pt_DFT_matrix;
        Eigen::MatrixXcd m_twid;
        Eigen::MatrixXcd m_N1_pt_DFT_matrix;
        Eigen::VectorXd m_Mask;

        int m_ibMax_N1;
        int m_ibMin_N1;

        Eigen::MatrixXcd m_N1_pt_DFT_matrix_reduced_up;
        Eigen::MatrixXcd m_N1_pt_DFT_matrix_reduced_low;

        int m_bInv;
        Eigen::MatrixXcd m_N1_pt_IDFT_matrix_reduced_up;
        Eigen::MatrixXcd m_N1_pt_IDFT_matrix_reduced_low;
    };

    class CStates {
        public:
            int m_iCurFrameIdx;

            Eigen::MatrixXcd m_matrix_sig_out_reduced;
            Eigen::MatrixXcd m_matrix_sigref_out_reduced;
            Eigen::MatrixXd m_matrix_ifft_out_reduced;

            Eigen::VectorXcd m_gcc_phat;
            double m_dDelay;
    };

public:
    ClongGCC(int _iBuffer);
    ~ClongGCC() {};

    double process(Eigen::VectorXd& _x, Eigen::VectorXd& _r);

    CParams m_cParam;
    CStates m_cStates;

protected:
    void partial_mic_fft(Eigen::VectorXd& _y);
    void partial_ref_fft(Eigen::VectorXd& _y);
    void partial_ifft(Eigen::VectorXcd& _y);
};

