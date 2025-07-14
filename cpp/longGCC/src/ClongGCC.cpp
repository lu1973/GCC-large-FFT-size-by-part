#include "ClongGCC.h"
#include <algorithm>
#include <cmath>
#include <format>
#include <Eigen/Core>
#include <iostream>

constexpr double weakEpsilon = Eigen::NumTraits<double>::dummy_precision();

ClongGCC::ClongGCC(int _iBuffer)
{
    m_cParam.m_iFS = 16000;
    m_cParam.m_iN = 42 * _iBuffer;
    m_cParam.m_iN1 = _iBuffer;
    m_cParam.m_dLambda = 0.1;
    m_cParam.m_iMaxlag = 1000;
    m_cParam.rho = 0.7;
    m_cParam.m = 4;
    m_cParam.m_sGccWeight = "PHAT";

    m_cParam.m_iN2 = m_cParam.m_iN / m_cParam.m_iN1;

    m_cParam.m_N2_pt_DFT_matrix.resize(m_cParam.m_iN2, m_cParam.m_iN2);
    std::complex<double> z(0, -2 * double(EIGEN_PI) / m_cParam.m_iN2);
    for (int i = 0; i < m_cParam.m_iN2; i++)
    {
        for (int j = 0; j < m_cParam.m_iN2; j++)
        {
            m_cParam.m_N2_pt_DFT_matrix(i, j) = pow(exp(z), i * j);
        }
    }

    m_cParam.m_twid.resize(m_cParam.m_iN1, m_cParam.m_iN2);
    z = std::complex<double>(0, -2 * double(EIGEN_PI) / m_cParam.m_iN);
    for (int i = 0; i < m_cParam.m_iN1; i++)
    {
        for (int j = 0; j < m_cParam.m_iN2; j++)
        {
            m_cParam.m_twid(i, j) = pow(exp(z), i * j);
        }
    }

    m_cParam.m_N1_pt_DFT_matrix.resize(m_cParam.m_iN1, m_cParam.m_iN1);
    z = std::complex<double>(0, -2 * double(EIGEN_PI) / m_cParam.m_iN1);
    for (int i = 0; i < m_cParam.m_iN1; i++)
    {
        for (int j = 0; j < m_cParam.m_iN1; j++)
        {
            m_cParam.m_N1_pt_DFT_matrix(i, j) = pow(exp(z), i * j);
        }
    }

    // Definition of the frequency mask
    int f0 = 300;
    int f1 = 6000;
    int bMin = std::max(0, int(f0 * m_cParam.m_iN / m_cParam.m_iFS));
    int bMax = int(f1 * m_cParam.m_iN / m_cParam.m_iFS);

    m_cParam.m_ibMax_N1 = int(ceil(bMax / m_cParam.m_iN2));
    int iWinLength = bMax - bMin + 1;
    Eigen::VectorXd windowMask(iWinLength);
    for (int n = 0; n < iWinLength; ++n) {
        windowMask(n) = 0.5 - 0.5 * std::cos(2 * double(EIGEN_PI) * n / (iWinLength - 1));
    }

    m_cParam.m_Mask.resize(m_cParam.m_iN);
    m_cParam.m_Mask.setZero();
    m_cParam.m_Mask(Eigen::seq(bMin, iWinLength + bMin - 1)) = windowMask;
    m_cParam.m_Mask(Eigen::seqN(Eigen::last, m_cParam.m_iN / 2, Eigen::fix<-1>)) = m_cParam.m_Mask(Eigen::seq(0, m_cParam.m_iN / 2 - 1));

    m_cParam.m_ibMin_N1 = std::max(int(floor(bMin / m_cParam.m_iN2)), 0);

    m_cParam.m_N1_pt_DFT_matrix_reduced_up.resize(m_cParam.m_ibMax_N1 - m_cParam.m_ibMin_N1 + 1, m_cParam.m_iN1);
    m_cParam.m_N1_pt_DFT_matrix_reduced_up = m_cParam.m_N1_pt_DFT_matrix(Eigen::seq(m_cParam.m_ibMin_N1-1, m_cParam.m_ibMax_N1-1), Eigen::all);

    m_cParam.m_N1_pt_DFT_matrix_reduced_low.resize(m_cParam.m_ibMax_N1 - m_cParam.m_ibMin_N1 + 1, m_cParam.m_iN1);
    m_cParam.m_N1_pt_DFT_matrix_reduced_low = m_cParam.m_N1_pt_DFT_matrix(Eigen::seq(m_cParam.m_iN1 - m_cParam.m_ibMax_N1, m_cParam.m_iN1 - m_cParam.m_ibMin_N1), Eigen::all);

    // for IFFT only
    m_cParam.m_bInv = int(ceil(m_cParam.m_iMaxlag / double(m_cParam.m_iN2)));
    m_cParam.m_N1_pt_IDFT_matrix_reduced_up.resize(m_cParam.m_bInv, m_cParam.m_iN1);
    m_cParam.m_N1_pt_IDFT_matrix_reduced_up = m_cParam.m_N1_pt_DFT_matrix(Eigen::seq(0, m_cParam.m_bInv - 1), Eigen::all);

    m_cParam.m_N1_pt_IDFT_matrix_reduced_low.resize(m_cParam.m_bInv, m_cParam.m_iN1);
    m_cParam.m_N1_pt_IDFT_matrix_reduced_low = m_cParam.m_N1_pt_DFT_matrix(Eigen::seq(m_cParam.m_iN1 - m_cParam.m_bInv-1, Eigen::last), Eigen::all);

    m_cStates.m_iCurFrameIdx = 0;

    m_cStates.m_matrix_sig_out_reduced.resize(m_cParam.m_iN1, m_cParam.m_iN2);
    m_cStates.m_matrix_sig_out_reduced.setZero();

    m_cStates.m_matrix_sigref_out_reduced.resize(m_cParam.m_iN1, m_cParam.m_iN2);
    m_cStates.m_matrix_sigref_out_reduced.setZero();

    m_cStates.m_matrix_ifft_out_reduced.resize(m_cParam.m_iN1, m_cParam.m_iN2);
    m_cStates.m_matrix_ifft_out_reduced.setZero();

    m_cStates.m_gcc_phat.resize(m_cParam.m_iN);       //  GCC in frequency domain
    m_cStates.m_gcc_phat.setZero();

    m_cStates.m_dDelay = 0;
}

/// \brief compute GCC
/// \param[in] _x : input signal
/// \param[in] _r : reference signal
/// \param[out] delay estimated between _x and _r
double ClongGCC::process(Eigen::VectorXd &_x, Eigen::VectorXd &_r)
{
    // compute partial FFT for the signal and the reference
    partial_mic_fft(_x);
    partial_ref_fft(_r);

    // compute partial IFFT of the GCC
    Eigen::VectorXcd l_gcc_part = m_cStates.m_gcc_phat(Eigen::seq(m_cStates.m_iCurFrameIdx * m_cParam.m_iN1, (m_cStates.m_iCurFrameIdx+1) * m_cParam.m_iN1-1));
    partial_ifft(l_gcc_part);

    // compute full FFT and GCC in the frequency domain
    if (m_cStates.m_iCurFrameIdx == (m_cParam.m_iN2-1))
    {
        // Flatten the matrix into a column vector (column-major order)
        Eigen::VectorXd temp = Eigen::VectorXd{ m_cStates.m_matrix_ifft_out_reduced.transpose().reshaped() };

        Eigen::VectorXd gcc_phat_time = temp / m_cParam.m_iN;
        Eigen::VectorXd gcc_phat_time_flipped(gcc_phat_time.size());
        gcc_phat_time_flipped(0) = gcc_phat_time(0);
        gcc_phat_time_flipped(Eigen::seq(1, Eigen::last)) = gcc_phat_time(Eigen::lastN(m_cParam.m_iN - 1).reverse());

        Eigen::VectorXd gcc_time_abs = gcc_phat_time_flipped.cwiseAbs().array();

        Eigen::VectorXd gcc_phat_time_abs(2 * m_cParam.m_iMaxlag + 1);
        gcc_phat_time_abs << gcc_time_abs(Eigen::lastN(m_cParam.m_iMaxlag + 1)), gcc_time_abs(Eigen::seq(0, m_cParam.m_iMaxlag-1));
        gcc_phat_time_abs = gcc_phat_time_abs.cwiseAbs();
        
        // Find index of the max
        Eigen::Index idxMax;
        double gccMax = gcc_phat_time_abs.maxCoeff(&idxMax);

        if (gccMax > 0)
        {
            // Lagrange polynomial peak extraction
            double g1 = gcc_phat_time_abs(idxMax - 1);
            double g2 = gccMax;
            double g3 = gcc_phat_time_abs(idxMax + 1);
            double f = (g1 - g3) / (2 * g1 - 4 * g2 + 2 * g3);
            // bias correction
            if (abs(f) > weakEpsilon)
            {
                f = (sqrt(32 * f * f + 1) - 1) / (8 * f);
            }
            // compensate offset
            m_cStates.m_dDelay = (int)idxMax + f - m_cParam.m_iMaxlag - 1;
        }
        m_cStates.m_matrix_ifft_out_reduced = m_cParam.m_dLambda * m_cStates.m_matrix_ifft_out_reduced;

        // Flatten the matrix into a column vector (column-major order)
        Eigen::VectorXcd xFFT = Eigen::VectorXcd{ m_cStates.m_matrix_sig_out_reduced.transpose().reshaped()};

        // introduce here freq masking
        Eigen::VectorXcd xFFT_std = xFFT.cwiseProduct(m_cParam.m_Mask);

        Eigen::VectorXcd rFFT = Eigen::VectorXcd{ m_cStates.m_matrix_sigref_out_reduced.transpose().reshaped() };
        // introduce here freq masking
        rFFT = rFFT.cwiseProduct(m_cParam.m_Mask);
        Eigen::VectorXcd est_corr = xFFT_std.cwiseProduct(rFFT.conjugate());

        Eigen::VectorXd abs_corr_1 = xFFT_std.cwiseAbs().array().max(weakEpsilon);
        Eigen::VectorXd abs_corr_2 = rFFT.cwiseAbs().array().max(weakEpsilon);
        Eigen::VectorXd abs_corr   = est_corr.cwiseAbs().array().max(weakEpsilon);

        if (m_cParam.m_sGccWeight == "rhoCSP")
        {
            m_cStates.m_gcc_phat = est_corr.array() / (abs_corr.array().pow(m_cParam.rho));
        }
        else if (m_cParam.m_sGccWeight == "CSPm")
        {
            m_cStates.m_gcc_phat = est_corr.array() / ((abs_corr_1.cwiseProduct(abs_corr_2)).array().pow(1 / m_cParam.m));
        }
        else
        {
            m_cStates.m_gcc_phat = est_corr.array() / abs_corr.array();
        }

        m_cStates.m_iCurFrameIdx = -1;

        m_cStates.m_matrix_sig_out_reduced    = m_cParam.m_dLambda * m_cStates.m_matrix_sig_out_reduced;
        m_cStates.m_matrix_sigref_out_reduced = m_cParam.m_dLambda * m_cStates.m_matrix_sigref_out_reduced;
    }
    m_cStates.m_iCurFrameIdx++;

    return m_cStates.m_dDelay;
}

/// \brief Compute partial FFT for selected mic signal
/// \param[in] _y : input signal
void ClongGCC::partial_mic_fft(Eigen::VectorXd& _y)
{
    // step 1: performing N2 - point DFTs of each one the Ni rows
    // optimization remark : in Matlab, it is better to do this large
    // multiplication than using the fact the N2 matrix is
    // symmetrical and doing a memory operation(flip + conj)
    Eigen::MatrixXcd matrix1 = _y * (m_cParam.m_N2_pt_DFT_matrix(Eigen::all, m_cStates.m_iCurFrameIdx).transpose());
    // step 2: multiplying elements of matrix with twiddle factor W_N* { n1k2 }
    Eigen::MatrixXcd matrix2 = matrix1.cwiseProduct(m_cParam.m_twid);
    // step 3: performing Nl - point DFTs of each column of array
    Eigen::MatrixXcd prod_up = m_cParam.m_N1_pt_DFT_matrix_reduced_up * matrix2;
    m_cStates.m_matrix_sig_out_reduced(Eigen::seq(m_cParam.m_ibMin_N1-1, m_cParam.m_ibMax_N1-1), Eigen::all ) = m_cStates.m_matrix_sig_out_reduced(Eigen::seq(m_cParam.m_ibMin_N1-1, m_cParam.m_ibMax_N1-1), Eigen::all) + (1 - m_cParam.m_dLambda) * prod_up;
    Eigen::MatrixXcd prod_down = m_cParam.m_N1_pt_DFT_matrix_reduced_low * matrix2;
    m_cStates.m_matrix_sig_out_reduced(Eigen::seq(m_cParam.m_iN1 - m_cParam.m_ibMax_N1, m_cParam.m_iN1 - m_cParam.m_ibMin_N1), Eigen::all) = m_cStates.m_matrix_sig_out_reduced(Eigen::seq(m_cParam.m_iN1 - m_cParam.m_ibMax_N1, m_cParam.m_iN1 - m_cParam.m_ibMin_N1), Eigen::all) + (1 - m_cParam.m_dLambda) * prod_down;

    //saveCplxMatrixBin(m_cStates.m_matrix_sig_out_reduced, "E:\\devC\\longGCC\\data\\m_matrix_sig_out_reduced.bin");
}

/// \brief Compute partial FFT for selected ref signal
/// \param[in] _y : input signal
void ClongGCC::partial_ref_fft(Eigen::VectorXd& _y)
{
    // step 1: performing N2 - point DFTs of each one the Ni rows
    // optimization remark : in Matlab, it is better to do this large
    // multiplication than using the fact the N2 matrix is
    // symmetrical and doing a memory operation(flip + conj)
    Eigen::MatrixXcd matrix1 = _y * (m_cParam.m_N2_pt_DFT_matrix(Eigen::all, m_cStates.m_iCurFrameIdx).transpose());
    // step 2: multiplying elements of matrix with twiddle factor W_N* { n1k2 }
    Eigen::MatrixXcd matrix2 = matrix1.cwiseProduct(m_cParam.m_twid);
    // step 3: performing Nl - point DFTs of each column of array
    Eigen::MatrixXcd prod_up = m_cParam.m_N1_pt_DFT_matrix_reduced_up * matrix2;
    m_cStates.m_matrix_sigref_out_reduced(Eigen::seq(m_cParam.m_ibMin_N1-1, m_cParam.m_ibMax_N1-1), Eigen::all) = m_cStates.m_matrix_sigref_out_reduced(Eigen::seq(m_cParam.m_ibMin_N1-1, m_cParam.m_ibMax_N1-1), Eigen::all) + (1 - m_cParam.m_dLambda) * prod_up;
    Eigen::MatrixXcd prod_down = m_cParam.m_N1_pt_DFT_matrix_reduced_low * matrix2;
    m_cStates.m_matrix_sigref_out_reduced(Eigen::seq(m_cParam.m_iN1 - m_cParam.m_ibMax_N1, m_cParam.m_iN1 - m_cParam.m_ibMin_N1), Eigen::all) = m_cStates.m_matrix_sigref_out_reduced(Eigen::seq(m_cParam.m_iN1 - m_cParam.m_ibMax_N1, m_cParam.m_iN1 - m_cParam.m_ibMin_N1), Eigen::all) + (1 - m_cParam.m_dLambda) * prod_down;
}

/// \brief Compute partial IFFT for selected ref channel
/// \param[in] _y : input signal
void ClongGCC::partial_ifft(Eigen::VectorXcd& _y)
{
    // step 1: performing N2 - point DFTs of each one the Ni rows
    Eigen::MatrixXcd matrix1 = _y * (m_cParam.m_N2_pt_DFT_matrix(Eigen::all, m_cStates.m_iCurFrameIdx).transpose());
    // step 2: multiplying elements of matrix with twiddle factor W_N* { n1k2 }
    Eigen::MatrixXcd matrix2 = matrix1.cwiseProduct(m_cParam.m_twid);
    // step 3: performing Nl - point DFTs of each column of array
    Eigen::MatrixXcd prod_up = m_cParam.m_N1_pt_IDFT_matrix_reduced_up * matrix2;
    m_cStates.m_matrix_ifft_out_reduced(Eigen::seq(0, m_cParam.m_bInv-1), Eigen::all) = m_cStates.m_matrix_ifft_out_reduced(Eigen::seq(0, m_cParam.m_bInv-1), Eigen::all) + (1 - m_cParam.m_dLambda) * prod_up.real();
    Eigen::MatrixXcd prod_down = m_cParam.m_N1_pt_IDFT_matrix_reduced_low * matrix2;
    m_cStates.m_matrix_ifft_out_reduced(Eigen::seq(m_cParam.m_iN1 - m_cParam.m_bInv-1, Eigen::last), Eigen::all) = m_cStates.m_matrix_ifft_out_reduced(Eigen::seq(m_cParam.m_iN1 - m_cParam.m_bInv-1, Eigen::last), Eigen::all) + (1 - m_cParam.m_dLambda) * prod_down.real();
}