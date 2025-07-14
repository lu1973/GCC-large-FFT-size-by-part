#include <iostream>
#include <string>
#include "CWaveTools.h"
#include "ClongGCC.h"

using namespace std;

int main()
{
    std::cout<<"Computation of GCC on a very large length buffer by using small frames." << std::endl;

    std::string file_name_x ("..\\..\\data\\226596.wav");
    unsigned long l_ulNbSamples_x;
    SWaveHeader sWaveHeader_x;
    HANDLE l_hWavFileIn_x = CWaveTools::GetInstance()->OpenWave(file_name_x, 0, &l_ulNbSamples_x, &sWaveHeader_x);

    std::string file_name_r("..\\..\\data\\226596_shifted.wav");
    unsigned long l_ulNbSamples_r;
    SWaveHeader sWaveHeader_r;
    HANDLE l_hWavFileIn_r = CWaveTools::GetInstance()->OpenWave(file_name_r, 0, &l_ulNbSamples_r, &sWaveHeader_r);

    unsigned long l_ulNbSamples = min(l_ulNbSamples_x, l_ulNbSamples_r);

    int l_iBuffer = 512;
    unsigned long l_ulPos = 0;

    Eigen::VectorXd l_Data_x;
    l_Data_x.resize(l_iBuffer);
    Eigen::VectorXd l_Data_r;
    l_Data_r.resize(l_iBuffer);

    ClongGCC* l_clongGCC;
    l_clongGCC = new ClongGCC(l_iBuffer);

    int l_iNframe = int(floor((l_ulNbSamples - l_iBuffer) / double(l_iBuffer)));
    Eigen::VectorXd resDelay;
    resDelay.resize(l_iNframe);

    for (int iF = 0 ; iF < l_iNframe ; iF++)
    {
        CWaveTools::GetInstance()->ReadWave16bits(l_hWavFileIn_x, l_Data_x, l_iBuffer, l_ulPos);
        CWaveTools::GetInstance()->ReadWave16bits(l_hWavFileIn_r, l_Data_r, l_iBuffer, l_ulPos);

        resDelay(iF) = l_clongGCC->process(l_Data_x, l_Data_r);

        if (iF > 1 && resDelay(iF) != resDelay(iF - 1))
        {
            std::cout << resDelay(iF) << std::endl;
        }


        l_ulPos += l_iBuffer;
    };

    delete l_clongGCC;

    if (l_hWavFileIn_x != nullptr)
        CloseHandle(l_hWavFileIn_x);
    if (l_hWavFileIn_r != nullptr)
        CloseHandle(l_hWavFileIn_r);

    return 0;
}