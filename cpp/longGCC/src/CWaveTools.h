#pragma once
#include <windows.h>
#include <string>
#include <vector>
#include <Eigen/Dense>

// Definition of WAVE format header
typedef struct
{
  char ChunkID[4];
  long ChunkSize;
  char Format[4];
  char Subchunk1ID[4];
  long Subchunk1Size;
  short AudioFormat;
  short NumChannels;
  long SampleRate;
  long ByteRate;
  short BlockAlign;
  short BitsPerSample;
} SWaveHeader;

class CWaveTools
{
public:
	CWaveTools(void);
	~CWaveTools(void) {};

	static CWaveTools * GetInstance();	///< design pattern singleton... 
	static void DestroyInstance();	///< design pattern singleton... 

	HANDLE OpenWave(std::string sFileName, unsigned long _ulOffset, unsigned long *ulSize_, SWaveHeader *sWaveHeader_);

	int ReadWave16bits(HANDLE hFile, Eigen::VectorXd &dData, int _iSize, unsigned long _ulPos);

protected:
	static CWaveTools * m_hInstance;	///< l'instance en cours de CWaveTools.
};

