#include "CWaveTools.h"

// current instance of CWaveTools
CWaveTools *CWaveTools::m_hInstance = NULL;

CWaveTools::CWaveTools(void)
{
}

CWaveTools * CWaveTools::GetInstance()
{
	if (m_hInstance == NULL) {
		m_hInstance = new CWaveTools();
	}
	return m_hInstance;
}

void CWaveTools::DestroyInstance()
{
	if (m_hInstance)
		delete m_hInstance;
}



/// \brief  Open a WAVE file and return a handle on the file with the pointer directly positionned on the data
/// \param[in] sFileName : filename
/// \param[in] _ulOffset : offset to start reading data
/// \param[out] ulSize_ : number of samples to read
/// \param[out] sWaveHeader_ : WAVE header
HANDLE CWaveTools::OpenWave(std::string sFileName, unsigned long _ulOffset, unsigned long *ulSize_, SWaveHeader *sWaveHeader_)
{
	// file opening
	HANDLE hFile = CreateFile(sFileName.c_str(), GENERIC_READ, 0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
	if (hFile == INVALID_HANDLE_VALUE) {
		return NULL;
	}

	DWORD dwRead;
	BYTE bHdr[512];
	int Offset = -1;

	// Search header to jump in wav file
	ReadFile(hFile, &bHdr, 512*sizeof(BYTE), &dwRead, NULL);
	int nRead, k;
	nRead = int(dwRead / sizeof(BYTE));
	BYTE bData[4];
	bData[0] = 'd';
	bData[1] = 'a';
	bData[2] = 't';
	bData[3] = 'a';
	BOOL bFindHeader = FALSE;
	for (k = 0 ; k < nRead ; k++) {
		if (memcmp(&bHdr[k], &bData, 4*sizeof(BYTE)) == 0) {
			bFindHeader = TRUE;
			Offset = k;
			break;
		}
	}

	if (bFindHeader == FALSE) {
		CloseHandle(hFile);
		return NULL;
	}

	SetFilePointer(hFile, 0, NULL, FILE_BEGIN);

	// store the header
	BYTE *WaveHdr;
	WaveHdr = new BYTE[Offset+8];
	ReadFile(hFile, WaveHdr, Offset+8, &dwRead, NULL);

	memcpy(sWaveHeader_, WaveHdr, sizeof(SWaveHeader));

	// get the Bytes number of the useful data only
	SetFilePointer(hFile, Offset+4, NULL, FILE_BEGIN);
	unsigned long ChunkSize;
	ReadFile(hFile, &ChunkSize, sizeof(unsigned long), &dwRead, NULL);

	(*ulSize_) = ChunkSize / (sWaveHeader_->BitsPerSample/8);

	// jump the beginning of the file to reach the useful data
	SetFilePointer(hFile, Offset+8 + (_ulOffset*sWaveHeader_->SampleRate*sWaveHeader_->BitsPerSample/8) , NULL, FILE_BEGIN);

	delete[] WaveHdr;

	return hFile;
}

/// \brief Read WAVE data (16 bits format) from a file handle given by the OpenWave function
/// \details return the data at the position _ulPos (0 is the beginning of the data)
/// The position pointer still remain at the position _ulPos+iSize after return
/// \param[in] hFile : file handle
/// \param[out] dData : vector with sample values
/// \param[in] _iSize : sample number to read
/// \param[in] _ulPos : offset in sample
int CWaveTools::ReadWave16bits(HANDLE hFile, Eigen::VectorXd &dData_, int _iSize, unsigned long _ulPos)
{
	DWORD dwRead;
	int iRead = 0;			// number of sample read

	// Get current position
	static DWORD dwBeginPos = SetFilePointer(hFile, 0, NULL, FILE_CURRENT);

	// Set new position
	SetFilePointer(hFile, dwBeginPos + _ulPos * sizeof(int16_t), NULL, FILE_BEGIN);

	int16_t* uData = NULL;
	uData = new int16_t[_iSize];

	// read buffer in 16 bits format
	ReadFile(hFile, uData, _iSize * sizeof(int16_t), &dwRead, NULL);

	iRead = dwRead / sizeof(int16_t);
	for (int SmpI = 0; SmpI < iRead ; SmpI ++) {
		dData_(SmpI) = double(float((uData[SmpI]) / 32768.0));
	}

	delete[] uData;

	return iRead;
}
