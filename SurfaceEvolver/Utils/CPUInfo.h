
/*
	Check SSE/AVX support.
	This application can detect the instruction support of
	SSE, SSE2, SSE3, SSSE3, SSE4.1, SSE4.2, SSE4a, SSE5, and AVX.
*/

#include <iostream>
#ifdef _MSC_VER
#include <intrin.h>
#include <array>
#include <string>
#endif

#ifdef __GNUC__

void __cpuid(int* cpuinfo, int info)
{
	__asm__ __volatile__(
		"xchg %%ebx, %%edi;"
		"cpuid;"
		"xchg %%ebx, %%edi;"
		:"=a" (cpuinfo[0]), "=D" (cpuinfo[1]), "=c" (cpuinfo[2]), "=d" (cpuinfo[3])
		: "0" (info)
	);
}

unsigned long long _xgetbv(unsigned int index)
{
	unsigned int eax, edx;
	__asm__ __volatile__(
		"xgetbv;"
		: "=a" (eax), "=d"(edx)
		: "c" (index)
	);
	return ((unsigned long long)edx << 32) | eax;
}

#endif

#ifndef CPUINFO_H_
#define CPUINFO_H_

std::string GetCPUInfo()
{
    // 4 is essentially hardcoded due to the __cpuid function requirements.
    // NOTE: Results are limited to whatever the sizeof(int) * 4 is...
    std::array<int, 4> integerBuffer = {};
    constexpr size_t sizeofIntegerBuffer = sizeof(int) * integerBuffer.size();

    std::array<char, 64> charBuffer = {};

    // The information you wanna query __cpuid for.
    // https://docs.microsoft.com/en-us/cpp/intrinsics/cpuid-cpuidex?view=vs-2019
    constexpr std::array<int, 3> functionIds = {
        // Manufacturer
        //  EX: "Intel(R) Core(TM"
        0x8000'0002,
        // Model
        //  EX: ") i7-8700K CPU @"
        0x8000'0003,
        // Clockspeed
        //  EX: " 3.70GHz"
        0x8000'0004
    };

    std::string cpu;

    for (int id : functionIds)
    {
        // Get the data for the current ID.
        __cpuid(integerBuffer.data(), id);

        // Copy the raw data from the integer buffer into the character buffer
        std::memcpy(charBuffer.data(), integerBuffer.data(), sizeofIntegerBuffer);

        // Copy that data into a std::string
        cpu += std::string(charBuffer.data());
    }

    return cpu;
}

/**/
bool isAVXSupported() {
	int cpuinfo[4];
	__cpuid(cpuinfo, 1);

	bool avxSupported = cpuinfo[2] & (1 << 28) || false;
	bool osxsaveSupported = cpuinfo[2] & (1 << 27) || false;
	if (osxsaveSupported && avxSupported)
	{
		// _XCR_XFEATURE_ENABLED_MASK = 0
		unsigned long long xcrFeatureMask = _xgetbv(0);
		avxSupported = (xcrFeatureMask & 0x6) == 0x6;
	}

	return avxSupported;
}

static std::string FilterSpecialSymbolsFromString(const std::string& inputString)
{
    std::string final;
    for (const auto& symbol : inputString)
    {
        if ((symbol >= '0' && symbol <= '9') || (symbol >= 'a' && symbol <= 'z') || (symbol >= 'A' && symbol <= 'Z') || symbol == '_')
        {
            final += symbol;
        }
    }
    return final;
}

std::string GetParsedCPUName()
{
	const auto unparsedResult = GetCPUInfo();
    return FilterSpecialSymbolsFromString(unparsedResult);
}

#endif