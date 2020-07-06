#ifndef KAIROSFFT_H
#define KAIROSFFT_H

#ifdef __cplusplus
extern "C" {
#endif

    #include <math.h>
    #include <stdlib.h>
    #include <string.h>

    typedef double real_t;
    
    typedef struct{
        real_t *Buffer;
        int WindowLength;
        real_t SampleTime;
        real_t OMEGA;
        real_t BaseFrequency;
        real_t *armonics;
        int NA;
        real_t *reH, *imH;
        real_t t, t_1;
        int Identifier;
        char *name;
        double latency;
        double NSAMPL;
    }RFFT_KairosConfig_t;

    
    #define MAX_ARMONICS_TO_CHECK        30
    typedef struct{
        real_t Min[MAX_ARMONICS_TO_CHECK];
        real_t Max[MAX_ARMONICS_TO_CHECK];
    }Kairos_Thresholds_t;

    typedef struct{
        Kairos_Thresholds_t f1,f2;
    }Kairos_PhaseThresholds_t;

    typedef union{
        struct{
            unsigned A:1;
            unsigned B:1;
            unsigned C:1;
            unsigned X:1;
            unsigned Y:1;
        };
        unsigned char reg;
    }Kairos_DampingPhaseState_t;
    extern Kairos_DampingPhaseState_t AlarmSignal;
    
    
    int FFT_KairosRecursive_Init(RFFT_KairosConfig_t *obj, real_t *armonics, const int NA, const real_t SampleTime, const real_t BaseFrequency, real_t *WorkSpace, real_t *DelayBuffer, int Identifier, char *name);
    void FFT_KairosRecursive_Iterate(RFFT_KairosConfig_t *obj, const real_t new_data);
    void FFT_KairosRecursive_End(RFFT_KairosConfig_t *obj);
    void FFT_KairosRecursive_Flush(RFFT_KairosConfig_t *obj);
    Kairos_DampingPhaseState_t DecisionLogic_Kairos(real_t *Af1Phase_a, real_t *Af2Phase_a, real_t *Af1Phase_b, real_t *Af2Phase_b, real_t *Af1Phase_c, real_t *Af2Phase_c, Kairos_PhaseThresholds_t *obj, const int nCheck);
    real_t TuneFrequency(real_t SampleTime, real_t BaseFrequency, real_t TrueSampleTime);
    

    
#ifdef __cplusplus
}
#endif

#endif /* KAIROSFFT_H */

