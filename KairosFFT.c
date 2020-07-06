#include "KairosFFT.h"

inline static real_t _sin(real_t x);
static real_t M_2PI = 2.0*M_PI;
#define _cos(_x_)  _sin((_x_)+M_PI_2)  // cos(x) = sin(x+pi/2) 
Kairos_DampingPhaseState_t AlarmSignal = {0};
/*================================================================================================================================================================*/
real_t FFT_KairosIO_delay(const real_t *DelayBuffer, real_t t){
   real_t y, y1, y2;
   int t1;
   t1 = (int)t;
   if((real_t)t1 == t){
       return DelayBuffer[t1-1];
   }
   else{
        y1 = DelayBuffer[t1-1];
        y2 = DelayBuffer[t1];  
        y = ( (t-t1)*(y2-y1) ) + y1;    
        return y;
    } 
}
/*================================================================================================================================================================*/
int FFT_KairosRecursive_Init(RFFT_KairosConfig_t *obj, real_t *armonics, const int NA, const real_t SampleTime, const real_t BaseFrequency, real_t *WorkSpace, real_t *DelayBuffer, int Identifier, char *name){
    memset(obj, 0, sizeof(RFFT_KairosConfig_t));
    obj->armonics = armonics;
    memset(WorkSpace, 0 , (2*NA)*sizeof(real_t));    
    obj->reH = WorkSpace;
    obj->imH = WorkSpace+NA;
    obj->NA = NA;
    obj->BaseFrequency = BaseFrequency;
    obj->SampleTime = SampleTime;
    obj->OMEGA = M_2PI*BaseFrequency;
    obj->NSAMPL = 1.0/(BaseFrequency*SampleTime);
    obj->WindowLength = (int)obj->NSAMPL+3;
    //obj->Buffer = (real_t *)calloc(obj->WindowLength, sizeof(real_t)); // heap slow, faster using stack
    obj->Buffer = DelayBuffer;
    obj->Identifier = Identifier;
    obj->name = name;
    if(obj->Buffer == NULL) return -1;
    return 0;
}
/*================================================================================================================================================================*/
void FFT_KairosRecursive_Iterate(RFFT_KairosConfig_t *obj, const real_t new_data){
    int i, j;
    real_t f1,f2,f3,f4, jOM, jOMt, OMt, jOMt_1, StjOM, sin_jOMt, cos_jOMt, sin_jOMt_1, cos_jOMt_1, a, b, c, n1_iMPi, c_StjOM;
    real_t absri, firstA = 0.0;

    for (i=obj->WindowLength-1;i>=1;i--) obj->Buffer[i]=obj->Buffer[i-1]; //samples shift
    obj->Buffer[0]=new_data;   
    
    obj->t += obj->SampleTime;
    f4 = obj->Buffer[0];
    f3 = obj->Buffer[1];
    f2 = FFT_KairosIO_delay(obj->Buffer, obj->NSAMPL+1); // last delay values using linear interpolation
    f1 = FFT_KairosIO_delay(obj->Buffer, obj->NSAMPL+2); // last delay values using linear interpolation        
    //f2 = obj->Buffer[obj->WindowLength-1];
    //f1 = obj->Buffer[obj->WindowLength-2];
    a = f4-f2;
    b = f3-f1;
    c = a-b;
    OMt = obj->OMEGA*obj->t;
    for(i=0,j=1;i<obj->NA;i++,j++){
        n1_iMPi = M_1_PI/j;
        jOM = j*obj->OMEGA;
        jOMt = j*OMt;
        sin_jOMt = _sin(jOMt);
        cos_jOMt = _cos(jOMt);
        jOMt_1 = jOM*obj->t_1;
        sin_jOMt_1 = _sin(jOMt_1);
        cos_jOMt_1 = _cos(jOMt_1);
        StjOM = obj->SampleTime*jOM;
        c_StjOM = c/StjOM;
        obj->reH[i]+= n1_iMPi * ( a*sin_jOMt -  b*sin_jOMt_1  +  c_StjOM*(cos_jOMt-cos_jOMt_1) ); 
        obj->imH[i]+= n1_iMPi * ( a*cos_jOMt -  b*cos_jOMt_1  -  c_StjOM*(sin_jOMt-sin_jOMt_1) );
        absri =  sqrt( (obj->reH[i]*obj->reH[i])  + (obj->imH[i]*obj->imH[i]) );
        if(!i) firstA = 100.0/absri;
        obj->armonics[i] = absri*firstA;            
    }
    obj->t_1 = obj->t;
}
/*================================================================================================================================================================*/
void FFT_KairosRecursive_Flush(RFFT_KairosConfig_t *obj){
    obj->t = obj->t_1 = 0.0;
    memset(obj->reH, 0, obj->NA);
    memset(obj->imH, 0, obj->NA);
    memset(obj->Buffer, 0, obj->WindowLength);
}
/*================================================================================================================================================================*/
void FFT_KairosRecursive_End(RFFT_KairosConfig_t *obj){
    //free(obj->Buffer);
    memset(obj, 0, sizeof(RFFT_KairosConfig_t));
}
/*================================================================================================================================================================*/
inline real_t _sin(real_t x) { //super-optimized version of sin(x) [accurate aprox] LITTLE-ENDIAN
    int k;
    real_t y, z;
    z = x*M_1_PI + 6755399441055744.0;
    k  = *((int *) &z);
    z = k*M_PI;
    x -= z;
    y  = x*x;
    z = (0.0073524681968701*y - 0.1652891139701474)*y + 0.9996919862959676;
    x *= z;
    k &= 1;
    k += k;
    return  (x - k*x);
}
/*================================================================================================================================================================*/
Kairos_DampingPhaseState_t DecisionLogic_Kairos(real_t *Af1Phase_a, real_t *Af2Phase_a, real_t *Af1Phase_b, real_t *Af2Phase_b, real_t *Af1Phase_c, real_t *Af2Phase_c, Kairos_PhaseThresholds_t *obj, const int nCheck){
    int cntAf1=0, cntBf1=0, cntCf1=0;
    int cntAf2=0, cntBf2=0, cntCf2=0;
    int i;
    Kairos_DampingPhaseState_t out;
    out.reg = 0;
    for(i=0;i<nCheck;i++){
        if(Af1Phase_a[i]>obj->f1.Min[i] && Af1Phase_a[i]<obj->f1.Max[i]) cntAf1++;
        if(Af2Phase_a[i]>obj->f2.Min[i] && Af2Phase_a[i]<obj->f2.Max[i]) cntAf2++;
        
        if(Af1Phase_b[i]>obj->f1.Min[i] && Af1Phase_b[i]<obj->f1.Max[i]) cntBf1++;
        if(Af2Phase_b[i]>obj->f2.Min[i] && Af2Phase_b[i]<obj->f2.Max[i]) cntBf2++;

        if(Af1Phase_c[i]>obj->f1.Min[i] && Af1Phase_c[i]<obj->f1.Max[i]) cntCf1++;
        if(Af2Phase_c[i]>obj->f2.Min[i] && Af2Phase_c[i]<obj->f2.Max[i]) cntCf2++;        
    }
    
    out.X = cntAf1>=nCheck | cntBf1>=nCheck | cntCf1>=nCheck;
    out.Y = cntAf2>=nCheck | cntBf2>=nCheck | cntCf2>=nCheck;
       
    out.A  = (cntAf1>=nCheck | cntAf2>=nCheck);
    out.B  = (cntBf1>=nCheck | cntBf2>=nCheck);
    out.C  = (cntCf1>=nCheck | cntCf2>=nCheck);
    return out;
}
/*================================================================================================================================================================*/
real_t TuneFrequency(real_t SampleTime, real_t BaseFrequency, real_t TrueSampleTime){
    int N1, N2;
    real_t freq = BaseFrequency;
    N1 = (int)(1.0/(BaseFrequency*TrueSampleTime))+2;
    N2 = (int)(1.0/(freq*SampleTime))+2;
    while(N1 != N2){
        freq = freq + 0.01;
        //printf("\r\n %g",freq);
        N2 = (int)(1.0/(freq*SampleTime))+2;
        //printf(" N1=%d N2=%d \r\n", N1, N2);
    }
    return freq;    
}
/*================================================================================================================================================================*/