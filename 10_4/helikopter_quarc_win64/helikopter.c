/*
 * helikopter.c
 *
 * Code generation for model "helikopter".
 *
 * Model version              : 1.198
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Thu Mar 08 13:18:19 2018
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "helikopter.h"
#include "helikopter_private.h"
#include "helikopter_dt.h"

/* Block signals (auto storage) */
B_helikopter_T helikopter_B;

/* Continuous states */
X_helikopter_T helikopter_X;

/* Block states (auto storage) */
DW_helikopter_T helikopter_DW;

/* Real-time model */
RT_MODEL_helikopter_T helikopter_M_;
RT_MODEL_helikopter_T *const helikopter_M = &helikopter_M_;

/*
 * Writes out MAT-file header.  Returns success or failure.
 * Returns:
 *      0 - success
 *      1 - failure
 */
int_T rt_WriteMat4FileHeader(FILE *fp, int32_T m, int32_T n, const char *name)
{
  typedef enum { ELITTLE_ENDIAN, EBIG_ENDIAN } ByteOrder;

  int16_T one = 1;
  ByteOrder byteOrder = (*((int8_T *)&one)==1) ? ELITTLE_ENDIAN : EBIG_ENDIAN;
  int32_T type = (byteOrder == ELITTLE_ENDIAN) ? 0: 1000;
  int32_T imagf = 0;
  int32_T name_len = (int32_T)strlen(name) + 1;
  if ((fwrite(&type, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&m, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&n, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&imagf, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&name_len, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(name, sizeof(char), name_len, fp) == 0)) {
    return(1);
  } else {
    return(0);
  }
}                                      /* end rt_WriteMat4FileHeader */

/*
 * This function updates continuous states using the ODE1 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE1_IntgData *id = (ODE1_IntgData *)rtsiGetSolverData(si);
  real_T *f0 = id->f[0];
  int_T i;
  int_T nXc = 4;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);
  rtsiSetdX(si, f0);
  helikopter_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helikopter_output(void)
{
  /* local block i/o variables */
  real_T rtb_Sum4[2];
  real_T rtb_Backgain;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T rtb_TmpSignalConversionAtToFile[6];
  real_T *lastU;
  real_T rtb_Derivative;
  int32_T i;
  if (rtmIsMajorTimeStep(helikopter_M)) {
    /* set solver stop time */
    if (!(helikopter_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helikopter_M->solverInfo,
                            ((helikopter_M->Timing.clockTickH0 + 1) *
        helikopter_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helikopter_M->solverInfo,
                            ((helikopter_M->Timing.clockTick0 + 1) *
        helikopter_M->Timing.stepSize0 + helikopter_M->Timing.clockTickH0 *
        helikopter_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helikopter_M)) {
    helikopter_M->Timing.t[0] = rtsiGetT(&helikopter_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helikopter_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

    /* S-Function Block: helikopter/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder(helikopter_DW.HILReadEncoderTimebase_Task,
        1, &helikopter_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          helikopter_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          helikopter_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          helikopter_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }

    /* Gain: '<S4>/Elevation: Count to rad' */
    helikopter_B.ElevationCounttorad = helikopter_P.ElevationCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o3;

    /* Gain: '<S6>/Gain' */
    helikopter_B.Gain = helikopter_P.Gain_Gain *
      helikopter_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    helikopter_B.Sum = helikopter_B.Gain +
      helikopter_P.elavation_offsetdeg_Value;

    /* Gain: '<S4>/Pitch: Count to rad' */
    helikopter_B.PitchCounttorad = helikopter_P.PitchCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S8>/Gain' */
    helikopter_B.Gain_i = helikopter_P.Gain_Gain_a *
      helikopter_B.PitchCounttorad;

    /* Gain: '<S4>/Travel: Count to rad' */
    helikopter_B.TravelCounttorad = helikopter_P.TravelCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o1;

    /* Gain: '<S11>/Gain' */
    helikopter_B.Gain_p = helikopter_P.Gain_Gain_ar *
      helikopter_B.TravelCounttorad;
  }

  /* TransferFcn: '<S4>/Travel: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helikopter_P.TravelTransferFcn_C *
    helikopter_X.TravelTransferFcn_CSTATE;
  rtb_Backgain += helikopter_P.TravelTransferFcn_D *
    helikopter_B.TravelCounttorad;

  /* Gain: '<S12>/Gain' */
  helikopter_B.Gain_d = helikopter_P.Gain_Gain_l * rtb_Backgain;

  /* TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helikopter_P.PitchTransferFcn_C *
    helikopter_X.PitchTransferFcn_CSTATE;
  rtb_Backgain += helikopter_P.PitchTransferFcn_D * helikopter_B.PitchCounttorad;

  /* Gain: '<S9>/Gain' */
  helikopter_B.Gain_b = helikopter_P.Gain_Gain_ae * rtb_Backgain;

  /* TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helikopter_P.ElevationTransferFcn_C *
    helikopter_X.ElevationTransferFcn_CSTATE;
  rtb_Backgain += helikopter_P.ElevationTransferFcn_D *
    helikopter_B.ElevationCounttorad;

  /* Gain: '<S7>/Gain' */
  helikopter_B.Gain_dg = helikopter_P.Gain_Gain_n * rtb_Backgain;

  /* Gain: '<S2>/Gain1' */
  helikopter_B.Gain1[0] = helikopter_P.Gain1_Gain * helikopter_B.Gain_p;
  helikopter_B.Gain1[1] = helikopter_P.Gain1_Gain * helikopter_B.Gain_d;
  helikopter_B.Gain1[2] = helikopter_P.Gain1_Gain * helikopter_B.Gain_i;
  helikopter_B.Gain1[3] = helikopter_P.Gain1_Gain * helikopter_B.Gain_b;
  helikopter_B.Gain1[4] = helikopter_P.Gain1_Gain * helikopter_B.Sum;
  helikopter_B.Gain1[5] = helikopter_P.Gain1_Gain * helikopter_B.Gain_dg;

  /* Sum: '<Root>/Sum5' incorporates:
   *  Constant: '<Root>/Constant1'
   */
  helikopter_B.Sum5 = helikopter_P.Constant1_Value + helikopter_B.Gain1[0];
  if (rtmIsMajorTimeStep(helikopter_M)) {
    /* SignalConversion: '<Root>/TmpSignal ConversionAtTo FileInport1' */
    rtb_TmpSignalConversionAtToFile[0] = helikopter_B.Sum5;
    for (i = 0; i < 5; i++) {
      rtb_TmpSignalConversionAtToFile[i + 1] = helikopter_B.Gain1[i + 1];
    }

    /* End of SignalConversion: '<Root>/TmpSignal ConversionAtTo FileInport1' */

    /* ToFile: '<Root>/To File' */
    {
      if (!(++helikopter_DW.ToFile_IWORK.Decimation % 1) &&
          (helikopter_DW.ToFile_IWORK.Count*7)+1 < 100000000 ) {
        FILE *fp = (FILE *) helikopter_DW.ToFile_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[7];
          helikopter_DW.ToFile_IWORK.Decimation = 0;
          u[0] = helikopter_M->Timing.t[1];
          u[1] = rtb_TmpSignalConversionAtToFile[0];
          u[2] = rtb_TmpSignalConversionAtToFile[1];
          u[3] = rtb_TmpSignalConversionAtToFile[2];
          u[4] = rtb_TmpSignalConversionAtToFile[3];
          u[5] = rtb_TmpSignalConversionAtToFile[4];
          u[6] = rtb_TmpSignalConversionAtToFile[5];
          if (fwrite(u, sizeof(real_T), 7, fp) != 7) {
            rtmSetErrorStatus(helikopter_M, "Error writing to MAT-file Data.mat");
            return;
          }

          if (((++helikopter_DW.ToFile_IWORK.Count)*7)+1 >= 100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file Data.mat.\n");
          }
        }
      }
    }
  }

  /* FromWorkspace: '<Root>/optimal input' */
  {
    real_T *pDataValues = (real_T *) helikopter_DW.optimalinput_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *) helikopter_DW.optimalinput_PWORK.TimePtr;
    int_T currTimeIndex = helikopter_DW.optimalinput_IWORK.PrevIndex;
    real_T t = helikopter_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[80]) {
      currTimeIndex = 79;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    helikopter_DW.optimalinput_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 2; ++elIdx) {
              (&rtb_Sum4[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 81;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 2; ++elIdx) {
              (&rtb_Sum4[0])[elIdx] = pDataValues[currTimeIndex + 1];
              pDataValues += 81;
            }
          }
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;

        {
          int_T elIdx;
          for (elIdx = 0; elIdx < 2; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&rtb_Sum4[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1, f2);
            pDataValues += 81;
          }
        }
      }
    }
  }

  /* Sum: '<Root>/Sum4' */
  rtb_Sum4[0] -= 0.0;
  rtb_Sum4[1] -= 0.0;

  /* Sum: '<Root>/Sum1' incorporates:
   *  Constant: '<Root>/Vd_bias'
   *  Gain: '<S5>/K_pd'
   *  Gain: '<S5>/K_pp'
   *  Sum: '<S5>/Sum2'
   *  Sum: '<S5>/Sum3'
   */
  helikopter_B.Sum1 = ((rtb_Sum4[0] - helikopter_B.Gain1[2]) * helikopter_P.K_pp
                       - helikopter_P.K_pd * helikopter_B.Gain1[3]) +
    helikopter_P.Vd_ff;
  if (rtmIsMajorTimeStep(helikopter_M)) {
  }

  /* Integrator: '<S3>/Integrator'
   *
   * Regarding '<S3>/Integrator':
   *  Limited Integrator
   */
  if (helikopter_X.Integrator_CSTATE >= helikopter_P.Integrator_UpperSat ) {
    helikopter_X.Integrator_CSTATE = helikopter_P.Integrator_UpperSat;
  } else if (helikopter_X.Integrator_CSTATE <= (helikopter_P.Integrator_LowerSat)
             ) {
    helikopter_X.Integrator_CSTATE = (helikopter_P.Integrator_LowerSat);
  }

  rtb_Backgain = helikopter_X.Integrator_CSTATE;

  /* Sum: '<S3>/Sum' */
  rtb_Derivative = rtb_Sum4[1] - helikopter_B.Gain1[4];

  /* Sum: '<Root>/Sum2' incorporates:
   *  Constant: '<Root>/Vs_bias'
   *  Gain: '<S3>/K_ed'
   *  Gain: '<S3>/K_ep'
   *  Sum: '<S3>/Sum1'
   */
  helikopter_B.Sum2 = ((helikopter_P.K_ep * rtb_Derivative + rtb_Backgain) -
                       helikopter_P.K_ed * helikopter_B.Gain1[5]) +
    helikopter_P.Vs_ff;
  if (rtmIsMajorTimeStep(helikopter_M)) {
  }

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  rtb_Backgain = (helikopter_B.Sum2 - helikopter_B.Sum1) *
    helikopter_P.Backgain_Gain;

  /* Gain: '<S3>/K_ei' */
  helikopter_B.K_ei = helikopter_P.K_ei * rtb_Derivative;
  if (rtmIsMajorTimeStep(helikopter_M)) {
  }

  /* Derivative: '<S4>/Derivative' */
  if ((helikopter_DW.TimeStampA >= helikopter_M->Timing.t[0]) &&
      (helikopter_DW.TimeStampB >= helikopter_M->Timing.t[0])) {
    rtb_Derivative = 0.0;
  } else {
    rtb_Derivative = helikopter_DW.TimeStampA;
    lastU = &helikopter_DW.LastUAtTimeA;
    if (helikopter_DW.TimeStampA < helikopter_DW.TimeStampB) {
      if (helikopter_DW.TimeStampB < helikopter_M->Timing.t[0]) {
        rtb_Derivative = helikopter_DW.TimeStampB;
        lastU = &helikopter_DW.LastUAtTimeB;
      }
    } else {
      if (helikopter_DW.TimeStampA >= helikopter_M->Timing.t[0]) {
        rtb_Derivative = helikopter_DW.TimeStampB;
        lastU = &helikopter_DW.LastUAtTimeB;
      }
    }

    rtb_Derivative = (helikopter_B.PitchCounttorad - *lastU) /
      (helikopter_M->Timing.t[0] - rtb_Derivative);
  }

  /* End of Derivative: '<S4>/Derivative' */

  /* Gain: '<S10>/Gain' */
  helikopter_B.Gain_l = helikopter_P.Gain_Gain_a1 * rtb_Derivative;
  if (rtmIsMajorTimeStep(helikopter_M)) {
  }

  /* Saturate: '<S4>/Back motor: Saturation' */
  if (rtb_Backgain > helikopter_P.BackmotorSaturation_UpperSat) {
    helikopter_B.BackmotorSaturation = helikopter_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Backgain < helikopter_P.BackmotorSaturation_LowerSat) {
    helikopter_B.BackmotorSaturation = helikopter_P.BackmotorSaturation_LowerSat;
  } else {
    helikopter_B.BackmotorSaturation = rtb_Backgain;
  }

  /* End of Saturate: '<S4>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(helikopter_M)) {
  }

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  rtb_Derivative = (helikopter_B.Sum1 + helikopter_B.Sum2) *
    helikopter_P.Frontgain_Gain;

  /* Saturate: '<S4>/Front motor: Saturation' */
  if (rtb_Derivative > helikopter_P.FrontmotorSaturation_UpperSat) {
    helikopter_B.FrontmotorSaturation =
      helikopter_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Derivative < helikopter_P.FrontmotorSaturation_LowerSat) {
    helikopter_B.FrontmotorSaturation =
      helikopter_P.FrontmotorSaturation_LowerSat;
  } else {
    helikopter_B.FrontmotorSaturation = rtb_Derivative;
  }

  /* End of Saturate: '<S4>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(helikopter_M)) {
    /* S-Function (hil_write_analog_block): '<S4>/HIL Write Analog' */

    /* S-Function Block: helikopter/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helikopter_DW.HILWriteAnalog_Buffer[0] = helikopter_B.FrontmotorSaturation;
      helikopter_DW.HILWriteAnalog_Buffer[1] = helikopter_B.BackmotorSaturation;
      result = hil_write_analog(helikopter_DW.HILInitialize_Card,
        helikopter_P.HILWriteAnalog_channels, 2,
        &helikopter_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_M, _rt_error_message);
      }
    }
  }
}

/* Model update function */
void helikopter_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S4>/Derivative' */
  if (helikopter_DW.TimeStampA == (rtInf)) {
    helikopter_DW.TimeStampA = helikopter_M->Timing.t[0];
    lastU = &helikopter_DW.LastUAtTimeA;
  } else if (helikopter_DW.TimeStampB == (rtInf)) {
    helikopter_DW.TimeStampB = helikopter_M->Timing.t[0];
    lastU = &helikopter_DW.LastUAtTimeB;
  } else if (helikopter_DW.TimeStampA < helikopter_DW.TimeStampB) {
    helikopter_DW.TimeStampA = helikopter_M->Timing.t[0];
    lastU = &helikopter_DW.LastUAtTimeA;
  } else {
    helikopter_DW.TimeStampB = helikopter_M->Timing.t[0];
    lastU = &helikopter_DW.LastUAtTimeB;
  }

  *lastU = helikopter_B.PitchCounttorad;

  /* End of Update for Derivative: '<S4>/Derivative' */
  if (rtmIsMajorTimeStep(helikopter_M)) {
    rt_ertODEUpdateContinuousStates(&helikopter_M->solverInfo);
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++helikopter_M->Timing.clockTick0)) {
    ++helikopter_M->Timing.clockTickH0;
  }

  helikopter_M->Timing.t[0] = rtsiGetSolverStopTime(&helikopter_M->solverInfo);

  {
    /* Update absolute timer for sample time: [0.002s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++helikopter_M->Timing.clockTick1)) {
      ++helikopter_M->Timing.clockTickH1;
    }

    helikopter_M->Timing.t[1] = helikopter_M->Timing.clockTick1 *
      helikopter_M->Timing.stepSize1 + helikopter_M->Timing.clockTickH1 *
      helikopter_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void helikopter_derivatives(void)
{
  XDot_helikopter_T *_rtXdot;
  _rtXdot = ((XDot_helikopter_T *) helikopter_M->ModelData.derivs);

  /* Derivatives for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += helikopter_P.TravelTransferFcn_A *
    helikopter_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += helikopter_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += helikopter_P.PitchTransferFcn_A *
    helikopter_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += helikopter_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE += helikopter_P.ElevationTransferFcn_A *
    helikopter_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE += helikopter_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S3>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helikopter_X.Integrator_CSTATE <= (helikopter_P.Integrator_LowerSat)
            );
    usat = ( helikopter_X.Integrator_CSTATE >= helikopter_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helikopter_B.K_ei > 0)) ||
        (usat && (helikopter_B.K_ei < 0)) ) {
      ((XDot_helikopter_T *) helikopter_M->ModelData.derivs)->Integrator_CSTATE =
        helikopter_B.K_ei;
    } else {
      /* in saturation */
      ((XDot_helikopter_T *) helikopter_M->ModelData.derivs)->Integrator_CSTATE =
        0.0;
    }
  }
}

/* Model initialize function */
void helikopter_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopter/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &helikopter_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopter_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options(helikopter_DW.HILInitialize_Card,
      "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopter_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(helikopter_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopter_M, _rt_error_message);
      return;
    }

    if ((helikopter_P.HILInitialize_set_analog_input_ && !is_switching) ||
        (helikopter_P.HILInitialize_set_analog_inpu_m && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &helikopter_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] = helikopter_P.HILInitialize_analog_input_mini;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &helikopter_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] = helikopter_P.HILInitialize_analog_input_maxi;
        }
      }

      result = hil_set_analog_input_ranges(helikopter_DW.HILInitialize_Card,
        helikopter_P.HILInitialize_analog_input_chan, 8U,
        &helikopter_DW.HILInitialize_AIMinimums[0],
        &helikopter_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter_P.HILInitialize_set_analog_output && !is_switching) ||
        (helikopter_P.HILInitialize_set_analog_outp_b && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &helikopter_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] = helikopter_P.HILInitialize_analog_output_min;
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &helikopter_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] = helikopter_P.HILInitialize_analog_output_max;
        }
      }

      result = hil_set_analog_output_ranges(helikopter_DW.HILInitialize_Card,
        helikopter_P.HILInitialize_analog_output_cha, 8U,
        &helikopter_DW.HILInitialize_AOMinimums[0],
        &helikopter_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter_P.HILInitialize_set_analog_outp_e && !is_switching) ||
        (helikopter_P.HILInitialize_set_analog_outp_j && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helikopter_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helikopter_P.HILInitialize_initial_analog_ou;
        }
      }

      result = hil_write_analog(helikopter_DW.HILInitialize_Card,
        helikopter_P.HILInitialize_analog_output_cha, 8U,
        &helikopter_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_M, _rt_error_message);
        return;
      }
    }

    if (helikopter_P.HILInitialize_set_analog_outp_p) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helikopter_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helikopter_P.HILInitialize_watchdog_analog_o;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (helikopter_DW.HILInitialize_Card,
         helikopter_P.HILInitialize_analog_output_cha, 8U,
         &helikopter_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter_P.HILInitialize_set_encoder_param && !is_switching) ||
        (helikopter_P.HILInitialize_set_encoder_par_m && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &helikopter_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = helikopter_P.HILInitialize_quadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode(helikopter_DW.HILInitialize_Card,
        helikopter_P.HILInitialize_encoder_channels, 8U,
        (t_encoder_quadrature_mode *)
        &helikopter_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter_P.HILInitialize_set_encoder_count && !is_switching) ||
        (helikopter_P.HILInitialize_set_encoder_cou_k && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &helikopter_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] = helikopter_P.HILInitialize_initial_encoder_c;
        }
      }

      result = hil_set_encoder_counts(helikopter_DW.HILInitialize_Card,
        helikopter_P.HILInitialize_encoder_channels, 8U,
        &helikopter_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter_P.HILInitialize_set_pwm_params_at && !is_switching) ||
        (helikopter_P.HILInitialize_set_pwm_params__f && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues = &helikopter_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helikopter_P.HILInitialize_pwm_modes;
        }
      }

      result = hil_set_pwm_mode(helikopter_DW.HILInitialize_Card,
        helikopter_P.HILInitialize_pwm_channels, 8U, (t_pwm_mode *)
        &helikopter_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_pwm_channels =
          helikopter_P.HILInitialize_pwm_channels;
        int32_T *dw_POModeValues = &helikopter_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE) {
            helikopter_DW.HILInitialize_POSortedChans[num_duty_cycle_modes] =
              p_HILInitialize_pwm_channels[i1];
            helikopter_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes] =
              helikopter_P.HILInitialize_pwm_frequency;
            num_duty_cycle_modes++;
          } else {
            helikopter_DW.HILInitialize_POSortedChans[7U - num_frequency_modes] =
              p_HILInitialize_pwm_channels[i1];
            helikopter_DW.HILInitialize_POSortedFreqs[7U - num_frequency_modes] =
              helikopter_P.HILInitialize_pwm_frequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(helikopter_DW.HILInitialize_Card,
          &helikopter_DW.HILInitialize_POSortedChans[0], num_duty_cycle_modes,
          &helikopter_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helikopter_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(helikopter_DW.HILInitialize_Card,
          &helikopter_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &helikopter_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helikopter_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues = &helikopter_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helikopter_P.HILInitialize_pwm_configuration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues = &helikopter_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = helikopter_P.HILInitialize_pwm_alignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &helikopter_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = helikopter_P.HILInitialize_pwm_polarity;
        }
      }

      result = hil_set_pwm_configuration(helikopter_DW.HILInitialize_Card,
        helikopter_P.HILInitialize_pwm_channels, 8U,
        (t_pwm_configuration *) &helikopter_DW.HILInitialize_POModeValues[0],
        (t_pwm_alignment *) &helikopter_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &helikopter_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs = &helikopter_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] = helikopter_P.HILInitialize_pwm_leading_deadb;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &helikopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helikopter_P.HILInitialize_pwm_trailing_dead;
        }
      }

      result = hil_set_pwm_deadband(helikopter_DW.HILInitialize_Card,
        helikopter_P.HILInitialize_pwm_channels, 8U,
        &helikopter_DW.HILInitialize_POSortedFreqs[0],
        &helikopter_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter_P.HILInitialize_set_pwm_outputs_a && !is_switching) ||
        (helikopter_P.HILInitialize_set_pwm_outputs_g && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helikopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helikopter_P.HILInitialize_initial_pwm_outpu;
        }
      }

      result = hil_write_pwm(helikopter_DW.HILInitialize_Card,
        helikopter_P.HILInitialize_pwm_channels, 8U,
        &helikopter_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_M, _rt_error_message);
        return;
      }
    }

    if (helikopter_P.HILInitialize_set_pwm_outputs_o) {
      {
        int_T i1;
        real_T *dw_POValues = &helikopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helikopter_P.HILInitialize_watchdog_pwm_outp;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (helikopter_DW.HILInitialize_Card,
         helikopter_P.HILInitialize_pwm_channels, 8U,
         &helikopter_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

  /* S-Function Block: helikopter/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader(helikopter_DW.HILInitialize_Card,
      helikopter_P.HILReadEncoderTimebase_samples_,
      helikopter_P.HILReadEncoderTimebase_channels, 3,
      &helikopter_DW.HILReadEncoderTimebase_Task);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopter_M, _rt_error_message);
    }
  }

  /* Start for ToFile: '<Root>/To File' */
  {
    char fileName[509] = "Data.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helikopter_M, "Error creating .mat file Data.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,7,0,"data")) {
      rtmSetErrorStatus(helikopter_M,
                        "Error writing mat file header to file Data.mat");
      return;
    }

    helikopter_DW.ToFile_IWORK.Count = 0;
    helikopter_DW.ToFile_IWORK.Decimation = -1;
    helikopter_DW.ToFile_PWORK.FilePtr = fp;
  }

  /* Start for FromWorkspace: '<Root>/optimal input' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.36413097631047159, 0.11143373894461261,
      -0.10676846607775466, -0.27223145316445352, -0.39175280559153636,
      -0.47193556578045476, -0.51906282030961814, -0.52359877559829882,
      -0.52359877559829882, -0.51834176496353113, -0.48689170426237471,
      -0.44654783884500204, -0.40051832597848441, -0.35151044016300975,
      -0.30175382422624336, -0.25305699429645295, -0.20684649438962749,
      -0.16421859828264262, -0.12597879018700969, -0.0926848533874275,
      -0.064674673413825143, -0.04208362922741074, -0.024848719922293306,
      -0.012694176936080418, -0.005096808064546589, -0.0012423498188599217,
      9.713834819321761E-8, 4.8333836730330497E-8, -4.2855744028049967E-7,
      -4.2855744028049967E-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.285523935767025, 0.30300414148504223, 0.31965353732631946,
      0.33487458946708876, 0.34792148358856795, 0.35786287360314778,
      0.36354331040947785, 0.36354225106439253, 0.35612369421219825,
      0.33916472097483097, 0.31009349644316275, 0.26579435981327482,
      0.20251033768720772, 0.11572034032074495, -9.9448766722938861E-7,
      -6.2935927519132127E-7, -9.3995220413000108E-7, -6.7938544080423865E-7,
      -1.0741310272446335E-6, -3.1497521746137116E-7, -7.0181139797765786E-7,
      -1.0415079990316022E-6, -3.7369774717399189E-7, -1.0645258920441089E-6,
      1.338829506002681E-7, -9.1777052191565832E-7, -2.0525472902114879E-7,
      -1.0951212058933846E-6, -9.6927894770118314E-7, -2.2752111696882431E-8,
      -8.8669822594102093E-7, 6.9136524297052164E-8, 1.1193973532717962E-6,
      -3.4601507864477696E-8, -3.6194314343330012E-7, -8.2139973491827665E-7,
      -7.3312230652968257E-7, -7.2901473185744489E-7, 7.4310718590747462E-7, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helikopter_DW.optimalinput_PWORK.TimePtr = (void *) pTimeValues0;
    helikopter_DW.optimalinput_PWORK.DataPtr = (void *) pDataValues0;
    helikopter_DW.optimalinput_IWORK.PrevIndex = 0;
  }

  /* InitializeConditions for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  helikopter_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  helikopter_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  helikopter_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S3>/Integrator' */
  helikopter_X.Integrator_CSTATE = helikopter_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S4>/Derivative' */
  helikopter_DW.TimeStampA = (rtInf);
  helikopter_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void helikopter_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopter/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(helikopter_DW.HILInitialize_Card);
    hil_monitor_stop_all(helikopter_DW.HILInitialize_Card);
    is_switching = false;
    if ((helikopter_P.HILInitialize_set_analog_out_ex && !is_switching) ||
        (helikopter_P.HILInitialize_set_analog_outp_c && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helikopter_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helikopter_P.HILInitialize_final_analog_outp;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((helikopter_P.HILInitialize_set_pwm_output_ap && !is_switching) ||
        (helikopter_P.HILInitialize_set_pwm_outputs_p && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helikopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helikopter_P.HILInitialize_final_pwm_outputs;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(helikopter_DW.HILInitialize_Card
                         , helikopter_P.HILInitialize_analog_output_cha,
                         num_final_analog_outputs
                         , helikopter_P.HILInitialize_pwm_channels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &helikopter_DW.HILInitialize_AOVoltages[0]
                         , &helikopter_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog(helikopter_DW.HILInitialize_Card,
            helikopter_P.HILInitialize_analog_output_cha,
            num_final_analog_outputs, &helikopter_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(helikopter_DW.HILInitialize_Card,
            helikopter_P.HILInitialize_pwm_channels, num_final_pwm_outputs,
            &helikopter_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helikopter_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(helikopter_DW.HILInitialize_Card);
    hil_monitor_delete_all(helikopter_DW.HILInitialize_Card);
    hil_close(helikopter_DW.HILInitialize_Card);
    helikopter_DW.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<Root>/To File' */
  {
    FILE *fp = (FILE *) helikopter_DW.ToFile_PWORK.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "Data.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopter_M, "Error closing MAT-file Data.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helikopter_M, "Error reopening MAT-file Data.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 7, helikopter_DW.ToFile_IWORK.Count, "data"))
      {
        rtmSetErrorStatus(helikopter_M,
                          "Error writing header for data to MAT-file Data.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopter_M, "Error closing MAT-file Data.mat");
        return;
      }

      helikopter_DW.ToFile_PWORK.FilePtr = (NULL);
    }
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/

/* Solver interface called by GRT_Main */
#ifndef USE_GENERATED_SOLVER

void rt_ODECreateIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEDestroyIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEUpdateContinuousStates(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

#endif

void MdlOutputs(int_T tid)
{
  helikopter_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  helikopter_update();
  UNUSED_PARAMETER(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  helikopter_initialize();
}

void MdlTerminate(void)
{
  helikopter_terminate();
}

/* Registration function */
RT_MODEL_helikopter_T *helikopter(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helikopter_P.Integrator_UpperSat = rtInf;
  helikopter_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helikopter_M, 0,
                sizeof(RT_MODEL_helikopter_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helikopter_M->solverInfo,
                          &helikopter_M->Timing.simTimeStep);
    rtsiSetTPtr(&helikopter_M->solverInfo, &rtmGetTPtr(helikopter_M));
    rtsiSetStepSizePtr(&helikopter_M->solverInfo,
                       &helikopter_M->Timing.stepSize0);
    rtsiSetdXPtr(&helikopter_M->solverInfo, &helikopter_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helikopter_M->solverInfo, (real_T **)
                         &helikopter_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helikopter_M->solverInfo,
      &helikopter_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helikopter_M->solverInfo, (&rtmGetErrorStatus
      (helikopter_M)));
    rtsiSetRTModelPtr(&helikopter_M->solverInfo, helikopter_M);
  }

  rtsiSetSimTimeStep(&helikopter_M->solverInfo, MAJOR_TIME_STEP);
  helikopter_M->ModelData.intgData.f[0] = helikopter_M->ModelData.odeF[0];
  helikopter_M->ModelData.contStates = ((real_T *) &helikopter_X);
  rtsiSetSolverData(&helikopter_M->solverInfo, (void *)
                    &helikopter_M->ModelData.intgData);
  rtsiSetSolverName(&helikopter_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helikopter_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helikopter_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helikopter_M->Timing.sampleTimes = (&helikopter_M->Timing.sampleTimesArray[0]);
    helikopter_M->Timing.offsetTimes = (&helikopter_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helikopter_M->Timing.sampleTimes[0] = (0.0);
    helikopter_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    helikopter_M->Timing.offsetTimes[0] = (0.0);
    helikopter_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helikopter_M, &helikopter_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helikopter_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helikopter_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helikopter_M, -1);
  helikopter_M->Timing.stepSize0 = 0.002;
  helikopter_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  helikopter_M->Sizes.checksums[0] = (535846574U);
  helikopter_M->Sizes.checksums[1] = (1792341678U);
  helikopter_M->Sizes.checksums[2] = (1917164069U);
  helikopter_M->Sizes.checksums[3] = (3798372910U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helikopter_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helikopter_M->extModeInfo,
      &helikopter_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helikopter_M->extModeInfo, helikopter_M->Sizes.checksums);
    rteiSetTPtr(helikopter_M->extModeInfo, rtmGetTPtr(helikopter_M));
  }

  helikopter_M->solverInfoPtr = (&helikopter_M->solverInfo);
  helikopter_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&helikopter_M->solverInfo, 0.002);
  rtsiSetSolverMode(&helikopter_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helikopter_M->ModelData.blockIO = ((void *) &helikopter_B);

  {
    int_T i;
    for (i = 0; i < 6; i++) {
      helikopter_B.Gain1[i] = 0.0;
    }

    helikopter_B.ElevationCounttorad = 0.0;
    helikopter_B.Gain = 0.0;
    helikopter_B.Sum = 0.0;
    helikopter_B.PitchCounttorad = 0.0;
    helikopter_B.Gain_i = 0.0;
    helikopter_B.TravelCounttorad = 0.0;
    helikopter_B.Gain_p = 0.0;
    helikopter_B.Gain_d = 0.0;
    helikopter_B.Gain_b = 0.0;
    helikopter_B.Gain_dg = 0.0;
    helikopter_B.Sum5 = 0.0;
    helikopter_B.Sum1 = 0.0;
    helikopter_B.Sum2 = 0.0;
    helikopter_B.K_ei = 0.0;
    helikopter_B.Gain_l = 0.0;
    helikopter_B.BackmotorSaturation = 0.0;
    helikopter_B.FrontmotorSaturation = 0.0;
  }

  /* parameters */
  helikopter_M->ModelData.defaultParam = ((real_T *)&helikopter_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helikopter_X;
    helikopter_M->ModelData.contStates = (x);
    (void) memset((void *)&helikopter_X, 0,
                  sizeof(X_helikopter_T));
  }

  /* states (dwork) */
  helikopter_M->ModelData.dwork = ((void *) &helikopter_DW);
  (void) memset((void *)&helikopter_DW, 0,
                sizeof(DW_helikopter_T));

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  helikopter_DW.TimeStampA = 0.0;
  helikopter_DW.LastUAtTimeA = 0.0;
  helikopter_DW.TimeStampB = 0.0;
  helikopter_DW.LastUAtTimeB = 0.0;
  helikopter_DW.HILWriteAnalog_Buffer[0] = 0.0;
  helikopter_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helikopter_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.B = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.P = &rtPTransTable;
  }

  /* Initialize Sizes */
  helikopter_M->Sizes.numContStates = (4);/* Number of continuous states */
  helikopter_M->Sizes.numY = (0);      /* Number of model outputs */
  helikopter_M->Sizes.numU = (0);      /* Number of model inputs */
  helikopter_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helikopter_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helikopter_M->Sizes.numBlocks = (60);/* Number of blocks */
  helikopter_M->Sizes.numBlockIO = (18);/* Number of block outputs */
  helikopter_M->Sizes.numBlockPrms = (141);/* Sum of parameter "widths" */
  return helikopter_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
