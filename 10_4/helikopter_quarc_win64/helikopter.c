/*
 * helikopter.c
 *
 * Code generation for model "helikopter".
 *
 * Model version              : 1.198
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Thu Mar 01 17:33:40 2018
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
  real_T rtb_optimalinput[2];
  real_T rtb_Sum3_o[6];
  real_T rtb_Backgain;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T rtb_TmpSignalConversionAtToFile[6];
  real_T *lastU;
  real_T rtb_Sum4[2];
  real_T rtb_Derivative;
  int32_T i;
  real_T tmp[6];
  int32_T i_0;
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

  /* FromWorkspace: '<Root>/x_optimal' */
  {
    real_T *pDataValues = (real_T *) helikopter_DW.x_optimal_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *) helikopter_DW.x_optimal_PWORK.TimePtr;
    int_T currTimeIndex = helikopter_DW.x_optimal_IWORK.PrevIndex;
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

    helikopter_DW.x_optimal_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 6; ++elIdx) {
              (&rtb_Sum3_o[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 81;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 6; ++elIdx) {
              (&rtb_Sum3_o[0])[elIdx] = pDataValues[currTimeIndex + 1];
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
          for (elIdx = 0; elIdx < 6; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&rtb_Sum3_o[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1, f2);
            pDataValues += 81;
          }
        }
      }
    }
  }

  /* Sum: '<Root>/Sum3' */
  tmp[0] = helikopter_B.Sum5;
  for (i = 0; i < 5; i++) {
    tmp[i + 1] = helikopter_B.Gain1[1 + i];
  }

  for (i = 0; i < 6; i++) {
    rtb_Sum3_o[i] = tmp[i] - rtb_Sum3_o[i];
  }

  /* End of Sum: '<Root>/Sum3' */

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
              (&rtb_optimalinput[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 81;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 2; ++elIdx) {
              (&rtb_optimalinput[0])[elIdx] = pDataValues[currTimeIndex + 1];
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
            (&rtb_optimalinput[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1,
              f2);
            pDataValues += 81;
          }
        }
      }
    }
  }

  /* Sum: '<Root>/Sum4' incorporates:
   *  Gain: '<Root>/Gain'
   */
  for (i = 0; i < 2; i++) {
    rtb_Derivative = 0.0;
    for (i_0 = 0; i_0 < 6; i_0++) {
      rtb_Derivative += helikopter_P.K_LQR[(i_0 << 1) + i] * rtb_Sum3_o[i_0];
    }

    rtb_Sum4[i] = rtb_optimalinput[i] - rtb_Derivative;
  }

  /* End of Sum: '<Root>/Sum4' */

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

  /* Start for FromWorkspace: '<Root>/x_optimal' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0 } ;

    static real_T pDataValues0[] = { 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1378421413625261, 3.1262155534579983,
      3.1033093000299643, 3.0666274151911783, 3.0144539223941584,
      2.9456562771175667, 2.8595077632935446, 2.7555515879651526,
      2.633505110490284, 2.4931956060320961, 2.334518576064299,
      2.1595167327017788, 1.9722003072130425, 1.7776737868169057,
      1.581278110675107, 1.3879726346119721, 1.2019716210539835,
      1.0265825031459077, 0.86418301570973988, 0.716283826635781,
      0.58363727410735511, 0.46636537815322782, 0.36409007523568188,
      0.27605559646428313, 0.20123766867525986, 0.1384372861860389,
      0.0863587389123823, 0.0436727667361332, 0.0090662899949057532,
      -0.018719519150506729, -0.04086102576863871, -0.058432058874766213,
      -0.072390451694978916, -0.083569984518432866, -0.092676830949857972,
      -0.10029008251402773, -0.10686635097076151, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.015002048909068425,
      -0.046506351618112118, -0.091625013712135411, -0.14672753935514374,
      -0.20869397118807928, -0.27519058110636674, -0.34459405529608839,
      -0.41582470131356869, -0.48818590989947463, -0.56123801783275185,
      -0.63470811987118858, -0.70000737345008046, -0.74926570195494546,
      -0.77810608158454686, -0.78558270456719481, -0.77322190425254,
      -0.74400405423195481, -0.70155647163230273, -0.64959794974467133,
      -0.59159675629583586, -0.53058621011370344, -0.46908758381650933,
      -0.40910121167018376, -0.35213791508559517, -0.29927171115609308,
      -0.25120152995688383, -0.20831418909462643, -0.1707438887049964,
      -0.13842590696490978, -0.11114323658164992, -0.088566026472527926,
      -0.07028413242451001, -0.055833571280850813, -0.044718131293815785,
      -0.03642738572570043, -0.030453006256678981, -0.02630507382693514,
      -0.023530284097554005, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.10602875205865553, 0.22266037932317662, 0.31888147181640647,
      0.3894436063114417, 0.43795507377677845, 0.46997264230390068,
      0.49051724877547076, 0.5034310014147434, 0.5114213858602934,
      0.51630439857701838, 0.51925862127063693, 0.46151018499522611,
      0.34813905297356113, 0.20383278841172603, 0.052841915745523889,
      -0.087361415747475551, -0.20650060497049738, -0.30000330209743947,
      -0.36722298853643576, -0.40993028329452247, -0.43119923907181512,
      -0.4346488029159491, -0.42396076814256589, -0.40259482465492485,
      -0.37363813854354172, -0.339741681598217, -0.3031113455436471,
      -0.26553253418426931, -0.22841115195168241, -0.19282349438389573,
      -0.15956709829077051, -0.1292094448517746, -0.10213104715871936,
      -0.078559684583869735, -0.058595823247298472, -0.042224632332873289,
      -0.029316002221770918, -0.019611154050671972, -0.012704237341955495, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.42411500823462212,
      0.46652650905808429, 0.38488436997291953, 0.282248537980141,
      0.19404586986134692, 0.1280702741084889, 0.082178425886280368,
      0.051655010557090465, 0.0319615377821998, 0.019532050866899828,
      0.011816890774474422, -0.23099374510164331, -0.45348452808665979,
      -0.57722505824734027, -0.60396349066480859, -0.56081332597199773,
      -0.47655675689208726, -0.37401078850776853, -0.26887874575598519,
      -0.17082917903234685, -0.08507582310917057, -0.013798255376535827,
      0.042752139093532683, 0.085463773950564412, 0.11582674444553241,
      0.13558582778129902, 0.14652134421827959, 0.15031524543751129,
      0.14848552893034761, 0.14235063027114675, 0.13302558437250084,
      0.12143061375598369, 0.10831359077222094, 0.094285450299398513,
      0.079855445346285039, 0.065484763657700715, 0.051634520444409478,
      0.038819392684395777, 0.027627666834865915, 0.01864642292532577, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0035450632474740718,
      0.010347403542814175, 0.020099773004038656, 0.032469386227554721,
      0.047088804302346764, 0.063543791804318558, 0.08135908521455848,
      0.099982397734347875, 0.11876718822738686, 0.13695514797417116,
      0.15365885920783459, 0.16784586149636782, 0.17832498227089058,
      0.18373616300200904, 0.18254550511092404, 0.17625220797667304,
      0.16617708590640209, 0.15346311861611667, 0.13908124060519769,
      0.12383978909895761, 0.10839681176843041, 0.093273628233094014,
      0.0788697913805047, 0.065477745775895241, 0.05329748604638674,
      0.042450637528193683, 0.032993233778697104, 0.024927624364284081,
      0.018212816531987138, 0.012773783845570388, 0.00850920012853432,
      0.0052979494230561662, 0.003004429002196992, 0.0014828156476273632,
      0.00058031071467695347, 0.00013947166240217007, -2.141606652279273E-22,
      2.3435266578412981E-23, 1.9843329920874446E-16, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.014180252989896287, 0.027209361181360418,
      0.03900947784489793, 0.049478452894064252, 0.058477672299168178,
      0.065819950007887151, 0.071261173640959671, 0.074493250079157566,
      0.075139161972155935, 0.072751838987137218, 0.066814844934653742,
      0.056748009154132889, 0.041916483098091033, 0.02164472292447378,
      -0.00476263156433996, -0.025173188537003982, -0.040300488281083874,
      -0.050855869161141626, -0.057527512043675927, -0.060965806024960413,
      -0.061771909322108777, -0.060492734141345587, -0.057615347410357239,
      -0.053568182418437887, -0.048721038918033976, -0.043387394072772215,
      -0.037829614997986323, -0.0322624376576521, -0.026859231329187772,
      -0.021756130745666995, -0.017058334868144277, -0.012845002821912617,
      -0.0091740816834366967, -0.0060864534182785163, -0.0036100197318016387,
      -0.0017633562090991336, -0.00055788664960868027, 1.1782127005121038E-21,
      7.9373276830372174E-16, 5.95789800716932E-16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helikopter_DW.x_optimal_PWORK.TimePtr = (void *) pTimeValues0;
    helikopter_DW.x_optimal_PWORK.DataPtr = (void *) pDataValues0;
    helikopter_DW.x_optimal_IWORK.PrevIndex = 0;
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
      0.52359877559829882, 0.22966806008431381, -0.012081074502814175,
      -0.19512759949696437, -0.326232533248753, -0.41396594850866508,
      -0.46626354634692968, -0.49026896243705548, -0.49227550238433476,
      -0.47770861624157313, -0.45116484681536895, -0.41646197843726168,
      -0.37671533211973329, -0.33441601274161431, -0.29151587054938571,
      -0.24950451189351031, -0.20948658121002364, -0.17225635935663269,
      -0.13835331543673796, -0.10812256211252616, -0.081756187327748878,
      -0.059332068064471513, -0.0408380496529649, -0.026179702602174091,
      -0.015184256644712976, -0.0075729943084719645, -0.0029340404341634033,
      -0.00067408833292906072, -2.8917680013454423E-7, -1.7019002532349467E-6,
      -1.7019002532349467E-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.056721011959585156, 0.066296685755752829, 0.076078734533552114,
      0.085756617247872674, 0.094937685190227411, 0.10313237111193911,
      0.10975280471943587, 0.11410394138163046, 0.11537825161369707,
      0.1126585065039433, 0.10491580061963469, 0.091021771474257065,
      0.0697597678801153, 0.039846109508065955, -3.47803001773024E-5,
      9.2472665949756757E-5, 0.00025441950812684456, 0.00045210456924198552,
      0.00068861821235043617, 0.00096502633354263515, 0.00128494944058539,
      0.0016447536536840278, 0.0020466964766689259, 0.0024836667716690583,
      0.0029498626023076496, 0.0034384484744345014, 0.0039345387343626779,
      0.0044235483381617688, 0.0048825869427918467, 0.0052883295517736923,
      0.0056090863933285764, 0.0058084977117739628, 0.0058445513309626166,
      0.0056705429517424641, 0.00523367405789404, 0.0044746983408144947,
      0.0033317144576182922, 0.0017393189160492164, 3.1749274270027468E-15, 0.0,
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
  helikopter_M->Sizes.checksums[0] = (3617501102U);
  helikopter_M->Sizes.checksums[1] = (5272714U);
  helikopter_M->Sizes.checksums[2] = (3821682694U);
  helikopter_M->Sizes.checksums[3] = (1259314972U);

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
  helikopter_M->Sizes.numBlocks = (63);/* Number of blocks */
  helikopter_M->Sizes.numBlockIO = (18);/* Number of block outputs */
  helikopter_M->Sizes.numBlockPrms = (153);/* Sum of parameter "widths" */
  return helikopter_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
