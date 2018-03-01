/*
 * helikopter.c
 *
 * Code generation for model "helikopter".
 *
 * Model version              : 1.194
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Thu Mar 01 15:43:01 2018
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
  real_T rtb_Sum3_o[4];
  real_T rtb_Backgain;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T rtb_TmpSignalConversionAtToFile[4];
  real_T *lastU;
  real_T rtb_Derivative;
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

  /* Gain: '<S12>/Gain' incorporates:
   *  TransferFcn: '<S4>/Travel: Transfer Fcn'
   */
  helikopter_B.Gain_d = (helikopter_P.TravelTransferFcn_C *
    helikopter_X.TravelTransferFcn_CSTATE + helikopter_P.TravelTransferFcn_D *
    helikopter_B.TravelCounttorad) * helikopter_P.Gain_Gain_l;

  /* Gain: '<S9>/Gain' incorporates:
   *  TransferFcn: '<S4>/Pitch: Transfer Fcn'
   */
  helikopter_B.Gain_b = (helikopter_P.PitchTransferFcn_C *
    helikopter_X.PitchTransferFcn_CSTATE + helikopter_P.PitchTransferFcn_D *
    helikopter_B.PitchCounttorad) * helikopter_P.Gain_Gain_ae;

  /* Gain: '<S7>/Gain' incorporates:
   *  TransferFcn: '<S4>/Elevation: Transfer Fcn'
   */
  helikopter_B.Gain_dg = (helikopter_P.ElevationTransferFcn_C *
    helikopter_X.ElevationTransferFcn_CSTATE +
    helikopter_P.ElevationTransferFcn_D * helikopter_B.ElevationCounttorad) *
    helikopter_P.Gain_Gain_n;

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
    rtb_TmpSignalConversionAtToFile[1] = helikopter_B.Gain1[1];
    rtb_TmpSignalConversionAtToFile[2] = helikopter_B.Gain1[2];
    rtb_TmpSignalConversionAtToFile[3] = helikopter_B.Gain1[3];

    /* ToFile: '<Root>/To File' */
    {
      if (!(++helikopter_DW.ToFile_IWORK.Decimation % 1) &&
          (helikopter_DW.ToFile_IWORK.Count*5)+1 < 100000000 ) {
        FILE *fp = (FILE *) helikopter_DW.ToFile_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[5];
          helikopter_DW.ToFile_IWORK.Decimation = 0;
          u[0] = helikopter_M->Timing.t[1];
          u[1] = rtb_TmpSignalConversionAtToFile[0];
          u[2] = rtb_TmpSignalConversionAtToFile[1];
          u[3] = rtb_TmpSignalConversionAtToFile[2];
          u[4] = rtb_TmpSignalConversionAtToFile[3];
          if (fwrite(u, sizeof(real_T), 5, fp) != 5) {
            rtmSetErrorStatus(helikopter_M, "Error writing to MAT-file Data.mat");
            return;
          }

          if (((++helikopter_DW.ToFile_IWORK.Count)*5)+1 >= 100000000) {
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

  /* FromWorkspace: '<Root>/From Workspace1' */
  {
    real_T *pDataValues = (real_T *) helikopter_DW.FromWorkspace1_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *) helikopter_DW.FromWorkspace1_PWORK.TimePtr;
    int_T currTimeIndex = helikopter_DW.FromWorkspace1_IWORK.PrevIndex;
    real_T t = helikopter_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
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

    helikopter_DW.FromWorkspace1_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&rtb_Sum3_o[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 141;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&rtb_Sum3_o[0])[elIdx] = pDataValues[currTimeIndex + 1];
              pDataValues += 141;
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
          for (elIdx = 0; elIdx < 4; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&rtb_Sum3_o[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1, f2);
            pDataValues += 141;
          }
        }
      }
    }
  }

  /* Sum: '<Root>/Sum3' */
  rtb_Sum3_o[0] = helikopter_B.Sum5 - rtb_Sum3_o[0];
  rtb_Sum3_o[1] = helikopter_B.Gain1[1] - rtb_Sum3_o[1];
  rtb_Sum3_o[2] = helikopter_B.Gain1[2] - rtb_Sum3_o[2];
  rtb_Sum3_o[3] = helikopter_B.Gain1[3] - rtb_Sum3_o[3];

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *) helikopter_DW.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *) helikopter_DW.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = helikopter_DW.FromWorkspace_IWORK.PrevIndex;
    real_T t = helikopter_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
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

    helikopter_DW.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          rtb_Backgain = pDataValues[currTimeIndex];
        } else {
          rtb_Backgain = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Backgain = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 141;
      }
    }
  }

  /* Sum: '<Root>/Sum1' incorporates:
   *  Constant: '<Root>/Constant'
   *  Constant: '<Root>/Vd_bias'
   *  DotProduct: '<Root>/Dot Product'
   *  Gain: '<S5>/K_pd'
   *  Gain: '<S5>/K_pp'
   *  Sum: '<Root>/Sum4'
   *  Sum: '<S5>/Sum2'
   *  Sum: '<S5>/Sum3'
   */
  helikopter_B.Sum1 = (((rtb_Backgain - (((helikopter_P.K_LQR_T[0] * rtb_Sum3_o
    [0] + helikopter_P.K_LQR_T[1] * rtb_Sum3_o[1]) + helikopter_P.K_LQR_T[2] *
    rtb_Sum3_o[2]) + helikopter_P.K_LQR_T[3] * rtb_Sum3_o[3])) -
                        helikopter_B.Gain1[2]) * helikopter_P.K_pp -
                       helikopter_P.K_pd * helikopter_B.Gain1[3]) +
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

  /* Sum: '<S3>/Sum' incorporates:
   *  Constant: '<Root>/e_ref'
   */
  rtb_Derivative = helikopter_P.e_ref_Value - helikopter_B.Gain1[4];

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

    if (rt_WriteMat4FileHeader(fp,5,0,"data")) {
      rtmSetErrorStatus(helikopter_M,
                        "Error writing mat file header to file Data.mat");
      return;
    }

    helikopter_DW.ToFile_IWORK.Count = 0;
    helikopter_DW.ToFile_IWORK.Decimation = -1;
    helikopter_DW.ToFile_PWORK.FilePtr = fp;
  }

  /* Start for FromWorkspace: '<Root>/From Workspace1' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues0[] = { 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1378421413592839, 3.1262155534547391,
      3.10330930002669, 3.0666274151878978, 3.0144539223908811,
      2.9456562771143004, 2.8595077632902974, 2.7555515879619352,
      2.6335051104871052, 2.49319560602897, 2.3345185760613303,
      2.1574113214684232, 1.9632257910248956, 1.7564180656089003,
      1.5436868085878352, 1.3320083929785957, 1.1270304834580145,
      0.93285515833409871, 0.75228628518120955, 0.58716723964309747,
      0.43867756498675531, 0.30755713767402187, 0.19426300997478957,
      0.099074126756229092, 0.021920806530826845, -0.037969106941276469,
      -0.082000864322075512, -0.11201438528174587, -0.13007758510389633,
      -0.13831437417159287, -0.13877705966936779, -0.13335998629863968,
      -0.12374720045059151, -0.11138646920390421, -0.0974828231134029,
      -0.08300595876779826, -0.068706969555291431, -0.055140858260406343,
      -0.042692107362439932, -0.03160125836584543, -0.021991002714563102,
      -0.013890734332292116, -0.0072588742757507844, -0.0020025647104054407,
      0.002005446697567539, 0.0049127405901952543, 0.006875861716385657,
      0.0080514538155846988, 0.00858963546312763, 0.0086293289706758283,
      0.008295244663195039, 0.0076962282956943686, 0.006924695862142346,
      0.0060569040432640791, 0.0051538330982572452, 0.00426248975147297,
      0.0034174686993073819, 0.0026426413668783565, 0.0019528684821971452,
      0.0013556582521045022, 0.00085271404542970272, 0.00044133436649277324,
      0.00011564356255156858, -0.00013235569250781389, -0.00031190731390423447,
      -0.00043279602967538069, -0.00050479773747891992, -0.0005372695804601276,
      -0.00053886160230484355, -0.00051733131808199152, -0.0004794429335224835,
      -0.000430934011108923, -0.00037653390199007076, -0.000320020058761219,
      -0.00026430026991701319, -0.00021151079805938937, -0.00016312227477006936,
      -0.00012004694419703517, -8.2742414656699576E-5, -5.1308450118581978E-5,
      -2.5574502475939959E-5, -5.1766529779689472E-6, 1.037659292575191E-5,
      2.1649613248022512E-5, 2.9235592659704489E-5, 3.3722267264324943E-5,
      3.5666504277670586E-5, 3.5576630136040989E-5, 3.39013932788114E-5,
      3.1024474013788195E-5, 2.7263515820455925E-5, 2.287273386368144E-5,
      1.8048242210739208E-5, 1.2935317111365279E-5, 7.6368669265227728E-6,
      2.2224002741056454E-6, -3.2632282955983968E-6, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.015002048907599654,
      -0.046506351616669328, -0.091625013710685488, -0.1467275393536602,
      -0.20869397118655766, -0.27519058110481381, -0.344594055294504,
      -0.41582470131193777, -0.488185909897812, -0.56123801783103233,
      -0.63470811986905074, -0.7084290183701184, -0.77674212177260027,
      -0.82723090166247226, -0.85092502808275117, -0.84671366243544866,
      -0.81991163808081657, -0.77670130049415376, -0.72227549261004731,
      -0.66047618215093951, -0.59395869862385942, -0.52448170924942483,
      -0.45317651079542004, -0.38075553287273273, -0.30861328090009987,
      -0.23955965388690414, -0.17612702952168705, -0.12005408383717231,
      -0.072252799287092717, -0.032947156269277134, -0.001850741989590518,
      0.021668293484421532, 0.038451143393701794, 0.049442924988258335,
      0.055614584363514338, 0.0579074573839277, 0.057195956851536472,
      0.054264445181049496, 0.049795003593374777, 0.04436339598788714,
      0.038441022606638453, 0.032401073530593079, 0.026527440227674455,
      0.02102523826289051, 0.016032045633401052, 0.011629175572019991,
      0.00785248450627074, 0.004702368398305296, 0.0021527265916808612,
      0.00015877403170192275, -0.0013363372284140261, -0.0023960654684935479,
      -0.00308612973269896, -0.0034711672740039348, -0.0036122837785182061,
      -0.0035653733856279676, -0.0033800842071532232, -0.00309930932820697,
      -0.002759091537215712, -0.0023888409188614406, -0.0020117768251900663,
      -0.0016455187142385861, -0.0013027632142556869, -0.000991997018728398,
      -0.00071820648407655053, -0.0004835548615754531, -0.00028800682970502484,
      -0.00012988737041569928, -6.3680858697319173E-6, 8.6121138400539629E-5,
      0.00015155353974716389, 0.00019403569116337402, 0.00021760043798454054,
      0.00022605537442453905, 0.00022287915688595488, 0.00021115788893962708,
      0.00019355409466641191, 0.00017230132380126852, 0.00014921811967047414,
      0.00012573585966160216, 0.00010293579207969986, 8.1591399501015831E-5,
      6.22129851240152E-5, 4.5092082798214182E-5, 3.0343919155859688E-5,
      1.79466999276136E-5, 7.7769495625143489E-6, -3.5949505738660423E-7,
      -6.700945919786586E-6, -1.1507675550961046E-5, -1.5043831264197291E-5,
      -1.7563126317966172E-5, -1.9297965102637152E-5, -2.0451698888363934E-5,
      -2.119379923023825E-5, -2.1657865100536735E-5, -2.1942512769684394E-5,
      -2.211524529118034E-5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.10602875205901979, 0.22266037932343921, 0.31888147181643506,
      0.38944360631128311, 0.43795507377658849, 0.46997264230375863,
      0.49051724877532715, 0.5034310014144936, 0.51142138586014851,
      0.51630439857669563, 0.5192586212677589, 0.521031154883602,
      0.48281092448924684, 0.35683541342492997, 0.16746103620184635,
      -0.029764323976778712, -0.18942647182334763, -0.30539416303172051,
      -0.38466082364566145, -0.43677392373401597, -0.47012016899865228,
      -0.49103682602378734, -0.5039579095627138, -0.51184381270528834,
      -0.50987388414221224, -0.48804466254646223, -0.44831756262191058,
      -0.39630216453879508, -0.33784122277996015, -0.27779727311258295,
      -0.21977757968645903, -0.16622356026477, -0.11861477339884366,
      -0.077685714293762656, -0.043618931364285546, -0.016205150806783978,
      0.00502861402445465, 0.020718804874740641, 0.031588306158296511,
      0.038388528098666441, 0.0418570731669063, 0.042688053272842821,
      0.041512596908236457, 0.038887462068548644, 0.035289978489165229,
      0.031117804035663815, 0.026692210046828733, 0.022263817548730479,
      0.018019894521923478, 0.014092495157252622, 0.010566875368890345,
      0.0074897544661524295, 0.0048771106679256267, 0.0027212982872497544,
      0.00099735755831783569, -0.00033154474070146313, -0.0013095531471430959,
      -0.0019844096098707168, -0.00240452940957717, -0.0026167899631451903,
      -0.0026649450044621283, -0.0025885724456534949, -0.0024224649675201326,
      -0.002196376781085647, -0.0019350469318810317, -0.0016584280488584019,
      -0.0013820588048495948, -0.0011175279486868573, -0.00087298712821052378,
      -0.00065367851328205844, -0.00046245122248631817, -0.00030024762121122,
      -0.00016654663052929449, -5.975626154244617E-5, 2.2448292570540459E-5,
      8.2841445299411563E-5, 0.00012441689470452769, 0.00015020646762743203,
      0.00016314327086410893, 0.00016596364540166062, 0.00016114216987295877,
      0.00015085401490627851, 0.00013695923182948101, 0.0001210039988478461,
      0.00010423438803721486, 8.76188107985282E-5, 7.18759115958594E-5,
      5.75052831554501E-5, 4.4818952827206746E-5, 3.3972129310560177E-5,
      2.4992198120239118E-5, 1.780541534180854E-5, 1.2261178033879124E-5,
      8.1541498388266747E-6, 5.2448818700426888E-6, 3.2798404127619175E-6,
      2.011781091002761E-6, 1.2208074913746717E-6, 7.34324476533056E-7, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.42411500823727227,
      0.4665265090591868, 0.38488436997349251, 0.28224853798090122,
      0.19404586986273076, 0.12807027411018959, 0.082178425887783166,
      0.051655010558174862, 0.031961537784128705, 0.019532050867697347,
      0.011816890765762315, 0.0070901344648815921, -0.15288092157591165,
      -0.50390204425575835, -0.75749750889082534, -0.78890144071299106,
      -0.63864859138476648, -0.46387076483198236, -0.31706664245425464,
      -0.20845240035190887, -0.133384981057036, -0.083666628099031179,
      -0.051684334154196723, -0.031543612568789285, 0.00787971425381348,
      0.087316886384509157, 0.15890839969971565, 0.20806159233397115,
      0.23384376703684881, 0.24017579867101801, 0.23207877370600491,
      0.21421607768826528, 0.19043514746521448, 0.16371623642183314,
      0.13626713171941759, 0.10965512223151538, 0.084935059326463655,
      0.0627607634026531, 0.04347800513573262, 0.027200887762988842,
      0.0138741802744686, 0.0033239204252551949, -0.0047018254569163248,
      -0.010500539357242131, -0.014389934316024547, -0.01668869781249651,
      -0.017702375953831197, -0.017713569990883877, -0.016975692105718879,
      -0.01570959745717429, -0.014102479151939974, -0.012308483609442531,
      -0.01045057519139808, -0.008623249521194358, -0.0068957629142185417,
      -0.0053156091945680637, -0.0039120336242573989, -0.0026994258494013523,
      -0.0016804791973166822, -0.00084904221276294913, -0.00019262016375862128,
      0.00030549023674366569, 0.00066442991404258025, 0.00090435274724707355,
      0.0010453193983275932, 0.001106475533599651, 0.0011054769775443605,
      0.0010581234261600814, 0.000978163283414466, 0.00087723446122299291,
      0.00076490916469209279, 0.00064881440660952453, 0.000534803964236834,
      0.000427161477456525, 0.00032881821796107831, 0.00024157261242461618,
      0.00016630179912959636, 0.00010315829320074912, 5.1747214455839342E-5,
      1.128149965933849E-5, -1.9285900605675681E-5, -4.1152618357589208E-5,
      -5.5579130798058159E-5, -6.3820930417407936E-5, -6.70784417333932E-5,
      -6.64623074456149E-5, -6.2971595301543416E-5, -5.74825122525054E-5,
      -5.0745319803841662E-5, -4.3387292557454495E-5, -3.5919723252152476E-5,
      -2.8747129604590533E-5, -2.2176947722585884E-5, -1.6428111271078027E-5,
      -1.1637070366004165E-5, -7.86016431999131E-6, -5.0722357779048506E-6,
      -3.163892889380581E-6, -1.9459305502346861E-6, -1.183447123036953E-6, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0 } ;

    helikopter_DW.FromWorkspace1_PWORK.TimePtr = (void *) pTimeValues0;
    helikopter_DW.FromWorkspace1_PWORK.DataPtr = (void *) pDataValues0;
    helikopter_DW.FromWorkspace1_IWORK.PrevIndex = 0;
  }

  /* Start for FromWorkspace: '<Root>/From Workspace' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5235987755997864,
      0.52359877559773482, 0.52359877559755752, 0.52359877559775247,
      0.52359877559826518, 0.523598775598392, 0.52359877559807744,
      0.52359877559763313, 0.52359877559872592, 0.52359877559711854,
      0.523598775585623, 0.523598775597618, 0.32964141753472442,
      -0.082195946473836629, -0.39016105633372555, -0.523598771544535,
      -0.52359877529793819, -0.52359877545191935, -0.52359877549299516,
      -0.52359877551053735, -0.5235987755261593, -0.523598775579307,
      -0.52359877559822243, -0.52359877559752532, -0.49033559411909139,
      -0.40501799176977038, -0.32447053715467344, -0.25079682012344823,
      -0.1853081702736471, -0.12865843371660385, -0.080975551375583765,
      -0.041984680129488451, -0.011118888704940903, 0.012384738098135506,
      0.029404372373842805, 0.040867877136132894, 0.047701497156645171,
      0.050791463352041342, 0.050956921177675737, 0.048932492835305284,
      0.045358789487693424, 0.040779271797987564, 0.035641989832585243,
      0.030304896713696263, 0.025043608782732076, 0.020060666534903869,
      0.015495526422572927, 0.011434677619788537, 0.0079214259628920536,
      0.0049650172426062636, 0.0025488829076262425, 0.0006378831873677435,
      -0.00081550343581588357, -0.0018660381306534559, -0.0025715781758488478,
      -0.0029898534467063293, -0.0031760655357053868, -0.0031812021337061921,
      -0.00305095637397614, -0.0028251433066478759, -0.0025375120562938919,
      -0.0022158612572670306, -0.0018823760056843567, -0.0015541159586354744,
      -0.0012435956807983219, -0.00095940938305024836, -0.00070686245925612233,
      -0.00048858147225664089, -0.00030508233640485512, -0.00015528334186793423,
      -3.6955378566604206E-5, 5.2893695090900226E-5, 0.0001177000395107964,
      0.00016108704505374662, 0.00018666580015861561, 0.00019788668817396089,
      0.00019793537470500478, 0.00018966627484209167, 0.00017556675122001636,
      0.00015774569831327888, 0.00013794073761740753, 0.00011753891743222823,
      9.7606527501949118E-5, 7.8924362195388731E-5, 6.2025466335813149E-5,
      4.72330548595817E-5, 3.4696899253016711E-5, 2.44270156373311E-5,
      1.6323973402267159E-5, 1.0205577062500813E-5, 5.83006967334796E-6,
      2.9163779061232931E-6, 1.1622769840101979E-6, 2.6168959237886567E-7,
      -7.741529589242848E-8, -1.1308305204924042E-7, -4.6748628876222309E-8,
      7.9664764827920109E-14, 7.9255567490106865E-14, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0 } ;

    helikopter_DW.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    helikopter_DW.FromWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    helikopter_DW.FromWorkspace_IWORK.PrevIndex = 0;
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

      if (rt_WriteMat4FileHeader(fp, 5, helikopter_DW.ToFile_IWORK.Count, "data"))
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
  helikopter_M->Sizes.checksums[0] = (2482973530U);
  helikopter_M->Sizes.checksums[1] = (1605493004U);
  helikopter_M->Sizes.checksums[2] = (1993769543U);
  helikopter_M->Sizes.checksums[3] = (3388951067U);

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
  helikopter_M->Sizes.numBlocks = (65);/* Number of blocks */
  helikopter_M->Sizes.numBlockIO = (18);/* Number of block outputs */
  helikopter_M->Sizes.numBlockPrms = (146);/* Sum of parameter "widths" */
  return helikopter_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
