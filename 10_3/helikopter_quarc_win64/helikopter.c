/*
 * helikopter.c
 *
 * Code generation for model "helikopter".
 *
 * Model version              : 1.195
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Thu Mar 08 13:30:21 2018
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
      3.1415926535897931, 3.1378421413513737, 3.1262155534475076,
      3.1033093000209022, 3.0666274151846253, 3.0144539223915232,
      2.9456562771206443, 2.8595077633046269, 2.755551587987362,
      2.6335051105279459, 2.493195606092121, 2.3345185761615079,
      2.1583785868490466, 1.9678013793910305, 1.7674143111840392,
      1.5625874980938674, 1.3586949815519545, 1.160609174858467,
      0.97237558041863781, 0.79697339345539675, 0.63645596478416488,
      0.49218023370945119, 0.36487408959506529, 0.2546521932177549,
      0.16110387797438525, 0.083403928095003715, 0.020423202350655527,
      -0.0291710392388966, -0.066828956411080664, -0.094046912933813162,
      -0.11231060497587365, -0.12305126577742863, -0.12761341876914015,
      -0.12723258646661925, -0.1230213984208462, -0.11596263112850513,
      -0.10690783149026063, -0.096580309256344463, -0.0855814242641799,
      -0.074399235247363554, -0.063418714305358712, -0.052932861619315749,
      -0.043154176408169348, -0.034226050896825379, -0.026233753290835122,
      -0.019214753017630821, -0.013168216831196757, -0.0080635681554018512,
      -0.003848054892789372, -0.00045331367007891537, 0.0021990479275836609,
      0.0041928041520926133, 0.0056133336218349488, 0.0065446340674039935,
      0.0070670446640132946, 0.007255588399225652, 0.0071788473489172782,
      0.0068982867461100288, 0.0064679486571853628, 0.0059344423489134,
      0.005337165546398237, 0.0047086983385413684, 0.0040753191533885406,
      0.0034575997404622406, 0.0028710432632298626, 0.0023267362793077061,
      0.0018319914728920563, 0.001390963446717242, 0.0010052246554293502,
      0.00067429267040593406, 0.00039610342959854491, 0.00016742798170867935,
      -1.5767470818722666E-5, -0.00015801663695189167, -0.00026410037795433934,
      -0.00033884224288085915, -0.00038694780217422506, -0.00041288294587281139,
      -0.00042078620319196063, -0.00041441021093319131, -0.00039708766735298051,
      -0.00037171741885080645, -0.000340766707165047, -0.00030628602753693873,
      -0.00026993349082189347, -0.00023300602612681614, -0.00019647519006521153,
      -0.00016102575186731615, -0.00012709559024953537, -9.4915759362575907E-5,
      -6.45498488866812E-5, -3.5931968247110289E-5, -8.9028158630426813E-6,
      1.6756662629216596E-5, 4.1294584513793962E-5, 6.4961146184283771E-5,
      8.7986382401692933E-5, 0.00011056407427959987, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.015002048902967588,
      -0.046506351610164379, -0.091625013701120417, -0.14672753933980562,
      -0.20869397116710756, -0.27519058107821587, -0.34459405525876891,
      -0.41582470126375731, -0.48818590983236321, -0.56123801773799753,
      -0.6347081197171518, -0.70455995724454534, -0.762308829826763,
      -0.80154827282266483, -0.81930725235538615, -0.81557006616235084,
      -0.79234322676864843, -0.75293437775401684, -0.70160874784766292,
      -0.64206971467962648, -0.57710292429355392, -0.50922457645224273,
      -0.44088758550394069, -0.37419326096817773, -0.31079979951222531,
      -0.25192290297209191, -0.19837696635290764, -0.15063166868343544,
      -0.10887182608562911, -0.073054768162941153, -0.042962643200918971,
      -0.01824861196154514, 0.0015233292153843778, 0.016844752188393097,
      0.028235069174665145, 0.036219198558278869, 0.041310088940965578,
      0.043995539973959032, 0.044728756072566281, 0.043922083773320249,
      0.041943410749472734, 0.039114740849886441, 0.035712502050676755,
      0.031969190429261891, 0.028076001098118069, 0.024186144751037128,
      0.020418594708480484, 0.016862053055750784, 0.013578964896142692,
      0.010609446395951169, 0.0079750249033366735, 0.0056821178842702113,
      0.0037252017875770427, 0.0020896423917380709, 0.00075417494615029215,
      -0.0003069641959326313, -0.0011222424059281304, -0.0017213523503977974,
      -0.00213402522778699, -0.0023891072047597822, -0.0025138688261266118,
      -0.0025335167353104449, -0.0024708776464043365, -0.0023462259036286458,
      -0.0021772279303877611, -0.0019789792203617346, -0.0017641120993983916,
      -0.0015429551598507021, -0.0013237279347927992, -0.0011127569579286916,
      -0.00091470178625859715, -0.00073278180480874293, -0.000568996659231811,
      -0.00042433495870892557, -0.00029896745440521411, -0.00019242223187259855,
      -0.0001037405694934802, -3.1613023975731892E-5, 2.5503974335942629E-5,
      6.9290179621708263E-5, 0.00010148099930956146, 0.00012380285204390272,
      0.00013792272381329824, 0.00014541015216104608, 0.0001477098640811745,
      0.00014612334954728354, 0.00014179775809244658, 0.00013572065177198828,
      0.00012871932884870289, 0.00012146364720444394, 0.00011447152785914877,
      0.00010811661483713556, 0.00010263791926990224, 9.81516928391746E-5,
      9.4666251982824354E-5, 9.2100950170501838E-5, 9.0310772812492854E-5,
      8.911758147801024E-5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.10602875205327919, 0.22266037931040097, 0.31888147179500592,
      0.3894436062811652, 0.43795507373723974, 0.46997264225343838,
      0.49051724871094826, 0.50343100132673324, 0.51142138573830032,
      0.51630439838192654, 0.519258620851929, 0.49368610960604914,
      0.40814697577673936, 0.27732939664137823, 0.12551368477112451,
      -0.026413004693958346, -0.1641584299654065, -0.2785266936811146,
      -0.3627499497286249, -0.42079914709287869, -0.45916046212351075,
      -0.47973823822577494, -0.48297975253954012, -0.47137001370813908,
      -0.44804077413625315, -0.41611941829220406, -0.37844222958824852,
      -0.33744552889031765, -0.29514261842798328, -0.25314128603163183,
      -0.21267964635040587, -0.17466933393753273, -0.13974052887506527,
      -0.10828596596594105, -0.080502410232294153, -0.056428777158800063,
      -0.035980468895367546, -0.018979742264847783, -0.0051820913524342777,
      0.0057012517244764685, 0.013984505231683689, 0.019991958515872175,
      0.02404572443913636, 0.026456296882856247, 0.027515575293444915,
      0.02749201903503241, 0.026627604786270836, 0.025136278075975065,
      0.023203613224826756, 0.020987422631600576, 0.018619084963370149,
      0.016205391096635824, 0.013830735580032194, 0.011559509151976703,
      0.0094385738596814948, 0.0074997261828547886, 0.0057620750148566838,
      0.0042342802737974534, 0.0029166142883153418, 0.0018028219917943973,
      0.000881767490739685, 0.00013886391839638794, -0.00044270915678906568,
      -0.000880990910267597, -0.001194413130341, -0.0014011461664031652,
      -0.0015185987479513054, -0.0015630527834777795, -0.0015494142984667136,
      -0.0014910622894857626, -0.0013997783112202266, -0.0012857409491471728,
      -0.0011575708553120448, -0.0010224136493575192, -0.00088604963938170586,
      -0.00075302094048982878, -0.0006267681198331036, -0.00050976994425835107,
      -0.00040368112948215853, -0.00030946417571426694, -0.00022751241885508294,
      -0.00015776232961517539, -9.9793860681817684E-5, -5.2918283661677681E-5,
      -1.6253485260561804E-5, 1.1212879076374181E-5, 3.0571628636577061E-5,
      4.2950666762612977E-5, 4.9482676731625373E-5, 5.1280387025967905E-5,
      4.9417629357450518E-5, 4.4914098428376625E-5, 3.8721328105773214E-5,
      3.1706935295947127E-5, 2.4633765079484385E-5, 1.8130573738751466E-5,
      1.2652290133681452E-5, 8.4330206446613281E-6, 5.4490049989459693E-6, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.42411500821730741,
      0.46652650903378806, 0.38488436994372066, 0.28224853794993804,
      0.19404586982959893, 0.12807027407009547, 0.082178425835340491,
      0.0516550104684409, 0.031961537651569179, 0.019532050579805747,
      0.011816889885310452, -0.10229004497821843, -0.34215653531193824,
      -0.52327031653614364, -0.60726284747571391, -0.60770675785503059,
      -0.55098170108049183, -0.45747305485753154, -0.33689302418474021,
      -0.23219678945171413, -0.15344526011722753, -0.082311104403755678,
      -0.012966057249759872, 0.046438955330904967, 0.093316958292844712,
      0.12768542338149727, 0.150708754821123, 0.16398680279702443,
      0.16921164185463825, 0.16800532959070666, 0.16184655873020473,
      0.15204124965679336, 0.13971522025517075, 0.12581825164179777,
      0.11113422293988841, 0.09629453229927723, 0.081793233059030968,
      0.0680029065273799, 0.055190603654954887, 0.043533372312943852,
      0.033133014034129743, 0.024029813142054813, 0.016215063698357603,
      0.00964228978018042, 0.0042371136476555341, -9.4225028349154662E-5,
      -0.0034576569897454344, -0.0059653068358822186, -0.0077306593992923734,
      -0.00886476236760385, -0.0094733506676208477, -0.00965477546163643,
      -0.0094986220611136575, -0.0090849057069210953, -0.0084837411638799716,
      -0.0077553907020059551, -0.006950604666691557, -0.0061111789589360557,
      -0.00527066393662758, -0.0044551691807829131, -0.003684217998917984,
      -0.0029716142840723234, -0.0023262922954409493, -0.00175312700861326,
      -0.0012536888749927475, -0.00082693213894779528, -0.00046981032089169524,
      -0.00017781613680503219, 5.4553945345129154E-5, 0.00023340804122466918,
      0.00036513591836300863, 0.00045614945359308045, 0.00051268038064137687,
      0.00054062882911896828, 0.00054545604520411773, 0.00053211480086837377,
      0.00050501128792776538, 0.00046799270759987575, 0.00042435526440563515,
      0.00037686782037243136, 0.0003278070327376011, 0.00027900036226049532,
      0.00023187388103429598, 0.00018750231338142512, 0.00014665919890532862,
      0.00010986546264860906, 7.7435003541676656E-5, 4.9516157805008771E-5,
      2.6128045176914721E-5, 7.1908464782352542E-6, -7.4510253732044127E-6,
      -1.8014118415430436E-5, -2.477107598954853E-5, -2.8057565938439239E-5,
      -2.8292675564985835E-5, -2.6012760062066553E-5, -2.1913129119414933E-5,
      -1.6877072655215379E-5, -1.1936057281996314E-5, -8.0243524503753083E-6,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0 } ;

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
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877557065593,
      0.5235987755648287, 0.52359877555398415, 0.52359877554595047,
      0.52359877553554213, 0.52359877551865264, 0.52359877549472145,
      0.52359877543879152, 0.5235987753775585, 0.52359877516608888,
      0.52359877440812641, 0.38856151693560531, 0.10947154576688319,
      -0.11008507929484046, -0.27695896968853256, -0.39795513806507843,
      -0.47968511189524515, -0.52317212167451765, -0.52359758885893459,
      -0.5235976653329798, -0.52352227216915193, -0.50347393523708062,
      -0.46500607547431955, -0.42080544680465126, -0.37350683555807906,
      -0.32525430068660577, -0.27774421552118833, -0.23227259445252346,
      -0.18978425234113133, -0.15092186790208051, -0.11607345081518917,
      -0.085417096181011981, -0.05896224476433104, -0.03658695170110076,
      -0.018070901997364233, -0.003124101767771162, 0.00858867597058542,
      0.0174274984173802, 0.023760658261973329, 0.027951507047099845,
      0.030348373830361489, 0.031277192136662961, 0.031036458786475471,
      0.029894160251148567, 0.028086322790823657, 0.025816869231145775,
      0.023258495674085292, 0.020554313921849133, 0.017820038464696847,
      0.015146529393943722, 0.012602533674664736, 0.010237496212285615,
      0.0080843386427271943, 0.0061621275094993957, 0.00447857434798449,
      0.0030323281761881413, 0.0018150360803460127, 0.00081316013712926091,
      9.5490327520515231E-6, -0.00061522934582098167, -0.0010817821136728937,
      -0.0014109958957998043, -0.0016233323509942598, -0.0017382938394307028,
      -0.0017740377727177073, -0.0017471183870815467, -0.00167233549054803,
      -0.0015626709929298828, -0.0014292955945502514, -0.0012816297687624986,
      -0.0011274450275163888, -0.00097299333024378611, -0.00082315432293468321,
      -0.00068159183027797306, -0.00055091263614977075, -0.00043282205481778667,
      -0.0003282721047272006, -0.0002375992441486923, -0.00016064961504931477,
      -9.6890574923673E-5, -4.5507986359221096E-5, -5.4892939661724343E-6,
      2.4307136746729174E-5, 4.50955969369442E-5, 5.8118439066833538E-5,
      6.4612010149800811E-5, 6.5781641290101316E-5, 6.2784055224599484E-5,
      5.6715362742299572E-5, 4.860256452775517E-5, 3.9396115260259257E-5,
      2.9960608447093088E-5, 2.1060013408334354E-5, 1.3333279773730962E-5,
      7.2559901527849288E-6, 3.085333720210612E-6, 7.9196522271045243E-7,
      2.6461258544509908E-13, 2.7287480603657809E-13, 0.0, 0.0, 0.0, 0.0, 0.0,
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
  helikopter_M->Sizes.checksums[0] = (1572777922U);
  helikopter_M->Sizes.checksums[1] = (2568827870U);
  helikopter_M->Sizes.checksums[2] = (1422397060U);
  helikopter_M->Sizes.checksums[3] = (3002865278U);

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
