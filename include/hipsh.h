/**
 * @file  hipsh.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
 *    $Revision: 1.6 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#ifndef HIPSH_H
#define HIPSH_H

#include "hips.h"

/**************** hipsh.h *********************/
/******************* seqord.c ***********************/
int h_seqord(struct header *hdi,struct header *hdo);
int h_seqord_i(struct header *hdi,struct header *hdo);
int h_seqord_f(struct header *hdi,struct header *hdo);
int h_seqord_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo);
int h_seqord_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo);
int h_invseqord(struct header *hdi,struct header *hdo);
int h_invseqord_i(struct header *hdi,struct header *hdo);
int h_invseqord_f(struct header *hdi,struct header *hdo);
int h_invseqord_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo);
int h_invseqord_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo);
int invorder(int index,int len);
/******************* abdou.c ***********************/
int h_abdou(struct header *hdi,struct header *hdo,int size);
int h_abdou_b(struct header *hdi,struct header *hdo,int size);
int h_abdou_B(unsigned char *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,int size);
/******************* hconvolv.c ***********************/
int h_hconvolve(struct header *hdi,struct header *hdo,int *mask,int nmask,int offset);
int h_hconvolve_i(struct header *hdi,struct header *hdo,int *mask,int nmask,int offset);
int h_hconvolve_f(struct header *hdi,struct header *hdo,float *mask,int nmask,int offset);
int h_hconvolve_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo,int *mask,int nmask,int offset);
int h_hconvolve_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,float *mask,int nmask,int offset);
/******************* nonisot.c ***********************/
int h_nonisot(struct header *hdi,struct header *hdc,int nmasks,struct header **hdm ,struct header *hdo);
int h_nonisot_i(struct header *hdi,struct header *hdc,int nmasks,struct header **hdm ,struct header *hdo);
int h_nonisot_I(int *imagei,int *imagec,int nmasks,int **imagem ,int *imageo,int nr,int nc,int nlpi,int nlpc,int *nrm,int *ncm,int *nlpm,int nlpo);
/******************* abs.c ***********************/
int h_abs(struct header *hdi,struct header *hdo);
int h_abs_i(struct header *hdi,struct header *hdo);
int h_abs_f(struct header *hdi,struct header *hdo);
int h_abs_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo);
int h_abs_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo);
/******************* absdiff.c ***********************/
int h_absdiff(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_absdiff_b(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_absdiff_s(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_absdiff_i(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_absdiff_f(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_absdiff_d(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_absdiff_B(unsigned char *imagei1,unsigned char *imagei2,unsigned char *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_absdiff_S(short *imagei1,short *imagei2,short *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_absdiff_I(int *imagei1,int *imagei2,int *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_absdiff_F(float *imagei1,float *imagei2,float *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_absdiff_D(double *imagei1,double *imagei2,double *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
/******************* add.c ***********************/
int h_add(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_add_bii(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_add_bsb(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_add_s(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_add_i(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_add_f(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_add_d(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_add_c(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_add_dc(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_add_ip(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_add_fp(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_add_BII(unsigned char *imagei1,int *imagei2,int *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_add_BSB(unsigned char *imagei1,short *imagei2,unsigned char *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_add_S(short *imagei1,short *imagei2,short *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_add_I(int *imagei1,int *imagei2,int *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_add_F(float *imagei1,float *imagei2,float *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_add_D(double *imagei1,double *imagei2,double *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_add_C(float *imagei1,float *imagei2,float *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_add_DC(double *imagei1,double *imagei2,double *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
/******************* addcos.c ***********************/
int h_addcos(struct header *hdi,struct header *hdo,double xf,double yf,double phase,double amplitude);
int h_addcos_f(struct header *hdi,struct header *hdo,double xf,double yf,double phase,double amplitude);
int h_addcos_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double xf,double yf,double phase,double amplitude);
/******************* addgabor.c ***********************/
int h_addgabor(struct header *hdi,struct header *hdo,double xm,double ym,double xf,double yf,double xs,double ys,double phase,double amplitude);
int h_addgabor_f(struct header *hdi,struct header *hdo,double xm,double ym,double xf,double yf,double xs,double ys,double phase,double amplitude);
int h_addgabor_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double xm,double ym,double xf,double yf,double xs,double ys,double phase,double amplitude);
/******************* affine.c ***********************/
int h_affine(struct header *hdi,struct header *hdo,double A,double B,double C,double a,double b,double c);
int h_affine_b(struct header *hdi,struct header *hdo,double A,double B,double C,double a,double b,double c);
int h_affine_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nor,int noc,int nlpo,double A,double B,double C,double a,double b,double c);
/******************* and.c ***********************/
int h_and(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_and_mp(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_and_lp(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_and_b(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_and_MP(unsigned char *imagei1,unsigned char *imagei2,unsigned char *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_and_LP(unsigned char *imagei1,unsigned char *imagei2,unsigned char *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_and_B(unsigned char *imagei1,unsigned char *imagei2,unsigned char *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
/******************* applylut.c ***********************/
int h_applylut(struct header *hdi,struct header *hdo,int count,unsigned char *lut);
int h_applylut_b(struct header *hdi,struct header *hdo,int count,unsigned char *lut);
int h_applylut_s(struct header *hdi,struct header *hdo,int count,short *lut);
int h_applylut_i(struct header *hdi,struct header *hdo,int count,int *lut);
int h_applylut_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,int count,unsigned char *lut);
int h_applylut_S(short *imagei,short *imageo,int nr,int nc,int nlpi,int nlpo,int count,short *lut);
int h_applylut_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo,int count,int *lut);
/******************* avg.c ***********************/
int h_avg(struct header *hdi1,struct header *hdi2,struct header *hdo,double wt1,double wt2);
int h_avg_b(struct header *hdi1,struct header *hdi2,struct header *hdo,double wt1,double wt2);
int h_avg_s(struct header *hdi1,struct header *hdi2,struct header *hdo,double wt1,double wt2);
int h_avg_i(struct header *hdi1,struct header *hdi2,struct header *hdo,double wt1,double wt2);
int h_avg_f(struct header *hdi1,struct header *hdi2,struct header *hdo,double wt1,double wt2);
int h_avg_B(unsigned char *imagei1,unsigned char *imagei2,unsigned char *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo,double wt1,double wt2);
int h_avg_S(short *imagei1,short *imagei2,short *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo,double wt1,double wt2);
int h_avg_I(int *imagei1,int *imagei2,int *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo,double wt1,double wt2);
int h_avg_F(float *imagei1,float *imagei2,float *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo,double wt1,double wt2);
/******************* bclean.c ***********************/
int h_bclean(struct header *hd,int size);
int h_bclean_b(struct header *hd,int size);
int h_bclean_B(unsigned char *image,int nr,int nc,int nlp,int size);
int addneigh(int xx,int yy,int comp);
int inlist(int xx,int yy);
/******************* bnoise.c ***********************/
int h_bnoise(struct header *hdi,struct header *hdo,int n,double p,union pixelval *addc,union pixelval *mulc);
int h_bnoise_b(struct header *hdi,struct header *hdo,int n,double p,union pixelval *addc,union pixelval *mulc);
int h_bnoise_i(struct header *hdi,struct header *hdo,int n,double p,union pixelval *addc,union pixelval *mulc);
int h_bnoise_f(struct header *hdi,struct header *hdo,int n,double p,union pixelval *addc,union pixelval *mulc);
int h_bnoise_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,int n,double p,int addc,int mulc);
int h_bnoise_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo,int n,double p,int addc,int mulc);
int h_bnoise_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,int n,double p,double addc,double mulc);
int compute_blookup(int n,double p);
/******************* checkers.c ***********************/
int h_checkers(struct header *hd,int highflag);
int h_checkers_b(struct header *hd,int highflag);
int h_checkers_B(unsigned char *image,int nr,int nc,int nlp,int highflag);
/******************* clearhis.c ***********************/
int h_clearhisto(struct hips_histo *histogram);
int h_Clearhisto(int nbins,int *histo);
/******************* colorkey.c ***********************/
int h_colorkey(struct header *hdc,int nimage,struct header **hdi ,struct header *hdo,int bflag);
int h_colorkey_bb(struct header *hdc,int nimage,struct header **hdi ,struct header *hdo,int bflag);
int h_colorkey_bi(struct header *hdc,int nimage,struct header **hdi ,struct header *hdo,int bflag);
int h_colorkey_bf(struct header *hdc,int nimage,struct header **hdi ,struct header *hdo,int bflag);
int h_colorkey_ib(struct header *hdc,int nimage,struct header **hdi ,struct header *hdo,int bflag);
int h_colorkey_ii(struct header *hdc,int nimage,struct header **hdi ,struct header *hdo,int bflag);
int h_colorkey_if(struct header *hdc,int nimage,struct header **hdi ,struct header *hdo,int bflag);
int h_colorkey_ipip(struct header *hdc,int nimage,struct header **hdi ,struct header *hdo,int bflag);
int h_colorkey_ipfp(struct header *hdc,int nimage,struct header **hdi ,struct header *hdo,int bflag);
int h_colorkey_BB(unsigned char *imagec,int nimage,unsigned char **imagei ,unsigned char *imageo,int nr,int nc,int nlpc,int *nlpi,int nlpo,int bflag);
int h_colorkey_BI(unsigned char *imagec,int nimage,int **imagei ,int *imageo,int nr,int nc,int nlpc,int *nlpi,int nlpo,int bflag);
int h_colorkey_BF(unsigned char *imagec,int nimage,float **imagei ,float *imageo,int nr,int nc,int nlpc,int *nlpi,int nlpo,int bflag);
int h_colorkey_IB(int *imagec,int nimage,unsigned char **imagei ,unsigned char *imageo,int nr,int nc,int nlpc,int *nlpi,int nlpo,int bflag);
int h_colorkey_II(int *imagec,int nimage,int **imagei ,int *imageo,int nr,int nc,int nlpc,int *nlpi,int nlpo,int bflag);
int h_colorkey_IF(int *imagec,int nimage,float **imagei ,float *imageo,int nr,int nc,int nlpc,int *nlpi,int nlpo,int bflag);
int h_colorkargs(int nimage,struct header **hdi );
int h_colorkargsp(int nimage,struct header **hdi );
/******************* combine.c ***********************/
int h_combine(struct header *hdi1,struct header *hdi2,struct header *hdo,int phasemagflag);
int h_combine_f(struct header *hdi1,struct header *hdi2,struct header *hdo,int phasemagflag);
int h_combine_d(struct header *hdi1,struct header *hdi2,struct header *hdo,int phasemagflag);
int h_combine_F(float *imagei1,float *imagei2,float *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo,int phasemagflag);
int h_combine_D(double *imagei1,double *imagei2,double *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo,int phasemagflag);
/******************* convolve.c ***********************/
int h_convolve(struct header *hdi,struct header *hdk,struct header *hdo,int firsti,int firstk,int numf);
int h_convolve_i(struct header *hdi,struct header *hdk,struct header *hdo,int firsti,int firstk,int numf);
int h_convolve_f(struct header *hdi,struct header *hdk,struct header *hdo,int firsti,int firstk,int numf);
int h_convolve_I(int *imagei,int *imagek,int *imageo,int nri,int nci,int tnri,int nlpi,int nrk,int nck,int tnrk,int nlpk,int nfrk,int nlpo,int firsti,int firstk,int numf);
int h_convolve_F(float *imagei,float *imagek,float *imageo,int nri,int nci,int tnri,int nlpi,int nrk,int nck,int tnrk,int nlpk,int nfrk,int nlpo,int firsti,int firstk,int numf);
/******************* copy.c ***********************/
int h_copy(struct header *hdi,struct header *hdo);
int h_copy_mp(struct header *hdi,struct header *hdo);
int h_copy_lp(struct header *hdi,struct header *hdo);
int h_copy_b(struct header *hdi,struct header *hdo);
int h_copy_s(struct header *hdi,struct header *hdo);
int h_copy_i(struct header *hdi,struct header *hdo);
int h_copy_f(struct header *hdi,struct header *hdo);
int h_copy_d(struct header *hdi,struct header *hdo);
int h_copy_c(struct header *hdi,struct header *hdo);
int h_copy_dc(struct header *hdi,struct header *hdo);
int h_copy_MP(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo);
int h_copy_LP(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo);
int h_copy_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo);
int h_copy_S(short *imagei,short *imageo,int nr,int nc,int nlpi,int nlpo);
int h_copy_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo);
int h_copy_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo);
int h_copy_D(double *imagei,double *imageo,int nr,int nc,int nlpi,int nlpo);
int h_copy_C(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo);
int h_copy_DC(double *imagei,double *imageo,int nr,int nc,int nlpi,int nlpo);
/******************* correl.c ***********************/
int h_correl(struct header *hdi1,struct header *hdi2,struct header *hdo,int dr0,int dc0);
int h_correl_f(struct header *hdi1,struct header *hdi2,struct header *hdo,int dr0,int dc0);
int h_correl_F(float *imagei1,float *imagei2,float *imageo,int nr1,int nc1,int nr2,int nc2,int nro,int nco,int nlpi1,int nlpi2,int nlpo,int dr0,int dc0);
/******************* dct.c ***********************/
int h_dct(struct header *ihd,struct header *ohd);
int h_dct_f(struct header *ihd,struct header *ohd);
int h_dct_d(struct header *ihd,struct header *ohd);
int h_dct_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo);
int h_dct_D(double *imagei,double *imageo,int nr,int nc,int nlpi,int nlpo);
int h_dctinv(struct header *ihd,struct header *ohd);
int h_dctinv_f(struct header *ihd,struct header *ohd);
int h_dctinv_d(struct header *ihd,struct header *ohd);
int h_dctinv_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo);
int h_dctinv_D(double *imagei,double *imageo,int nr,int nc,int nlpi,int nlpo);
int h_dct2d_f(float *veci,float *veco,int logrows,int logcols,int nlpi,int nlpo);
int h_dctinv2d_f(float *veci,float *veco,int logrows,int logcols,int nlpi,int nlpo);
int h_dct2d_d(double *veci,double *veco,int logrows,int logcols,int nlpi,int nlpo);
int h_dctinv2d_d(double *veci,double *veco,int logrows,int logcols,int nlpi,int nlpo);
int h_dctfree(void);
/******************* diff.c ***********************/
int h_diff(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_diff_s(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_diff_i(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_diff_ibi(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_diff_f(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_diff_d(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_diff_c(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_diff_dc(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_diff_ip(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_diff_fp(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_diff_S(short *imagei1,short *imagei2,short *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_diff_I(int *imagei1,int *imagei2,int *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_diff_IBI(int *imagei1,unsigned char *imagei2,int *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_diff_F(float *imagei1,float *imagei2,float *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_diff_D(double *imagei1,double *imagei2,double *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_diff_C(float *imagei1,float *imagei2,float *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_diff_DC(double *imagei1,double *imagei2,double *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
/******************* discedge.c ***********************/
int h_discedge2(struct header *hdi,struct header *hdo,int size,double varcrit,int edgethresh);
int h_discedge2_b(struct header *hdi,struct header *hdo,int size,double varcrit,int edgethresh);
int h_discedge2_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,int size,double varcrit,int edgethresh);
/******************* disphist.c ***********************/
int h_disphist(struct hips_histo *histo,struct header *hd,int barwidth,int barheight,int maxcnt,int borderw,int borderg);
int h_disphist_b(struct hips_histo *histo,struct header *hd,int barwidth,int barheight,int maxcnt,int borderw,int borderg);
int h_disphist_B(int *histo,int nbins,unsigned char *image,int nlp,int barwidth,int barheight,int maxcnt,int borderw,int borderg);
/******************* dither.c ***********************/
int h_dither(struct header *hdi,struct header *hdo);
int h_dither_b(struct header *hdi,struct header *hdo);
int h_dither_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo);
/******************* div.c ***********************/
int h_div(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_div_b(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_div_s(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_div_i(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_div_f(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_div_fc(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_div_d(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_div_ddc(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_div_c(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_div_cf(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_div_dc(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_div_dcd(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_div_ip(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_div_fp(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_div_B(unsigned char *imagei1,unsigned char *imagei2,unsigned char *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_div_S(short *imagei1,short *imagei2,short *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_div_I(int *imagei1,int *imagei2,int *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_div_F(float *imagei1,float *imagei2,float *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_div_FC(float *imagei1,float *imagei2,float *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_div_D(double *imagei1,double *imagei2,double *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_div_DDC(double *imagei1,double *imagei2,double *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_div_C(float *imagei1,float *imagei2,float *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_div_CF(float *imagei1,float *imagei2,float *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_div_DC(double *imagei1,double *imagei2,double *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_div_DCD(double *imagei1,double *imagei2,double *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
/******************* divscale.c ***********************/
int h_divscale(struct header *hdi,struct header *hdo,union pixelval *b);
int h_divscale_s(struct header *hdi,struct header *hdo,union pixelval *b);
int h_divscale_ib(struct header *hdi,struct header *hdo,union pixelval *b);
int h_divscale_i(struct header *hdi,struct header *hdo,union pixelval *b);
int h_divscale_if(struct header *hdi,struct header *hdo,union pixelval *b);
int h_divscale_f(struct header *hdi,struct header *hdo,union pixelval *b);
int h_divscale_c(struct header *hdi,struct header *hdo,union pixelval *b);
int h_divscale_dc(struct header *hdi,struct header *hdo,union pixelval *b);
int h_divscale_S(short *imagei,short *imageo,int nr,int nc,int nlpi,int nlpo,short b);
int h_divscale_Ss(short *imagei,short *imageo,int nr,int nc,int nlpi,int nlpo,short b);
int h_divscale_Sf(short *imagei,short *imageo,int nr,int nc,int nlpi,int nlpo,short b);
int h_divscale_IB(int *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,int b);
int h_divscale_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo,int b);
int h_divscale_IF(int *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double b);
int h_divscale_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double b);
int h_divscale_C(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double b);
int h_divscale_DC(double *imagei,double *imageo,int nr,int nc,int nlpi,int nlpo,double b);
/******************* dotdiff.c ***********************/
int h_dotdiff(struct header *hdi,struct header *hdt,struct header *hdo);
int h_dotdiff_b(struct header *hdi,struct header *hdt,struct header *hdo);
int h_dotdiff_S(unsigned char *imagei,short *imaget,unsigned char *imageo,int nr,int nc,int nlpi,int nlpt,int nlpo);
int weight(int x,int y);
/******************* enlarge.c ***********************/
int h_enlarge(struct header *hdi,struct header *hdo,int xf,int yf);
int h_enlarge_b(struct header *hdi,struct header *hdo,int xf,int yf);
int h_enlarge_i(struct header *hdi,struct header *hdo,int xf,int yf);
int h_enlarge_f(struct header *hdi,struct header *hdo,int xf,int yf);
int h_enlarge_c(struct header *hdi,struct header *hdo,int xf,int yf);
int h_enlarge_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,int xf,int yf);
int h_enlarge_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo,int xf,int yf);
int h_enlarge_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,int xf,int yf);
int h_enlarge_C(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,int xf,int yf);
/******************* entropy.c ***********************/
double h_entropy(int *table,int count,int pairflag);
/******************* entropyc.c ***********************/
int h_entropycnt(struct header *hd,int *table,int pairflag);
int h_entropycnt_b(struct header *hd,int *table,int pairflag);
int h_entropycnt_B(unsigned char *image,int nr,int nc,int nlp,int *table,int pairflag);
/******************* exp.c ***********************/
int h_exp(struct header *hdi,struct header *hdo);
int h_exp_b(struct header *hdi,struct header *hdo);
int h_exp_s(struct header *hdi,struct header *hdo);
int h_exp_i(struct header *hdi,struct header *hdo);
int h_exp_f(struct header *hdi,struct header *hdo);
int h_exp_B(unsigned char *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo);
int h_exp_S(short *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo);
int h_exp_I(int *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo);
int h_exp_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo);
/******************* extract.c ***********************/
int h_extract(struct header *hdi,struct header *hdo,int frow,int fcol,int nrows,int ncols);
/******************* extremum.c ***********************/
int h_extremum(struct header *hdi,struct header *hdo,int size);
int h_extremum_b(struct header *hdi,struct header *hdo,int size);
int h_extremum_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,int size);
/******************* fastaddc.c ***********************/
int h_fastaddcos(struct header *hdi,struct header *hdo,int xf,int yf,double phase,double amplitude);
int h_fastaddcos_f(struct header *hdi,struct header *hdo,int xf,int yf,double phase,double amplitude);
int h_fastaddcos_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,int xf,int yf,double phase,double amplitude);
/******************* fft2.c ***********************/
int h_fft_ri_c(float *rvec,float *ivec,int loglen);
int h_fft2d_ri_c(float *rvec,float *ivec,int loglen);
int h_fft2dgen_ri_c(float *rvec,float *ivec,int logrows,int logcols);
int h_fftn_ri_c(float *rvec,float *ivec,int loglen,int nskip);
int h_fft_ri_dc(double *rvec,double *ivec,int loglen);
int h_fft2d_ri_dc(double *rvec,double *ivec,int loglen);
int h_fft2dgen_ri_dc(double *rvec,double *ivec,int logrows,int logcols);
int h_fftn_ri_dc(double *rvec,double *ivec,int loglen,int nskip);
/******************* filter.c ***********************/
int h_filter(struct header *hdi,struct header *hdo,struct hips_filter *filter);
int h_filter_f(struct header *hdi,struct header *hdo,struct hips_filter *filter);
int h_filter_d(struct header *hdi,struct header *hdo,struct hips_filter *filter);
int h_filter_c(struct header *hdi,struct header *hdo,struct hips_filter *filter);
int h_filter_dc(struct header *hdi,struct header *hdo,struct hips_filter *filter);
int h_filter_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,int method,int disttype,int ftype,double dmetric,double lowcut,int loworder,double highcut,int highorder);
int h_filter_D(double *imagei,double *imageo,int nr,int nc,int nlpi,int nlpo,int method,int disttype,int ftype,double dmetric,double lowcut,int loworder,double highcut,int highorder);
int h_filter_C(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,int method,int disttype,int ftype,double dmetric,double lowcut,int loworder,double highcut,int highorder);
int h_filter_DC(double *imagei,double *imageo,int nr,int nc,int nlpi,int nlpo,int method,int disttype,int ftype,double dmetric,double lowcut,int loworder,double highcut,int highorder);
int calcdist(int disttype,double dmetric,int nr2,int nc2);
double modulate(void);
/******************* flipquad.c ***********************/
int h_flipquad(struct header *hdi,struct header *hdo);
int h_flipquad_b(struct header *hdi,struct header *hdo);
int h_flipquad_f(struct header *hdi,struct header *hdo);
int h_flipquad_d(struct header *hdi,struct header *hdo);
int h_flipquad_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo);
int h_flipquad_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo);
int h_flipquad_D(double *imagei,double *imageo,int nr,int nc,int nlpi,int nlpo);
/******************* fourtr.c ***********************/
int h_fourtr(struct header *hd);
int h_fourtr_c(struct header *hd);
int h_fourtr_dc(struct header *hd);
int h_fourtr_C(float (*image)[2],int nr,int nc,int nlp);
int h_fourtr_DC(double (*image)[2],int nr,int nc,int nlp);
int h_fft_c(float (*vec)[2],int loglen);
int h_fft2d_c(float (*vec)[2],int loglen,int nlp);
int h_fft2dgen_c(float (*vec)[2],int logrows,int logcols,int nlp);
int h_fftn_c(float (*vec)[2],int loglen,int nskip);
int h_fft_dc(double (*vec)[2],int loglen);
int h_fft2d_dc(double (*vec)[2],int loglen,int nlp);
int h_fft2dgen_dc(double (*vec)[2],int logrows,int logcols,int nlp);
int h_fftn_dc(double (*vec)[2],int loglen,int nskip);
/******************* fourtr3d.c ***********************/
int h_fourtr3d(struct header *hd);
int h_fourtr3d_c(struct header *hd);
int h_fourtr3d_C(float (*image)[2],int nr,int nc,int nf);
/******************* gaussmas.c ***********************/
double h_gaussmask(double sigma,int nmask,float *maskarr,int precision);
/******************* gnoise.c ***********************/
int h_gnoise(struct header *hdi,struct header *hdo,double p,int fastflag);
int h_gnoise_s(struct header *hdi,struct header *hdo,double p,int fastflag);
int h_gnoise_f(struct header *hdi,struct header *hdo,double p,int fastflag);
int h_gnoise_S(short *imagei,short *imageo,int nr,int nc,int nlpi,int nlpo,double p,int fastflag);
int h_gnoise_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double p,int fastflag);
double rand_gauss(void);
int srand_g(unsigned int x);
double rand_g(void);
/******************* greybar.c ***********************/
int h_greybar(struct header *hd,int width,double low,double step);
int h_greybar_b(struct header *hd,int width,double low,double step);
int h_greybar_B(unsigned char *image,int nr,int nc,int nlp,int width,double low,double step);
/******************* gridwarp.c ***********************/
int h_gridwarp(struct header *hdi,struct header *hdo,int nx,int ny,float *ogridx,float *ogridy,float *igridx,float *igridy);
int h_gridwarp_b(struct header *hdi,struct header *hdo,int nx,int ny,float *ogridx,float *ogridy,float *igridx,float *igridy);
int h_gridwarp_B(unsigned char *imagei,unsigned char *imageo,int nri,int nci,int nlpi,int nro,int nco,int nlpo,int nx,int ny,float *ogridx,float *ogridy,float *igridx,float *igridy);
int truncup(double x);
int truncdown(double x);
/******************* halftone.c ***********************/
int h_halftone(struct header *hdi,struct header *hdo);
int h_halftone2(struct header *hdi,struct header *hdo,int lower,int upper,int rflag,int alpha,int beta,int gamma,int delta);
int h_halftone_b(struct header *hdi,struct header *hdo,int lower,int upper,int rflag,int alpha,int beta,int gamma,int delta);
int h_halftone_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,int lower,int upper,int rflag,int alpha,int beta,int gamma,int delta);
/******************* hardthre.c ***********************/
int h_hardthresh(struct header *hdi,struct header *hdo,union pixelval *thresh);
int h_hardthresh_b(struct header *hdi,struct header *hdo,union pixelval *thresh);
int h_hardthresh_i(struct header *hdi,struct header *hdo,union pixelval *thresh);
int h_hardthresh_f(struct header *hdi,struct header *hdo,union pixelval *thresh);
int h_hardthresh_c(struct header *hdi,struct header *hdo,union pixelval *thresh);
int h_hardthresh_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,unsigned char thresh);
int h_hardthresh_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo,int thresh);
int h_hardthresh_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double thresh);
int h_hardthresh_C(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double thresh);
/******************* histo.c ***********************/
int h_histo(struct header *hd,struct hips_histo *histogram,int nzflag,int *count);
int h_histo_b(struct header *hd,struct hips_histo *histogram,int nzflag,int *count);
int h_histo_sb(struct header *hd,struct hips_histo *histogram,int nzflag,int *count);
int h_histo_s(struct header *hd,struct hips_histo *histogram,int nzflag,int *count);
int h_histo_us(struct header *hd,struct hips_histo *histogram,int nzflag,int *count);
int h_histo_i(struct header *hd,struct hips_histo *histogram,int nzflag,int *count);
int h_histo_ui(struct header *hd,struct hips_histo *histogram,int nzflag,int *count);
int h_histo_f(struct header *hd,struct hips_histo *histogram,int nzflag,int *count);
int h_histo_d(struct header *hd,struct hips_histo *histogram,int nzflag,int *count);
int h_histo_c(struct header *hd,struct hips_histo *histogram,int nzflag,int *count);
int h_histo_dc(struct header *hd,struct hips_histo *histogram,int nzflag,int *count);
int h_histo_B(unsigned char *image,int nr,int nc,int nlp,int nbins,int *histo,unsigned char min,unsigned char width,int nzflag,int *count);
int h_histo_SB(char *image,int nr,int nc,int nlp,int nbins,int *histo,char min,char width,int nzflag,int *count);
int h_histo_S(short *image,int nr,int nc,int nlp,int nbins,int *histo,short min,short width,int nzflag,int *count);
int h_histo_US(unsigned short *image,int nr,int nc,int nlp,int nbins,int *histo,unsigned short min,unsigned short width,int nzflag,int *count);
int h_histo_I(int *image,int nr,int nc,int nlp,int nbins,int *histo,int min,int width,int nzflag,int *count);
int h_histo_UI(unsigned int *image,int nr,int nc,int nlp,int nbins,int *histo,unsigned int min,unsigned int width,int nzflag,int *count);
int h_histo_F(float *image,int nr,int nc,int nlp,int nbins,int *histo,double min,double width,int nzflag,int *count);
int h_histo_D(double *image,int nr,int nc,int nlp,int nbins,int *histo,double min,double width,int nzflag,int *count);
int h_histo_C(float *image,int nr,int nc,int nlp,int nbins,int *histo,double min,double width,int nzflag,int *count);
int h_histo_DC(double *image,int nr,int nc,int nlp,int nbins,int *histo,double min,double width,int nzflag,int *count);
/******************* histoeq.c ***********************/
int h_histoeq(struct hips_histo *histogram,int count,unsigned char *map);
int h_histoeq_b(struct hips_histo *histogram,int count,unsigned char *map);
int h_histoeq_B(int nbins,int *histo,int count,unsigned char *map);
/******************* ienlarge.c ***********************/
int h_ienlarge3(struct header *hdi1,struct header *hdi2,struct header *hdo,int xf,int yf,int tf,int t);
int h_ienlarge3_b(struct header *hdi1,struct header *hdi2,struct header *hdo,int xf,int yf,int tf,int t);
int h_ienlarge3_i(struct header *hdi1,struct header *hdi2,struct header *hdo,int xf,int yf,int tf,int t);
int h_ienlarge3_f(struct header *hdi1,struct header *hdi2,struct header *hdo,int xf,int yf,int tf,int t);
int h_ienlarge3_c(struct header *hdi1,struct header *hdi2,struct header *hdo,int xf,int yf,int tf,int t);
int h_ienlarge3_B(unsigned char *imagei1,unsigned char *imagei2,unsigned char *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo,int xf,int yf,int tf,int t);
int h_ienlarge3_I(int *imagei1,int *imagei2,int *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo,int xf,int yf,int tf,int t);
int h_ienlarge3_F(float *imagei1,float *imagei2,float *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo,int xf,int yf,int tf,int t);
int h_ienlarge3_C(float *imagei1,float *imagei2,float *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo,int xf,int yf,int tf,int t);
int ienl3_alloc(int xf,int yf,int tf);
/******************* invert.c ***********************/
int h_invert(struct header *hdi,struct header *hdo);
int h_invert_lp(struct header *hdi,struct header *hdo);
int h_invert_mp(struct header *hdi,struct header *hdo);
int h_invert_b(struct header *hdi,struct header *hdo);
int h_invert_i(struct header *hdi,struct header *hdo);
int h_invert_f(struct header *hdi,struct header *hdo);
int h_invert_MP(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo);
int h_invert_LP(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo);
int h_invert_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo);
int h_invert_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo);
int h_invert_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo);
/******************* invfourt.c ***********************/
int h_invfourtr(struct header *hd);
int h_invfourtr_c(struct header *hd);
int h_invfourtr_dc(struct header *hd);
int h_invfourtr_C(float (*image)[2],int nr,int nc,int nlp);
int h_invfourtr_DC(double (*image)[2],int nr,int nc,int nlp);
/******************* linscale.c ***********************/
int h_linscale(struct header *hdi,struct header *hdo,float b,float c);
int h_linscale_b(struct header *hdi,struct header *hdo,float b,float c);
int h_linscale_s(struct header *hdi,struct header *hdo,float b,float c);
int h_linscale_i(struct header *hdi,struct header *hdo,float b,float c);
int h_linscale_f(struct header *hdi,struct header *hdo,float b,float c);
int h_linscale_B(unsigned char *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,float b,float c);
int h_linscale_Bs(unsigned char *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,float b,float c);
int h_linscale_Bf(unsigned char *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,float b,float c);
int h_linscale_S(short *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,float b,float c);
int h_linscale_Ss(short *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,float b,float c);
int h_linscale_Sf(short *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,float b,float c);
int h_linscale_I(int *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,float b,float c);
int h_linscale_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,float b,float c);
/******************* log.c ***********************/
int h_log(struct header *hdi,struct header *hdo,double offset);
int h_log_b(struct header *hdi,struct header *hdo,double offset);
int h_log_s(struct header *hdi,struct header *hdo,double offset);
int h_log_i(struct header *hdi,struct header *hdo,double offset);
int h_log_f(struct header *hdi,struct header *hdo,double offset);
int h_log_B(unsigned char *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double offset);
int h_log_Bs(unsigned char *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double offset);
int h_log_Bf(unsigned char *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double offset);
int h_log_S(short *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double offset);
int h_log_Ss(short *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double offset);
int h_log_Sf(short *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double offset);
int h_log_I(int *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double offset);
int h_log_Is(int *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double offset);
int h_log_If(int *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double offset);
int h_log_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double offset);
/******************* mask.c ***********************/
int h_mask(struct header *hdi,struct hips_mask *mask,struct header *hdo);
int h_mask_bif(struct header *hdi,struct hips_mask *mask,struct header *hdo);
int h_mask_bff(struct header *hdi,struct hips_mask *mask,struct header *hdo);
int h_mask_iif(struct header *hdi,struct hips_mask *mask,struct header *hdo);
int h_mask_iff(struct header *hdi,struct hips_mask *mask,struct header *hdo);
int h_mask_fff(struct header *hdi,struct hips_mask *mask,struct header *hdo);
int h_mask_BIF(unsigned char *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,int nmasks,int maskfunc,int **masks ,int *mrows,int *mcols,int *mrowoff,int *mcoloff);
int h_mask_BFF(unsigned char *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,int nmasks,int maskfunc,float **masks ,int *mrows,int *mcols,int *mrowoff,int *mcoloff);
int h_mask_IIF(int *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,int nmasks,int maskfunc,int **masks ,int *mrows,int *mcols,int *mrowoff,int *mcoloff);
int h_mask_IFF(int *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,int nmasks,int maskfunc,float **masks ,int *mrows,int *mcols,int *mrowoff,int *mcoloff);
int h_mask_FFF(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,int nmasks,int maskfunc,float **masks ,int *mrows,int *mcols,int *mrowoff,int *mcoloff);
int h_mask_alloc(int nr,int nc,int nlpi,int *mrows,int *mcols,int *mrowoff,int *mcoloff);
int h_mask_valloc(int format);
float h_mask_value_i(void);
float h_mask_value_f(void);
/******************* max.c ***********************/
int h_max(struct header *hd,union pixelval *maxval,int nzflag);
int h_max_b(struct header *hd,union pixelval *maxval,int nzflag);
int h_max_sb(struct header *hd,union pixelval *maxval,int nzflag);
int h_max_s(struct header *hd,union pixelval *maxval,int nzflag);
int h_max_us(struct header *hd,union pixelval *maxval,int nzflag);
int h_max_i(struct header *hd,union pixelval *maxval,int nzflag);
int h_max_ui(struct header *hd,union pixelval *maxval,int nzflag);
int h_max_f(struct header *hd,union pixelval *maxval,int nzflag);
int h_max_d(struct header *hd,union pixelval *maxval,int nzflag);
unsigned char h_max_B(unsigned char *image,int nr,int nc,int nlp,int nzflag);
char h_max_SB(char *image,int nr,int nc,int nlp,int nzflag);
short h_max_S(short *image,int nr,int nc,int nlp,int nzflag);
unsigned short h_max_US(unsigned short *image,int nr,int nc,int nlp,int nzflag);
int h_max_I(int *image,int nr,int nc,int nlp,int nzflag);
unsigned int h_max_UI(unsigned int *image,int nr,int nc,int nlp,int nzflag);
float h_max_F(float *image,int nr,int nc,int nlp,int nzflag);
double h_max_D(double *image,int nr,int nc,int nlp,int nzflag);
/******************* maxabsp.c ***********************/
int h_maxabsp(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_maxabsp_s(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_maxabsp_i(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_maxabsp_f(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_maxabsp_d(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_maxabsp_ip(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_maxabsp_fp(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_maxabsp_S(short *imagei1,short *imagei2,short *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_maxabsp_I(int *imagei1,int *imagei2,int *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_maxabsp_F(float *imagei1,float *imagei2,float *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_maxabsp_D(double *imagei1,double *imagei2,double *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
/******************* maxhisto.c ***********************/
int h_maxhisto(struct hips_histo *histo);
int h_Maxhisto(int *histo,int nbins);
/******************* maxp.c ***********************/
int h_maxp(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_maxp_b(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_maxp_s(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_maxp_i(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_maxp_f(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_maxp_d(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_maxp_ip(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_maxp_fp(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_maxp_B(unsigned char *imagei1,unsigned char *imagei2,unsigned char *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_maxp_S(short *imagei1,short *imagei2,short *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_maxp_I(int *imagei1,int *imagei2,int *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_maxp_F(float *imagei1,float *imagei2,float *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_maxp_D(double *imagei1,double *imagei2,double *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
/******************* mean.c ***********************/
int h_mean(struct header *hd,float *mean,int nzflag);
int h_mean_f(struct header *hd,float *mean,int nzflag);
float h_mean_F(float *image,int nr,int nc,int nlp,int nzflag);
/******************* median.c ***********************/
int h_median(struct header *hdi,struct header *hdo,int size);
int h_median_b(struct header *hdi,struct header *hdo,int size);
int h_median_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,int size);
int sselect(int k,int *lo,int *hi);
/******************* minabsp.c ***********************/
int h_minabsp(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_minabsp_s(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_minabsp_i(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_minabsp_f(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_minabsp_d(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_minabsp_ip(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_minabsp_fp(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_minabsp_S(short *imagei1,short *imagei2,short *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_minabsp_I(int *imagei1,int *imagei2,int *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_minabsp_F(float *imagei1,float *imagei2,float *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_minabsp_D(double *imagei1,double *imagei2,double *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
/******************* minmax.c ***********************/
int h_minmax(struct header *hd,union pixelval *minval,union pixelval *maxval,int nzflag);
int h_minmax_b(struct header *hd,union pixelval *minval,union pixelval *maxval,int nzflag);
int h_minmax_sb(struct header *hd,union pixelval *minval,union pixelval *maxval,int nzflag);
int h_minmax_s(struct header *hd,union pixelval *minval,union pixelval *maxval,int nzflag);
int h_minmax_us(struct header *hd,union pixelval *minval,union pixelval *maxval,int nzflag);
int h_minmax_i(struct header *hd,union pixelval *minval,union pixelval *maxval,int nzflag);
int h_minmax_ui(struct header *hd,union pixelval *minval,union pixelval *maxval,int nzflag);
int h_minmax_f(struct header *hd,union pixelval *minval,union pixelval *maxval,int nzflag);
int h_minmax_d(struct header *hd,union pixelval *minval,union pixelval *maxval,int nzflag);
int h_minmax_c(struct header *hd,union pixelval *minval,union pixelval *maxval,int nzflag);
int h_minmax_dc(struct header *hd,union pixelval *minval,union pixelval *maxval,int nzflag);
int h_minmax_B(unsigned char *image,int nr,int nc,int nlp,unsigned char *rmin,unsigned char *rmax,int nzflag);
int h_minmax_SB(char *image,int nr,int nc,int nlp,char *rmin,char *rmax,int nzflag);
int h_minmax_S(short *image,int nr,int nc,int nlp,short *rmin,short *rmax,int nzflag);
int h_minmax_US(unsigned short *image,int nr,int nc,int nlp,unsigned short *rmin,unsigned short *rmax,int nzflag);
int h_minmax_I(int *image,int nr,int nc,int nlp,int *rmin,int *rmax,int nzflag);
int h_minmax_UI(unsigned int *image,int nr,int nc,int nlp,unsigned int *rmin,unsigned int *rmax,int nzflag);
int h_minmax_F(float *image,int nr,int nc,int nlp,float *rmin,float *rmax,int nzflag);
int h_minmax_D(double *image,int nr,int nc,int nlp,double *rmin,double *rmax,int nzflag);
int h_minmax_C(float *image,int nr,int nc,int nlp,float *rmin,float *rmax,int nzflag);
int h_minmax_DC(double *image,int nr,int nc,int nlp,double *rmin,double *rmax,int nzflag);
/******************* minp.c ***********************/
int h_minp(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_minp_b(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_minp_s(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_minp_i(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_minp_f(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_minp_d(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_minp_ip(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_minp_fp(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_minp_B(unsigned char *imagei1,unsigned char *imagei2,unsigned char *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_minp_S(short *imagei1,short *imagei2,short *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_minp_I(int *imagei1,int *imagei2,int *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_minp_F(float *imagei1,float *imagei2,float *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_minp_D(double *imagei1,double *imagei2,double *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
/******************* minroi.c ***********************/
int h_minroi(struct header *hd,int gray);
int h_minroi_b(struct header *hd,int gray);
int h_minroi_B(unsigned char *image,int nr,int nc,int nlp,int *frow,int *fcol,int *lrow,int *lcol,int gray);
/******************* morphdil.c ***********************/
int h_morphdil(struct header *hdi,struct header *hde,struct header *hdo,int centerr,int centerc,int gray);
int h_morphdil_b(struct header *hdi,struct header *hde,struct header *hdo,int centerr,int centerc,int gray);
int h_morphdil_B(unsigned char *imagei,unsigned char *imagee,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,int nre,int nce,int nlpe,int centerr,int centerc,int gray);
/******************* morphero.c ***********************/
int h_morphero(struct header *hdi,struct header *hde,struct header *hdt,struct header *hdo,int centerr,int centerc,int gray);
int h_morphero_b(struct header *hdi,struct header *hde,struct header *hdt,struct header *hdo,int centerr,int centerc,int gray);
int h_morphero_B(unsigned char *imagei,unsigned char *imagee,int *imaget,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,int nlpt,int nre,int nce,int nlpe,int centerr,int centerc,int gray);
/******************* mul.c ***********************/
int h_mul(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_mul_b(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_mul_s(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_mul_i(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_mul_f(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_mul_fc(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_mul_d(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_mul_ddc(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_mul_c(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_mul_dc(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_mul_ip(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_mul_fp(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_mul_B(unsigned char *imagei1,unsigned char *imagei2,unsigned char *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_mul_S(short *imagei1,short *imagei2,short *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_mul_I(int *imagei1,int *imagei2,int *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_mul_F(float *imagei1,float *imagei2,float *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_mul_FC(float *imagei1,float *imagei2,float *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_mul_D(double *imagei1,double *imagei2,double *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_mul_DDC(double *imagei1,double *imagei2,double *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_mul_C(float *imagei1,float *imagei2,float *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_mul_DC(double *imagei1,double *imagei2,double *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
/******************* mulscale.c ***********************/
int h_mulscale(struct header *hdi,struct header *hdo,union pixelval *b);
int h_mulscale_b(struct header *hdi,struct header *hdo,union pixelval *b);
int h_mulscale_s(struct header *hdi,struct header *hdo,union pixelval *b);
int h_mulscale_i(struct header *hdi,struct header *hdo,union pixelval *b);
int h_mulscale_f(struct header *hdi,struct header *hdo,union pixelval *b);
int h_mulscale_d(struct header *hdi,struct header *hdo,union pixelval *b);
int h_mulscale_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,unsigned char b);
int h_mulscale_Bs(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,unsigned char b);
int h_mulscale_Bf(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,unsigned char b);
int h_mulscale_S(short *imagei,short *imageo,int nr,int nc,int nlpi,int nlpo,short b);
int h_mulscale_Ss(short *imagei,short *imageo,int nr,int nc,int nlpi,int nlpo,short b);
int h_mulscale_Sf(short *imagei,short *imageo,int nr,int nc,int nlpi,int nlpo,short b);
int h_mulscale_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo,int b);
int h_mulscale_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double b);
int h_mulscale_D(double *imagei,double *imageo,int nr,int nc,int nlpi,int nlpo,double b);
/******************* neg.c ***********************/
int h_neg(struct header *hdi,struct header *hdo);
int h_neg_mp(struct header *hdi,struct header *hdo);
int h_neg_lp(struct header *hdi,struct header *hdo);
int h_neg_b(struct header *hdi,struct header *hdo);
int h_neg_s(struct header *hdi,struct header *hdo);
int h_neg_i(struct header *hdi,struct header *hdo);
int h_neg_f(struct header *hdi,struct header *hdo);
int h_neg_MP(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo);
int h_neg_LP(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo);
int h_neg_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo);
int h_neg_S(short *imagei,short *imageo,int nr,int nc,int nlpi,int nlpo);
int h_neg_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo);
int h_neg_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo);
/******************* noise.c ***********************/
int h_noise(struct header *hdi,struct header *hdo,double p,int *counter,int bpp);
int h_noise_b(struct header *hdi,struct header *hdo,double p,int *counter,int bpp);
int h_noise_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,double p,int *counter,int bpp);
/******************* or.c ***********************/
int h_or(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_or_mp(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_or_lp(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_or_b(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_or_MP(unsigned char *imagei1,unsigned char *imagei2,unsigned char *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_or_LP(unsigned char *imagei1,unsigned char *imagei2,unsigned char *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_or_B(unsigned char *imagei1,unsigned char *imagei2,unsigned char *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
/******************* pixmap.c ***********************/
int h_pixmap(struct header *hdi,struct header *hdo,unsigned char *map);
int h_pixmap_b(struct header *hdi,struct header *hdo,unsigned char *map);
int h_pixmap_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,unsigned char *map);
/******************* power.c ***********************/
int h_power(struct header *hdi,struct header *hdo,double power);
int h_power_b(struct header *hdi,struct header *hdo,double power);
int h_power_s(struct header *hdi,struct header *hdo,double power);
int h_power_i(struct header *hdi,struct header *hdo,double power);
int h_power_f(struct header *hdi,struct header *hdo,double power);
int h_power_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,double power);
int h_power_Bs(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,double power);
int h_power_Bf(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,double power);
int h_power_S(short *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double power);
int h_power_Ss(short *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double power);
int h_power_Sf(short *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double power);
int h_power_I(int *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double power);
int h_power_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double power);
/******************* pyrdisp.c ***********************/
#if 0
int h_pyrdisp(struct FIMAGEtag *pyr,int botlev,int toplev,struct header *hdo,int cflag,int margin);
int h_pyrdisp_i(struct IIMAGEtag *pyr,int botlev,int toplev,struct header *hdo,int cflag,int margin);
int h_pyrdisp_f(struct FIMAGEtag *pyr,int botlev,int toplev,struct header *hdo,int cflag,int margin);
int h_pyrdisp_I(struct IIMAGEtag *pyr,int botlev,int toplev,int *imageo,int nlpo,int cflag,int margin);
int h_pyrdisp_F(struct FIMAGEtag *pyr,int botlev,int toplev,float *imageo,int nlpo,int cflag,int margin);
#endif
/******************* quadscal.c ***********************/
int h_quadscale(struct header *hdi,struct header *hdo,double a,double b,double c);
int h_quadscale_b(struct header *hdi,struct header *hdo,double a,double b,double c);
int h_quadscale_s(struct header *hdi,struct header *hdo,double a,double b,double c);
int h_quadscale_i(struct header *hdi,struct header *hdo,double a,double b,double c);
int h_quadscale_f(struct header *hdi,struct header *hdo,double a,double b,double c);
int h_quadscale_B(unsigned char *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double a,double b,double c);
int h_quadscale_Bs(unsigned char *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double a,double b,double c);
int h_quadscale_Bf(unsigned char *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double a,double b,double c);
int h_quadscale_S(short *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double a,double b,double c);
int h_quadscale_Ss(short *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double a,double b,double c);
int h_quadscale_Sf(short *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double a,double b,double c);
int h_quadscale_I(int *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double a,double b,double c);
int h_quadscale_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double a,double b,double c);
/******************* reduce.c ***********************/
int h_reduce(struct header *hdi,struct header *hdo,int xf,int yf);
int h_reduce_bi(struct header *hdi,struct header *hdo,int xf,int yf);
int h_reduce_s(struct header *hdi,struct header *hdo,int xf,int yf);
int h_reduce_i(struct header *hdi,struct header *hdo,int xf,int yf);
int h_reduce_f(struct header *hdi,struct header *hdo,int xf,int yf);
int h_reduce_c(struct header *hdi,struct header *hdo,int xf,int yf);
int h_reduce_BI(unsigned char *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo,int xf,int yf);
int h_reduce_S(short *imagei,short *imageo,int nr,int nc,int nlpi,int nlpo,int xf,int yf);
int h_reduce_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo,int xf,int yf);
int h_reduce_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,int xf,int yf);
int h_reduce_C(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,int xf,int yf);
/******************* reflect.c ***********************/
int h_reflect(struct header *hdi,struct header *hdo);
int h_reflect_b(struct header *hdi,struct header *hdo);
int h_reflect_i(struct header *hdi,struct header *hdo);
int h_reflect_f(struct header *hdi,struct header *hdo);
int h_reflect_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo);
int h_reflect_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo);
int h_reflect_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo);
/******************* rot180.c ***********************/
int h_rot180(struct header *hdi,struct header *hdo);
int h_rot180_b(struct header *hdi,struct header *hdo);
int h_rot180_i(struct header *hdi,struct header *hdo);
int h_rot180_f(struct header *hdi,struct header *hdo);
int h_rot180_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo);
int h_rot180_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo);
int h_rot180_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo);
/******************* rot90.c ***********************/
int h_rot90(struct header *hdi,struct header *hdo,int dirflag);
int h_rot90_b(struct header *hdi,struct header *hdo,int dirflag);
int h_rot90_i(struct header *hdi,struct header *hdo,int dirflag);
int h_rot90_f(struct header *hdi,struct header *hdo,int dirflag);
int h_rot90_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,int dirflag);
int h_rot90_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo,int dirflag);
int h_rot90_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,int dirflag);
/******************* sample.c ***********************/
int h_sample(struct header *hdi,struct header *hdo,int ratex,int ratey,int offsetx,int offsety);
int h_sample_b(struct header *hdi,struct header *hdo,int ratex,int ratey,int offsetx,int offsety);
int h_sample_s(struct header *hdi,struct header *hdo,int ratex,int ratey,int offsetx,int offsety);
int h_sample_i(struct header *hdi,struct header *hdo,int ratex,int ratey,int offsetx,int offsety);
int h_sample_f(struct header *hdi,struct header *hdo,int ratex,int ratey,int offsetx,int offsety);
int h_sample_d(struct header *hdi,struct header *hdo,int ratex,int ratey,int offsetx,int offsety);
int h_sample_c(struct header *hdi,struct header *hdo,int ratex,int ratey,int offsetx,int offsety);
int h_sample_dc(struct header *hdi,struct header *hdo,int ratex,int ratey,int offsetx,int offsety);
int h_sample_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,int ratex,int ratey,int offsetx,int offsety);
int h_sample_S(short *imagei,short *imageo,int nr,int nc,int nlpi,int nlpo,int ratex,int ratey,int offsetx,int offsety);
int h_sample_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo,int ratex,int ratey,int offsetx,int offsety);
int h_sample_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,int ratex,int ratey,int offsetx,int offsety);
int h_sample_D(double *imagei,double *imageo,int nr,int nc,int nlpi,int nlpo,int ratex,int ratey,int offsetx,int offsety);
int h_sample_C(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,int ratex,int ratey,int offsetx,int offsety);
int h_sample_DC(double *imagei,double *imageo,int nr,int nc,int nlpi,int nlpo,int ratex,int ratey,int offsetx,int offsety);
/******************* scaleadd.c ***********************/
int h_scaleadd(struct header *hdi,struct header *hdo,double s);
int h_scaleadd_f(struct header *hdi,struct header *hdo,double s);
int h_scaleadd_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double s);
/******************* sepconv.c ***********************/
int h_sepconv(struct header *hdi,struct header *hdt,struct header *hdo,int *maskh,int nmaskh,int offseth,int *maskv,int nmaskv,int offsetv);
int h_sepconv_i(struct header *hdi,struct header *hdt,struct header *hdo,int *maskh,int nmaskh,int offseth,int *maskv,int nmaskv,int offsetv);
int h_sepconv_f(struct header *hdi,struct header *hdt,struct header *hdo,float *maskh,int nmaskh,int offseth,float *maskv,int nmaskv,int offsetv);
/******************* setimage.c ***********************/
int h_setimage(struct header *hd,union pixelval *val);
int h_setimage_mp(struct header *hd,union pixelval *val);
int h_setimage_lp(struct header *hd,union pixelval *val);
int h_setimage_b(struct header *hd,union pixelval *val);
int h_setimage_sb(struct header *hd,union pixelval *val);
int h_setimage_s(struct header *hd,union pixelval *val);
int h_setimage_us(struct header *hd,union pixelval *val);
int h_setimage_i(struct header *hd,union pixelval *val);
int h_setimage_ui(struct header *hd,union pixelval *val);
int h_setimage_f(struct header *hd,union pixelval *val);
int h_setimage_d(struct header *hd,union pixelval *val);
int h_setimage_c(struct header *hd,union pixelval *val);
int h_setimage_dc(struct header *hd,union pixelval *val);
int h_setimage_MP(unsigned char *image,int nr,int nc,int nlp,unsigned char val);
int h_setimage_LP(unsigned char *image,int nr,int nc,int nlp,unsigned char val);
int h_setimage_B(unsigned char *image,int nr,int nc,int nlp,unsigned char val);
int h_setimage_SB(char *image,int nr,int nc,int nlp,char val);
int h_setimage_S(short *image,int nr,int nc,int nlp,short val);
int h_setimage_US(unsigned short *image,int nr,int nc,int nlp,unsigned short val);
int h_setimage_I(int *image,int nr,int nc,int nlp,int val);
int h_setimage_UI(unsigned int *image,int nr,int nc,int nlp,unsigned int val);
int h_setimage_F(float *image,int nr,int nc,int nlp,double val);
int h_setimage_D(double *image,int nr,int nc,int nlp,double val);
int h_setimage_C(float *image,int nr,int nc,int nlp,float *val);
int h_setimage_DC(double *image,int nr,int nc,int nlp,double *val);
/******************* shift.c ***********************/
int h_shift(struct header *hdi,struct header *hdo,int shift);
int h_shift_b(struct header *hdi,struct header *hdo,int shift);
int h_shift_i(struct header *hdi,struct header *hdo,int shift);
int h_shift_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,int shift);
int h_shift_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo,int shift);
/******************* shufflea.c ***********************/
int h_shuffleadd(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_shuffleadd_bsb(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_shuffleadd_f(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_shuffleadd_BSB(unsigned char *imagei1,short *imagei2,unsigned char *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_shuffleadd_F(float *imagei1,float *imagei2,float *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int shuffle(int *p,int n);
/******************* slice.c ***********************/
int h_slice(struct header *hdi,struct header *hdo,int rowcol,int vflag);
int h_slice_b(struct header *hdi,struct header *hdo,int rowcol,int vflag);
int h_slice_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,int rowcol,int vflag);
/******************* softthre.c ***********************/
int h_softthresh(struct header *hdi,struct header *hdo,union pixelval *thresh);
int h_softthresh_b(struct header *hdi,struct header *hdo,union pixelval *thresh);
int h_softthresh_i(struct header *hdi,struct header *hdo,union pixelval *thresh);
int h_softthresh_f(struct header *hdi,struct header *hdo,union pixelval *thresh);
int h_softthresh_c(struct header *hdi,struct header *hdo,union pixelval *thresh);
int h_softthresh_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,unsigned char thresh);
int h_softthresh_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo,int thresh);
int h_softthresh_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double thresh);
int h_softthresh_C(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double thresh);
/******************* stats.c ***********************/
int h_stats(struct header *hd,struct hips_stats *stats,int nzflag);
int h_stats_b(struct header *hd,struct hips_stats *stats,int nzflag);
int h_stats_f(struct header *hd,struct hips_stats *stats,int nzflag);
int h_stats_B(unsigned char *image,int nr,int nc,int nlp,int *rnelem,unsigned char *rmin,unsigned char *rmax,double *rsum,double *rssq,int nzflag);
int h_stats_F(float *image,int nr,int nc,int nlp,int *rnelem,float *rmin,float *rmax,double *rsum,double *rssq,int nzflag);
/******************* stretch.c ***********************/
int h_stretch(struct header *hdi,struct header *hdo,double d,double expt1,double expt2,union pixelval *mval);
int h_stretch_b(struct header *hdi,struct header *hdo,double d,double expt1,double expt2,union pixelval *mval);
int h_stretch_s(struct header *hdi,struct header *hdo,double d,double expt1,double expt2,union pixelval *mval);
int h_stretch_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,double d,double expt1,double expt2,unsigned char mval);
int h_stretch_S(short *imagei,short *imageo,int nr,int nc,int nlpi,int nlpo,double d,double expt1,double expt2,short mval);
int h_stretch_Sf(short *imagei,short *imageo,int nr,int nc,int nlpi,int nlpo,double d,double expt1,double expt2,short mval);
int h_stretch_Ss(short *imagei,short *imageo,int nr,int nc,int nlpi,int nlpo,double d,double expt1,double expt2,short mval);
/******************* stretchi.c ***********************/
int h_stretchimg(struct header *hdi,struct header *hdo);
int h_stretchimg_b(struct header *hdi,struct header *hdo);
int h_stretchimg_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nor,int noc,int nlpo);
/******************* thicken.c ***********************/
int h_thicken(struct header *hdi,struct header *hdo);
int h_thicken_b(struct header *hdi,struct header *hdo);
int h_thicken_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo);
/******************* thin.c ***********************/
int h_thin(struct header *hd,int *passflag,int sflag,int vflag,int f);
int h_thin_b(struct header *hd,int *passflag,int sflag,int vflag,int f);
int h_thin_B(unsigned char *image,int nr,int nc,int nlp,int *passflag,int sflag,int vflag,int f);
int addmult(int row,int col,int label);
int neighbor(int r,int c,int d);
/******************* translat.c ***********************/
int h_translate(struct header *hdi,struct header *hdo,int shiftx,int shifty,int palinflag);
int h_translate_b(struct header *hdi,struct header *hdo,int shiftx,int shifty,int palinflag);
int h_translate_s(struct header *hdi,struct header *hdo,int shiftx,int shifty,int palinflag);
int h_translate_i(struct header *hdi,struct header *hdo,int shiftx,int shifty,int palinflag);
int h_translate_f(struct header *hdi,struct header *hdo,int shiftx,int shifty,int palinflag);
int h_translate_d(struct header *hdi,struct header *hdo,int shiftx,int shifty,int palinflag);
int h_translate_c(struct header *hdi,struct header *hdo,int shiftx,int shifty,int palinflag);
int h_translate_dc(struct header *hdi,struct header *hdo,int shiftx,int shifty,int palinflag);
int h_translate_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo,int shiftx,int shifty,int palinflag);
int h_translate_S(short *imagei,short *imageo,int nr,int nc,int nlpi,int nlpo,int shiftx,int shifty,int palinflag);
int h_translate_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo,int shiftx,int shifty,int palinflag);
int h_translate_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,int shiftx,int shifty,int palinflag);
int h_translate_D(double *imagei,double *imageo,int nr,int nc,int nlpi,int nlpo,int shiftx,int shifty,int palinflag);
int h_translate_C(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,int shiftx,int shifty,int palinflag);
int h_translate_DC(double *imagei,double *imageo,int nr,int nc,int nlpi,int nlpo,int shiftx,int shifty,int palinflag);
/******************* transpos.c ***********************/
int h_transpose(struct header *hdi,struct header *hdo);
int h_transpose_b(struct header *hdi,struct header *hdo);
int h_transpose_i(struct header *hdi,struct header *hdo);
int h_transpose_f(struct header *hdi,struct header *hdo);
int h_transpose_B(unsigned char *imagei,unsigned char *imageo,int nr,int nc,int nlpi,int nlpo);
int h_transpose_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo);
int h_transpose_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo);
/******************* vconvolv.c ***********************/
int h_vconvolve(struct header *hdi,struct header *hdo,int *mask,int nmask,int offset);
int h_vconvolve_i(struct header *hdi,struct header *hdo,int *mask,int nmask,int offset);
int h_vconvolve_f(struct header *hdi,struct header *hdo,float *mask,int nmask,int offset);
int h_vconvolve_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo,int *mask,int nmask,int offset);
int h_vconvolve_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,float *mask,int nmask,int offset);
/******************* walshtr.c ***********************/
int h_walshtr(struct header *hd);
int h_walshtr_i(struct header *hd);
int h_walshtr_f(struct header *hd);
int h_walshtr_I(int *image,int nr,int nc);
int h_walshtr_F(float *image,int nr,int nc);
int h_fwt_i(int *vec,int loglen);
int h_fwt_f(float *vec,int loglen);
/******************* wgauss.c ***********************/
int h_wgauss(struct header *hdi,struct header *hdo,double rowmu,double colmu,double rowsigma,double colsigma,double factor);
int h_wgauss_f(struct header *hdi,struct header *hdo,double rowmu,double colmu,double rowsigma,double colsigma,double factor);
int h_wgauss_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double rowmu,double colmu,double rowsigma,double colsigma,double factor);
/******************* wtsum3.c ***********************/
int h_wtsum3(struct header *hdi1,struct header *hdi2,struct header *hdi3,struct header *hdo,double wt1,double wt2,double wt3);
int h_wtsum3_b(struct header *hdi1,struct header *hdi2,struct header *hdi3,struct header *hdo,double wt1,double wt2,double wt3);
int h_wtsum3_s(struct header *hdi1,struct header *hdi2,struct header *hdi3,struct header *hdo,double wt1,double wt2,double wt3);
int h_wtsum3_i(struct header *hdi1,struct header *hdi2,struct header *hdi3,struct header *hdo,double wt1,double wt2,double wt3);
int h_wtsum3_f(struct header *hdi1,struct header *hdi2,struct header *hdi3,struct header *hdo,double wt1,double wt2,double wt3);
int h_wtsum3_B(unsigned char *imagei1,unsigned char *imagei2,unsigned char *imagei3,unsigned char *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpi3,int nlpo,double wt1,double wt2,double wt3);
int h_wtsum3_S(short *imagei1,short *imagei2,short *imagei3,short *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpi3,int nlpo,double wt1,double wt2,double wt3);
int h_wtsum3_I(int *imagei1,int *imagei2,int *imagei3,int *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpi3,int nlpo,double wt1,double wt2,double wt3);
int h_wtsum3_F(float *imagei1,float *imagei2,float *imagei3,float *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpi3,int nlpo,double wt1,double wt2,double wt3);
/******************* xor.c ***********************/
int h_xor(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_xor_mp(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_xor_lp(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_xor_b(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_xor_MP(unsigned char *imagei1,unsigned char *imagei2,unsigned char *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_xor_LP(unsigned char *imagei1,unsigned char *imagei2,unsigned char *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
int h_xor_B(unsigned char *imagei1,unsigned char *imagei2,unsigned char *imageo,int nr,int nc,int nlpi1,int nlpi2,int nlpo);
/******************* zc.c ***********************/
int h_zc(struct header *hdi,struct header *hdo,double error,int nflag,int zflag);
int h_zc_f(struct header *hdi,struct header *hdo,double error,int nflag,int zflag);
int h_zc_F(float *imagei,float *imageo,int nr,int nc,int nlpi,int nlpo,double error,int nflag,int zflag);
/******************* invft3d.c ***********************/
int h_invfourtr3d(struct header *hd);
int h_invfourtr3d_c(struct header *hd);
int h_invfourtr3d_C(float (*image)[2],int nr,int nc,int nf);

int h_canny(struct header *Isrc, struct header *Idst, double sigma,
            int mask_size,double lfrac,double hfrac,int dothin);

int h_logmap(struct header *phd,struct header *phd_LM, LMtable *ptable);
int h_LMtable_uniform(LMtable *ptable,struct header *phd_w,int Irows,int Icols,
                      float a,int R,int Fx,int Fy,float disk_size) ;

#endif
