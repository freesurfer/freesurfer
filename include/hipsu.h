/**
 * @file  hipsu.h
 * @brief Header file for hips utility functions
 *
 */
/*
 * Original Author: Bruce Fischl
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


#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>

/**************** hipsu.h *********************/
/* Header file for hips utility functions */
/******************* add_head.c ***********************/
int history_set(struct header *hda,...);
int history_append(struct header *hda,...);
int history_indentadd(struct header *hd,char *s);
/******************* addplot.c ***********************/
int addpoint(char *buf,int index,int limit,double b,double x,double y,double z);
int addvec(char *buf,int index,int limit,double b,double x1,double y1,double z1,double x2,double y2,double z2);
int addend(char *buf,int index,int limit,double x,double y,double z);
/******************* getplot.c ***********************/
int getplot(char *buf,int index,int *op,double *b,double *x1,double *y1,double *z1,double *x2,double *y2,double *z2);
/******************* trans_fr.c ***********************/
int trans_frame(char *buf,int nbuf,double *shift_v,double (*rot_m)[3],int *flags);
/******************* alloc_hi.c ***********************/
int alloc_histo(struct hips_histo *histo,union pixelval *min,union pixelval *max,int nbins,int format);
int alloc_histobins(struct hips_histo *histo);
/******************* alloc_im.c ***********************/
int alloc_image(struct header *hd);
int alloc_imagez(struct header *hd);
int free_image(struct header *hd);
/******************* cmapio.c ***********************/
int readcmap(char *filename,int maxcount,int *count,unsigned char *r,unsigned char *g,unsigned char *b);
/******************* col1toco.c ***********************/
int h_col1torgb(struct header *hdi,struct header *hdo,int fr);
int h_col1torgbz(struct header *hdi,struct header *hdo,int fr);
int h_col1tozrgb(struct header *hdi,struct header *hdo,int fr);
int h_col1tobgr(struct header *hdi,struct header *hdo,int fr);
int h_col1tobgrz(struct header *hdi,struct header *hdo,int fr);
int h_col1tozbgr(struct header *hdi,struct header *hdo,int fr);
int h_btorgb_1(struct header *hdi,struct header *hdo,int fr);
int h_itorgb_1(struct header *hdi,struct header *hdo,int fr);
int h_btorgbz_1(struct header *hdi,struct header *hdo,int fr);
int h_itorgbz_1(struct header *hdi,struct header *hdo,int fr);
int h_btozrgb_1(struct header *hdi,struct header *hdo,int fr);
int h_itozrgb_1(struct header *hdi,struct header *hdo,int fr);
int h_btobgr_1(struct header *hdi,struct header *hdo,int fr);
int h_itobgr_1(struct header *hdi,struct header *hdo,int fr);
int h_btobgrz_1(struct header *hdi,struct header *hdo,int fr);
int h_itobgrz_1(struct header *hdi,struct header *hdo,int fr);
int h_btozbgr_1(struct header *hdi,struct header *hdo,int fr);
int h_itozbgr_1(struct header *hdi,struct header *hdo,int fr);
int h_1to3_b(unsigned char *pi,unsigned char *po,int np);
int h_1to3_i(int *pi,unsigned char *po,int np);
int h_1to4_b(unsigned char *pi,unsigned char *po,int np);
int h_1to4_i(int *pi,unsigned char *po,int np);
/******************* col3tob.c ***********************/
int h_col3tob(struct header *hdi,struct header *hdo,int fr);
int h_rgbtob_1(struct header *hdi,struct header *hdo,int fr);
int h_rgbztob_1(struct header *hdi,struct header *hdo,int fr);
int h_zrgbtob_1(struct header *hdi,struct header *hdo,int fr);
int h_bgrtob_1(struct header *hdi,struct header *hdo,int fr);
int h_bgrztob_1(struct header *hdi,struct header *hdo,int fr);
int h_zbgrtob_1(struct header *hdi,struct header *hdo,int fr);
int h_3to1_b(unsigned char *pi,unsigned char *po,int np);
int h_4to1_b(unsigned char *pi,unsigned char *po,int np);
/******************* col3toco.c ***********************/
int h_torgb(struct header *hdi,struct header *hdo);
int h_torgbz(struct header *hdi,struct header *hdo);
int h_tozrgb(struct header *hdi,struct header *hdo);
int h_tobgr(struct header *hdi,struct header *hdo);
int h_tobgrz(struct header *hdi,struct header *hdo);
int h_tozbgr(struct header *hdi,struct header *hdo);
int h_rgbztorgb(struct header *hdi,struct header *hdo);
int h_zrgbtorgb(struct header *hdi,struct header *hdo);
int h_bgrtorgb(struct header *hdi,struct header *hdo);
int h_bgrztorgb(struct header *hdi,struct header *hdo);
int h_zbgrtorgb(struct header *hdi,struct header *hdo);
int h_rgbtorgbz(struct header *hdi,struct header *hdo);
int h_zrgbtorgbz(struct header *hdi,struct header *hdo);
int h_bgrtorgbz(struct header *hdi,struct header *hdo);
int h_bgrztorgbz(struct header *hdi,struct header *hdo);
int h_zbgrtorgbz(struct header *hdi,struct header *hdo);
int h_rgbtozrgb(struct header *hdi,struct header *hdo);
int h_rgbztozrgb(struct header *hdi,struct header *hdo);
int h_bgrtozrgb(struct header *hdi,struct header *hdo);
int h_bgrztozrgb(struct header *hdi,struct header *hdo);
int h_zbgrtozrgb(struct header *hdi,struct header *hdo);
int h_rgbtobgr(struct header *hdi,struct header *hdo);
int h_rgbztobgr(struct header *hdi,struct header *hdo);
int h_zrgbtobgr(struct header *hdi,struct header *hdo);
int h_bgrztobgr(struct header *hdi,struct header *hdo);
int h_zbgrtobgr(struct header *hdi,struct header *hdo);
int h_rgbtobgrz(struct header *hdi,struct header *hdo);
int h_rgbztobgrz(struct header *hdi,struct header *hdo);
int h_zrgbtobgrz(struct header *hdi,struct header *hdo);
int h_bgrtobgrz(struct header *hdi,struct header *hdo);
int h_zbgrtobgrz(struct header *hdi,struct header *hdo);
int h_rgbtozbgr(struct header *hdi,struct header *hdo);
int h_rgbztozbgr(struct header *hdi,struct header *hdo);
int h_zrgbtozbgr(struct header *hdi,struct header *hdo);
int h_bgrtozbgr(struct header *hdi,struct header *hdo);
int h_bgrztozbgr(struct header *hdi,struct header *hdo);
int h_30to3(struct header *hdi,struct header *hdo);
int h_03to3(struct header *hdi,struct header *hdo);
int h_3to3_flip(struct header *hdi,struct header *hdo);
int h_30to3_flip(struct header *hdi,struct header *hdo);
int h_03to3_flip(struct header *hdi,struct header *hdo);
int h_3to30(struct header *hdi,struct header *hdo);
int h_03to30(struct header *hdi,struct header *hdo);
int h_3to30_flip(struct header *hdi,struct header *hdo);
int h_30to30_flip(struct header *hdi,struct header *hdo);
int h_03to30_flip(struct header *hdi,struct header *hdo);
int h_3to03(struct header *hdi,struct header *hdo);
int h_30to03(struct header *hdi,struct header *hdo);
int h_3to03_flip(struct header *hdi,struct header *hdo);
int h_30to03_flip(struct header *hdi,struct header *hdo);
int h_03to03_flip(struct header *hdi,struct header *hdo);
/******************* col3toi.c ***********************/
int h_col3toi(struct header *hdi,struct header *hdo,int fr);
int h_rgbtoi_1(struct header *hdi,struct header *hdo,int fr);
int h_rgbztoi_1(struct header *hdi,struct header *hdo,int fr);
int h_zrgbtoi_1(struct header *hdi,struct header *hdo,int fr);
int h_bgrtoi_1(struct header *hdi,struct header *hdo,int fr);
int h_bgrztoi_1(struct header *hdi,struct header *hdo,int fr);
int h_zbgrtoi_1(struct header *hdi,struct header *hdo,int fr);
int h_3to1_i(unsigned char *pi,int *po,int np);
int h_4to1_i(unsigned char *pi,int *po,int np);
/******************* conversi.c ***********************/
int find_closest(struct header *hd,int *typeslist);
int ffind_closest(struct header *hd,int *typeslist,char *fname);
int pfind_closest(int pfmt,int *typeslist,char *fname);
int in_typeslist(int type,int *typeslist);
int find_method(int typein,int typeout);
int ffind_method(int typein,int typeout,char *fname);
int set_conversion(struct header *hd,struct header *hdp,int *typeslist);
int fset_conversion(struct header *hd,struct header *hdp,int *typeslist,char *fname);
int pset_conversion(struct header *hd,struct header *hdp,int ptype,char *fname);
int hconvert(struct header *hd,struct header *hdp,int method,int fr);
int fhconvert(struct header *hd,struct header *hdp,int method,int fr,char *fname);
int hconvertback(struct header *hd,struct header *hdp,int method,int fr);
int setupconvback(struct header *hd,struct header *hdp,struct header *hdcb);
int read_imagec(struct header *hd,struct header *hdp,int method,int fr);
int fread_imagec(FILE *fp,struct header *hd,struct header *hdp,int method,int fr,char *fname);
int write_imagec(struct header *hd,struct header *hdp,int method,int flag,int fr);
int fwrite_imagec(FILE *fp,struct header *hd,struct header *hdp,int method,int flag,int fr,char *fname);
/******************* cut_fram.c ***********************/
int cut_frame(char *inbuf,int nbytes,char *outbuf,int limit,double x0,double y0,double xn,double yn);
/******************* desc_set.c ***********************/
int desc_set(struct header *hd,char *s);
int desc_append(struct header *hd,char *s);
int desc_set2(struct header *hda,...);
int desc_append2(struct header *hda,...);
int desc_indentadd(struct header *hd,char *s);
/******************* dup_head.c ***********************/
int dup_header(struct header *hd1,struct header *hd2);
int dup_headern(struct header *hd1,struct header *hd2);
/******************* formathe.c ***********************/
char *formatheader(struct header *h);
char *formatheadera(struct header *h);
char *formatheaderc(struct header *h,int aflag);
int adddec(int i);
int addstr(char *s);
/******************* fread_ol.c ***********************/
int fread_oldhdr(FILE *fp,struct header *hd,char *firsts,char *fname);
//int getline(FILE *fp,char **s ,int *l);
int swallownl(FILE *fp);
int hfgets(char *s,int n,FILE *fp);
/******************* free_hea.c ***********************/
int free_header(struct header *hd);
int free_hdrcon(struct header *hd);
/******************* futils.c ***********************/
FILE *ffopen(char *name,char *mode);
FILE *ffreopen(char *name,char *mode,FILE *stream1);
int ffread(char *ptr,int size,int nelem,FILE *stream);
int ffwrite(char *ptr,int size,int nelem,FILE *stream);
FILE *hfopenr(char *filename);
/******************* halloc.c ***********************/
unsigned char *halloc(int i,int j);
unsigned char *hmalloc(unsigned long i);
char *memalloc(int nelem,unsigned long elsize);
/******************* hdepth.c ***********************/
int hgetdepth(struct header *hd);
int hsetdepth(struct header *hd,int depth);
/******************* herrs.c ***********************/
/******************* hformatn.c ***********************/
char *hformatname(int pfmt);
char *hformatname_f(int pfmtf,int pfmtt);
char *hformatname_t(int pfmtf,int pfmtt);
/******************* hips_ext.c ***********************/
/******************* hsizepix.c ***********************/
unsigned long hsizepix(int pfmt);
unsigned long hbitsperpixel(int pfmt);
/******************* htypes.c ***********************/
/******************* init_hea.c ***********************/
int init_header(struct header *hd,char *onm,char *snm,int nfr,char *odt,int rw,int cl,int pfmt,int nc,char *desc);
int init_hdr_alloc(struct header *hd,char *onm,char *snm,int nfr,char *odt,int rw,int cl,int pfmt,int nc,char *desc);
int init_header_d(struct header *hd,char *onm,int onmd,char *snm,int snmd,int nfr,char *odt,int odtd,int rw,int cl,int pfmt,int nc,char *desc,int descd);
int init_hdr_alloc_d(struct header *hd,char *onm,int onmd,char *snm,int snmd,int nfr,char *odt,int odtd,int rw,int cl,int pfmt,int nc,char *desc,int descd);
/******************* maxforma.c ***********************/
int hformatlevel(int pfmt);
int maxformat(int pfmt1,int pfmt2,int *typeslist,char *fname1,char *fname2);
/******************* pa_defau.c ***********************/
int set_defaults(void);
int set_flag_defaults(void);
int set_filename_defaults(void);
/******************* pa_exter.c ***********************/
/******************* pa_main.c ***********************/
int parseargs(char *va_alist);
int parseargsu(char *va_alist);
/******************* pa_parse.c ***********************/
int parse_command_line(int argc,char **argv );
int accept_flag(int *arg_count,char ***arg_list  );
int accept_filter_flag(int *arg_count,char ***arg_list  ,int *is_filter_flag);
int accept_filter_parameters(struct flag_key *flag_ptr,int *arg_count,char ***arg_list  );
int accept_parameter_value(union generic_ptr *parameter_ptr,int parameter_type,char *value,int *is_valid);
int accept_parameter_list(union generic_ptr *parameter_ptr,char ***arglist  ,int *argcount,int *is_valid);
void disable_mutex_flags(char **flag_list );
int accept_standard_flag(int *arg_count,char ***arg_list  ,int *is_standard_flag);
int accept_filenames(int *arg_count,char ***arg_list  );
int accept_file_value(char **filename_ptr ,char *value,int *accept_stdin,int *is_valid_filename);
/******************* pa_table.c ***********************/
int get_flag_data(char **varguments );
int count_flags(struct flag_format *flag_list);
void add_flag_to_table(struct flag_format *flag_ptr,struct flag_key *table_ptr,char **varguments );
struct flag_key *find_flag(char *value);
void unlock_all_flags(void);
void lock_flags(char **flag_list );
void get_filename_data(char **varguments );
/******************* pa_types.c ***********************/
int isboolean(char *string);
int ischaracter(char *string);
int isstring(char *string);
int isinteger(char *string);
int isreal(char *string);
int isfilename(char *string);
int isflag(char *string);
/******************* pa_usage.c ***********************/
void build_usage_message(void);
void add_group_to_usage(struct flag_key *flag_ptr);
void add_flag_to_usage(struct flag_format *flag_ptr,char *groupusage);
char *get_ptype_text(int ptype);
void add_filenames_to_usage(void);
void print_usage(void);
/******************* perr.c ***********************/
int install_perr_printf(int (*new_printf)(const char *,...));
int perr(int err,...);
/******************* pix_code.c ***********************/
int pix_code(char *buf,int nbytes,unsigned char *pic,int rows,int nlp);
int codepoint(double b,double x,double y,unsigned char *pic,int rows,int nlp);
int codevec(double b,double x1,double y1,double x2,double y2,unsigned char *pic,int rows,int nlp);
int codeendp(double x,double y,unsigned char *pic,int rows,int nlp);
int fillvect(int ib,double fx0,double fy0,double fx1,double fy1,unsigned char *pic,int rows,int nlp);
/******************* pyralloc.c ***********************/
#if 0
int def_ipyr(struct IIMAGEtag *pyr,int lev,int nr,int nc);
int alloc_ipyr(struct IIMAGEtag *pyr,int botlev,int toplev);
int free_ipyr(struct IIMAGEtag *pyr,int botlev,int toplev);
int def_fpyr(struct FIMAGEtag *pyr,int lev,int nr,int nc);
int alloc_fpyr(struct FIMAGEtag *pyr,int botlev,int toplev);
int free_fpyr(struct FIMAGEtag *pyr,int botlev,int toplev);
int alloc_iimage(struct IIMAGEtag *img);
int free_iimage(struct IIMAGEtag img);
int free_fimage(struct FIMAGEtag img);
int alloc_fimage(struct FIMAGEtag *img);
#endif
int **_alloc_iimage(int nr,int nc) ;
float **_alloc_fimage(int nr,int nc) ;
int _free_iimage(int **img );
int _free_fimage(float **img );
/******************* pyrcopy.c ***********************/
#if 0
int copy_itoii(struct header *hd,struct IIMAGEtag img);
int copy_ftoff(struct header *hd,struct FIMAGEtag img);
#endif
/******************* pyrfilti.c ***********************/
#if 0
int getpyrfilters(char *filename,struct FILTERtag *rf,struct FILTERtag *ef);
int read_1dfilter(struct FILTERtag *f,FILE *stream,char *filename);
int default_1dfilter(struct FILTERtag *f);
#endif
/******************* pyrio.c ***********************/
#if 0
int read_iimage(FILE *stream,struct IIMAGEtag img,int fr,char *fname);
int read_fimage(FILE *stream,struct FIMAGEtag img,int fr,char *fname);
int write_iimage(FILE *stream,struct IIMAGEtag img,int fr);
int write_fimage(FILE *stream,struct FIMAGEtag img,int fr);
int read_ipyr(FILE *stream,struct IIMAGEtag *pyr,int botlev,int toplev,int fr,char *fname);
int read_fpyr(FILE *stream,struct FIMAGEtag *pyr,int botlev,int toplev,int fr,char *fname);
int write_ipyr(FILE *stream,struct IIMAGEtag *pyr,int botlev,int toplev,int fr);
int write_fpyr(FILE *stream,struct FIMAGEtag *pyr,int botlev,int toplev,int fr);
int _read_iimage(FILE *stream,int **img ,int nr,int nc,int fr,char *fname);
int _read_fimage(FILE *stream,float **img ,int nr,int nc,int fr,char *fname);
int _write_iimage(FILE *stream,int **p ,int nr,int nc,int fr);
int _write_fimage(FILE *stream,float **p ,int nr,int nc,int fr);
/******************* pyrnumpi.c ***********************/
int pyrnumpix(int toplev,int nr,int nc);
/******************* pyrredex.c ***********************/
int freduce(struct FIMAGEtag *pyr,int botlev,int toplev,struct FILTERtag rf,int rtype);
int fexpand(struct FIMAGEtag *pyr,int botlev,int toplev,struct FILTERtag ef,int mode,int rtype);
int _freduce_odd(struct FIMAGEtag in,struct FIMAGEtag out,struct FILTERtag rf,int rtype);
int _fexpand_odd(struct FIMAGEtag in,struct FIMAGEtag out,struct FILTERtag ef,int mode,int rtype);
int ireduce(struct IIMAGEtag *pyr,int botlev,int toplev,struct FILTERtag rf,int rtype);
int iexpand(struct IIMAGEtag *pyr,int botlev,int toplev,struct FILTERtag ef,int mode,int rtype);
int _ireduce_odd(struct IIMAGEtag in,struct IIMAGEtag out,struct FILTERtag rf,int rtype);
int _iexpand_odd(struct IIMAGEtag in,struct IIMAGEtag out,struct FILTERtag ef,int mode,int rtype);
/******************* pyrrefle.c ***********************/
int hor_reflectf(struct FIMAGEtag h,int border,int rtype);
int ver_reflectf(struct FIMAGEtag v,int border,int rtype);
int reflectf(struct FIMAGEtag image,int border,int rtype);
int hor_reflecti(struct IIMAGEtag h,int border,int rtype);
int ver_reflecti(struct IIMAGEtag v,int border,int rtype);
int reflecti(struct IIMAGEtag image,int border,int rtype);
#endif
/******************* read_fra.c ***********************/
int read_frame(FILE *fp,char *buf,int buf_limit,int *flags,double *shift_vector,double (*rot_matrix)[3],int fr,char *fname);
/******************* read_hea.c ***********************/
int read_header(struct header *hd);
int fread_header(FILE *fp,struct header *hd,const char *fname);
/******************* read_his.c ***********************/
int read_histo(struct hips_histo *histo,int fr);
int fread_histo(FILE *fp,struct hips_histo *histo,int fr,char *fname);
int hdr_to_histo(struct header *hd,struct hips_histo *histo);
/******************* read_hut.c ***********************/
int read_hdr_a(struct header *hd);
int fread_hdr_a(FILE *fp,struct header *hd,char *fname);
int read_hdr_cpf(struct header *hd,int *typelist);
int fread_hdr_cpf(FILE *fp,struct header *hd,int *typelist,char *fname);
int read_hdr_cpfa(struct header *hd,int *typelist);
int fread_hdr_cpfa(FILE *fp,struct header *hd,int *typelist,char *fname);
int fread_hdr_cpfac(FILE *fp,struct header *hd,int *typelist,char *fname,int flagc,int flaga);
int read_hdr_cc(struct header *hd,struct header *chd,int mask);
int fread_hdr_cc(FILE *fp,struct header *hd,struct header *chd,int mask,char *fname);
int read_hdr_cca(struct header *hd,struct header *chd,int mask);
int fread_hdr_cca(FILE *fp,struct header *hd,struct header *chd,int mask,char *fname);
int fread_hdr_ccac(FILE *fp,struct header *hd,struct header *chd,int mask,char *fname,int flaga);
/******************* read_ima.c ***********************/
int read_image(struct header *hd,int fr);
int fread_image(FILE *fp,struct header *hd,int fr,const char *fname);
/******************* read_roi.c ***********************/
int read_roi(struct header *hd,int fr);
int fread_roi(FILE *fp,struct header *hd,int fr,const char *fname);
/******************* rgb.c ***********************/
int h_btorgb(struct header *hdr,struct header *hdg,struct header *hdb,struct header *hdo);
int h_btorgbz(struct header *hdr,struct header *hdg,struct header *hdb,struct header *hdo);
int h_rgbtob(struct header *hdi,struct header *hdr,struct header *hdg,struct header *hdb);
int h_rgbztob(struct header *hdi,struct header *hdr,struct header *hdg,struct header *hdb);
int h_rgbtob2(struct header *hdi,struct header *hdo,char *color);
int h_rgbztob2(struct header *hdi,struct header *hdo,char *color);
/******************* setforma.c ***********************/
int setformat(struct header *hd,int pfmt);
int setpyrformat(struct header *hd,int pfmt,int toplev);
/******************* setroi.c ***********************/
int setroi(struct header *hd,int fr,int fc,int nr,int nc);
int setroi2(struct header *hd,struct hips_roi *roi);
int getroi(struct header *hd,struct hips_roi *roi);
int clearroi(struct header *hd);
/******************* setsize.c ***********************/
int setsize(struct header *hd,int r,int c);
/******************* strsave.c ***********************/
char *strsave(char *s);
/******************* tob.c ***********************/
int h_tob(struct header *hdi,struct header *hdo);
int h_mptob(struct header *hdi,struct header *hdo);
int h_lptob(struct header *hdi,struct header *hdo);
int h_stob(struct header *hdi,struct header *hdo);
int h_itob(struct header *hdi,struct header *hdo);
int h_ftob(struct header *hdi,struct header *hdo);
/******************* toc.c ***********************/
int h_toc(struct header *hdi,struct header *hdo);
int h_itoc(struct header *hdi,struct header *hdo);
int h_ftoc(struct header *hdi,struct header *hdo);
int h_dtoc(struct header *hdi,struct header *hdo);
int h_dctoc(struct header *hdi,struct header *hdo);
/******************* tod.c ***********************/
int h_tod(struct header *hdi,struct header *hdo);
int h_itod(struct header *hdi,struct header *hdo);
int h_ftod(struct header *hdi,struct header *hdo);
int h_ctod(struct header *hdi,struct header *hdo);
int h_dctod(struct header *hdi,struct header *hdo);
/******************* todc.c ***********************/
int h_todc(struct header *hdi,struct header *hdo);
int h_itodc(struct header *hdi,struct header *hdo);
int h_ftodc(struct header *hdi,struct header *hdo);
int h_dtodc(struct header *hdi,struct header *hdo);
int h_ctodc(struct header *hdi,struct header *hdo);
/******************* tof.c ***********************/
int h_tof(struct header *hdi,struct header *hdo);
int h_btof(struct header *hdi,struct header *hdo);
int h_stof(struct header *hdi,struct header *hdo);
int h_itof(struct header *hdi,struct header *hdo);
int h_dtof(struct header *hdi,struct header *hdo);
int h_ctof(struct header *hdi,struct header *hdo);
int h_dctof(struct header *hdi,struct header *hdo);
/******************* toi.c ***********************/
int h_toi(struct header *hdi,struct header *hdo);
int h_mptoi(struct header *hdi,struct header *hdo);
int h_lptoi(struct header *hdi,struct header *hdo);
int h_btoi(struct header *hdi,struct header *hdo);
int h_sbtoi(struct header *hdi,struct header *hdo);
int h_ustoi(struct header *hdi,struct header *hdo);
int h_stoi(struct header *hdi,struct header *hdo);
int h_uitoi(struct header *hdi,struct header *hdo);
int h_ftoi(struct header *hdi,struct header *hdo);
int h_dtoi(struct header *hdi,struct header *hdo);
int h_ctoi(struct header *hdi,struct header *hdo);
int h_dctoi(struct header *hdi,struct header *hdo);
/******************* tolp.c ***********************/
int h_tolp(struct header *hdi,struct header *hdo);
int h_btolp(struct header *hdi,struct header *hdo);
int h_itolp(struct header *hdi,struct header *hdo);
/******************* tomp.c ***********************/
int h_tomp(struct header *hdi,struct header *hdo);
int h_btomp(struct header *hdi,struct header *hdo);
int h_itomp(struct header *hdi,struct header *hdo);
/******************* tos.c ***********************/
int h_tos(struct header *hdi,struct header *hdo);
int h_btos(struct header *hdi,struct header *hdo);
int h_sbtos(struct header *hdi,struct header *hdo);
int h_itos(struct header *hdi,struct header *hdo);
int h_ftos(struct header *hdi,struct header *hdo);
/******************* tosb.c ***********************/
int h_tosb(struct header *hdi,struct header *hdo);
int h_stosb(struct header *hdi,struct header *hdo);
int h_itosb(struct header *hdi,struct header *hdo);
/******************* toui.c ***********************/
int h_toui(struct header *hdi,struct header *hdo);
int h_itoui(struct header *hdi,struct header *hdo);
/******************* tous.c ***********************/
int h_tous(struct header *hdi,struct header *hdo);
int h_itous(struct header *hdi,struct header *hdo);
/******************* type_is_.c ***********************/
int type_is_col3(struct header *hd);
int ptype_is_col3(int pfmt);
/******************* update_h.c ***********************/
int update_header(struct header *hd,int argc,char **argv );
int update_headern(struct header *hd,int argc,char **argv );
int update_headerc(struct header *hd,int argc,char **argv ,int pflag);
/******************* view_fra.c ***********************/
int view_frame(char *inbuf,int nbytes,char *outbuf,int limit,double dist);
/******************* write_fr.c ***********************/
int write_frame(FILE *fp,char *buf,int nbytes,double *shift_v,double (*rot_m)[3],int fr);
/******************* write_he.c ***********************/
int write_headeru(struct header *hd,int argc,char **argv );
int write_headeru2(struct header *hd,struct header *hdp,int argc,char **argv ,int flag);
int write_headerun(struct header *hd,int argc,char **argv );
int write_header(struct header *hd);
int fwrite_header(FILE *fp,struct header *hd,const char *fname);
/******************* write_hi.c ***********************/
int write_histo(struct hips_histo *histo,int fr);
int fwrite_histo(FILE *fp,struct hips_histo *histo,int fr,char *fname);
int histo_to_hdr(struct header *hd,struct hips_histo *histo);
/******************* write_im.c ***********************/
int write_image(struct header *hd,int fr);
int fwrite_image(FILE *fp,struct header *hd,int fr,const char *fname);
/******************* write_ro.c ***********************/
int write_roi(struct header *hd,int fr);
int fwrite_roi(FILE *fp,struct header *hd,int fr,char *fname);
/******************* wsubs.c ***********************/
int wnocr(FILE *fp,char *s);
int dfprintf(FILE *fp,int i,char *fname);
/******************* xparam.c ***********************/
int setparam(struct header *hda,...);
int setparamd(struct header *hda,...);
int getparam(struct header *hda,...);
int clearparam(struct header *hd,char *name);
struct extpar *findparam(struct header *hd,char *name);
int checkname(char *name);
int mergeparam(struct header *hd1,struct header *hd2);
struct extpar *grepparam(struct header *hd,char *name);
char *strstr1(char *s1,char *s2);
