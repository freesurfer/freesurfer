#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "cma.h"
#include "diag.h"
#include "gcamorph.h"
#include "macros.h"
#include "mri.h"

static int non_hippo_labels[] = {
    entorhinal_cortex,
    Amygdala,
    Cerebral_White_Matter,
    Cerebral_Cortex,
    lateral_ventricle,
    Inf_Lat_Vent,
    Left_Cerebral_Cortex,
    Right_Cerebral_Cortex,
    Left_Cerebral_White_Matter,
    Right_Cerebral_White_Matter,
    Left_Inf_Lat_Vent,
    Right_Inf_Lat_Vent,
    Right_Lateral_Ventricle,
    Left_Lateral_Ventricle,
    Left_Thalamus,
    Right_Thalamus,
    Left_choroid_plexus,
    Right_choroid_plexus,
    Left_Amygdala,
    Right_Amygdala,
};
#define NUM_NON_HIPPO_LABELS (sizeof(non_hippo_labels) / sizeof(non_hippo_labels[0]))

MRI *HIPPOremoveNonHippoLabels(MRI *mri_src, MRI *mri_dst)
{
  int label;

  mri_dst = MRIcopy(mri_src, mri_dst);

  for (unsigned i = 0; i < NUM_NON_HIPPO_LABELS; i++) {
    label = non_hippo_labels[i];
    MRIreplaceValues(mri_dst, mri_dst, label, 0);
  }
  return (mri_dst);
}
MRI *HIPPOestimateIntensityImage(MRI *mri_hippo_labels, MRI *mri_aseg, MRI *mri_intensity, MRI *mri_dst)
{
  GCAM_LABEL_TRANSLATION_TABLE gcam_ltt;
  int i, x, y, z, in_label, out_label;

  if (mri_dst == NULL) {
    mri_dst = MRIalloc(mri_hippo_labels->width, mri_hippo_labels->height, mri_hippo_labels->depth, MRI_FLOAT);
    MRIcopyHeader(mri_hippo_labels, mri_dst);
  }

  memset(&gcam_ltt, 0, sizeof(gcam_ltt));
  gcam_ltt.nlabels = 0;

  memset(gcam_ltt.means, 0, sizeof(gcam_ltt.means));
  memset(gcam_ltt.scales, 0, sizeof(gcam_ltt.scales));

  /* don't use inf_lat_vent label as it may be too small to
           give reliable estimate of density */
  gcam_ltt.input_labels[gcam_ltt.nlabels] = alveus;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_White_Matter;
  gcam_ltt.scales[gcam_ltt.nlabels] = 1.0 /* .95 */;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) {
    gcam_ltt.second_labels[gcam_ltt.nlabels] = Left_Hippocampus;
    gcam_ltt.second_label_pct[gcam_ltt.nlabels] = .25;
    gcam_ltt.nlabels++;
  }
  gcam_ltt.input_labels[gcam_ltt.nlabels] = left_alveus;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Cerebral_White_Matter;
  gcam_ltt.scales[gcam_ltt.nlabels] = 1.0 /* .95 */;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) {
    gcam_ltt.second_labels[gcam_ltt.nlabels] = Left_Hippocampus;
    gcam_ltt.second_label_pct[gcam_ltt.nlabels] = .25;
    gcam_ltt.nlabels++;
  }
  gcam_ltt.input_labels[gcam_ltt.nlabels] = right_alveus;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_White_Matter;
  gcam_ltt.scales[gcam_ltt.nlabels] = 1.0 /* .95 */;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) {
    gcam_ltt.second_labels[gcam_ltt.nlabels] = Right_Hippocampus;
    gcam_ltt.second_label_pct[gcam_ltt.nlabels] = .25;
    gcam_ltt.nlabels++;
  }

  gcam_ltt.input_labels[gcam_ltt.nlabels] = perforant_pathway;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_White_Matter;
  gcam_ltt.scales[gcam_ltt.nlabels] = 1.0 /* .95 */;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = parasubiculum;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;

  gcam_ltt.input_labels[gcam_ltt.nlabels] = presubiculum;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = subiculum;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;

  gcam_ltt.input_labels[gcam_ltt.nlabels] = right_presubiculum;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = left_presubiculum;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;

  gcam_ltt.input_labels[gcam_ltt.nlabels] = right_subiculum;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = left_subiculum;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;

  gcam_ltt.input_labels[gcam_ltt.nlabels] = CA1;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;

  gcam_ltt.input_labels[gcam_ltt.nlabels] = right_CA1;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = left_CA1;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;

  gcam_ltt.input_labels[gcam_ltt.nlabels] = right_CA2_3;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = left_CA2_3;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = right_CA4_DG;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = left_CA4_DG;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;

  gcam_ltt.input_labels[gcam_ltt.nlabels] = CA2;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = CA3;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = CA4;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = GC_DG;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = HATA;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;

  gcam_ltt.input_labels[gcam_ltt.nlabels] = left_fimbria;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Cerebral_White_Matter;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) {
    gcam_ltt.second_labels[gcam_ltt.nlabels] = Left_Hippocampus;
    gcam_ltt.second_label_pct[gcam_ltt.nlabels] = .25;
    gcam_ltt.nlabels++;
  }
  gcam_ltt.input_labels[gcam_ltt.nlabels] = right_fimbria;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_White_Matter;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) {
    gcam_ltt.second_labels[gcam_ltt.nlabels] = Right_Hippocampus;
    gcam_ltt.second_label_pct[gcam_ltt.nlabels] = .25;
    gcam_ltt.nlabels++;
  }
  gcam_ltt.input_labels[gcam_ltt.nlabels] = fimbria;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_White_Matter;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) {
    gcam_ltt.second_labels[gcam_ltt.nlabels] = Right_Hippocampus;
    gcam_ltt.second_label_pct[gcam_ltt.nlabels] = .25;
    gcam_ltt.nlabels++;
  }

  gcam_ltt.input_labels[gcam_ltt.nlabels] = lateral_ventricle;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Lateral_Ventricle;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = molecular_layer_HP;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;

  gcam_ltt.input_labels[gcam_ltt.nlabels] = hippocampal_fissure;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Inf_Lat_Vent;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = right_hippocampal_fissure;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Inf_Lat_Vent;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = left_hippocampal_fissure;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Inf_Lat_Vent;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;

  gcam_ltt.input_labels[gcam_ltt.nlabels] = entorhinal_cortex;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_Cortex;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = molecular_layer_subiculum;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = Amygdala;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = Cerebral_White_Matter;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_White_Matter;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = Cerebral_Cortex;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_Cortex;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;

  gcam_ltt.input_labels[gcam_ltt.nlabels] = Inf_Lat_Vent;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Inf_Lat_Vent;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;

#if 0
	gcam_ltt.input_labels[gcam_ltt.nlabels] = Left_hippocampal_fissure ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Inf_Lat_Vent ;
	if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

	gcam_ltt.input_labels[gcam_ltt.nlabels] = Left_CADG_head ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Hippocampus ;
	if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

	gcam_ltt.input_labels[gcam_ltt.nlabels] = Left_subiculum ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Hippocampus ;
	if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

	gcam_ltt.input_labels[gcam_ltt.nlabels] = Left_fimbria ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Cerebral_White_Matter ;
	gcam_ltt.scales[gcam_ltt.nlabels] = 1.0 /* .95 */ ;
	if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels]))
  {
    gcam_ltt.second_labels[gcam_ltt.nlabels] = Left_Hippocampus ;
    gcam_ltt.second_label_pct[gcam_ltt.nlabels] = .25 ;
		gcam_ltt.nlabels++ ;
  }

	gcam_ltt.input_labels[gcam_ltt.nlabels] = Right_hippocampal_fissure ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Inf_Lat_Vent ;
	if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

	gcam_ltt.input_labels[gcam_ltt.nlabels] = Right_CADG_head ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

	gcam_ltt.input_labels[gcam_ltt.nlabels] = Right_subiculum ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

	gcam_ltt.input_labels[gcam_ltt.nlabels] = Right_fimbria ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_White_Matter ;
	gcam_ltt.scales[gcam_ltt.nlabels] = 1.0 /* .95 */ ;
	if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels]))
  {
    gcam_ltt.second_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
    gcam_ltt.second_label_pct[gcam_ltt.nlabels] = .25 ;
		gcam_ltt.nlabels++ ;
  }
#endif

  gcam_ltt.input_labels[gcam_ltt.nlabels] = Right_choroid_plexus;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) {
    gcam_ltt.second_labels[gcam_ltt.nlabels] = Right_Inf_Lat_Vent;
    gcam_ltt.second_label_pct[gcam_ltt.nlabels] = .25;
    gcam_ltt.nlabels++;
  }

  gcam_ltt.input_labels[gcam_ltt.nlabels] = Left_choroid_plexus;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Hippocampus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) {
    gcam_ltt.second_labels[gcam_ltt.nlabels] = Left_Inf_Lat_Vent;
    gcam_ltt.second_label_pct[gcam_ltt.nlabels] = .25;
    gcam_ltt.nlabels++;
  }

  gcam_ltt.input_labels[gcam_ltt.nlabels] = Left_Thalamus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;

  gcam_ltt.input_labels[gcam_ltt.nlabels] = Right_Thalamus;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;

  gcam_ltt.input_labels[gcam_ltt.nlabels] = Conn_Tissue;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Lateral_Ventricle;
  if (MRIlabelInVolume(mri_hippo_labels, gcam_ltt.input_labels[gcam_ltt.nlabels])) gcam_ltt.nlabels++;

  for (i = 0; i < gcam_ltt.nlabels; i++)
    if (FZERO(gcam_ltt.scales[i])) gcam_ltt.scales[i] = 1;
  for (i = 0; i < gcam_ltt.nlabels; i++)
    if (FZERO(gcam_ltt.scales[i])) {
      gcam_ltt.scales[i] = 1;
      gcam_ltt.done[i] = 0;
      gcam_ltt.means[i] = 0;
    }

  for (x = 0; x < mri_hippo_labels->width; x++) {
    for (y = 0; y < mri_hippo_labels->height; y++) {
      for (z = 0; z < mri_hippo_labels->depth; z++) {
        if (x == Gx && y == Gy && z == Gz) DiagBreak();

        in_label = MRIgetVoxVal(mri_hippo_labels, x, y, z, 0);
        if (IS_UNKNOWN(in_label)) continue;
        out_label = 0;
        for (i = 0; i < gcam_ltt.nlabels; i++)
          if (gcam_ltt.input_labels[i] == in_label) break;

        if (i >= gcam_ltt.nlabels) {
          if (MRIlabelInVolume(mri_aseg, in_label) == 0) {
            printf("could not find translation for label %s (%d)!!!!!!!!!!\n", cma_label_to_name(in_label), in_label);
            continue;
          }
          printf("allocating new translation entry for %s (%d)\n", cma_label_to_name(in_label), in_label);
          i = gcam_ltt.nlabels;
          gcam_ltt.input_labels[gcam_ltt.nlabels] = in_label;
          gcam_ltt.output_labels[gcam_ltt.nlabels] = in_label;
          gcam_ltt.scales[i] = 1;
          gcam_ltt.done[i] = 0;
          gcam_ltt.means[i] = 0;
          gcam_ltt.nlabels++;
        }
        out_label = gcam_ltt.output_labels[i];

        if (gcam_ltt.done[i] == 0)  // estimate mean intensity
        {
          HISTOGRAM *h;
          char fname[STRLEN];
          int peak, j;

          h = MRIhistogramLabel(mri_intensity, mri_aseg, out_label, -1);
          sprintf(fname, "histo%d.dat", out_label);
          HISTOplot(h, fname);
          peak = HISTOfindHighestPeakInRegion(h, 0, h->nbins);
          gcam_ltt.means[i] = h->bins[peak];
          for (j = 0; j < gcam_ltt.nlabels; j++)
            if (gcam_ltt.output_labels[j] == out_label) {
              gcam_ltt.done[j] = 1;
              gcam_ltt.means[j] = gcam_ltt.means[i];
            }
          printf("%s: %2.0f\n", cma_label_to_name(out_label), h->bins[peak]);
          HISTOfree(&h);
        }

        MRIsetVoxVal(mri_dst, x, y, z, 0, gcam_ltt.means[i]);
      }
    }
  }

  return (mri_dst);
}
