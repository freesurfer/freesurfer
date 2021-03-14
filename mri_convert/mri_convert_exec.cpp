/**
 * @brief performs all kinds of conversion and reformatting of MRI volume files
 *
 */
/*
 * Original Author: Bruce Fischl (Apr 16, 1997)
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

#include "mri_convert.hpp"

#include "DICOMRead.h"
#include "cma.h"
#include "diag.h"
#include "fio.h"
#include "fmriutils.h"
#include "gcamorph.h"
#include "mri2.h"
#include "mri2020.hpp"
#include "mri_conform.h"
#include "mri_identify.h"
#include "mriio.hpp"

#include <vector>

#include <absl/strings/str_join.h>
#include <boost/algorithm/string.hpp>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

auto main(int argc, char const *argv[]) -> int {
  std::vector<char const *> args(argv, argv + argc);
  static auto               err_logger = spdlog::stderr_color_mt("mri-convert");
  spdlog::set_level(spdlog::level::warn);
  auto env     = ENV();
  auto cmdargs = CMDARGS(args);
  if (!good_cmdline_args(&cmdargs, &env)) {
    spdlog::get("mri-convert")
        ->critical("Error while parsing command line arguments.");
    return 1;
  }

  MRI *               mri{};
  MRI *               mri2{};
  MRI *               mri_template{};
  MRI *               mri_in_like{};
  int                 i{};
  int                 err{};
  float               invert_val{-1.0};
  int                 conform_width{};
  int                 in_volume_type{};
  int                 out_volume_type{MRI_VOLUME_TYPE_UNKNOWN};
  bool                sizes_good_flag{};
  std::string         tmpstr{};
  char *              stem{};
  char *              ext{};
  std::array<char, 4> ostr{};
  float               i_dot_j{};
  float               i_dot_k{};
  float               j_dot_k{};
  float               fov_x{};
  float               fov_y{};
  float               fov_z{};
  int                 nframes{};
  std::string         in_name_only;
  LTA *               lta_transform{};
  MRI *               mri_transformed{};
  MRI *               mritmp{};
  int                 transform_type{-1};
  MATRIX *            inverse_transform_matrix{};
  bool                read_parcellation_volume_flag{};
  int                 read_otl_flags{};
  int                 temp_type{};
  FILE *              fptmp{};
  int                 j{};
  VOL_GEOM            vgtmp;
  LT *                lt{};
  MATRIX *            T{};
  std::string         cmdline{};
  float               v{};
  MATRIX *            cras{};
  MATRIX *            vmid{};
  int                 c{};
  int                 r{};
  int                 s{};
  int                 f{};
  int                 c1{};
  int                 c2{};
  int                 r1{};
  int                 r2{};
  int                 s1{};
  int                 s2{};
  STAT_TABLE *        StatTable{};
  STAT_TABLE *        OutStatTable{};

  DiagInit(nullptr, nullptr, nullptr);

  // TODO(aboualiaa): Implement a safe version of make_cmd_version_string (in
  // utils/version.cpp) then delete this loop
  for (auto arg : cmdargs.raw) {
    cmdline = absl::StrJoin(cmdargs.raw, " ");
  }

  Progname = GET_PROGRAM_NAME();

  sizes_good_flag = true;

  if (!sizes_good_flag) {
    spdlog::get("mri-convert")->critical("sizes_good_flag is not set");
    fs::util::cli::usage_message(stdout);
    return 1;
  }

  /* ----- copy file name (only -- strip '@' and '#') ----- */
  in_name_only = fs::mri::io::getVolumeName(cmdargs.in_name);

  /* If input type is spm and N_Zero_Pad_Input < 0, set to 3*/
  if (cmdargs.in_type_string == "spm" && N_Zero_Pad_Input < 0) {
    N_Zero_Pad_Input = 3;
  }

  /* ----- catch unknown volume types ----- */
  if ((cmdargs.force_in_type_flag) &&
      cmdargs.forced_in_type == MRI_VOLUME_TYPE_UNKNOWN) {
    fmt::fprintf(stderr,
                 "\n%s: unknown input "
                 "volume type %s\n",
                 Progname, cmdargs.in_type_string.data());
    fs::util::cli::usage_message(stdout);
    return 1;
  }

  if ((cmdargs.force_template_type_flag) &&
      cmdargs.forced_template_type == MRI_VOLUME_TYPE_UNKNOWN) {
    fmt::fprintf(stderr, "\n%s: unknown template volume type %s\n", Progname,
                 cmdargs.template_type_string.data());
    fs::util::cli::usage_message(stdout);
    return 1;
  }

  if ((cmdargs.force_out_type_flag) &&
      cmdargs.forced_out_type == MRI_VOLUME_TYPE_UNKNOWN &&
      (cmdargs.ascii_flag == 0)) {
    fmt::fprintf(stderr, "\n%s: unknown output volume type %s\n", Progname,
                 cmdargs.out_type_string.data());
    fs::util::cli::usage_message(stdout);
    return 1;
  }

  /* ----- warn if read only is desired and an
     output volume is specified or the output info flag is set ----- */
  if ((cmdargs.read_only_flag) && cmdargs.out_name[0] != '\0') {
    fmt::fprintf(stderr,
                 "%s: warning: read only flag is set; "
                 "nothing will be written to %s\n",
                 Progname, cmdargs.out_name.data());
  }
  if ((cmdargs.read_only_flag) &&
      ((cmdargs.out_info_flag) || (cmdargs.out_matrix_flag))) {
    fmt::fprintf(stderr,
                 "%s: warning: read only flag is set; "
                 "no output information will be printed\n",
                 Progname);
  }

  /* ----- get the type of the output ----- */
  if (!cmdargs.read_only_flag) {
    if ((!cmdargs.force_out_type_flag) && (!cmdargs.out_stats_table_flag)) {
      // if(!read_only_flag && !no_write_flag)
      // because conform_flag value changes depending on type below
      {
        out_volume_type = mri_identify(cmdargs.out_name.data());
        if (out_volume_type == MRI_VOLUME_TYPE_UNKNOWN) {
          fmt::fprintf(stderr, "%s: can't determine type of output volume\n",
                       Progname);
          return 1;
        }
      }
    } else {
      out_volume_type = cmdargs.forced_out_type;
    }

    // if output type is COR, then it is always conformed
    if (out_volume_type == MRI_CORONAL_SLICE_DIRECTORY) {
      cmdargs.conform_flag = TRUE;
    }
  }

  /* ----- catch the parse-only flag ----- */
  if (cmdargs.parse_only_flag) {
    fmt::printf("input volume name: %s\n", cmdargs.in_name.data());
    fmt::printf("input name only: %s\n", in_name_only.data());
    fmt::printf("output volume name: %s\n", cmdargs.out_name.data());
    fmt::printf("parse_only_flag = %d\n", cmdargs.parse_only_flag);
    fmt::printf("conform_flag = %d\n", cmdargs.conform_flag);
    fmt::printf("conform_size = %f\n", cmdargs.conform_size);
    fmt::printf("in_info_flag = %d\n", cmdargs.in_info_flag);
    fmt::printf("out_info_flag = %d\n", cmdargs.out_info_flag);
    fmt::printf("in_matrix_flag = %d\n", cmdargs.in_matrix_flag);
    fmt::printf("out_matrix_flag = %d\n", cmdargs.out_matrix_flag);

    if (cmdargs.force_in_type_flag) {
      fmt::printf("input type is %d\n", cmdargs.forced_in_type);
    }
    if (cmdargs.force_out_type_flag) {
      fmt::printf("output type is %d\n", cmdargs.forced_out_type);
    }

    if (cmdargs.subject_name[0] != '\0') {
      fmt::printf("subject name is %s\n", cmdargs.subject_name.data());
    }

    if (invert_val >= 0) {
      fmt::printf("inversion, value is %g\n", invert_val);
    }

    if (!cmdargs.reorder_vals.empty()) {
      fmt::printf("reordering, values are %d %d %d\n", cmdargs.reorder_vals[0],
                  cmdargs.reorder_vals[1], cmdargs.reorder_vals[2]);
    }

    if (!cmdargs.reorder4_vals.empty()) {
      fmt::printf("reordering, values are %d %d %d %d\n",
                  cmdargs.reorder4_vals[0], cmdargs.reorder4_vals[1],
                  cmdargs.reorder4_vals[2], cmdargs.reorder4_vals[3]);
    }

    fmt::printf("translation of otl labels is %s\n",
                cmdargs.translate_labels_flag ? "on" : "off");

    return 0;
  }

  /* ----- check for a gdf image stem if the output type is gdf ----- */
  if (out_volume_type == GDF_FILE && cmdargs.gdf_image_stem.empty()) {
    fmt::fprintf(stderr,
                 "%s: GDF output type, "
                 "but no GDF image file stem\n",
                 Progname);
    return 1;
  }

  /* ----- read the in_like volume ----- */
  if (cmdargs.in_like_flag) {
    fmt::printf("reading info from %s...\n", cmdargs.in_like_name.data());
    mri_in_like = MRIreadInfo(cmdargs.in_like_name.data());
    if (mri_in_like == nullptr) {
      return 1;
    }
  }

  /* ----- read the volume ----- */
  in_volume_type = MRI_VOLUME_TYPE_UNKNOWN;
  if (!cmdargs.in_stats_table_flag) {
    if (cmdargs.force_in_type_flag) {
      in_volume_type = cmdargs.forced_in_type;
    } else {
      in_volume_type = mri_identify(in_name_only.data());
    }
    if (in_volume_type == MRI_VOLUME_TYPE_UNKNOWN) {
      errno = 0;
      if (fio_FileExistsReadable(in_name_only.data()) == 0) {
        fmt::printf("ERROR: file %s does not exist\n", in_name_only.data());
      } else {
        fmt::printf("ERROR: cannot determine file type for %s \n",
                    in_name_only.data());
      }
      if (cmdargs.in_like_flag) {
        MRIfree(&mri_in_like);
      }
      return 1;
    }
  }

  if ((cmdargs.roi_flag) && in_volume_type != GENESIS_FILE) {
    errno = 0;
    ErrorPrintf(ERROR_BADPARM, "rois must be in GE format");
    if (cmdargs.in_like_flag) {
      MRIfree(&mri_in_like);
    }
    return 1;
  }

  if (((cmdargs.zero_ge_z_offset_flag)) && in_volume_type != DICOM_FILE) // E/
  {
    cmdargs.zero_ge_z_offset_flag = FALSE;
    fmt::fprintf(stderr, "Not a GE dicom volume: -zgez "
                         "= --zero_ge_z_offset option ignored.\n");
  }

  fmt::printf("reading from %s...\n", in_name_only.data());

  if (in_volume_type == MGH_MORPH) {
    GCA_MORPH *gcam;
    GCA_MORPH *gcam_out;
    gcam = GCAMread(in_name_only.data());
    if (gcam == nullptr) {
      ErrorExit(ERROR_NOFILE, "%s: could not read input morph from %s",
                Progname, in_name_only.data());
    }
    if (cmdargs.downsample2_flag) {
      gcam_out = GCAMdownsample2(gcam);
    } else {
      gcam_out = gcam;
    }

    GCAMwrite(gcam_out, cmdargs.out_name.data());
    return 0;
  }

  if (in_volume_type == OTL_FILE) {

    if ((!cmdargs.in_like_flag) && (!cmdargs.in_n_k_flag)) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM,
                  "parcellation read: must specify"
                  "a volume depth with either in_like or in_k_count");
      return 1;
    }

    if (!cmdargs.color_file_flag) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM, "parcellation read: must specify a"
                                 "color file name");
      if (cmdargs.in_like_flag) {
        MRIfree(&mri_in_like);
      }
      return 1;
    }

    read_parcellation_volume_flag = TRUE;
    if ((cmdargs.read_only_flag) &&
        ((cmdargs.in_info_flag) || (cmdargs.in_matrix_flag)) &&
        (!cmdargs.in_stats_flag)) {
      read_parcellation_volume_flag = FALSE;
    }

    read_otl_flags = 0x00;

    if (read_parcellation_volume_flag) {
      read_otl_flags |= READ_OTL_READ_VOLUME_FLAG;
    }

    if (cmdargs.fill_parcellation_flag) {
      read_otl_flags |= READ_OTL_FILL_FLAG;
    } else {
      fmt::printf("notice: unfilled parcellations "
                  "currently unimplemented\n");
      fmt::printf("notice: filling outlines\n");
      read_otl_flags |= READ_OTL_FILL_FLAG;
    }

    if (cmdargs.translate_labels_flag) {
      read_otl_flags |= READ_OTL_TRANSLATE_LABELS_FLAG;
    }

    if (cmdargs.zero_outlines_flag) {
      read_otl_flags |= READ_OTL_ZERO_OUTLINES_FLAG;
    }

    if (cmdargs.in_like_flag) {
      mri = MRIreadOtl(cmdargs.in_name.data(), mri_in_like->width,
                       mri_in_like->height, mri_in_like->depth,
                       cmdargs.color_file_name.data(), read_otl_flags);
    } else {
      mri = MRIreadOtl(cmdargs.in_name.data(), 0, 0, cmdargs.in_n_k,
                       cmdargs.color_file_name.data(), read_otl_flags);
    }

    if (mri == nullptr) {
      if (cmdargs.in_like_flag) {
        MRIfree(&mri_in_like);
      }
      return 1;
    }

#if 0
    if (out_volume_type == MGH_MORPH) {
      GCA_MORPH *gcam;

      if (mri->nframes != 3)
        ErrorExit(ERROR_UNSUPPORTED,
                  "%s: input volume must have 3 frames not %d to write to warp",
                  Progname, mri->nframes);
      gcam = GCAMalloc(mri->width, mri->height, mri->depth);
      GCAMinit(gcam, mri, NULL, NULL, 0);
      GCAMreadWarpFromMRI(gcam, mri);
      GCAMwrite(gcam, out_name);
      return 0;
    }
#endif

    /* ----- smooth the parcellation if requested ----- */
    if (cmdargs.smooth_parcellation_count != 0) {
      fmt::printf("smoothing parcellation...\n");
      mri2 = MRIsmoothParcellation(mri, cmdargs.smooth_parcellation_count);
      if (mri2 == nullptr) {
        if (cmdargs.in_like_flag) {
          MRIfree(&mri_in_like);
        }
        return 1;
      }
      MRIfree(&mri);
      mri = mri2;
    }

    cmdargs.resample_type_val = SAMPLE_NEAREST;
    cmdargs.no_scale_flag     = TRUE;

  } else if (cmdargs.roi_flag) {
    if ((!cmdargs.in_like_flag) && (!cmdargs.in_n_k_flag)) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM, "roi read: must specify a volume"
                                 "depth with either in_like or in_k_count");
      if (cmdargs.in_like_flag) {
        MRIfree(&mri_in_like);
      }
      return 1;
    }

    if (cmdargs.in_like_flag) {
      mri = MRIreadGeRoi(cmdargs.in_name.data(), mri_in_like->depth);
    } else {
      mri = MRIreadGeRoi(cmdargs.in_name.data(), cmdargs.in_n_k);
    }

    if (mri == nullptr) {
      if (cmdargs.in_like_flag) {
        MRIfree(&mri_in_like);
      }
      return 1;
    }

    cmdargs.resample_type_val = SAMPLE_NEAREST;
    cmdargs.no_scale_flag     = TRUE;
  } else {
    if ((cmdargs.read_only_flag) &&
        ((cmdargs.in_info_flag) || (cmdargs.in_matrix_flag)) &&
        (!cmdargs.in_stats_flag)) {
      if (cmdargs.force_in_type_flag) {
        mri = MRIreadHeader(cmdargs.in_name.data(), in_volume_type);
      } else {
        mri = MRIreadInfo(cmdargs.in_name.data());
      }
    } else {
      if (cmdargs.force_in_type_flag) {
        // fmt::printf("MRIreadType()\n");
        mri = MRIreadType(cmdargs.in_name.data(), in_volume_type);
      } else {
        if (cmdargs.nthframe < 0) {
          if (!cmdargs.in_stats_table_flag) {
            mri = MRIread(cmdargs.in_name.data());
          } else {
            fmt::printf("Loading in stat table %s\n", cmdargs.in_name.data());
            StatTable = LoadStatTable(cmdargs.in_name.data());
            if (StatTable == nullptr) {
              fmt::printf("ERROR: loading y %s as a stat table\n",
                          cmdargs.in_name.data());
              return 1;
            }
            mri = StatTable->mri;
          }
        } else {
          mri = MRIreadEx(cmdargs.in_name.data(), cmdargs.nthframe);
        }
      }
    }
  }

  if (mri == nullptr) {
    if (cmdargs.in_like_flag) {
      MRIfree(&mri_in_like);
    }
    return 1;
  }

  if (cmdargs.outside_val > 0) {
    mri->outside_val = cmdargs.outside_val;
  }
  if (cmdargs.upsample_flag) {
    fmt::printf("UpsampleFactor = %d\n", cmdargs.upsample_factor);
    mritmp = MRIupsampleN(mri, nullptr, cmdargs.upsample_factor);
    MRIfree(&mri);
    mri = mritmp;
  }

  if (cmdargs.slice_crop_flag) {
    fmt::printf("Cropping slices from %d to %d\n", cmdargs.slice_crop_start,
                cmdargs.slice_crop_stop);
    mri2 = MRIcrop(mri, 0, 0, cmdargs.slice_crop_start, mri->width - 1,
                   mri->height - 1, cmdargs.slice_crop_stop);
    if (mri2 == nullptr) {
      return 1;
    }
    MRIfree(&mri);
    mri = mri2;
  }

  if (cmdargs.left_right_reverse) {
    // Performs a left-right reversal of the geometry by finding the
    // dimension that is most left-right oriented and negating its
    // column in the vox2ras matrix. The pixel data itself is not
    // changed. It would also have worked to have negated the first
    // row in the vox2ras, but negating a column is the header
    // equivalent to reversing the pixel data in a dimension. This
    // also adjusts the CRAS so that if the pixel data are
    // subsequently reversed then a header registration will bring
    // the new and the original volume into perfect registration.
    fmt::printf("WARNING: applying left-right reversal to the input geometry\n"
                "without reversing the pixel data. This will likely make \n"
                "the volume geometry WRONG, so make sure you know what you  \n"
                "are doing.\n");

    fmt::printf("  CRAS Before: %g %g %g\n", mri->c_r, mri->c_a, mri->c_s);
    T                = MRIxfmCRS2XYZ(mri, 0);
    vmid             = MatrixAlloc(4, 1, MATRIX_REAL);
    vmid->rptr[4][1] = 1;

    MRIdircosToOrientationString(mri, ostr.data());
    if (ostr[0] == 'L' || ostr[0] == 'R') {
      fmt::printf("  Reversing geometry for the columns\n");
      vmid->rptr[1][1] = mri->width / 2.0 - 1;
      vmid->rptr[2][1] = mri->height / 2.0;
      vmid->rptr[3][1] = mri->depth / 2.0;
      cras             = MatrixMultiply(T, vmid, NULL);
      mri->x_r *= -1.0;
      mri->x_a *= -1.0;
      mri->x_s *= -1.0;
    }
    if (ostr[1] == 'L' || ostr[1] == 'R') {
      fmt::printf("  Reversing geometry for the rows\n");
      vmid->rptr[1][1] = mri->width / 2.0;
      vmid->rptr[2][1] = mri->height / 2.0 - 1;
      vmid->rptr[3][1] = mri->depth / 2.0;
      cras             = MatrixMultiply(T, vmid, NULL);
      mri->y_r *= -1.0;
      mri->y_a *= -1.0;
      mri->y_s *= -1.0;
    }
    if (ostr[2] == 'L' || ostr[2] == 'R') {
      fmt::printf("  Reversing geometry for the slices\n");
      vmid->rptr[1][1] = mri->width / 2.0;
      vmid->rptr[2][1] = mri->height / 2.0;
      vmid->rptr[3][1] = mri->depth / 2.0 - 1;
      cras             = MatrixMultiply(T, vmid, NULL);
      mri->z_r *= -1.0;
      mri->z_a *= -1.0;
      mri->z_s *= -1.0;
    }
    mri->c_r = cras->rptr[1][1];
    mri->c_a = cras->rptr[2][1];
    mri->c_s = cras->rptr[3][1];
    fmt::printf("  CRAS After %g %g %g\n", mri->c_r, mri->c_a, mri->c_s);
    MatrixFree(&T);
    MatrixFree(&vmid);
    MatrixFree(&cras);
  }

  if (cmdargs.left_right_swap_label) {
    fmt::printf("Performing left-right swap of labels\n");
    // Good for aseg, aparc+aseg, wmparc, etc. Does not change geometry
    mri2 = MRIlrswapAseg(mri);
    MRIfree(&mri);
    mri = mri2;
  }

  if (cmdargs.left_right_reverse_pix) {
    // Performs a left-right reversal of the pixels by finding the
    // dimension that is most left-right oriented and reversing
    // the order of the pixels. The geometry itself is not
    // changed.
    fmt::printf("WARNING: applying left-right reversal to the input pixels\n"
                "without changing geometry. This will likely make \n"
                "the volume geometry WRONG, so make sure you know what you  \n"
                "are doing.\n");

    MRIdircosToOrientationString(mri, ostr.data());
    mri2 = MRIcopy(mri, nullptr);
    if (ostr[0] == 'L' || ostr[0] == 'R') {
      fmt::printf("  Reversing pixels for the columns\n");
      for (c = 0; c < mri->width; c++) {
        c2 = mri->width - c - 1;
        for (r = 0; r < mri->height; r++) {
          for (s = 0; s < mri->depth; s++) {
            for (f = 0; f < mri->nframes; f++) {
              v = MRIgetVoxVal(mri, c, r, s, f);
              MRIsetVoxVal(mri2, c2, r, s, f, v);
            }
          }
        }
      }
    }
    if (ostr[1] == 'L' || ostr[1] == 'R') {
      fmt::printf("  Reversing pixels for the rows\n");
      for (c = 0; c < mri->width; c++) {
        for (r = 0; r < mri->height; r++) {
          r2 = mri->height - r - 1;
          for (s = 0; s < mri->depth; s++) {
            for (f = 0; f < mri->nframes; f++) {
              v = MRIgetVoxVal(mri, c, r, s, f);
              MRIsetVoxVal(mri2, c, r2, s, f, v);
            }
          }
        }
      }
    }
    if (ostr[2] == 'L' || ostr[2] == 'R') {
      fmt::printf("  Reversing pixels for the slices\n");
      for (c = 0; c < mri->width; c++) {
        for (r = 0; r < mri->height; r++) {
          for (s = 0; s < mri->depth; s++) {
            s2 = mri->depth - s - 1;
            for (f = 0; f < mri->nframes; f++) {
              v = MRIgetVoxVal(mri, c, r, s, f);
              MRIsetVoxVal(mri2, c, r, s2, f, v);
            }
          }
        }
      }
    }
    MRIfree(&mri);
    mri = mri2;
  }

  if ((cmdargs.left_right_mirror_flag) ||
      (cmdargs.left_right_keep_flag)) //(mr, 2012-03-08)
  {
    // Mirror one half of the image into the other half (left-right)
    // the split point is at the middle of the dimension that is
    // most left-right oriented. The header geometry is not changed
    // If LeftRightKeep, then don't mirror but fill with zero and
    // only keep the specified hemisphere
    fmt::printf("WARNING: Mirroring %s half into the other half,\n"
                "or masking one half of the image is only meaningful if\n"
                "the image is upright and centerd, see make_upright.\n",
                cmdargs.left_right_mirror_hemi.data());

    MRIdircosToOrientationString(mri, ostr.data());
    fmt::printf("  Orientation string: %s\n", ostr.data());
    mri2 = MRIcopy(mri, nullptr);
    v    = 0;
    if (ostr[0] == 'L' || ostr[0] == 'R') {
      fmt::printf("  Mirror or keep pixels for the columns\n");
      for (c = 0; c < mri->width / 2; c++) {
        c1 = c;
        if (ostr[0] == toupper(cmdargs.left_right_mirror_hemi[0])) {
          c1 = mri->width - c - 1;
        }
        c2 = mri->width - c1 - 1;
        for (r = 0; r < mri->height; r++) {
          for (s = 0; s < mri->depth; s++) {
            for (f = 0; f < mri->nframes; f++) {
              if (cmdargs.left_right_mirror_flag) {
                v = MRIgetVoxVal(mri, c1, r, s, f);
              }
              MRIsetVoxVal(mri2, c2, r, s, f, v);
            }
          }
        }
      }
    }
    if (ostr[1] == 'L' || ostr[1] == 'R') {
      fmt::printf("  Mirror or keep pixels for the rows\n");
      for (c = 0; c < mri->width; c++) {
        for (r = 0; r < mri->height / 2; r++) {
          r1 = r;
          if (ostr[1] == toupper(cmdargs.left_right_mirror_hemi[0])) {
            r1 = mri->height - r - 1;
          }
          r2 = mri->height - r1 - 1;
          for (s = 0; s < mri->depth; s++) {
            for (f = 0; f < mri->nframes; f++) {
              if (cmdargs.left_right_mirror_flag) {
                v = MRIgetVoxVal(mri, c, r1, s, f);
              }
              MRIsetVoxVal(mri2, c, r2, s, f, v);
            }
          }
        }
      }
    }
    if (ostr[2] == 'L' || ostr[2] == 'R') {
      fmt::printf("  Mirror or keep pixels for the slices\n");
      for (c = 0; c < mri->width; c++) {
        for (r = 0; r < mri->height; r++) {
          for (s = 0; s < mri->depth / 2; s++) {
            s1 = s;
            if (ostr[2] == toupper(cmdargs.left_right_mirror_hemi[0])) {
              s1 = mri->depth - s - 1;
            }
            s2 = mri->depth - s1 - 1;
            for (f = 0; f < mri->nframes; f++) {
              if (cmdargs.left_right_mirror_flag) {
                v = MRIgetVoxVal(mri, c, r, s1, f);
              }
              MRIsetVoxVal(mri2, c, r, s2, f, v);
            }
          }
        }
      }
    }
    MRIfree(&mri);
    mri = mri2;
  }

  if (cmdargs.flip_cols) {
    // Reverses the columns WITHOUT changing the geometry in the
    // header.  Know what you are going here.
    fmt::printf("WARNING: flipping cols without changing geometry\n");
    mri2 = MRIcopy(mri, nullptr);
    fmt::printf("type %d %d\n", mri->type, mri2->type);
    for (c = 0; c < mri->width; c++) {
      c2 = mri->width - c - 1;
      for (r = 0; r < mri->height; r++) {
        for (s = 0; s < mri->depth; s++) {
          for (f = 0; f < mri->nframes; f++) {
            v = MRIgetVoxVal(mri, c, r, s, f);
            MRIsetVoxVal(mri2, c2, r, s, f, v);
          }
        }
      }
    }
    fmt::printf("type %d %d\n", mri->type, mri2->type);
    MRIfree(&mri);
    mri = mri2;
  }

  if (cmdargs.slice_reverse) {
    fmt::printf("Reversing slices, updating vox2ras\n");
    mri2 = MRIreverseSlices(mri, nullptr);
    if (mri2 == nullptr) {
      return 1;
    }
    MRIfree(&mri);
    mri = mri2;
  }

  if (cmdargs.SliceBias) {
    fmt::printf("Applying Half-Cosine Slice Bias, Alpha = %g\n",
                cmdargs.SliceBiasAlpha);
    MRIhalfCosBias(mri, cmdargs.SliceBiasAlpha, mri);
  }

  if (cmdargs.reduce > 0) {
    MRI *mri_tmp;
    int  r;

    if (mri->depth == 1) {
      for (r = 0; r < cmdargs.reduce; r++) {
        mri_tmp = MRIreduce2D(mri, nullptr);
        MRIfree(&mri);
        mri = mri_tmp;
      }
    } else {
      for (r = 0; r < cmdargs.reduce; r++) {
        mri_tmp = MRIreduce(mri, nullptr);
        MRIfree(&mri);
        mri = mri_tmp;
      }
    }
  }

  if (cmdargs.fwhm > 0) {
    fmt::printf("Smoothing input at fwhm = %g, gstd = %g\n", cmdargs.fwhm,
                cmdargs.gstd);
    MRIgaussianSmooth(mri, cmdargs.gstd, 1, mri);
  }

  if (!FEQUAL(cmdargs.rescale_factor, 1.0)) {
    // Rescale so that the global mean of input is rescale_factor
    double globalmean = 0;
    for (c = 0; c < mri->width; c++) {
      for (r = 0; r < mri->height; r++) {
        for (s = 0; s < mri->depth; s++) {
          for (f = 0; f < mri->nframes; f++) {
            v = MRIgetVoxVal(mri, c, r, s, f);
            globalmean += v;
          }
        }
      }
    }
    globalmean /= static_cast<double>(mri->width * mri->height * mri->depth *
                                      mri->nframes);
    fmt::printf("Global rescaling input mean from %g to %g\n", globalmean,
                cmdargs.rescale_factor);
    MRIscalarMul(mri, mri, cmdargs.rescale_factor / globalmean);
  }

  MRIaddCommandLine(mri, cmdline.data());
  if (!FEQUAL(cmdargs.scale_factor, 1.0)) {
    fmt::printf("scaling input volume by %2.3f\n", cmdargs.scale_factor);
    MRIscalarMul(mri, mri, cmdargs.scale_factor);
  }

  if (cmdargs.zero_ge_z_offset_flag) // E/
  {
    mri->c_s = 0.0;
  }

  fmt::printf("TR=%2.2f, TE=%2.2f, TI=%2.2f, flip angle=%2.2f\n", mri->tr,
              mri->te, mri->ti, DEGREES(mri->flip_angle));
  if (in_volume_type != OTL_FILE) {
    if (cmdargs.fill_parcellation_flag) {
      fmt::printf("fill_parcellation flag ignored on a "
                  "non-parcellation read\n");
    }
    if (cmdargs.smooth_parcellation_count != 0) {
      fmt::printf("smooth_parcellation flag ignored on a "
                  "non-parcellation read\n");
    }
  }

  /* ----- apply the in_like volume if it's been read ----- */
  if (cmdargs.in_like_flag) {
    if (mri->width != mri_in_like->width ||
        mri->height != mri_in_like->height ||
        mri->depth != mri_in_like->depth) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM, "volume sizes do not match\n");
      ErrorPrintf(ERROR_BADPARM,
                  "%s: (width, height, depth, frames) "
                  "= (%d, %d, %d, %d)\n",
                  cmdargs.in_name.c_str(), mri->width, mri->height, mri->depth,
                  mri->nframes);
      ErrorPrintf(ERROR_BADPARM,
                  "%s: (width, height, depth, frames) "
                  "= (%d, %d, %d, %d)\n",
                  cmdargs.in_like_name.c_str(), mri_in_like->width,
                  mri_in_like->height, mri_in_like->depth,
                  mri_in_like->nframes);
      MRIfree(&mri);
      MRIfree(&mri_in_like);
      return 1;
    }
    if (mri->nframes != mri_in_like->nframes) {
      fmt::printf("INFO: frames are not the same\n");
    }
    temp_type = mri->type;

    mritmp = MRIcopyHeader(mri_in_like, mri);
    if (mritmp == nullptr) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM,
                  "error copying information from "
                  "%s structure to %s structure\n",
                  cmdargs.in_like_name.c_str(), cmdargs.in_name.data());
      MRIfree(&mri);
      MRIfree(&mri_in_like);
      return 1;
    }

    mri->type = temp_type;
    MRIfree(&mri_in_like);
  }

  if (mri->ras_good_flag == 0) {
    fmt::printf("WARNING: it does not appear that there "
                "was sufficient information\n"
                "in the input to assign orientation to the volume... \n");
    if (cmdargs.force_ras_good) {
      fmt::printf("However, you have specified that the "
                  "default orientation should\n"
                  "be used with by adding --force_ras_good "
                  "on the command-line.\n");
      mri->ras_good_flag = 1;
    }
    if ((cmdargs.in_i_direction_flag) || (cmdargs.in_j_direction_flag) ||
        (cmdargs.in_k_direction_flag)) {
      fmt::printf("However, you have specified one or more "
                  "orientations on the \n"
                  "command-line using -i?d or --in-?-direction (?=i,j,k).\n");
      mri->ras_good_flag = 1;
    }
  }

  /* ----- apply command-line parameters ----- */
  if (cmdargs.in_i_size_flag) {
    mri->xsize = cmdargs.in_i_size;
  }
  if (cmdargs.in_j_size_flag) {
    mri->ysize = cmdargs.in_j_size;
  }
  if (cmdargs.in_k_size_flag) {
    mri->zsize = cmdargs.in_k_size;
  }
  if (cmdargs.in_i_direction_flag) {
    mri->x_r           = cmdargs.in_i_directions[0];
    mri->x_a           = cmdargs.in_i_directions[1];
    mri->x_s           = cmdargs.in_i_directions[2];
    mri->ras_good_flag = 1;
  }
  if (cmdargs.in_j_direction_flag) {
    mri->y_r           = cmdargs.in_j_directions[0];
    mri->y_a           = cmdargs.in_j_directions[1];
    mri->y_s           = cmdargs.in_j_directions[2];
    mri->ras_good_flag = 1;
  }
  if (cmdargs.in_k_direction_flag) {
    mri->z_r           = cmdargs.in_k_directions[0];
    mri->z_a           = cmdargs.in_k_directions[1];
    mri->z_s           = cmdargs.in_k_directions[2];
    mri->ras_good_flag = 1;
  }
  if (cmdargs.delete_cmds) {
    mri->ncmds = 0;
  }
  if (!cmdargs.new_transform_fname.empty()) {
    fmt::printf("Changing xform name to %s\n",
                cmdargs.new_transform_fname.data());
    strcpy(mri->transform_fname, cmdargs.new_transform_fname.data());
  }
  if (!cmdargs.in_orientation_string.empty()) {
    fmt::printf("Setting input orientation to %s\n",
                cmdargs.in_orientation_string.data());
    MRIorientationStringToDircos(mri, cmdargs.in_orientation_string.data());
    mri->ras_good_flag = 1;
  }
  if (!cmdargs.in_center.empty()) {
    mri->c_r = cmdargs.in_center[0];
    mri->c_a = cmdargs.in_center[1];
    mri->c_s = cmdargs.in_center[2];
  }
  if (!cmdargs.delta_in_center.empty()) {
    mri->c_r += cmdargs.delta_in_center[0];
    mri->c_a += cmdargs.delta_in_center[1];
    mri->c_s += cmdargs.delta_in_center[2];
  }
  if (cmdargs.subject_name[0] != '\0') {
    strcpy(mri->subject_name, cmdargs.subject_name.data());
  }

  if (cmdargs.in_tr_flag) {
    mri->tr = cmdargs.in_tr;
  }
  if (cmdargs.in_ti_flag) {
    mri->ti = cmdargs.in_ti;
  }
  if (cmdargs.in_te_flag) {
    mri->te = cmdargs.in_te;
  }
  if (cmdargs.in_flip_angle_flag) {
    mri->flip_angle = cmdargs.in_flip_angle;
  }

  /* ----- correct starts, ends, and fov if necessary ----- */
  if ((cmdargs.in_i_size_flag) || (cmdargs.in_j_size_flag) ||
      (cmdargs.in_k_size_flag)) {

    fov_x = mri->xsize * mri->width;
    fov_y = mri->ysize * mri->height;
    fov_z = mri->zsize * mri->depth;

    mri->xend   = fov_x / 2.0;
    mri->xstart = -mri->xend;
    mri->yend   = fov_y / 2.0;
    mri->ystart = -mri->yend;
    mri->zend   = fov_z / 2.0;
    mri->zstart = -mri->zend;

    mri->fov = (fov_x > fov_y ? (fov_x > fov_z ? fov_x : fov_z)
                              : (fov_y > fov_z ? fov_y : fov_z));
  }

  /* ----- give a warning for non-orthogonal directions ----- */
  i_dot_j = mri->x_r * mri->y_r + mri->x_a * mri->y_a + mri->x_s * mri->y_s;
  i_dot_k = mri->x_r * mri->z_r + mri->x_a * mri->z_a + mri->x_s * mri->z_s;
  j_dot_k = mri->y_r * mri->z_r + mri->y_a * mri->z_a + mri->y_s * mri->z_s;
  if (fabs(i_dot_j) > CLOSE_ENOUGH || fabs(i_dot_k) > CLOSE_ENOUGH ||
      fabs(i_dot_k) > CLOSE_ENOUGH) {
    fmt::printf("warning: input volume axes are not orthogonal:"
                "i_dot_j = %.6f, i_dot_k = %.6f, j_dot_k = %.6f\n",
                i_dot_j, i_dot_k, j_dot_k);
  }
  fmt::printf("i_ras = (%g, %g, %g)\n", mri->x_r, mri->x_a, mri->x_s);
  fmt::printf("j_ras = (%g, %g, %g)\n", mri->y_r, mri->y_a, mri->y_s);
  fmt::printf("k_ras = (%g, %g, %g)\n", mri->z_r, mri->z_a, mri->z_s);

  /* ----- catch the in info flag ----- */
  if (cmdargs.in_info_flag) {
    fmt::printf("input structure:\n");
    MRIdump(mri, stdout);
  }

  if (cmdargs.in_matrix_flag) {
    MATRIX *i_to_r;
    i_to_r = extract_i_to_r(mri);
    if (i_to_r != nullptr) {
      fmt::printf("input ijk -> ras:\n");
      MatrixPrint(stdout, i_to_r);
      MatrixFree(&i_to_r);
    } else {
      fmt::printf("error getting input matrix\n");
    }
  }

  /* ----- catch the in stats flag ----- */
  if (cmdargs.in_stats_flag) {
    MRIprintStats(mri, stdout);
  }

  // Load the transform
  if (cmdargs.transform_flag) {
    fmt::printf("INFO: Reading transformation from file %s...\n",
                cmdargs.transform_fname.data());
    transform_type = TransformFileNameType(cmdargs.transform_fname.data());
    if (transform_type == MNI_TRANSFORM_TYPE ||
        transform_type == TRANSFORM_ARRAY_TYPE ||
        transform_type == REGISTER_DAT || transform_type == FSLREG_TYPE) {
      fmt::printf("Reading transform with LTAreadEx()\n");
      // lta_transform = LTAread(transform_fname);
      lta_transform = LTAreadEx(cmdargs.transform_fname.data());
      if (lta_transform == nullptr) {
        fmt::fprintf(stderr, "ERROR: Reading transform from file %s\n",
                     cmdargs.transform_fname.data());
        return 1;
      }
      if (transform_type == FSLREG_TYPE) {
        MRI *tmp = nullptr;
        if (!cmdargs.out_like_flag) {
          fmt::printf("ERROR: fslmat does not have the information "
                      "on the dst volume\n");
          fmt::printf(
              "ERROR: you must give option '--like volume' to specify the"
              " dst volume info\n");
          MRIfree(&mri);
          return 1;
        }
        // now setup dst volume info
        tmp = MRIreadHeader(cmdargs.out_like_name.data(),
                            MRI_VOLUME_TYPE_UNKNOWN);
        // flsmat does not contain src and dst info
        LTAmodifySrcDstGeom(lta_transform, mri, tmp);
        // add src and dst information
        LTAchangeType(lta_transform, LINEAR_VOX_TO_VOX);
        MRIfree(&tmp);
      } // end FSLREG_TYPE

      if (cmdargs.DevXFM != 0) {
        fmt::printf("INFO: devolving XFM (%s)\n",
                    cmdargs.devolvexfm_subject.data());
        fmt::printf("-------- before ---------\n");
        MatrixPrint(stdout, lta_transform->xforms[0].m_L);
        T = DevolveXFM(cmdargs.devolvexfm_subject.data(),
                       lta_transform->xforms[0].m_L, nullptr);
        if (T == nullptr) {
          return 1;
        }
        fmt::printf("-------- after ---------\n");
        MatrixPrint(stdout, lta_transform->xforms[0].m_L);
        fmt::printf("-----------------\n");
      } // end DevXFM

      if (cmdargs.invert_transform_flag) {
        inverse_transform_matrix =
            MatrixInverse(lta_transform->xforms[0].m_L, nullptr);
        if (inverse_transform_matrix == nullptr) {
          fmt::fprintf(stderr, "ERROR: inverting transform\n");
          MatrixPrint(stdout, lta_transform->xforms[0].m_L);
          return 1;
        }

        MatrixFree(&(lta_transform->xforms[0].m_L));
        lta_transform->xforms[0].m_L = inverse_transform_matrix;
        // reverse src and dst target info.
        // since it affects the c_ras values of the result
        // in LTAtransform()
        // question is what to do when transform src info is invalid.
        lt = &lta_transform->xforms[0];
        if (lt->src.valid == 0) {
          char  buf[512];
          char *p;
          MRI * mriOrig;

          fmt::fprintf(stderr, "INFO: Trying to get the source "
                               "volume information from the transform name\n");
          // copy transform filename
          strcpy(buf, cmdargs.transform_fname.data());
          // reverse look for the first '/'
          p = strrchr(buf, '/');
          if (p != nullptr) {
            p++;
            *p = '\0';
            // set the terminator. i.e.
            // ".... mri/transforms" from "
            //..../mri/transforms/talairach.xfm"
          } else // no / present means only a filename is given
          {
            strcpy(buf, "./");
          }
          strcat(buf, "../orig"); // go to mri/orig from
          // mri/transforms/
          // check whether we can read header info or not
          mriOrig = MRIreadHeader(buf, MRI_VOLUME_TYPE_UNKNOWN);
          if (mriOrig != nullptr) {
            getVolGeom(mriOrig, &lt->src);
            fmt::fprintf(stderr, "INFO: Succeeded in retrieving "
                                 "the source volume info.\n");
          } else {
            fmt::printf("INFO: Failed to find %s as a source volume.  \n"
                        "      The inverse c_(ras) may not be valid.\n",
                        buf);
          }
        } // end not valid
        copyVolGeom(&lt->dst, &vgtmp);
        copyVolGeom(&lt->src, &lt->dst);
        copyVolGeom(&vgtmp, &lt->src);
      } // end invert_transform_flag

    } // end transform type
    // if type is non-linear, load transform below when applying it...
  } // end transform_flag

  if (cmdargs.reslice_like_flag) {

    if (cmdargs.force_template_type_flag) {
      fmt::printf("reading template info from (type %s) volume %s...\n",
                  cmdargs.template_type_string.data(),
                  cmdargs.reslice_like_name.data());
      mri_template = MRIreadHeader(cmdargs.reslice_like_name.data(),
                                   cmdargs.forced_template_type);
      if (mri_template == nullptr) {
        fmt::fprintf(stderr, "error reading from volume %s\n",
                     cmdargs.reslice_like_name.data());
        return 1;
      }
    } else {
      fmt::printf("reading template info from volume %s...\n",
                  cmdargs.reslice_like_name.data());
      mri_template = MRIreadInfo(cmdargs.reslice_like_name.data());
      if (mri_template == nullptr) {
        fmt::fprintf(stderr, "error reading from volume %s\n",
                     cmdargs.reslice_like_name.data());
        return 1;
      }
      // if we loaded a transform above, set the lta target geometry from
      // mri_template:
      if (lta_transform != nullptr) {
        lta_transform = LTAchangeType(lta_transform, LINEAR_RAS_TO_RAS);
        LTAmodifySrcDstGeom(lta_transform, nullptr, mri_template);
      }
    }
  } else {
    mri_template = MRIallocHeader(mri->width, mri->height, mri->depth,
                                  mri->type, mri->nframes);
    MRIcopyHeader(mri, mri_template);

    // if we loaded a transform above, get target geometry from lta:
    if ((lta_transform != nullptr) && lta_transform->xforms[0].dst.valid == 1) {
      useVolGeomToMRI(&lta_transform->xforms[0].dst, mri_template);
    }

    if (cmdargs.conform_flag) {
      conform_width = 256;
      if (cmdargs.conform_min_flag) {
        cmdargs.conform_size = MRIfindMinSize(mri, &conform_width);
      } else {
        if (cmdargs.conform_width_256_flag) {
          conform_width = 256; // force it
        } else {
          conform_width = MRIfindRightSize(mri, cmdargs.conform_size);
        }
      }
      mri_template = MRIconformedTemplate(
          mri, conform_width, cmdargs.conform_size, cmdargs.conf_keep_dc);
    } else if (!cmdargs.voxel_size.empty()) {
      mri_template = MRIallocHeader(mri->width, mri->height, mri->depth,
                                    mri->type, mri->nframes);
      MRIcopyHeader(mri, mri_template);

      mri_template->nframes = mri->nframes;

      mri_template->width = static_cast<int>(
          ceil(mri->xsize * mri->width / cmdargs.voxel_size[0]));
      mri_template->height = static_cast<int>(
          ceil(mri->ysize * mri->height / cmdargs.voxel_size[1]));
      mri_template->depth = static_cast<int>(
          ceil(mri->zsize * mri->depth / cmdargs.voxel_size[2]));

      mri_template->xsize = cmdargs.voxel_size[0];
      mri_template->ysize = cmdargs.voxel_size[1];
      mri_template->zsize = cmdargs.voxel_size[2];

      mri_template->xstart = -mri_template->width / 2;
      mri_template->xend   = mri_template->width / 2;
      mri_template->ystart = -mri_template->height / 2;
      mri_template->yend   = mri_template->height / 2;
      mri_template->zstart = -mri_template->depth / 2;
      mri_template->zend   = mri_template->depth / 2;

    } else if (!cmdargs.downsample_factor.empty()) {
      mri_template = MRIallocHeader(mri->width, mri->height, mri->depth,
                                    mri->type, mri->nframes);
      MRIcopyHeader(mri, mri_template);

      mri_template->nframes = mri->nframes;

      mri_template->width =
          static_cast<int>(ceil(mri->width / cmdargs.downsample_factor[0]));
      mri_template->height =
          static_cast<int>(ceil(mri->height / cmdargs.downsample_factor[1]));
      mri_template->depth =
          static_cast<int>(ceil(mri->depth / cmdargs.downsample_factor[2]));

      mri_template->xsize *= cmdargs.downsample_factor[0];
      mri_template->ysize *= cmdargs.downsample_factor[1];
      mri_template->zsize *= cmdargs.downsample_factor[2];

      mri_template->xstart = -mri_template->xsize * mri_template->width / 2;
      mri_template->xend   = mri_template->xsize * mri_template->width / 2;
      mri_template->ystart = -mri_template->ysize * mri_template->height / 2;
      mri_template->yend   = mri_template->ysize * mri_template->height / 2;
      mri_template->zstart = -mri_template->zsize * mri_template->depth / 2;
      mri_template->zend   = mri_template->zsize * mri_template->depth / 2;
    }
  }

  /* ----- apply command-line parameters ----- */
  if (cmdargs.out_i_size_flag) {
    float scale;
    scale               = mri_template->xsize / cmdargs.out_i_size;
    mri_template->xsize = cmdargs.out_i_size;
    mri_template->width = nint(mri_template->width * scale);
  }
  if (cmdargs.out_j_size_flag) {
    float scale;
    scale                = mri_template->ysize / cmdargs.out_j_size;
    mri_template->ysize  = cmdargs.out_j_size;
    mri_template->height = nint(mri_template->height * scale);
  }
  if (cmdargs.out_k_size_flag) {
    float scale;
    scale               = mri_template->zsize / cmdargs.out_k_size;
    mri_template->zsize = cmdargs.out_k_size;
    mri_template->depth = nint(mri_template->depth * scale);
  }
  if (cmdargs.out_n_i_flag) {
    mri_template->width = cmdargs.out_n_i;
  }
  if (cmdargs.out_n_j_flag) {
    mri_template->height = cmdargs.out_n_j;
  }
  if (cmdargs.out_n_k_flag) {
    mri_template->depth = cmdargs.out_n_k;
    mri_template->imnr1 = mri_template->imnr0 + cmdargs.out_n_k - 1;
  }
  if (cmdargs.out_i_direction_flag) {
    mri_template->x_r = cmdargs.out_i_directions[0];
    mri_template->x_a = cmdargs.out_i_directions[1];
    mri_template->x_s = cmdargs.out_i_directions[2];
  }
  if (cmdargs.out_j_direction_flag) {
    mri_template->y_r = cmdargs.out_j_directions[0];
    mri_template->y_a = cmdargs.out_j_directions[1];
    mri_template->y_s = cmdargs.out_j_directions[2];
  }
  if (cmdargs.out_k_direction_flag) {
    mri_template->z_r = cmdargs.out_k_directions[0];
    mri_template->z_a = cmdargs.out_k_directions[1];
    mri_template->z_s = cmdargs.out_k_directions[2];
  }
  if (!cmdargs.out_orientation_string.empty()) {
    fmt::printf("Setting output orientation to %s\n",
                cmdargs.out_orientation_string.data());
    MRIorientationStringToDircos(mri_template,
                                 cmdargs.out_orientation_string.data());
  }
  if (!cmdargs.out_center.empty()) {
    mri_template->c_r = cmdargs.out_center[0];
    mri_template->c_a = cmdargs.out_center[1];
    mri_template->c_s = cmdargs.out_center[2];
  }

  /* ----- correct starts, ends, and fov if necessary ----- */
  if ((cmdargs.out_i_size_flag) || (cmdargs.out_j_size_flag) ||
      (cmdargs.out_k_size_flag) || (cmdargs.out_n_i_flag) ||
      (cmdargs.out_n_j_flag) || (cmdargs.out_n_k_flag)) {

    fov_x                = mri_template->xsize * mri_template->width;
    fov_y                = mri_template->ysize * mri_template->height;
    fov_z                = mri_template->zsize * mri_template->depth;
    mri_template->xend   = fov_x / 2.0;
    mri_template->xstart = -mri_template->xend;
    mri_template->yend   = fov_y / 2.0;
    mri_template->ystart = -mri_template->yend;
    mri_template->zend   = fov_z / 2.0;
    mri_template->zstart = -mri_template->zend;

    mri_template->fov = (fov_x > fov_y ? (fov_x > fov_z ? fov_x : fov_z)
                                       : (fov_y > fov_z ? fov_y : fov_z));
  }

  /* ----- give a warning for non-orthogonal directions ----- */
  i_dot_j = mri_template->x_r * mri_template->y_r +
            mri_template->x_a * mri_template->y_a +
            mri_template->x_s * mri_template->y_s;
  i_dot_k = mri_template->x_r * mri_template->z_r +
            mri_template->x_a * mri_template->z_a +
            mri_template->x_s * mri_template->z_s;
  j_dot_k = mri_template->y_r * mri_template->z_r +
            mri_template->y_a * mri_template->z_a +
            mri_template->y_s * mri_template->z_s;
  if (fabs(i_dot_j) > CLOSE_ENOUGH || fabs(i_dot_k) > CLOSE_ENOUGH ||
      fabs(i_dot_k) > CLOSE_ENOUGH) {
    fmt::printf("warning: output volume axes are not orthogonal:"
                "i_dot_j = %.6f, i_dot_k = %.6f, j_dot_k = %.6f\n",
                i_dot_j, i_dot_k, j_dot_k);
    fmt::printf("i_ras = (%g, %g, %g)\n", mri_template->x_r, mri_template->x_a,
                mri_template->x_s);
    fmt::printf("j_ras = (%g, %g, %g)\n", mri_template->y_r, mri_template->y_a,
                mri_template->y_s);
    fmt::printf("k_ras = (%g, %g, %g)\n", mri_template->z_r, mri_template->z_a,
                mri_template->z_s);
  }
  if (cmdargs.out_data_type >= 0) {
    mri_template->type = cmdargs.out_data_type;
  }

  /* ----- catch the mri_template info flag ----- */
  if (cmdargs.template_info_flag) {
    fmt::printf("template structure:\n");
    MRIdump(mri_template, stdout);
  }

  /* ----- exit here if read only is desired ----- */
  if (cmdargs.read_only_flag) {
    return 0;
  }

  /* ----- apply a transformation if requested ----- */
  if (cmdargs.transform_flag) {
    fmt::printf("INFO: Applying transformation from file %s...\n",
                cmdargs.transform_fname.data());
    transform_type = TransformFileNameType(cmdargs.transform_fname.data());
    if (transform_type == MNI_TRANSFORM_TYPE ||
        transform_type == TRANSFORM_ARRAY_TYPE ||
        transform_type == REGISTER_DAT || transform_type == FSLREG_TYPE) {

      // in these cases the lta_transform was loaded above
      if (lta_transform == nullptr) {
        fmt::fprintf(stderr, "ERROR: lta should have been loaded \n");
        return 1;
      }

      /* Think about calling MRIlinearTransform() here; need vox2vox
         transform. Can create NN version. In theory, LTAtransform()
         can handle multiple transforms, but the inverse assumes only
         one. NN is good for ROI*/

      fmt::printf("---------------------------------\n");
      fmt::printf("INFO: Transform Matrix (%s)\n",
                  LTAtransformTypeName(lta_transform->type));
      MatrixPrint(stdout, lta_transform->xforms[0].m_L);
      fmt::printf("---------------------------------\n");

      /* LTAtransform() runs either MRIapplyRASlinearTransform()
         for RAS2RAS or MRIlinearTransform() for Vox2Vox. */
      /* MRIlinearTransform() calls MRIlinearTransformInterp() */
      if (cmdargs.out_like_flag) {
        MRI *tmp = nullptr;
        fmt::printf(
            "INFO: transform dst into the like-volume (resample_type %d)\n",
            cmdargs.resample_type_val);
        tmp = MRIreadHeader(cmdargs.out_like_name.data(),
                            MRI_VOLUME_TYPE_UNKNOWN);
        mri_transformed =
            MRIalloc(tmp->width, tmp->height, tmp->depth, mri->type);
        if (mri_transformed == nullptr) {
          ErrorExit(ERROR_NOMEMORY, "could not allocate memory");
        }
        MRIcopyHeader(tmp, mri_transformed);
        MRIfree(&tmp);
        tmp             = nullptr;
        mri_transformed = LTAtransformInterp(
            mri, mri_transformed, lta_transform, cmdargs.resample_type_val);
      } else {
        if (!cmdargs.out_center.empty()) {
          MATRIX *m;
          MATRIX *mtmp;
          m    = MRIgetResampleMatrix(mri_template, mri);
          mtmp = MRIrasXformToVoxelXform(mri, mri, lta_transform->xforms[0].m_L,
                                         nullptr);
          MatrixFree(&lta_transform->xforms[0].m_L);
          lta_transform->xforms[0].m_L = MatrixMultiply(m, mtmp, NULL);
          MatrixFree(&m);
          MatrixFree(&mtmp);
          lta_transform->type = LINEAR_VOX_TO_VOX;
        }

        fmt::printf("Applying LTAtransformInterp (resample_type %d)\n",
                    cmdargs.resample_type_val);
        if (Gdiag_no > 0) {
          fmt::printf(
              "Dumping LTA ---------++++++++++++-----------------------\n");
          LTAdump(stdout, lta_transform);
          fmt::printf(
              "========-------------++++++++++++-------------------=====\n");
        }
        mri_transformed = LTAtransformInterp(mri, nullptr, lta_transform,
                                             cmdargs.resample_type_val);
      } // end out_like_flag treatment

      if (mri_transformed == nullptr) {
        fmt::fprintf(stderr, "ERROR: applying transform to volume\n");
        return 1;
      }
      if (!cmdargs.out_center.empty()) {
        mri_transformed->c_r = mri_template->c_r;
        mri_transformed->c_a = mri_template->c_a;
        mri_transformed->c_s = mri_template->c_s;
      }
      LTAfree(&lta_transform);
      MRIfree(&mri);
      mri = mri_transformed;
    } else if (transform_type == MORPH_3D_TYPE) {
      fmt::printf("Applying morph_3d ...\n");
      // this is a non-linear vox-to-vox transform
      // note that in this case trilinear
      // interpolation is always used, and -rt
      // option has no effect! -xh
      TRANSFORM *tran = TransformRead(cmdargs.transform_fname.data());
      // check whether the volume to be morphed and the morph have the same
      // dimensions
      if (tran == nullptr) {
        ErrorExit(ERROR_NOFILE, "%s: could not read xform from %s\n", Progname,
                  cmdargs.transform_fname.c_str());
      }
      if (!cmdargs.invert_transform_flag) {
        fmt::printf("morphing to atlas with resample type %d\n",
                    cmdargs.resample_type_val);
        mri_transformed =
            GCAMmorphToAtlas(mri, static_cast<GCA_MORPH *>(tran->xform),
                             nullptr, 0, cmdargs.resample_type_val);
        useVolGeomToMRI(&(static_cast<GCA_MORPH *>(tran->xform))->atlas,
                        mri_template);
      } else // invert
      {
        auto *gcam = static_cast<GCA_MORPH *>(tran->xform);

        // check whether the volume to be morphed and the morph have the same
        // dimensions
#if 0
        if ((gcam->image.width == mri->width) &&
            (gcam->image.height == mri->height) &&
            (gcam->image.depth == mri->depth))
          mri_transformed = MRIclone(mri, NULL);
        else // when morphing atlas to subject using -ait <3d morph> must
             // resample atlas to 256 before applying morph
        {
          MRI *mri_tmp;
          mri_transformed = MRIalloc(gcam->image.width, gcam->image.height,
                                     gcam->image.depth, mri->type);
          useVolGeomToMRI(&gcam->image, mri_transformed);
          mri_tmp = MRIresample(mri, mri_transformed, resample_type_val);
          MRIfree(&mri);
          MRIfree(&mri_template);
          mri = mri_tmp;
        }
#endif
        fmt::printf("morphing from atlas with resample type %d\n",
                    cmdargs.resample_type_val);
        mri_transformed = GCAMmorphFromAtlas(mri, gcam, mri_transformed,
                                             cmdargs.resample_type_val);
        MRIwrite(mri_transformed, "t.mgz");
      }
      TransformFree(&tran);
      MRIfree(&mri);
      mri = mri_transformed;
    } else {
      fmt::fprintf(stderr, "unknown transform type in file %s\n",
                   cmdargs.transform_fname.data());
      return 1;
    }
  } else if (((cmdargs.out_like_flag)) &&
             (!cmdargs.out_stats_table_flag)) // flag set but no transform
  {
    // modify direction cosines to use the
    // out-like volume direction cosines
    //
    //   src --> RAS
    //    |       |
    //    |       |
    //    V       V
    //   dst --> RAS   where dst i_to_r is taken from out volume
    MRI *   tmp     = nullptr;
    MATRIX *src2dst = nullptr;

    fmt::printf("INFO: transform src into the like-volume: %s\n",
                cmdargs.out_like_name.data());
    tmp = MRIreadHeader(cmdargs.out_like_name.data(), MRI_VOLUME_TYPE_UNKNOWN);
    mri_transformed = MRIalloc(tmp->width, tmp->height, tmp->depth, mri->type);
    if (mri_transformed == nullptr) {
      ErrorExit(ERROR_NOMEMORY, "could not allocate memory");
    }
    MRIcopyHeader(tmp, mri_transformed);
    MRIfree(&tmp);
    tmp = nullptr;
    // just to make sure
    MATRIX *tmp2;
    if ((mri->i_to_r__) == nullptr) {
      tmp2 = extract_i_to_r(mri);
      AffineMatrixAlloc(&(mri->i_to_r__));
      SetAffineMatrix(mri->i_to_r__, tmp2);
    } else {
      tmp2 = MatrixAlloc(4, 4, MATRIX_REAL);
      GetAffineMatrix(tmp2, mri->i_to_r__);
    }

    if ((mri_transformed->r_to_i__) == nullptr) {
      mri_transformed->r_to_i__ = extract_r_to_i(mri_transformed);
    }
    // got the transform
    src2dst = MatrixMultiply(mri_transformed->r_to_i__, tmp2, NULL);
    // now get the values (tri-linear)
    MRIlinearTransform(mri, mri_transformed, src2dst);
    MatrixFree(&src2dst);
    MatrixFree(&tmp2);
    MRIfree(&mri);
    mri = mri_transformed;
  }

  /* ----- change type if necessary ----- */
  if (mri->type != mri_template->type && !cmdargs.nochange_flag) {
    fmt::printf("changing data type from %s to %s (noscale = %d)...\n",
                MRItype2str(mri->type), MRItype2str(mri_template->type),
                cmdargs.no_scale_flag);
    mri2 = MRISeqchangeType(mri, mri_template->type, 0.0, 0.999,
                            static_cast<int>(cmdargs.no_scale_flag));
    if (mri2 == nullptr) {
      fmt::printf("ERROR: MRISeqchangeType\n");
      return 1;
    }
    /* (mr) now always map 0 to 0 (was only in conform case before)
       should this be part of MRISeqchangeType? */
    /* (dng) this used to use MRImask, but that did not work for multiple
     * frame*/
    /* Neither mr nor dng can remember where a problem came up that required
     * this fix*/

    fs::mri::new_vox_getter vox_getter =
        fs::mri::get_typed_new_vox_getter_chunked(mri);
    fs::mri::new_vox_setter vox_setter =
        fs::mri::get_typed_new_vox_setter_chunked(mri2);

    if (mri->ischunked) {
      for (size_t index{0}; index < mri->vox_total; index++) {
        v = vox_getter(mri, index);
        if ((v == 0)) {
          vox_setter(mri2, index, 0);
        }
      }
    } else {
      for (c = 0; c < mri->width; c++) {
        for (r = 0; r < mri->height; r++) {
          for (s = 0; s < mri->depth; s++) {
            for (f = 0; f < mri->nframes; f++) {
              v = MRIgetVoxVal(mri, c, r, s, f);
              if (v == 0) {
                MRIsetVoxVal(mri2, c, r, s, f, 0);
              }
            }
          }
        }
      }
    }
    MRIfree(&mri);
    mri = mri2;
  }

#if 0
  printf("mri  vs.  mri_template\n");
  fmt::printf(" dim  ( %i , %i , %i ) vs ( %i , %i , %i)\n", mri->width,
              mri->height, mri->depth, mri_template->width,
              mri_template->height, mri_template->depth);
  fmt::printf(" size ( %.9g , %.9g , %.9g ) vs ( %.9g , %.9g , %.9g)\n",
              mri->xsize, mri->ysize, mri->zsize, mri_template->xsize,
              mri_template->ysize, mri_template->zsize);
  fmt::printf(" xras ( %.9g , %.9g , %.9g ) vs ( %.9g , %.9g , %.9g)\n",
              mri->x_r, mri->x_a, mri->x_s, mri_template->x_r,
              mri_template->x_a, mri_template->x_s);
  fmt::printf(" yras ( %.9g , %.9g , %.9g ) vs ( %.9g , %.9g , %.9g)\n",
              mri->y_r, mri->x_a, mri->y_s, mri_template->y_r,
              mri_template->y_a, mri_template->y_s);
  fmt::printf(" zras ( %.9g , %.9g , %.9g ) vs ( %.9g , %.9g , %.9g)\n",
              mri->z_r, mri->x_a, mri->z_s, mri_template->z_r,
              mri_template->z_a, mri_template->z_s);
  fmt::printf(" cras ( %.9g , %.9g , %.9g ) vs ( %.9g , %.9g , %.9g)\n",
              mri->c_r, mri->c_a, mri->c_s, mri_template->c_r,
              mri_template->c_a, mri_template->c_s);
#endif

  /* ----- reslice if necessary and not performed during transform ----- */
  constexpr auto eps =
      static_cast<const float>(1e-05); /* (mr) do eps testing to avoid reslicing
                        due to tiny differences, e.g. from IO */
  if ((!cmdargs.out_like_flag) && (transform_type != MORPH_3D_TYPE) &&
      (fabs(mri->xsize - mri_template->xsize) > eps ||
       fabs(mri->ysize - mri_template->ysize) > eps ||
       fabs(mri->zsize - mri_template->zsize) > eps ||
       mri->width != mri_template->width ||
       mri->height != mri_template->height ||
       mri->depth != mri_template->depth ||
       fabs(mri->x_r - mri_template->x_r) > eps ||
       fabs(mri->x_a - mri_template->x_a) > eps ||
       fabs(mri->x_s - mri_template->x_s) > eps ||
       fabs(mri->y_r - mri_template->y_r) > eps ||
       fabs(mri->y_a - mri_template->y_a) > eps ||
       fabs(mri->y_s - mri_template->y_s) > eps ||
       fabs(mri->z_r - mri_template->z_r) > eps ||
       fabs(mri->z_a - mri_template->z_a) > eps ||
       fabs(mri->z_s - mri_template->z_s) > eps ||
       fabs(mri->c_r - mri_template->c_r) > eps ||
       fabs(mri->c_a - mri_template->c_a) > eps ||
       fabs(mri->c_s - mri_template->c_s) > eps)) {
    fmt::printf("Reslicing using ");
    switch (cmdargs.resample_type_val) {
    case SAMPLE_TRILINEAR:
      fmt::printf("trilinear interpolation \n");
      break;
    case SAMPLE_NEAREST:
      fmt::printf("nearest \n");
      break;
      /*    case SAMPLE_SINC:
            fmt::printf("sinc \n");
            break;*/
    case SAMPLE_CUBIC:
      fmt::printf("cubic \n");
      break;
    case SAMPLE_WEIGHTED:
      fmt::printf("weighted \n");
      break;
    case SAMPLE_VOTE:
      fmt::printf("voting \n");
      break;
    }
    mri2 = MRIresample(mri, mri_template, cmdargs.resample_type_val);
    if (mri2 == nullptr) {
      return 1;
    }
    MRIfree(&mri);
    mri = mri2;
  }

  /* ----- invert contrast if necessary ----- */
  if (invert_val >= 0) {
    fmt::printf("inverting contrast...\n");
    mri2 = MRIinvertContrast(mri, nullptr, invert_val);
    if (mri2 == nullptr) {
      return 1;
    }
    MRIfree(&mri);
    mri = mri2;
  }

  /* ----- reorder if necessary ----- */
  if (!cmdargs.reorder_vals.empty()) {
    fmt::printf("reordering axes...\n");
    mri2 = MRIreorder(mri, nullptr, cmdargs.reorder_vals[0],
                      cmdargs.reorder_vals[1], cmdargs.reorder_vals[2]);
    if (mri2 == nullptr) {
      fmt::fprintf(stderr, "error reordering axes\n");
      return 1;
    }
    MRIfree(&mri);
    mri = mri2;
  }

  /* ----- reorder if necessary ----- */
  if (!cmdargs.reorder4_vals.empty()) {
    fmt::printf("reordering all axes...\n");
    fmt::printf("reordering, values are %d %d %d %d\n",
                cmdargs.reorder4_vals[0], cmdargs.reorder4_vals[1],
                cmdargs.reorder4_vals[2], cmdargs.reorder4_vals[3]);
    mri2 = MRIreorder4(mri, cmdargs.reorder4_vals.data());
    if (mri2 == nullptr) {
      fmt::fprintf(stderr, "error reordering all axes\n");
      return 1;
    }
    MRIfree(&mri);
    mri = mri2;
  }

  /* ----- store the gdf file stem ----- */
  strcpy(mri->gdf_image_stem, cmdargs.gdf_image_stem.data());

  /* ----- catch the out info flag ----- */
  if (cmdargs.out_info_flag) {
    fmt::printf("output structure:\n");
    MRIdump(mri, stdout);
  }

  /* ----- catch the out matrix flag ----- */
  if (cmdargs.out_matrix_flag) {
    MATRIX *i_to_r;
    i_to_r = extract_i_to_r(mri);
    if (i_to_r != nullptr) {
      fmt::printf("output ijk -> ras:\n");
      MatrixPrint(stdout, i_to_r);
      MatrixFree(&i_to_r);
    } else {
      fmt::printf("error getting output matrix\n");
    }
  }

  /* ----- catch the out stats flag ----- */
  if (cmdargs.out_stats_flag) {
    MRIprintStats(mri, stdout);
  }

  if (cmdargs.erode_seg_flag) {
    fmt::printf("Eroding segmentation %d\n", cmdargs.n_erode_seg);
    mri2 = MRIerodeSegmentation(mri, nullptr, cmdargs.n_erode_seg, 0);
    if (mri2 == nullptr) {
      return 1;
    }
    MRIfree(&mri);
    mri = mri2;
  }
  if (cmdargs.dil_seg_flag) {
    mritmp = nullptr;
    if (cmdargs.dil_seg_mask[0] == '\0') {
      fmt::printf("Dilating segmentation %d\n", cmdargs.n_dil_seg);
    } else {
      fmt::printf("Dilating segmentation, mask %s\n",
                  cmdargs.dil_seg_mask.data());
      mritmp = MRIread(cmdargs.dil_seg_mask.data());
      if (mritmp == nullptr) {
        fmt::printf("ERROR: could not read %s\n", cmdargs.dil_seg_mask.data());
        return 1;
      }
    }
    mri2 = MRIdilateSegmentation(mri, nullptr, cmdargs.n_dil_seg, mritmp, 0,
                                 0.5, &i);
    if (mri2 == nullptr) {
      return 1;
    }
    MRIfree(&mri);
    mri = mri2;
  }

  if (cmdargs.frame_flag) {
    if (cmdargs.mid_frame_flag) {
      nframes           = 1;
      cmdargs.frames[0] = nint(mri->nframes / 2);
    }
    if (nframes == 1) {
      fmt::printf("keeping frame %d\n", cmdargs.frames[0]);
    } else {
      fmt::printf("keeping frames");
      for (f = 0; f < nframes; f++) {
        fmt::printf(" %d", cmdargs.frames[f]);
      }
      fmt::printf("\n");
    }
    mri2 = MRIallocSequence(mri->width, mri->height, mri->depth, mri->type,
                            nframes);
    if (mri2 == nullptr) {
      ErrorExit(ERROR_NOMEMORY, "could not allocate memory");
    }
    MRIcopyHeader(mri, mri2);
    MRIcopyPulseParameters(mri, mri2);
    for (f = 0; f < nframes; f++) {
      if (cmdargs.frames[f] < 0 || cmdargs.frames[f] >= mri->nframes) {
        fmt::printf("   ERROR: valid frame numbers are between 0 and %d\n",
                    mri->nframes - 1);
        return 1;
      } else {
        for (s = 0; s < mri2->depth; s++) {
          for (r = 0; r < mri2->height; r++) {
            for (c = 0; c < mri2->width; c++) {
              MRIsetVoxVal(mri2, c, r, s, f,
                           MRIgetVoxVal(mri, c, r, s, cmdargs.frames[f]));
            }
          }
        }
      }
    }
    MRIfree(&mri);
    mri = mri2;
  }

  /* ----- drop the last ndrop frames (drop before skip) ------ */
  if (cmdargs.ndrop > 0) {
    fmt::printf("Dropping last %d frames\n", cmdargs.ndrop);
    if (mri->nframes <= cmdargs.ndrop) {
      fmt::printf("   ERROR: can't drop, volume only has %d frames\n",
                  mri->nframes);
      return 1;
    }
    mri2 = fMRIndrop(mri, cmdargs.ndrop, nullptr);
    if (mri2 == nullptr) {
      return 1;
    }
    MRIfree(&mri);
    mri = mri2;
  }

  /* ----- skip the first nskip frames ------ */
  if (cmdargs.nskip > 0) {
    fmt::printf("Skipping %d frames\n", cmdargs.nskip);
    if (mri->nframes <= cmdargs.nskip) {
      fmt::printf("   ERROR: can't skip, volume only has %d frames\n",
                  mri->nframes);
      return 1;
    }
    mri2 = fMRInskip(mri, cmdargs.nskip, nullptr);
    if (mri2 == nullptr) {
      return 1;
    }
    MRIfree(&mri);
    mri = mri2;
  }
  if (cmdargs.subsample_flag) {
    fmt::printf("SubSample: Start = %d  Delta = %d, End = %d\n",
                cmdargs.SubSampStart, cmdargs.SubSampDelta, cmdargs.SubSampEnd);
    mri2 = fMRIsubSample(mri, cmdargs.SubSampStart, cmdargs.SubSampDelta,
                         cmdargs.SubSampEnd, nullptr);
    if (mri2 == nullptr) {
      return 1;
    }
    MRIfree(&mri);
    mri = mri2;
  }

  if (cmdargs.cutends_flag) {
    fmt::printf("Cutting ends: n = %d\n", cmdargs.ncutends);
    mri2 = MRIcutEndSlices(mri, cmdargs.ncutends);
    if (mri2 == nullptr) {
      return 1;
    }
    MRIfree(&mri);
    mri = mri2;
  }

  /* ----- crops ---------*/
  if (cmdargs.crop_flag) {
    MRI *mri_tmp;
    int  x0;
    int  y0;
    int  z0;

    x0 = cmdargs.crop_center[0];
    y0 = cmdargs.crop_center[1];
    z0 = cmdargs.crop_center[2];

    x0 -= cmdargs.cropsize[0] / 2;
    y0 -= cmdargs.cropsize[1] / 2;
    z0 -= cmdargs.cropsize[2] / 2;
    mri_tmp = MRIextract(mri, nullptr, x0, y0, z0, cmdargs.cropsize[0],
                         cmdargs.cropsize[1], cmdargs.cropsize[2]);
    MRIfree(&mri);
    mri = mri_tmp;
  }

  if (!FEQUAL(cmdargs.out_scale_factor, 1.0)) {
    fmt::printf("scaling output intensities by %2.3f\n",
                cmdargs.out_scale_factor);
    MRIscalarMul(mri, mri, cmdargs.out_scale_factor);
  }

  if (cmdargs.sphinx_flag) {
    fmt::printf("Changing orientation information to Sphinx\n");
    MRIhfs2Sphinx(mri);
  }

  if (cmdargs.AutoAlign != nullptr) {
    mri->AutoAlign = cmdargs.AutoAlign;
  }

  // ----- modify color lookup table (specified by --ctab option) -----
  if (cmdargs.colortablefile == "remove") {
    // remove an embedded ctab
    if (mri->ct != nullptr) {
      std::cout << "removing color lookup table" << std::endl;
      CTABfree(&mri->ct);
    }
  } else if (!cmdargs.colortablefile.empty()) {
    // add a user-specified ctab
    std::cout << "embedding color lookup table" << std::endl;
    if (mri->ct != nullptr) {
      CTABfree(&mri->ct);
    }
    mri->ct = CTABreadASCII(cmdargs.colortablefile.data());
    if (mri->ct == nullptr) {
      fs::fatal() << "could not read lookup table from "
                  << cmdargs.colortablefile.data();
    }
  }

  /*------ Finally, write the output -----*/

  if (cmdargs.out_stats_table_flag) {
    fmt::printf("Writing as Stats-Table to %s using template %s\n",
                cmdargs.out_name.data(), cmdargs.out_like_name.data());

    OutStatTable = InitStatTableFromMRI(mri, cmdargs.out_like_name.data());
    WriteStatTable(cmdargs.out_name.data(), OutStatTable);

    fmt::printf("done\n");
    return 0;
  }

  if (cmdargs.ascii_flag != 0) {
    fmt::printf("Writing as ASCII to %s\n", cmdargs.out_name.data());
    fptmp = fopen(cmdargs.out_name.data(), "we");
    if (cmdargs.ascii_flag == 1 || cmdargs.ascii_flag == 2) {
      for (f = 0; f < mri->nframes; f++) {
        for (s = 0; s < mri->depth; s++) {
          for (r = 0; r < mri->height; r++) {
            for (c = 0; c < mri->width; c++) {
              if (cmdargs.ascii_flag == 1) {
                fmt::fprintf(fptmp, "%lf\n", MRIgetVoxVal(mri, c, r, s, f));
              }
              if (cmdargs.ascii_flag == 2) {
                fmt::fprintf(fptmp, "%3d %3d %3d %3d %lf\n", c, r, s, f,
                             MRIgetVoxVal(mri, c, r, s, f));
              }
            }
          }
        }
      }
    }
    if (cmdargs.ascii_flag == 3) {
      for (s = 0; s < mri->depth; s++) {
        for (r = 0; r < mri->height; r++) {
          for (c = 0; c < mri->width; c++) {
            for (f = 0; f < mri->nframes; f++) {
              fmt::fprintf(fptmp, "%lf", MRIgetVoxVal(mri, c, r, s, f));
              if (f != mri->nframes - 1) {
                fmt::fprintf(fptmp, " ");
              }
            }
            fmt::fprintf(fptmp, "\n");
          }
        }
      }
    }
    fclose(fptmp);
    fmt::printf("done\n");
    return 0;
  }
  if (!cmdargs.no_write_flag) {
    if (!cmdargs.split_frames_flag) {
      fmt::printf("writing to %s...\n", cmdargs.out_name.data());
      if (cmdargs.force_out_type_flag) {
        err = MRIwriteType(mri, cmdargs.out_name.data(), out_volume_type);
        if (err != NO_ERROR) {
          fmt::printf("ERROR: failure writing %s as volume type %d\n",
                      cmdargs.out_name.data(), out_volume_type);
          return 1;
        }
      } else {
        err = MRIwrite(mri, cmdargs.out_name.data());
        if (err != NO_ERROR) {
          fmt::printf("ERROR: failure writing %s\n", cmdargs.out_name.data());
          return 1;
        }
      }
    } else {
      stem = IDstemFromName(cmdargs.out_name.data());
      ext  = IDextensionFromName(cmdargs.out_name.data());

      fmt::printf("splitting frames, stem = %s, ext = %s\n", stem, ext);
      mri2 = nullptr;
      for (i = 0; i < mri->nframes; i++) {
        mri2   = MRIcopyFrame(mri, mri2, i, 0);
        tmpstr = fmt::sprintf("%s%04d.%s", stem, i, ext);
        fmt::printf("%2d %s\n", i, tmpstr.data());
        err = MRIwrite(mri2, tmpstr.data());
        if (err != NO_ERROR) {
          fmt::printf("ERROR: failure writing %s\n", tmpstr.data());
          return 1;
        }
      }
    }
  }
  // free memory
  // MRIfree(&mri); // This causes a seg fault with change of type
  MRIfree(&mri_template);

  return 0;
}
