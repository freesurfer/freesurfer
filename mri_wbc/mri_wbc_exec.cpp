/*
 * Original Author: Douglas N. Greve
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

#include "mri_wbc.hpp"

#include <boost/filesystem.hpp>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "diag.h"
#include "fmriutils.h"
#include "fsenv.h"
#include "matfile.h"
#include "mri2.h"
#include "mrimorph.h"
#include "romp_support.h"
#include "version.h"

namespace fsys = boost::filesystem; // sadly fs is already defined: freesurfer

int main(int argc, char *argv[]) {
  int  nargs, err, n;
  WBC *wbc;

  cmdargs             = (CMDARGS *)calloc(sizeof(CMDARGS), 1);
  cmdargs->distthresh = 10;
  cmdargs->DoDist     = 0;
  cmdargs->DoMat      = 0;
  cmdargs->nrholist   = 0;
  cmdargs->ForceFail  = 0;

  nargs = handleVersionOption(argc, argv, "mri_wbc");
  if (nargs && argc - nargs == 1)
    exit(0);
  argc -= nargs;
  cmdline = argv2cmdline(argc, argv);
  uname(&uts);
  getcwd(cwd, 2000);

  Progname = argv[0];
  argc--;
  argv++;
  ErrorInit(NULL, NULL, NULL);
  DiagInit(NULL, NULL, NULL);
  if (argc == 0)
    usage_exit();
  parse_commandline(argc, argv);
  check_options();
  if (checkoptsonly)
    return (0);
  dump_options(stdout);

  if (cmdargs->outdir) {
    err = mkdir(cmdargs->outdir, 0777);
    if (err != 0 && errno != EEXIST) {
      printf("ERROR: creating directory %s\n", cmdargs->outdir);
      perror(NULL);
      exit(1);
    }
  }

  // TODO(aboualiaa): do you really need this copying?
  // yes you do, cebause the workers below take a reference to wbc and access
  // its contents
  // TODO(aboualiaa): fix this dependency
  wbc.distthresh = cmdargs.distthresh;
  wbc.DoDist     = cmdargs.DoDist;
  wbc.DoMat      = cmdargs.DoMat;
  wbc.rholist    = cmdargs.rholist;
  wbc.DoTest     = cmdargs.DoTest;

  if (cmdargs.DoTest) {
    WBCtestSynth(&wbc);
  }

  err = WBCnframes(&wbc);
  if (err != 0) {
    err_logger->critical("WBCnframes failed");
    exit(1);
  }

  spdlog::debug("preapre, work, then finish");

  WBCprep(&wbc);
  WholeBrainCon(&wbc);
  WBCfinish(&wbc);

  if (wbc.DoMat) {
    spdlog::info("Writing matfile {}", cmdargs.matfile);
    err = MatlabWrite(wbc.M, cmdargs.matfile, "M");
    if (err != 0) {
      err_logger->critical("couldn't write matfile {}", cmdargs.matfile);
      exit(1);
    }
  }

  if (wbc.DoTest) {
    err = WBCtestCheck(&wbc);
    spdlog::error(err);
  }

  if (!cmdargs.volcon.empty()) {
    err = MRIwrite(wbc.volcon, cmdargs.volcon);
    if (err != 0) {
      err_logger->critical("failed to write mri {}", cmdargs.volcon);
      exit(1);
    }
    if (wbc.DoDist) {
      err = MRIwrite(wbc.volconS, cmdargs.volconS);
      if (err != 0) {
        err_logger->critical("failed to write mri {}", cmdargs.volconS);
        exit(1);
      }
      err = MRIwrite(wbc.volconL, cmdargs.volconL);
      if (err != 0) {
        err_logger->critical("failed to write mri {}", cmdargs.volconL);
        exit(1);
      }
    }
  }
  if (!cmdargs.lhcon.empty()) {
    err = MRIwrite(wbc.lhcon, cmdargs.lhcon);
    if (err != 0) {
      err_logger->critical("failed to write mri {}", cmdargs.lhcon);
      exit(1);
    }
    if (wbc.DoDist) {
      err = MRIwrite(wbc.lhconS, cmdargs.lhconS);
      if (err != 0) {
        err_logger->critical("failed to write mri {}", cmdargs.lhconS);
        exit(1);
      }
      err = MRIwrite(wbc.lhconL, cmdargs.lhconL);
      if (err != 0) {
        err_logger->critical("failed to write mri {}", cmdargs.lhconL);
        exit(1);
      }
    }
  }
  if (!cmdargs.rhcon.empty()) {
    err = MRIwrite(wbc.rhcon, cmdargs.rhcon);
    if (err != 0) {
      err_logger->critical("failed to write mri {}", cmdargs.rhcon);
      exit(1);
    }
    if (wbc.DoDist) {
      err = MRIwrite(wbc.rhconS, cmdargs.rhconS);
      err_logger->critical("failed to write mri {}", cmdargs.rhconS);
      if (err != 0) {
        exit(1);
      }
      err = MRIwrite(wbc.rhconL, cmdargs.rhconL);
      if (err != 0) {
        err_logger->critical("failed to write mri {}", cmdargs.rhconL);
        exit(1);
      }
    }
  }

  if (!cmdargs.volrhomean.empty()) {
    err = MRIwrite(wbc.volrhomean, cmdargs.volrhomean);
    if (err != 0) {
      err_logger->critical("failed to write mri {}", cmdargs.volrhomean);
      exit(1);
    }
  }
  if (!cmdargs.lhrhomean.empty()) {
    err = MRIwrite(wbc.lhrhomean, cmdargs.lhrhomean);
    if (err != 0) {
      err_logger->critical("failed to write mri {}", cmdargs.lhrhomean);
      exit(1);
    }
  }
  if (!cmdargs.rhrhomean.empty()) {
    err = MRIwrite(wbc.rhrhomean, cmdargs.rhrhomean);
    if (err != 0) {
      err_logger->critical("failed to write mri {}", cmdargs.rhrhomean);
      exit(1);
    }
  }

  spdlog::info("finished, check output directory {}", cmdargs.outdir);
  exit(0);
}

// MARK: - Helpers
namespace wbc {

void initArgDesc(podesc *desc, CMDARGS *cmdargs) {

  desc->add_options()                                               /**/
                                                                    /**/
      ("help,h",                                                    /**/
       "print out information on how to use this program and exit") /**/
                                                                    /**/
      ("version,v",                                                 /**/
       "print out version and exit")                                /**/
                                                                    /**/
      ("debug,d",                                                   /**/
       po::bool_switch(&cmdargs->debug),                            /**/
       "turn on debugging")                                         /**/
                                                                    /**/
      ("fvol",                                                      /**/
       po::value<std::string>(&(cmdargs->fvol)),                    /**/
       "functional volume")                                         /**/
                                                                    /**/
      ("outdir,o",                                                  /**/
       po::value<std::string>(&(cmdargs->outdir)),                  /**/
       "output directory")                                          /**/
                                                                    /**/
      ("volmask",                                                   /**/
       po::value<std::string>(&(cmdargs->volmask)),                 /**/
       "mask for functional volume")                                /**/
                                                                    /**/
      ("lh",                                                        /**/
       po::value<std::vector<std::string>>()->multitoken(),         /**/
       "functional surface of left hemisphere")                     /**/
                                                                    /**/
      ("lhmask",                                                    /**/
       po::value<std::string>(&(cmdargs->lhmask)),                  /**/
       "mask for functional surface of left hemisphere")            /**/
                                                                    /**/
      ("lhlabel",                                                   /**/
       po::value<std::string>(&(cmdargs->lhlabel)),                 /**/
       "label mask for functional surface of left hemisphere")      /**/
                                                                    /**/
      ("rh",                                                        /**/
       po::value<std::vector<std::string>>()->multitoken(),         /**/
       "functional surface of right hemisphere")                    /**/
                                                                    /**/
      ("rhmask",                                                    /**/
       po::value<std::string>(&(cmdargs->rhmask)),                  /**/
       "mask for functional surface of right hemisphere")           /**/
                                                                    /**/
      ("rhlabel",                                                   /**/
       po::value<std::string>(&(cmdargs->rhlabel)),                 /**/
       "label mask for functional surface of right hemisphere")     /**/
                                                                    /**/
      ("rho",                                                       /**/
       po::value<std::vector<double>>()->multitoken(),              /**/
       "rhothresh")                                                 /**/
                                                                    /**/
      ("dist",                                                      /**/
       po::value<double>(&(cmdargs->distthresh))                    /**/
           ->default_value(defaultdistthresh),                      /**/
       "distthresh")                                                /**/
                                                                    /**/
      ("volrhomean",                                                /**/
       po::value<std::string>(&(cmdargs->volrhomean)),              /**/
       "volume rhomean")                                            /**/
                                                                    /**/
      ("lhrhomean",                                                 /**/
       po::value<std::string>(&(cmdargs->lhrhomean)),               /**/
       "left hemisphere rhomean")                                   /**/
                                                                    /**/
      ("rhrhomean",                                                 /**/
       po::value<std::string>(&(cmdargs->rhrhomean)),               /**/
       "right hemisphere rhomean")                                  /**/
                                                                    /**/
      ("volcon",                                                    /**/
       po::value<std::string>(&(cmdargs->volcon)),                  /**/
       "volume con")                                                /**/
                                                                    /**/
      ("lhcon",                                                     /**/
       po::value<std::string>(&(cmdargs->lhcon)),                   /**/
       "left hemisphere con")                                       /**/
                                                                    /**/
      ("rhcon",                                                     /**/
       po::value<std::string>(&(cmdargs->rhcon)),                   /**/
       "right hemisphere con")                                      /**/
                                                                    /**/
      ("volconS",                                                   /**/
       po::value<std::string>(&(cmdargs->volconS)),                 /**/
       "volume conS")                                               /**/
                                                                    /**/
      ("lhconS",                                                    /**/
       po::value<std::string>(&(cmdargs->lhconS)),                  /**/
       "left hemisphere conS")                                      /**/
                                                                    /**/
      ("rhconS",                                                    /**/
       po::value<std::string>(&(cmdargs->rhconS)),                  /**/
       "right hemisphere conS")                                     /**/
                                                                    /**/
      ("volconL",                                                   /**/
       po::value<std::string>(&(cmdargs->volconL)),                 /**/
       "volume conL")                                               /**/
                                                                    /**/
      ("lhconL",                                                    /**/
       po::value<std::string>(&(cmdargs->lhconL)),                  /**/
       "left hemisphere conL")                                      /**/
                                                                    /**/
      ("rhconL",                                                    /**/
       po::value<std::string>(&(cmdargs->rhconL)),                  /**/
       "right hemisphere conL")                                     /**/
                                                                    /**/
      ("mat",                                                       /**/
       po::value<std::string>(&(cmdargs->matfile)),                 /**/
       "matrix")                                                    /**/
                                                                    /**/
      ("uppersub-test",                                             /**/
       po::value<int>(),                                            /**/
       "runs Index2UpperSubscriptTest(int) then exits")             /**/
                                                                    /**/
      ("test,t",                                                    /**/
       po::bool_switch(&cmdargs->DoTest),                           /**/
       "testing flag")                                              /**/
                                                                    /**/
      ("test-fail",                                                 /**/
       po::bool_switch(&cmdargs->ForceFail),                        /**/
       "force fail flag")                                           /**/
                                                                    /**/
      ("save-test",                                                 /**/
       po::bool_switch(&cmdargs->SaveTest),                         /**/
       "save test flag")                                            /**/
                                                                    /**/
      ("threads,j",                                                 /**/
       po::value<int>(&(cmdargs->nthreads)),                        /**/
       "number of threads to use")                                  /**/
                                                                    /**/
      ("diag",                                                      /**/
       po::value<int>(&Gdiag_no),                                   /**/
       "diagnostics")                                               /**/
                                                                    /**/
      ("diag-show",                                                 /**/
       "show diagnostics")                                          /**/
                                                                    /**/
      ("diag-verbose",                                              /**/
       "verbose diagnostics")                                       /**/
                                                                    /**/
      ("all-info",                                                  /**/
       "print all info about environment")                          /**/
      ;
}

static bool good_cmdline_args(CMDARGS *cmdargs, ENV *env, WBC *wbc) {
  podesc desc("\nUSAGE: mri_wbc <options> --lh <lhsurface> -o "
              "<outdir>\n\nAvailable Options");
  povm   vm;

  initArgDesc(&desc, cmdargs);
  auto args = cmdargs->raw;
  auto ac   = static_cast<int>(args.size());
  auto av   = args.data();
  if (ac == 1) {
    print_usage(desc, env);
    return false;
  }

  try {
    auto parsed_opts =
        po::command_line_parser(ac, av).options(desc).style(cl_style).run();
    po::store(parsed_opts, vm);
  } catch (std::exception const &e) {
    spdlog::get("stderr")->critical(e.what());
    return false;
  }

  if (vm.count("help") != 0U) {
    print_help(desc, env);
    return false;
  }

  if ((vm.count("version") != 0U) || (vm.count("all-info") != 0U)) {
    handle_version_option(vm.count("all-info") != 0U, args, env->vcid,
                          "$Name:  $");
    return false;
  }

  try {
    po::notify(vm);
  } catch (std::exception const &e) {
    spdlog::get("stderr")->critical(e.what());
    return false;
  }

  if (cmdargs->debug) {
    spdlog::set_level(spdlog::level::debug); // Set global log level to debug
  } else {
    spdlog::set_level(spdlog::level::warn);
  }

  if (vm.count("lh") != 0U) {
    auto opts = vm["lh"].as<std::vector<std::string>>();
    if (!vm["lh"].empty() && (opts.size() == 3 || opts.size() == 2)) {
      cmdargs->flh       = opts[0];
      cmdargs->lhsurface = opts[1];
      if (opts.size() == 3) {
        cmdargs->lhsurface2 = opts[2];
      }
    } else {
      spdlog::get("stderr")->critical(
          "option --lh requires either 2 or 3 arguments, for "
          "example: --lh lhffile lh.surface <lh.inflated>");
      return false;
    }
  }

  if (vm.count("rh") != 0U) {
    auto opts = vm["rh"].as<std::vector<std::string>>();
    if (!vm["rh"].empty() && (opts.size() == 3 || opts.size() == 2)) {
      cmdargs->frh       = opts[0];
      cmdargs->rhsurface = opts[1];
      if (opts.size() == 3) {
        cmdargs->rhsurface2 = opts[2];
      }
    } else {
      spdlog::get("stderr")->critical(
          "option --rh requires either 2 or 3 arguments, for "
          "example: --rh rhffile rh.surface <rh.inflated>");
      return false;
    }
  }

  if (vm.count("rho") != 0U) {
    cmdargs->rholist = vm["rho"].as<std::vector<double>>();
    if (cmdargs->rholist.empty()) {
      cmdargs->rholist.push_back(fifth);
    }
  }

  if (vm.count("dist") != 0U) {
    cmdargs->DoDist = true;
  }

  if (vm.count("mat") != 0U) {
    cmdargs->DoMat = true;
  }

  // TODO(aboualiaa): what does this do?
  if (vm.count("uppersub-test") != 0U) {
    auto n = Index2UpperSubscriptTest(vm["uppersub-test"].as<int>());
    exit(n);
  }

  // if enabled testing indirectly, i.e. not with the --test flag
  if (cmdargs->ForceFail || cmdargs->SaveTest) {
    cmdargs->DoTest = true;
  }

  // if testing the save function without output directory → die
  if (cmdargs->SaveTest) {
    if (vm.count("outdir") == 0U) {
      spdlog::get("stderr")->critical("no output directory specified");
      return false;
    }
    cmdargs->DoTest = true;
  }

  // TODO(aboualiaa): unuglify
  if (shouldSave(cmdargs, vm)) {
    auto vols = {"vol", "lh", "rh"};
    for (auto v : vols) {
      if (v == "vol") {
        if (vm.count("fvol") == 0U) {
          continue;
        }
      } else if (vm.count(v) == 0U) {
        continue;
      }
      int         volNameInt = v == "vol" ? 1 : v == "lh" ? 2 : 3;
      std::string prefix(cmdargs->outdir + "/" + v);
      std::string conName(prefix + ".con" + ".nii.gz");
      std::string conNameS(prefix + ".con" + "S" + ".nii.gz");
      std::string conNameL(prefix + ".con" + "L" + ".nii.gz");
      std::string rhoMean(prefix + ".rho.mean" + ".nii.gz");
      switch (volNameInt) {
      case 1:
        cmdargs->volcon     = conName;
        cmdargs->volrhomean = rhoMean;
        if (cmdargs->DoDist) {
          cmdargs->volconS = conNameS;
          cmdargs->volconL = conNameL;
        }
        break;
      case 2:
        cmdargs->lhcon     = conName;
        cmdargs->lhrhomean = rhoMean;
        if (cmdargs->DoDist) {
          cmdargs->lhconS = conNameS;
          cmdargs->lhconL = conNameL;
        }
        break;
      case 3:
        cmdargs->rhcon     = conName;
        cmdargs->rhrhomean = rhoMean;
        if (cmdargs->DoDist) {
          cmdargs->rhconS = conNameS;
          cmdargs->rhconL = conNameL;
        }
        break;
      default:
        break;
      }
      spdlog::debug(conName);
      spdlog::debug(conNameS);
      spdlog::debug(conNameL);
      spdlog::debug(rhoMean);
    }
  }

  if (!cmdargs->DoTest) {
    // if no input format specified and not testing → die
    if ((vm.count("fvol") == 0U) && (vm.count("lh") == 0U) &&
        (vm.count("rh") == 0U)) {
      spdlog::get("stderr")->critical("ERROR: no input files specified");
      return false;
    }

    // if no output directory specified and not testing → die
    if ((vm.count("outdir") == 0U)) {
      spdlog::get("stderr")->critical("ERROR: no output directory specified\n");
      return false;
    }

    spdlog::debug("checking for volumes to read because not testing");

    if (!cmdargs->fvol.empty()) {
      wbc->fvol = MRIread(cmdargs->fvol.c_str());
      if (wbc->fvol == nullptr) {
        spdlog::get("stderr")->critical("could not read volume {}",
                                        cmdargs->fvol);
        exit(1);
      }
    }
    if (!cmdargs->volmask.empty()) {
      wbc->volmask = MRIread(cmdargs->volmask.c_str());
      if (wbc->volmask == nullptr) {
        spdlog::get("stderr")->critical("could not read volume {}",
                                        cmdargs->volmask);
        exit(1);
      }
    }
    if (!cmdargs->flh.empty()) {
      wbc->flh = MRIread(cmdargs->flh.c_str());
      if (wbc->flh == nullptr) {
        spdlog::get("stderr")->critical("could not read volume {}",
                                        cmdargs->flh);
        exit(1);
      }
      wbc->lh = MRISread(cmdargs->lhsurface.c_str());
      if (wbc->lh == nullptr) {
        spdlog::get("stderr")->critical("could not read volume {}",
                                        cmdargs->lhsurface);
        exit(1);
      }
      if (!cmdargs->lhlabel.empty()) {
        wbc->lhlabel = LabelRead(nullptr, cmdargs->lhlabel.c_str());
        if (wbc->lhlabel == nullptr) {
          spdlog::get("stderr")->critical("could not read volume {}",
                                          cmdargs->lhlabel);
          exit(1);
        }
      }
      if (!cmdargs->lhmask.empty()) {
        wbc->lhmask = MRIread(cmdargs->lhmask.c_str());
        if (wbc->lhmask == nullptr) {
          spdlog::get("stderr")->critical("could not read volume {}",
                                          cmdargs->lhmask);
          exit(1);
        }
      }
      if (!cmdargs->lhsurface2.empty()) {
        wbc->lh2 = MRISread(cmdargs->lhsurface2.c_str());
        if (wbc->lh2 == nullptr) {
          spdlog::get("stderr")->critical("could not read volume {}",
                                          cmdargs->lhsurface2);
          exit(1);
        }
      }
    }
    if (!cmdargs->frh.empty()) {
      wbc->frh = MRIread(cmdargs->frh.c_str());
      if (wbc->frh == nullptr) {
        spdlog::get("stderr")->critical("could not read volume {}",
                                        cmdargs->frh);
        exit(1);
      }
      wbc->rh = MRISread(cmdargs->rhsurface.c_str());
      if (wbc->rh == nullptr) {
        spdlog::get("stderr")->critical("could not read volume {}",
                                        cmdargs->rhsurface);
        exit(1);
      }
      if (!cmdargs->rhlabel.empty()) {
        wbc->rhlabel = LabelRead(nullptr, cmdargs->rhlabel.c_str());
        if (wbc->rhlabel == nullptr) {
          spdlog::get("stderr")->critical("could not read volume {}",
                                          cmdargs->rhlabel);
          exit(1);
        }
      }
      if (!cmdargs->rhmask.empty()) {
        wbc->rhmask = MRIread(cmdargs->rhmask.c_str());
        if (wbc->rhmask == nullptr) {
          spdlog::get("stderr")->critical("could not read volume {}",
                                          cmdargs->rhmask);
          exit(1);
        }
      }
      if (!cmdargs->rhsurface2.empty()) {
        wbc->rh2 = MRISread(cmdargs->rhsurface2.c_str());
        if (wbc->rh2 == nullptr) {
          spdlog::get("stderr")->critical("could not read volume {}",
                                          cmdargs->rhsurface2);
          exit(1);
        }
      }
    }
  } else {
    spdlog::debug("doing test, creating wbcsynth");
    wbc->wbcsynth               = wbc::WBCSYNTH();
    wbc->wbcsynth.volres        = defaultRes;
    wbc->wbcsynth.voldim        = defaultDim;
    wbc->wbcsynth.nframes       = five;
    wbc->wbcsynth.nshorttarg[0] = defaultnshorttarg;
    wbc->wbcsynth.nlongtarg[0]  = defaultnlongtarg;
    wbc->wbcsynth.nshorttarg[1] = defaultnshorttarg1;
    wbc->wbcsynth.nlongtarg[1]  = defaultnlongtarg1;
    wbc->wbcsynth.ForceFail     = static_cast<int>(cmdargs->ForceFail);
    WBCtestSynth(wbc);
  }

  // pass number of desired threads to openMP
  if (vm.count("threads") != 0U) {
#ifdef HAVE_OPENMP
    omp_set_num_threads(cmdargs->nthreads);
#endif
  }

  // TODO(aboualiaa): what does this do?
  if (vm.count("diag-show") != 0U) {
    Gdiag = (Gdiag & DIAG_SHOW);
  }

  // TODO(aboualiaa): what does this do?
  if (vm.count("diag-verbose") != 0U) {
    Gdiag = (Gdiag & DIAG_VERBOSE);
  }
  return true;
}

// MARK: - Workers

MRI *WholeBrainCon(WBC *wbc) {
  std::vector<MRI *>              conth;
  std::vector<MRI *>              conSth;
  std::vector<MRI *>              conLth;
  std::vector<MRI *>              rhomean;
  int                             nthreads = 1;
  int                             nthrho;
  int64                           nperthread0;
  std::vector<int64>              nperthread;
  int64                           ntot;
  int                             threadno; // thread number
  int64                           npairs;
  int64                           ia;
  int64                           ib;
  int                             k1; // note: these are redefined in loop
  Timer                           timer;
  std::vector<double *>           pfold;
  std::vector<std::vector<float>> pf;
  float                           val;
  std::vector<int64>              k1a;
  std::vector<int64>              k1b;
  std::vector<int64>              k2a;
  std::vector<int64>              k2b;

  // Load time courses into a pointer structure
  pf.reserve(wbc->ntot);
  for (k1 = 0; k1 < wbc->ntot; k1++) {
    int t;
    pf[k1].reserve(wbc->f->nframes);
    for (t = 0; t < wbc->f->nframes; t++) {
      pf[k1][t] = MRIgetVoxVal(wbc->fnorm, k1, 0, 0, t);
    }
  }

  // This is for testing. In most apps, M will be huge
  if (wbc->DoMat) {
    wbc->M = MatrixAlloc(wbc->ntot, wbc->ntot, MATRIX_REAL);
  }

#ifdef HAVE_OPENMP
  nthreads = omp_get_max_threads();
#endif
  npairs = static_cast<int64>(wbc->ntot) * (wbc->ntot - 1) / 2;
  nperthread0 =
      static_cast<int64>(round(static_cast<double>(npairs) / nthreads -
                               1)); // don't use nint(), need long
  if (nperthread0 <= 0) {
    nperthread0 = 1;
  }
  fmt::printf("ntot = %d, nthreads = %d, npairs = %ld, nperthread0 = %ld\n",
              wbc->ntot, nthreads, npairs, nperthread0);

  // Compute number of pairs per thread
  nperthread.reserve(nthreads);
  ntot = 0;
  for (threadno = 0; threadno < nthreads; threadno++) {
    nperthread[threadno] = nperthread0; // number of pairs per thread
    ntot += nperthread0;
  }
  if (ntot != npairs) {
    nperthread[nthreads - 1] += npairs - ntot;
  }

  // Allocate con volumes for each thread
  rhomean.reserve(nthreads);
  conth.reserve(nthreads);
  if (wbc->DoDist) {
    conSth.reserve(nthreads);
    conLth.reserve(nthreads);
  }
  for (threadno = 0; threadno < nthreads; threadno++) {
    rhomean[threadno] = MRIallocSequence(wbc->ntot, 1, 1, MRI_FLOAT, 1);
    conth[threadno] =
        MRIallocSequence(wbc->ntot, 1, 1, MRI_FLOAT, wbc->rholist.size());
    if (wbc->DoDist) {
      conSth[threadno] =
          MRIallocSequence(wbc->ntot, 1, 1, MRI_FLOAT, wbc->rholist.size());
      conLth[threadno] =
          MRIallocSequence(wbc->ntot, 1, 1, MRI_FLOAT, wbc->rholist.size());
    }
  }

  // This is just a test
  ntot = 0;
  for (threadno = 0; threadno < nthreads; threadno++) {
    fmt::printf("thread %d %ld\n", threadno, nperthread[threadno]);
    ntot += nperthread[threadno];
  }
  fmt::printf("ntotpairs = %ld vs npairs %ld, diff %ld\n", ntot, npairs,
              ntot - npairs); // should be equal

  k1a.reserve(nthreads);
  k1b.reserve(nthreads);
  k2a.reserve(nthreads);
  k2b.reserve(nthreads);
  ia = 0;
  for (threadno = 0; threadno < nthreads; threadno++) {
    ib = ia + nperthread[threadno];
    Index2UpperSubscript(wbc->ntot, ia, &k1a[threadno], &k2a[threadno]);
    Index2UpperSubscript(wbc->ntot, ib - 1, &k1b[threadno], &k2b[threadno]);
    fmt::printf("thread %d %12ld %12ld   %7d %7d   %7d %7d\n", threadno, ia, ib,
                k1a[threadno], k2a[threadno], k1b[threadno], k2b[threadno]);
    ia = ib;
  }

  fmt::printf("rho thresholds (%d): ", wbc->rholist.size());
  for (auto rho : wbc->rholist) {
    fmt::printf("%lf ", rho);
  }
  fmt::printf("\n");

  fmt::printf("Starting WBC loop\n");

  timer.reset();
  ROMP_PF_begin
#ifdef HAVE_OPENMP
// this is definitly wrong :D
#pragma omp parallel for default(none)                                         \
    shared(nthreads, k1a, k2a, k1b, k2b, wbc, npairs, timer, pf, val, rhomean, \
           conth, conSth, conLth) if_ROMP(experimental)
#endif

      for (threadno = 0; threadno < nthreads; threadno++) {
    ROMP_PFLB_begin int k1;
    int                 k2;
    int                 t;
    float               n;
    int                thno; // thread number, used to save thread number of omp
    int                q;
    float              ct1;
    float              ct2;
    int                k1start;
    int                k1stop;
    int                k2start;
    int                k2stop;
    int                k2min;
    int                k2max;
    int                nthrho;
    int64              nthpair;
    float              rho;
    double             dx;
    double             dy;
    double             dz;
    double             dist;
    double *           pf1old;
    std::vector<float> pf1;
    double *           pf2old;
    std::vector<float> pf2;
    double             x1;
    double             y1;
    double             z1;
    double             x2;
    double             y2;
    double             z2;
    double             rhothresh;
    MRI *              conDth;
    thno = threadno;
#ifdef HAVE_OPENMP
    thno = omp_get_thread_num(); // actual thread number
#endif

    k1start = k1a[thno];
    k1stop  = k1b[thno];
    k2start = k2a[thno];
    k2stop  = k2b[thno];

    nthpair = 0;
    for (k1 = k1start; k1 <= k1stop; k1++) {
      if (k1 == k1start) {
        k2min = k2start;
      } else {
        k2min = k1 + 1;
      }
      if (k1 == k1stop) {
        k2max = k2stop;
      } else {
        k2max = wbc->ntot - 1;
      }
      for (k2 = k2min; k2 <= k2max; k2++) {
        nthpair++;
        if (thno == 0 && (nthpair == 1 ||
                          (nthpair % static_cast<int64>(tenToEighth) == 0))) {
          fmt::printf("%4.1f%%, t=%5.1f min\n",
                      nthreads * 100.0 * nthpair / npairs, timer.minutes());
        }

        // compute rho
        rho = 0;
        pf1 = pf[k1];
        pf2 = pf[k2];
        for (t = 0; t < wbc->f->nframes; t++) {
          rho += pf1[t] * pf2[t]; // is this right? below is the old code
                                  // rho += (*pf1) * (*pf2);
                                  // pf1++;
                                  // pf2++;
        }
        if (wbc->M != nullptr) {
          wbc->M->rptr[k1 + 1][k2 + 1] = rho;
          wbc->M->rptr[k2 + 1][k1 + 1] = rho;
        }

        val = MRIgetVoxVal(rhomean[thno], k1, 0, 0, 0);
        MRIsetVoxVal(rhomean[thno], k1, 0, 0, 0, val + rho);

        for (auto nthrho : wbc->rholist) {
          if (fabs(rho) < nthrho) {
            continue;
          }
          n = MRIgetVoxVal(conth[thno], k1, 0, 0, nthrho);
          MRIsetVoxVal(conth[thno], k1, 0, 0, nthrho, n + 1);
          n = MRIgetVoxVal(conth[thno], k2, 0, 0, nthrho);
          MRIsetVoxVal(conth[thno], k2, 0, 0, nthrho, n + 1);
        }

        if (!wbc->DoDist) {
          continue;
        }

        ct1 = MRIgetVoxVal(wbc->coordtype, k1, 0, 0, 0);
        ct2 = MRIgetVoxVal(wbc->coordtype, k2, 0, 0, 0);
        if (ct1 == ct2) {
          x1   = MRIgetVoxVal(wbc->xyz, k1, 0, 0, 0);
          y1   = MRIgetVoxVal(wbc->xyz, k1, 0, 0, 1);
          z1   = MRIgetVoxVal(wbc->xyz, k1, 0, 0, 2);
          x2   = MRIgetVoxVal(wbc->xyz, k2, 0, 0, 0);
          y2   = MRIgetVoxVal(wbc->xyz, k2, 0, 0, 1);
          z2   = MRIgetVoxVal(wbc->xyz, k2, 0, 0, 2);
          dx   = x1 - x2;
          dy   = y1 - y2;
          dz   = z1 - z2;
          dist = sqrt(dx * dx + dy * dy + dz * dz);
          if (dist < wbc->distthresh) {
            q = 1; // short dist
          } else {
            q = 2; // long dist
          }
          if (q == 1 && (((wbc->lh2 != nullptr) && ct1 == 1) ||
                         ((wbc->rh2 != nullptr) && ct1 == 2))) {
            /* Surface point that has been declared short using the
            first surface (prob ?h.white). Now check with the second
            surface (prob ?h.inflated). This can only cause a short to
            be recoded to long. The problem this is attempting to
            solve is that it is very expensive to compute the actual
            distance between two surface points. The actual distance is
            always >= the euclidean dist between the ?h.white points, so,
            if that distance is > thresh, it must be a long connection.
            If it is shorter, then it could be a short connection or
            it could be a long with a short euclidean (eg, jumping across
            a gyrus or sulcus) so this checks the inflated. The problem
            with the inflated is that it is not distortion-free.
            */
            x1   = MRIgetVoxVal(wbc->xyz2, k1, 0, 0, 0);
            y1   = MRIgetVoxVal(wbc->xyz2, k1, 0, 0, 1);
            z1   = MRIgetVoxVal(wbc->xyz2, k1, 0, 0, 2);
            x2   = MRIgetVoxVal(wbc->xyz2, k2, 0, 0, 0);
            y2   = MRIgetVoxVal(wbc->xyz2, k2, 0, 0, 1);
            z2   = MRIgetVoxVal(wbc->xyz2, k2, 0, 0, 2);
            dx   = x1 - x2;
            dy   = y1 - y2;
            dz   = z1 - z2;
            dist = sqrt(dx * dx + dy * dy + dz * dz);
            if (dist >= wbc->distthresh) {
              q = 2; // long dist
            }
          }
        } else {
          q = 2; // long dist
        }
        if (q == 1) {
          conDth = conSth[thno];
        }
        if (q == 2) {
          conDth = conLth[thno];
        }

        for (auto nthrho : wbc->rholist) {
          if (fabs(rho) < nthrho) {
            continue;
          }
          n = MRIgetVoxVal(conDth, k1, 0, 0, nthrho);
          MRIsetVoxVal(conDth, k1, 0, 0, nthrho, n + 1);
          n = MRIgetVoxVal(conDth, k2, 0, 0, nthrho);
          MRIsetVoxVal(conDth, k2, 0, 0, nthrho, n + 1);
        }

      } // k2
    }   // k1
    ROMP_PFLB_end
  } // thread
  ROMP_PF_end

      // Sum up the threads
      for (threadno = 0; threadno < nthreads; threadno++) {
    MRIadd(wbc->rhomean, rhomean[threadno], wbc->rhomean);
    MRIadd(wbc->con, conth[threadno], wbc->con);
    if (wbc->DoDist) {
      MRIadd(wbc->conS, conSth[threadno], wbc->conS);
      MRIadd(wbc->conL, conLth[threadno], wbc->conL);
    }
  }

  // Divide number of connections by total possible
  std::cout << "Scaling by ntot-1 " << wbc->ntot - 1 << std::endl;
  MRImultiplyConst(wbc->rhomean, 1.0 / (wbc->ntot - 1), wbc->rhomean);
  MRImultiplyConst(wbc->con, 1.0 / (wbc->ntot - 1), wbc->con);
  if (wbc->DoDist) {
    MRImultiplyConst(wbc->conS, 1.0 / (wbc->ntot - 1), wbc->conS);
    MRImultiplyConst(wbc->conL, 1.0 / (wbc->ntot - 1), wbc->conL);
  }

  // Clean up
  for (threadno = 0; threadno < nthreads; threadno++) {
    MRIfree(&rhomean[threadno]);
    MRIfree(&conth[threadno]);
    if (wbc->DoDist) {
      MRIfree(&conSth[threadno]);
      MRIfree(&conLth[threadno]);
    }
  }

  fmt::printf(" t = %6.4f, ntot = %d\n", timer.seconds(), wbc->ntot);
  fflush(stdout);
  return (wbc->con);
}

int WBCprep(WBC *wbc) {
  int     nthvox;
  int     c;
  int     r;
  int     s;
  int     t;
  int     k;
  int     vtxno;
  MATRIX *V;
  MATRIX *crs;
  MATRIX *xyz = nullptr;
  float   val;
  float   voxvolume;

  if (wbc->fvol != nullptr) {
    wbc->nframes = wbc->fvol->nframes;
  } else if (wbc->flh != nullptr) {
    wbc->nframes = wbc->flh->nframes;
  } else if (wbc->frh != nullptr) {
    wbc->nframes = wbc->frh->nframes;
  }

  if (wbc->lhlabel != nullptr) {
    wbc->lhmask = MRISlabel2Mask(wbc->lh, wbc->lhlabel, nullptr);
  }
  if (wbc->rhlabel != nullptr) {
    wbc->rhmask = MRISlabel2Mask(wbc->rh, wbc->rhlabel, nullptr);
  }

  wbc->nvolmask = 0;
  wbc->nlhmask  = 0;
  wbc->nrhmask  = 0;
  if (wbc->fvol != nullptr) {
    if (wbc->volmask != nullptr) {
      wbc->nvolmask = MRIcountAboveThreshold(wbc->volmask, half);
    } else {
      wbc->nvolmask = wbc->fvol->width * wbc->fvol->height * wbc->fvol->depth;
    }
  }
  if (wbc->flh != nullptr) {
    if (wbc->lhmask != nullptr) {
      wbc->nlhmask = MRIcountAboveThreshold(wbc->lhmask, half);
    } else {
      wbc->nlhmask = wbc->flh->width;
    }
  }
  if (wbc->frh != nullptr) {
    if (wbc->rhmask != nullptr) {
      wbc->nrhmask = MRIcountAboveThreshold(wbc->rhmask, half);
    } else {
      wbc->nrhmask = wbc->frh->width;
    }
  }

  wbc->ntot = wbc->nvolmask + wbc->nlhmask + wbc->nrhmask;
  fmt::printf("nmask %d %d %d   %d\n", wbc->nvolmask, wbc->nlhmask,
              wbc->nrhmask, wbc->ntot);

  wbc->coordtype = MRIalloc(wbc->ntot, 1, 1, MRI_INT);
  wbc->vvol =
      MRIallocSequence(wbc->ntot, 1, 1, MRI_FLOAT, 1); // vertex/voxel volume
  wbc->vertexno = MRIalloc(wbc->ntot, 1, 1, MRI_INT);
  wbc->xyz      = MRIallocSequence(wbc->ntot, 1, 1, MRI_FLOAT, 3);
  if ((wbc->lh2 != nullptr) || (wbc->rh2 != nullptr)) {
    wbc->xyz2 = MRIallocSequence(wbc->ntot, 1, 1, MRI_FLOAT, 3);
  }
  wbc->f = MRIallocSequence(wbc->ntot, 1, 1, MRI_FLOAT, wbc->nframes);

  nthvox = 0;
  if (wbc->fvol != nullptr) {
    // Not sure it is necessary/useful to map into a particular mm space as
    // long as it is a mm space so the distance threshold is right. This tkreg
    // space does not match that of the surface, but this should be ok becase
    // we never compute the distance from volume to surface, it is always
    // long distance. If this turns out not to be the case, then will need
    // a registraion matrix
    voxvolume       = wbc->fvol->xsize * wbc->fvol->ysize * wbc->fvol->zsize;
    V               = MRIxfmCRS2XYZtkreg(wbc->fvol);
    crs             = MatrixAlloc(4, 1, MATRIX_REAL);
    crs->rptr[4][1] = 1;
    for (s = 0; s < wbc->fvol->depth; s++) {
      for (c = 0; c < wbc->fvol->width; c++) {
        for (r = 0; r < wbc->fvol->height; r++) {
          if ((wbc->volmask != nullptr) &&
              MRIgetVoxVal(wbc->volmask, c, r, s, 0) < half) {
            continue;
          }
          MRIsetVoxVal(wbc->coordtype, nthvox, 0, 0, 0, 0);
          MRIsetVoxVal(wbc->vvol, nthvox, 0, 0, 0, voxvolume);
          crs->rptr[1][1] = static_cast<float>(c);
          crs->rptr[2][1] = static_cast<float>(r);
          crs->rptr[3][1] = static_cast<float>(s);
          xyz             = MatrixMultiplyD(V, crs, xyz);
          for (k = 0; k < 3; k++) {
            MRIsetVoxVal(wbc->xyz, nthvox, 0, 0, k, xyz->rptr[k + 1][1]);
          }
          for (t = 0; t < wbc->nframes; t++) {
            val = MRIgetVoxVal(wbc->fvol, c, r, s, t);
            MRIsetVoxVal(wbc->f, nthvox, 0, 0, t, val);
          } // time
          nthvox++;
        } // row
      }   // col
    }     // slice
    MatrixFree(&V);
    MatrixFree(&xyz);
    MatrixFree(&crs);
  }
  if (wbc->flh != nullptr) {
    for (vtxno = 0; vtxno < wbc->flh->width; vtxno++) {
      if ((wbc->lhmask != nullptr) &&
          MRIgetVoxVal(wbc->lhmask, vtxno, 0, 0, 0) < half) {
        continue;
      }
      MRIsetVoxVal(wbc->vertexno, nthvox, 0, 0, 0, static_cast<float>(vtxno));
      MRIsetVoxVal(wbc->coordtype, nthvox, 0, 0, 0, 1);
      MRIsetVoxVal(wbc->xyz, nthvox, 0, 0, 0, wbc->lh->vertices[vtxno].x);
      MRIsetVoxVal(wbc->xyz, nthvox, 0, 0, 1, wbc->lh->vertices[vtxno].y);
      MRIsetVoxVal(wbc->xyz, nthvox, 0, 0, 2, wbc->lh->vertices[vtxno].z);
      if (wbc->lh2 != nullptr) {
        MRIsetVoxVal(wbc->xyz2, nthvox, 0, 0, 0, wbc->lh->vertices[vtxno].x);
        MRIsetVoxVal(wbc->xyz2, nthvox, 0, 0, 1, wbc->lh->vertices[vtxno].y);
        MRIsetVoxVal(wbc->xyz2, nthvox, 0, 0, 2, wbc->lh->vertices[vtxno].z);
      }
      for (t = 0; t < wbc->flh->nframes; t++) {
        val = MRIgetVoxVal(wbc->flh, vtxno, 0, 0, t);
        MRIsetVoxVal(wbc->f, nthvox, 0, 0, t, val);
      } // time
      nthvox++;
    } // vertex
  }
  if (wbc->frh != nullptr) {
    for (vtxno = 0; vtxno < wbc->frh->width; vtxno++) {
      if ((wbc->rhmask != nullptr) &&
          MRIgetVoxVal(wbc->rhmask, vtxno, 0, 0, 0) < half) {
        continue;
      }
      MRIsetVoxVal(wbc->vertexno, nthvox, 0, 0, 0, static_cast<float>(vtxno));
      MRIsetVoxVal(wbc->coordtype, nthvox, 0, 0, 0, 2);
      MRIsetVoxVal(wbc->xyz, nthvox, 0, 0, 0, wbc->rh->vertices[vtxno].x);
      MRIsetVoxVal(wbc->xyz, nthvox, 0, 0, 1, wbc->rh->vertices[vtxno].y);
      MRIsetVoxVal(wbc->xyz, nthvox, 0, 0, 2, wbc->rh->vertices[vtxno].z);
      if (wbc->rh2 != nullptr) {
        MRIsetVoxVal(wbc->xyz2, nthvox, 0, 0, 0, wbc->rh->vertices[vtxno].x);
        MRIsetVoxVal(wbc->xyz2, nthvox, 0, 0, 1, wbc->rh->vertices[vtxno].y);
        MRIsetVoxVal(wbc->xyz2, nthvox, 0, 0, 2, wbc->rh->vertices[vtxno].z);
      }
      for (t = 0; t < wbc->frh->nframes; t++) {
        val = MRIgetVoxVal(wbc->frh, vtxno, 0, 0, t);
        MRIsetVoxVal(wbc->f, nthvox, 0, 0, t, val);
      } // time
      nthvox++;
    } // vertex
  }

  // normalize by stddev
  wbc->fnorm = MRIframeNorm(wbc->f, nullptr, nullptr);

  wbc->rhomean = MRIallocSequence(wbc->ntot, 1, 1, MRI_FLOAT, 1);
  wbc->con = MRIallocSequence(wbc->ntot, 1, 1, MRI_FLOAT, wbc->rholist.size());
  if (wbc->DoDist) {
    wbc->conS =
        MRIallocSequence(wbc->ntot, 1, 1, MRI_FLOAT, wbc->rholist.size());
    wbc->conL =
        MRIallocSequence(wbc->ntot, 1, 1, MRI_FLOAT, wbc->rholist.size());
  }

  return (0);
}

int WBCfinish(WBC *wbc) {
  int   nthvox;
  int   c;
  int   r;
  int   s;
  int   vtxno;
  int   nthrho;
  float val;

  nthvox = 0;
  if (wbc->fvol != nullptr) {
    wbc->volrhomean = MRIallocSequence(wbc->fvol->width, wbc->fvol->height,
                                       wbc->fvol->depth, MRI_FLOAT, 1);
    MRIcopyHeader(wbc->fvol, wbc->volrhomean);
    MRIcopyPulseParameters(wbc->fvol, wbc->volrhomean);
    wbc->volcon =
        MRIallocSequence(wbc->fvol->width, wbc->fvol->height, wbc->fvol->depth,
                         MRI_FLOAT, wbc->rholist.size());
    MRIcopyHeader(wbc->fvol, wbc->volcon);
    MRIcopyPulseParameters(wbc->fvol, wbc->volcon);
    if (wbc->DoDist) {
      wbc->volconS =
          MRIallocSequence(wbc->fvol->width, wbc->fvol->height,
                           wbc->fvol->depth, MRI_FLOAT, wbc->rholist.size());
      MRIcopyHeader(wbc->fvol, wbc->volconS);
      MRIcopyPulseParameters(wbc->fvol, wbc->volconS);
      wbc->volconL =
          MRIallocSequence(wbc->fvol->width, wbc->fvol->height,
                           wbc->fvol->depth, MRI_FLOAT, wbc->rholist.size());
      MRIcopyHeader(wbc->fvol, wbc->volconL);
      MRIcopyPulseParameters(wbc->fvol, wbc->volconL);
    }
    for (s = 0; s < wbc->fvol->depth; s++) {
      for (c = 0; c < wbc->fvol->width; c++) {
        for (r = 0; r < wbc->fvol->height; r++) {
          if ((wbc->volmask != nullptr) &&
              MRIgetVoxVal(wbc->volmask, c, r, s, 0) < half) {
            continue;
          }
          val = MRIgetVoxVal(wbc->rhomean, nthvox, 0, 0, 0);
          MRIsetVoxVal(wbc->volrhomean, c, r, s, 0, val);
          for (auto nthrho : wbc->rholist) {
            val = MRIgetVoxVal(wbc->con, nthvox, 0, 0, nthrho);
            MRIsetVoxVal(wbc->volcon, c, r, s, nthrho, val);
            if (wbc->DoDist) {
              val = MRIgetVoxVal(wbc->conS, nthvox, 0, 0, nthrho);
              MRIsetVoxVal(wbc->volconS, c, r, s, nthrho, val);
              val = MRIgetVoxVal(wbc->conL, nthvox, 0, 0, nthrho);
              MRIsetVoxVal(wbc->volconL, c, r, s, nthrho, val);
            }
          } // rho
          nthvox++;
        } // row
      }   // col
    }     // slice
  }
  if (wbc->flh != nullptr) {
    wbc->lhrhomean = MRIallocSequence(wbc->flh->width, wbc->flh->height,
                                      wbc->flh->depth, MRI_FLOAT, 1);
    MRIcopyHeader(wbc->flh, wbc->lhrhomean);
    MRIcopyPulseParameters(wbc->flh, wbc->lhrhomean);
    wbc->lhcon =
        MRIallocSequence(wbc->flh->width, wbc->flh->height, wbc->flh->depth,
                         MRI_FLOAT, wbc->rholist.size());
    MRIcopyHeader(wbc->flh, wbc->lhcon);
    MRIcopyPulseParameters(wbc->flh, wbc->lhcon);
    if (wbc->DoDist) {
      wbc->lhconS =
          MRIallocSequence(wbc->flh->width, wbc->flh->height, wbc->flh->depth,
                           MRI_FLOAT, wbc->rholist.size());
      MRIcopyHeader(wbc->flh, wbc->lhconS);
      MRIcopyPulseParameters(wbc->flh, wbc->lhconS);
      wbc->lhconL =
          MRIallocSequence(wbc->flh->width, wbc->flh->height, wbc->flh->depth,
                           MRI_FLOAT, wbc->rholist.size());
      MRIcopyHeader(wbc->flh, wbc->lhconL);
      MRIcopyPulseParameters(wbc->flh, wbc->lhconL);
    }
    for (vtxno = 0; vtxno < wbc->flh->width; vtxno++) {
      if ((wbc->lhmask != nullptr) &&
          MRIgetVoxVal(wbc->lhmask, vtxno, 0, 0, 0) < half) {
        continue;
      }
      val = MRIgetVoxVal(wbc->rhomean, nthvox, 0, 0, 0);
      MRIsetVoxVal(wbc->lhrhomean, vtxno, 0, 0, 0, val);
      for (auto nthrho : wbc->rholist) {
        val = MRIgetVoxVal(wbc->con, nthvox, 0, 0, nthrho);
        MRIsetVoxVal(wbc->lhcon, vtxno, 0, 0, nthrho, val);
        if (wbc->DoDist) {
          val = MRIgetVoxVal(wbc->conS, nthvox, 0, 0, nthrho);
          MRIsetVoxVal(wbc->lhconS, vtxno, 0, 0, nthrho, val);
          val = MRIgetVoxVal(wbc->conL, nthvox, 0, 0, nthrho);
          MRIsetVoxVal(wbc->lhconL, vtxno, 0, 0, nthrho, val);
        }
      } // rho
      nthvox++;
    } // vertex
  }
  if (wbc->frh != nullptr) {
    wbc->rhrhomean = MRIallocSequence(wbc->frh->width, wbc->frh->height,
                                      wbc->frh->depth, MRI_FLOAT, 1);
    MRIcopyHeader(wbc->frh, wbc->rhrhomean);
    MRIcopyPulseParameters(wbc->frh, wbc->rhrhomean);
    wbc->rhcon =
        MRIallocSequence(wbc->frh->width, wbc->frh->height, wbc->frh->depth,
                         MRI_FLOAT, wbc->rholist.size());
    MRIcopyHeader(wbc->frh, wbc->rhcon);
    MRIcopyPulseParameters(wbc->frh, wbc->rhcon);
    if (wbc->DoDist) {
      wbc->rhconS =
          MRIallocSequence(wbc->frh->width, wbc->frh->height, wbc->frh->depth,
                           MRI_FLOAT, wbc->rholist.size());
      MRIcopyHeader(wbc->frh, wbc->rhconS);
      MRIcopyPulseParameters(wbc->frh, wbc->rhconS);
      wbc->rhconL =
          MRIallocSequence(wbc->frh->width, wbc->frh->height, wbc->frh->depth,
                           MRI_FLOAT, wbc->rholist.size());
      MRIcopyHeader(wbc->frh, wbc->rhconL);
      MRIcopyPulseParameters(wbc->frh, wbc->rhconL);
    }
    for (vtxno = 0; vtxno < wbc->frh->width; vtxno++) {
      if ((wbc->rhmask != nullptr) &&
          MRIgetVoxVal(wbc->rhmask, vtxno, 0, 0, 0) < half) {
        continue;
      }
      val = MRIgetVoxVal(wbc->rhomean, nthvox, 0, 0, 0);
      MRIsetVoxVal(wbc->rhrhomean, vtxno, 0, 0, 0, val);
      for (auto nthrho : wbc->rholist) {
        val = MRIgetVoxVal(wbc->con, nthvox, 0, 0, nthrho);
        MRIsetVoxVal(wbc->rhcon, vtxno, 0, 0, nthrho, val);
        if (wbc->DoDist) {
          val = MRIgetVoxVal(wbc->conS, nthvox, 0, 0, nthrho);
          MRIsetVoxVal(wbc->rhconS, vtxno, 0, 0, nthrho, val);
          val = MRIgetVoxVal(wbc->conL, nthvox, 0, 0, nthrho);
          MRIsetVoxVal(wbc->rhconL, vtxno, 0, 0, nthrho, val);
        }
      } // rho
      nthvox++;
    } // vertex
  }
  return (0);
}

int Index2UpperSubscript(int N, int64 i, int64 *r, int64 *c) {
  int64 i1;
  int64 r1;
  int64 c1;
  int64 ir1;
  i1 = i + 1;
  r1 = ceil(
      (-(1 - two * static_cast<float>(N)) -
       sqrt(pow((1 - two * static_cast<float>(N)), 2) - eightPointZero * i1)) /
      two);
  ir1 = N * r1 - (r1 * (r1 + 1)) / 2;
  c1  = N - (ir1 - i1);
  *r  = r1 - 1;
  *c  = c1 - 1;
  return (0);
}

// TODO(aboualiaa): move to tests file
int Index2UpperSubscriptTest(int N) {
  int   r;
  int   c;
  int64 rt;
  int64 ct;
  int64 i;
  int64 err;

  fmt::printf("Index2UpperSubscriptTest(): N = %d\n", N);
  err = 0;
  i   = 0;
  for (r = 0; r < N - 1; r++) {
    for (c = r + 1; c < N; c++) {
      Index2UpperSubscript(N, i, &rt, &ct);
      if (r != rt || c != ct) {
        err++;
        fmt::printf("ERROR: %4d %4d %4ld  %4d %4d\n", r, c, i, rt, ct);
      }
      i = i + 1;
    }
  }
  fmt::printf("Found %ld errors\n", err);
  if (err != 0) {
    return (1);
  }

  return (0);
}

int WBCnframes(WBC *wbc) {
  if (wbc->fvol != nullptr) {
    wbc->nframes = wbc->fvol->nframes;
    if (wbc->flh != nullptr) {
      if (wbc->fvol->nframes != wbc->flh->nframes) {
        fmt::printf("ERROR: nframes mismatch fvol and flh %d %d\n",
                    wbc->fvol->nframes, wbc->flh->nframes);
        return (1);
      }
    }
    if (wbc->frh != nullptr) {
      if (wbc->fvol->nframes != wbc->frh->nframes) {
        fmt::printf("ERROR: nframes mismatch fvol and frh %d %d\n",
                    wbc->fvol->nframes, wbc->frh->nframes);
        return (1);
      }
    }
  }
  if ((wbc->flh != nullptr) && (wbc->frh != nullptr)) {
    wbc->nframes = wbc->flh->nframes;
    if (wbc->flh->nframes != wbc->frh->nframes) {
      fmt::printf("ERROR: nframes mismatch flh and frh %d %d\n",
                  wbc->flh->nframes, wbc->frh->nframes);
      return (1);
    }
  }
  std::cout << "nframes = " << wbc->nframes << std::endl;
  return (0);
}

// TODO(aboualiaa): move to tests file
WBC *WBCtestSynth(WBC *wbc) {
  FSENV *             fsenv;
  double              x0;
  double              y0;
  double              z0;
  double              x;
  double              y;
  double              z;
  double              dx;
  double              dy;
  double              dz;
  double              d;
  float               volres;
  std::vector<double> wf;
  double              dmax = 0;
  int                 nshort;
  int                 nlong;
  int                 t;
  int                 vno;
  int                 hemi;
  int                 voldim;
  int                 nframes;
  int                 c0;
  int                 r0;
  int                 s0;
  int                 c2;
  std::string         tmpstr;
  char *              hemistr;
  MRIS *              surf;
  MRIS *              surf2;
  MRI *               func;
  MRI *               mask;

  fsenv        = FSENVgetenv();
  volres       = wbc->wbcsynth.volres;
  voldim       = wbc->wbcsynth.voldim;
  nframes      = wbc->wbcsynth.nframes;
  wbc->nframes = nframes;

  if (wbc->wbcsynth.ForceFail != 0) {
    wbc->wbcsynth.distthresh = 2 * wbc->distthresh;
  } else {
    wbc->wbcsynth.distthresh = wbc->distthresh;
  }

  // if(fabs(wbc->distthresh - wbc->wbcsynth->distthresh) < .00001) dmax =
  // 10e10;

  spdlog::debug("create synthetic time series");

  // create synthetic time seriesa
  wf.reserve(wbc->nframes);
  for (t = 0; t < wbc->nframes; t++) {
    wf[t] = drand48() - half;
  }
  wbc->wbcsynth.wf = wf;

  spdlog::debug("mri alloc sequence");

  wbc->fvol = MRIallocSequence(voldim, voldim, voldim, MRI_FLOAT, wbc->nframes);
  wbc->fvol->xsize = volres;
  wbc->fvol->ysize = volres;
  wbc->fvol->zsize = volres;
  wbc->fvol->tr    = max_rows;

  wbc->volmask        = MRIallocSequence(voldim, voldim, voldim, MRI_FLOAT, 1);
  wbc->volmask->xsize = volres;
  wbc->volmask->ysize = volres;
  wbc->volmask->zsize = volres;
  wbc->volmask->tr    = max_rows;

  c0 = nint(static_cast<float>(wbc->fvol->width) / two);
  r0 = nint(static_cast<float>(wbc->fvol->height) / two);
  s0 = nint(static_cast<float>(wbc->fvol->depth) / two);
  // Column for "long" connections
  c2 = c0 + static_cast<int>(ceil(wbc->wbcsynth.distthresh / volres)) + 1;
  if (wbc->wbcsynth.ForceFail != 0) {
    c2 = c0 + static_cast<int>(floor(wbc->distthresh / volres)) - 1;
  }
  fmt::printf("synth vol: %d %d %d, c2=%d, res=%lf, dim=%d\n", c0, r0, s0, c2,
              volres, voldim);
  wbc->wbcsynth.c0 = c0;
  wbc->wbcsynth.r0 = r0;
  wbc->wbcsynth.s0 = s0;
  wbc->wbcsynth.c2 = c2;

  spdlog::debug("loop of setvoxvals");

  for (t = 0; t < wbc->nframes; t++) {
    // Set crs0 to have 2 short connections
    MRIsetVoxVal(wbc->fvol, c0, r0, s0, t, wf[t]);
    MRIsetVoxVal(wbc->fvol, c0, r0 + 1, s0, t, wf[t]);
    MRIsetVoxVal(wbc->fvol, c0, r0 - 1, s0, t, wf[t]);
    // and 2 longs plus the longs from the surfaces
    MRIsetVoxVal(wbc->fvol, c2, r0, s0, t, wf[t]);
    MRIsetVoxVal(wbc->fvol, c2, r0 + 1, s0, t, wf[t]);
    if (wbc->wbcsynth.ForceFail != 0) {
      MRIsetVoxVal(wbc->fvol, c2, r0 + 1, s0 + 1, t, wf[t]);
    }
  }
  MRIsetVoxVal(wbc->volmask, c0, r0, s0, 0, 1);
  MRIsetVoxVal(wbc->volmask, c0, r0 + 1, s0, 0, 1);
  MRIsetVoxVal(wbc->volmask, c0, r0 - 1, s0, 0, 1);
  MRIsetVoxVal(wbc->volmask, c2, r0, s0, 0, 1);
  MRIsetVoxVal(wbc->volmask, c2, r0 + 1, s0, 0, 1);
  MRIsetVoxVal(wbc->volmask, c2, r0 + 1, s0 + 1, 0,
               1); // mask but no signal, unless ForceFail
  MRIsetVoxVal(wbc->volmask, c2, r0 + 1, s0 - 1, 0, 1); // mask but no signal
  wbc->wbcsynth.nshortvol = 2; // number of short con to/from crs0
  wbc->wbcsynth.nlongvol  = 2; // number of long con to/from crs0

  wbc->wbcsynth.v0 = 0;
  for (hemi = 0; hemi < 2; hemi++) {
    if (hemi == 0) {
      hemistr = "lh";
    } else {
      hemistr = "rh";
    }

    spdlog::debug("reading volumes");

    // FIXME: does this work as intended
    // no it doesn't, crashes on dev build cebause freesurfer_home is null
    std::string root;
    if ((fsenv != nullptr) && (fsenv->FREESURFER_HOME != nullptr)) {
      root = fsenv->FREESURFER_HOME;
    } else {
      root = "/Users/aboualiaa";
    }
    tmpstr.append(root)
        .append("/subjects/fsaverage/surf/")
        .append(hemistr)
        .append(".white");
    spdlog::info("reading: {}", tmpstr);
    surf = MRISread(tmpstr.c_str());
    if (surf == nullptr) {
      spdlog::get("stderr")->critical("failed to read: {}", tmpstr);
      exit(1);
    }
    tmpstr.replace(tmpstr.find(".white"), sizeof(".white") - 1, ".inflated");
    spdlog::info("reading volume: {}", tmpstr);
    surf2 = MRISread(tmpstr.c_str());
    if (surf2 == nullptr) {
      spdlog::get("stderr")->critical("failed to read: {}", tmpstr);
      exit(1);
    }
    // fmt::sprintf(tmpstr,"%s/subjects/fsaverage/label/%s.aparc.annot",fsenv->FREESURFER_HOME,hemistr);
    // err = MRISreadAnnotation(surf, tmpstr);

    spdlog::debug("alloc sequence again");

    func     = MRIallocSequence(surf->nvertices, 1, 1, MRI_FLOAT, wbc->nframes);
    func->tr = max_rows;
    mask     = MRIallocSequence(surf->nvertices, 1, 1, MRI_FLOAT, 1);

    spdlog::debug("set the first vertex waveform");

    // set the first vertex waveform
    MRIsetVoxVal(mask, 0, 0, 0, 0, 1);
    for (t = 0; t < wbc->nframes; t++) {
      MRIsetVoxVal(func, 0, 0, 0, t, wf[t]);
    }

    x0     = surf->vertices[0].x;
    y0     = surf->vertices[0].y;
    z0     = surf->vertices[0].z;
    nshort = 0;
    nlong  = 0;
    for (vno = 0; vno < surf->nvertices; vno++) {
      x  = surf->vertices[vno].x;
      y  = surf->vertices[vno].y;
      z  = surf->vertices[vno].z;
      dx = (x - x0);
      dy = (y - y0);
      dz = (z - z0);
      d  = sqrt(dx * dx + dy * dy + dz * dz);
      if (d < zeroPointNine * wbc->wbcsynth.distthresh &&
          nshort < gsl::at(wbc->wbcsynth.nshorttarg, hemi)) {
        for (t = 0; t < wbc->nframes; t++) {
          MRIsetVoxVal(func, vno, 0, 0, t, wf[t]);
        }
        MRIsetVoxVal(mask, vno, 0, 0, 0, 1);
        // fmt::printf("%s short %d %d\n",hemistr,nshort,vno);
        nshort++;
      }
      if (d > onePointOne * wbc->wbcsynth.distthresh &&
          nlong < gsl::at(wbc->wbcsynth.nlongtarg, hemi)) {
        for (t = 0; t < wbc->nframes; t++) {
          MRIsetVoxVal(func, vno, 0, 0, t, wf[t]);
        }
        MRIsetVoxVal(mask, vno, 0, 0, 0, 1);
        // fmt::printf("%s long %d %d\n",hemistr,nlong,vno);
        nlong++;
      }
    }

    spdlog::debug("copying");

    gsl::at(wbc->wbcsynth.nshort, hemi) = nshort - 1; // remove self
    gsl::at(wbc->wbcsynth.nlong, hemi)  = nlong;
    fmt::printf("%s short %d long %d\n", hemistr, nshort, nlong);
    if (hemi == 0) {
      wbc->lh     = surf;
      wbc->lh2    = surf2;
      wbc->flh    = func;
      wbc->lhmask = mask;
    } else {
      wbc->rh     = surf;
      wbc->rh2    = surf2;
      wbc->frh    = func;
      wbc->rhmask = mask;
    }
  } // hemi

  return (wbc);
}

// TODO(aboualiaa): move to tests file
int WBCtestCheck(WBC *wbc) {
  int n;
  int nexp;
  int v0;
  int c0;
  int r0;
  int s0;
  int err;

  v0 = wbc->wbcsynth.v0;
  c0 = wbc->wbcsynth.c0;
  r0 = wbc->wbcsynth.r0;
  s0 = wbc->wbcsynth.s0;

  fmt::printf("WBCtestCheck(): ------------\n");
  fmt::printf("Vol short=%d, long=%d\n", wbc->wbcsynth.nshortvol,
              wbc->wbcsynth.nlongvol);
  fmt::printf("LH short=%d, long=%d\n", wbc->wbcsynth.nshort[0],
              wbc->wbcsynth.nlong[0]);
  fmt::printf("RH short=%d, long=%d\n", wbc->wbcsynth.nshort[1],
              wbc->wbcsynth.nlong[1]);
  fmt::printf("DistThresh: wbc = %lf, synth = %lf\n", wbc->distthresh,
              wbc->wbcsynth.distthresh);
  fmt::printf("ForceFail: %d\n", wbc->wbcsynth.ForceFail);

  err = 0;

  // expected number of cons from lhv0, rhv0, or crs0 to/from everyone else
  nexp = wbc->wbcsynth.nshortvol + wbc->wbcsynth.nlongvol +
         wbc->wbcsynth.nshort[0] + wbc->wbcsynth.nlong[0] +
         wbc->wbcsynth.nshort[1] + wbc->wbcsynth.nlong[1];
  // The above numbers exclude self from the shorts. There are 3 "selves"
  // so add 3, but then subtract one to account for the self for which
  // the comutation is being done, so +2.
  nexp += 2;
  fmt::printf("Expecting %d total connections\n", nexp);
  n = nint(MRIgetVoxVal(wbc->volcon, c0, r0, s0, 0) *
           static_cast<float>(wbc->ntot - 1));
  fmt::printf("  Vol %d ", n);
  if (n != nexp) {
    fmt::printf(" ERROR\n");
    err++;
  } else {
    fmt::printf(" PASS\n");
  }
  n = nint(MRIgetVoxVal(wbc->lhcon, v0, 0, 0, 0) *
           static_cast<float>(wbc->ntot - 1));
  fmt::printf("  lh %d ", n);
  if (n != nexp) {
    fmt::printf(" ERROR\n");
    err++;
  } else {
    fmt::printf(" PASS\n");
  }
  n = nint(MRIgetVoxVal(wbc->rhcon, v0, 0, 0, 0) *
           static_cast<float>(wbc->ntot - 1));
  fmt::printf("  rh %d ", n);
  if (n != nexp) {
    fmt::printf(" ERROR\n");
    err++;
  } else {
    fmt::printf(" PASS\n");
  }

  nexp = wbc->wbcsynth.nshortvol;
  fmt::printf("Volume short: expecting %3d  ", nexp);
  n = nint(MRIgetVoxVal(wbc->volconS, c0, r0, s0, 0) *
           static_cast<float>(wbc->ntot - 1));
  fmt::printf("found %4d ", n);
  if (n != nexp) {
    fmt::printf(" ERROR\n");
    err++;
  } else {
    fmt::printf(" PASS\n");
  }

  nexp = wbc->wbcsynth.nlongvol + wbc->wbcsynth.nshort[0] +
         wbc->wbcsynth.nlong[0] + wbc->wbcsynth.nshort[1] +
         wbc->wbcsynth.nlong[1];
  nexp += 2;
  fmt::printf("Volume long:  expecting %3d  ", nexp);
  n = nint(MRIgetVoxVal(wbc->volconL, c0, r0, s0, 0) *
           static_cast<float>(wbc->ntot - 1));
  fmt::printf("found %4d ", n);
  if (n != nexp) {
    fmt::printf(" ERROR\n");
    err++;
  } else {
    fmt::printf(" PASS\n");
  }

  nexp = wbc->wbcsynth.nshort[0];
  fmt::printf("lh short:       expecting %3d  ", nexp);
  n = nint(MRIgetVoxVal(wbc->lhconS, v0, 0, 0, 0) *
           static_cast<float>(wbc->ntot - 1));
  fmt::printf("found %4d ", n);
  if (n != nexp) {
    fmt::printf(" ERROR\n");
    err++;
  } else {
    fmt::printf(" PASS\n");
  }

  nexp = wbc->wbcsynth.nshortvol + wbc->wbcsynth.nlongvol +
         wbc->wbcsynth.nlong[0] + wbc->wbcsynth.nshort[1] +
         wbc->wbcsynth.nlong[1];
  nexp += 2;
  fmt::printf("lh long:        expecting %3d  ", nexp);
  n = nint(MRIgetVoxVal(wbc->lhconL, v0, 0, 0, 0) *
           static_cast<float>(wbc->ntot - 1));
  fmt::printf("found %4d ", n);
  if (n != nexp) {
    fmt::printf(" ERROR\n");
    err++;
  } else {
    fmt::printf(" PASS\n");
  }

  nexp = wbc->wbcsynth.nshort[1];
  fmt::printf("rh short:       expecting %3d  ", nexp);
  n = nint(MRIgetVoxVal(wbc->rhconS, v0, 0, 0, 0) *
           static_cast<float>(wbc->ntot - 1));
  fmt::printf("found %4d ", n);
  if (n != nexp) {
    fmt::printf(" ERROR\n");
    err++;
  } else {
    fmt::printf(" PASS\n");
  }

  nexp = wbc->wbcsynth.nshortvol + wbc->wbcsynth.nlongvol +
         wbc->wbcsynth.nshort[0] + wbc->wbcsynth.nlong[0] +
         wbc->wbcsynth.nlong[1];
  nexp += 2;
  fmt::printf("rh long:        expecting %3d  ", nexp);
  n = nint(MRIgetVoxVal(wbc->rhconL, v0, 0, 0, 0) *
           static_cast<float>(wbc->ntot - 1));
  fmt::printf("found %4d ", n);
  if (n != nexp) {
    fmt::printf(" ERROR\n");
    err++;
  } else {
    fmt::printf(" PASS\n");
  }

  if (err != 0) {
    fmt::printf("Overall FAILED with %d errors\n", err);
    if (wbc->wbcsynth.ForceFail != 0) {
      fmt::printf(" ... but, a failure was forced, so this is a good thing\n");
      if (err != expectedErrs) {
        fmt::printf(" ... however, expected 9 errors but got %d\n", err);
        return (err);
      }
      return (0);
    }
    return (err);
  }

  fmt::printf("Overall PASS\n");
  return (0);
}
} // namespace wbc
