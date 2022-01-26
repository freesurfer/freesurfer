#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>

#include "mri.h"
#include "mri2.h"


#include "mri_segcentroids.help.xml.h"
static void printHelp(int exit_val) {
  outputHelpXml(mri_segcentroids_help_xml, mri_segcentroids_help_xml_len);
  exit(exit_val);
}


// structure to represent the centroid of each label
struct Centroid {
  int id;
  std::string labelname;
  float mass, x, y, z;
};


// command line inpute parser
class InputParser {
public:
  std::string segfile, weightsfile, ltafile, outfile, pointset,ctabfile;
  bool include_zero;
  int precision;

  InputParser(int &argc, char **argv) {

    // check for any input
    if (argc == 1) printHelp(1);

    include_zero = false;
    int i = 1;
    std::string opt;
    while (i < argc) {
      std::string opt(argv[i]);

      // -----------------------
      //           --i
      // -----------------------
      if (opt  == "--i") {
        i++;
        if ((i >= argc) || (ISOPTION(*argv[i]))) {
          std::cerr << "ERROR: must specify valid path to input segmentation after '--i'\n";
          exit(1);
        }
        segfile = argv[i];
      }
      // -----------------------
      //           --o
      // -----------------------
      else if (opt == "--o") {
        i++;
        if ((i >= argc) || (ISOPTION(*argv[i]))) {
          std::cerr << "ERROR: must specify output filename with '--o'\n";
          exit(1);
        }
        outfile = argv[i];
      }
      // -----------------------
      //           --p
      // -----------------------
      else if (opt == "--p") {
        i++;
        if ((i >= argc) || (ISOPTION(*argv[i]))) {
          std::cerr << "ERROR: must specify pointset filename with '--p'\n";
          exit(1);
        }
        pointset = argv[i];
      }
      // -----------------------
      //       --weights
      // -----------------------
      else if (opt == "--weights") {
        i++;
        if ((i >= argc) || (ISOPTION(*argv[i]))) {
          std::cerr << "ERROR: must specify path to weight volume after '--weights'\n";
          exit(1);
        }
        weightsfile  = argv[i];
      }
      // -----------------------
      //         --reg
      // -----------------------
      else if (opt == "--reg") {
        i++;
        if ((i >= argc) || (ISOPTION(*argv[i]))) {
          std::cerr << "ERROR: must specify valid path to lta after '--reg'\n";
          exit(1);
        }
        ltafile = argv[i];
      }
      // -----------------------
      //         --ctab
      // -----------------------
      else if (opt == "--ctab") {
        i++;
        if ((i >= argc) || (ISOPTION(*argv[i]))) {
          std::cerr << "ERROR: must specify path to color lookup table after '--ctab'\n";
          exit(1);
        }
        ctabfile = argv[i];
      }
      // -----------------------
      //     --ctab-default
      // -----------------------
      else if (opt == "--ctab-default") {
        std::string fs_home = std::getenv("FREESURFER_HOME");
        if (fs_home.empty()) {
          std::cerr << "WARNING: FREESURFER_HOME is not set"
                    << " - cannot find default lookup table\n";
        }
        ctabfile = fs_home + "/FreeSurferColorLUT.txt";
      }
      // -----------------------
      //     --include-zero
      // -----------------------
      else if (opt == "--include-zero") {
        include_zero = true;
      }
      // -----------------------
      //          help
      // -----------------------
      else if ((opt == "-h") || (opt == "--help")) {
        printHelp(0);
      }
      // -----------------------
      //         unknown
      // -----------------------
      else {
        std::cerr << "ERROR: unknown option " << argv[i] << std::endl;
        exit(1);
      }

      i++;
    }
  }
};


int main(int argc, char **argv) {

  InputParser input(argc, argv);

  // load segmentation
  if (input.segfile.empty()) {
    std::cerr << "ERROR: must specify path to input segmentation with '--i'\n";
    exit(1);
  }
  MRI *seg = MRIread(input.segfile.c_str());
  if (!seg) {
    std::cerr << "ERROR: loading segmentation " << input.segfile << std::endl;
    exit(1);
  }

  // check for output file name
  if (input.outfile.empty()) {
    std::cerr << "ERROR: must specify output filename with '--o'\n";
    exit(1);
  }

  // load weights volume
  MRI *weights = NULL;
  if (!input.weightsfile.empty()) {
    std::cerr << "using weights from  " << input.weightsfile << std::endl;
    weights = MRIread(input.weightsfile.c_str());
    if (!weights) {
      std::cerr << "ERROR: loading volume " << input.weightsfile << std::endl;
      exit(1);
    }
  }

  // use lta
  LTA *lta = NULL;
  if (!input.ltafile.empty()) {
    std::cerr << "using linear transform: " << input.ltafile << std::endl;
    lta = LTAreadEx(input.ltafile.c_str());
    if (lta->type != LINEAR_RAS_TO_RAS) lta = LTAchangeType(lta, LINEAR_RAS_TO_RAS);
  }

  // load colortable
  COLOR_TABLE *ctab = NULL;
  if (!input.ctabfile.empty()) {
    ctab = CTABreadASCII(input.ctabfile.c_str());
    if (!ctab) {
      std::cerr << "ERROR: reading " << input.ctabfile << std::endl;
      exit(1);
    }
  }


  // -------------------- compute centroids --------------------


  std::map<int, Centroid> centroids;
  int numids, label_chars, id_chars, valid_id;
  char char_name[500];
  double x, y, z, wx, wy, wz, weight;
  float fx, fy, fz;
  int max_id_chars = 2, max_label_chars = 10;  // used for table formatting

  int *idlist = MRIsegIdList(seg, &numids, 0);
  std::vector<int> ids(idlist, idlist + numids);

  for (int i = 0 ; i < numids; i++) {
    // create centroid for each label
    Centroid centroid = Centroid();
    centroid.id = ids[i];

    if ((!input.include_zero) && (centroid.id == 0)) continue;

    // get label name from color table
    if (ctab) {
      CTABisEntryValid(ctab, centroid.id, &valid_id);
      if (!valid_id) {
        std::cerr << "WARNING: cannot find ID " << centroid.id
                  << " in " << input.ctabfile << "... ignoring this label\n";
        continue;
      }
      CTABcopyName(ctab, centroid.id, char_name, sizeof(char_name));
      centroid.labelname = char_name;
      // this is for table column formatting
      label_chars = centroid.labelname.length();
      if (label_chars > max_label_chars) max_label_chars = label_chars;
    }

    // this is also for table column formatting
    std::ostringstream ss;
    ss << centroid.id;
    id_chars = ss.str().length();
    // the line below is only c++11 compatible :(
    // id_chars = (std::to_string(centroid.id).length());
    if (id_chars > max_id_chars) max_id_chars = id_chars;

    centroid.mass = 0;
    centroids[centroid.id] = centroid;
  }

  std::cout << "computing centers of mass for " << centroids.size() << " structures\n";

  // iterate through voxels
  int voxid;
  for (int col = 0 ; col < seg->width ; col++) {
    for (int row = 0 ; row < seg->height ; row++) {
      for (int slice = 0 ; slice < seg->depth ; slice++) {
        voxid = MRIgetVoxVal(seg, col, row, slice, 0);
        if (centroids.find(voxid) != centroids.end()) {
          MRIvoxelToWorld(seg, col, row, slice, &x, &y, &z);
          // apply voxel weighting if provided
          if (weights) {
            MRIworldToVoxel(weights, x, y, z, &wx, &wy,& wz);
            MRIsampleVolume(weights, wx, wy, wz, &weight);
          } else {
            weight = 1;
          }
          // apply linear transform if provided
          if (lta) {
            // unfortunately, LTAworldToWorldEx only accepts floats
            LTAworldToWorldEx(lta, x, y, z, &fx, &fy, &fz);
            x = fx; y = fy; z = fz;
          }
          centroids[voxid].x += x * weight;
          centroids[voxid].y += y * weight;
          centroids[voxid].z += z * weight;
          centroids[voxid].mass += weight;
        }
      }
    }
  }

  MRIfree(&seg);
  if (lta) LTAfree(&lta);

  // compute centers of mass
  std::map<int, Centroid>::iterator it;
  for (it = centroids.begin(); it != centroids.end(); it++) {
    Centroid *c = &it->second;
    c->x /= c->mass;
    c->y /= c->mass;
    c->z /= c->mass;
  }


  // -------------------- write table --------------------


  std::ofstream tablefile(input.outfile);
  if (!tablefile.is_open()) {
    std::cerr << "could not open writable file at " << input.outfile << std::endl;
    exit(1);
  }

  std::cout << "writing results to " << input.outfile << std::endl;

  int precision = 4;
  int cwidth = precision + 8;

  // table header

  if(input.ltafile.length()) tablefile << "# Reg " << input.ltafile << std::endl;
  tablefile << std::setw(max_id_chars) << "# ID";
  if (ctab) {
    tablefile << std::left << "    " << std::setw(max_label_chars)
              << "StructName" << std::right;
  }
  tablefile << " " << std::setw(cwidth) << "R" + std::string(precision, ' ')
            << " " << std::setw(cwidth) << "A" + std::string(precision, ' ')
            << " " << std::setw(cwidth) << "S" + std::string(precision, ' ')
            << std::endl;

  int table_width = max_id_chars + (cwidth * 3) + 3;
  if (ctab) table_width += max_label_chars + 4;
  tablefile << "# ";
  tablefile << std::string(table_width, '-') << std::endl;

  // table body

  tablefile << std::fixed << std::setprecision(precision);
  for (it = centroids.begin(); it != centroids.end(); it++) {
    Centroid c = it->second;

    tablefile << std::setw(max_id_chars) << c.id;
    if (ctab) {
      tablefile << std::left << "    " << std::setw(max_label_chars) 
                << c.labelname << std::right;
    }
    tablefile << " " << std::setw(cwidth) << c.x
              << " " << std::setw(cwidth) << c.y
              << " " << std::setw(cwidth) << c.z
              << std::endl;
  }

  tablefile.close();

  // pointset
  if(!input.pointset.empty()) {
    FILE *fp = fopen(input.pointset.c_str(),"w");
    for (it = centroids.begin(); it != centroids.end(); it++) {
      Centroid c = it->second;
      fprintf(fp,"%7.4f %7.4f %7.4f\n",c.x,c.y,c.z);
    }
    fprintf(fp,"info\n");
    fprintf(fp,"numpoints %d\n",(int)centroids.size());
    fprintf(fp,"UseRealRAS 1\n");
    fclose(fp);
  }


  std::cout << "mri_segcentroids done" << std::endl;

  return 0;
}
