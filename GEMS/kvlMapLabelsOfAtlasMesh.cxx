/**
 * @file  kvlMapLabelsOfAtlasMesh.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:40 $
 *    $Revision: 1.3 $
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
#include "kvlAtlasMeshCollection.h"
#include "kvlCompressionLookupTable.h"



int main( int argc, char* argv[] )
{

  //
  if ( argc < 3 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " meshCollectionFileName compressionLookupTable [ mode=0 ]" << std::endl;
    std::cerr << "  where mode 0 will map to global tissue types " << std::endl;
    std::cerr << "  where mode 1 will map to right hippocampus " << std::endl;
    std::cerr << "  where mode 2 will map to left hippocampus " << std::endl;
    exit( -1 );
  }

  // Read the mesh collection
  const std::string  meshCollectionFileName( argv[ 1 ] );
  kvl::AtlasMeshCollection::Pointer  collection = kvl::AtlasMeshCollection::New();
  if ( !collection->Read( meshCollectionFileName.c_str() ) )
  {
    std::cerr << "Couldn't read mesh collection from file " << meshCollectionFileName << std::endl;
    exit( -1 );
  }

  // Read the compressionLookupTable
  const std::string  compressionLookupTableFileName( argv[ 2 ] );
  kvl::CompressionLookupTable::Pointer  compressor = kvl::CompressionLookupTable::New();
  if ( !compressor->Read( compressionLookupTableFileName.c_str() ) )
  {
    std::cerr << "Couldn't read compressionLookupTable from file " << compressionLookupTableFileName << std::endl;
    exit( -1 );
  }

  int mode = 0;
  if ( argc > 3 )
  {
    std::istringstream  modeStream( argv[ 3 ] );
    modeStream >> mode;
  }

  // Initialize new alphas with zeros
  kvl::AtlasMesh::PointDataContainer::ConstPointer  oldParameters = collection->GetPointParameters();
  kvl::AtlasMesh::PointDataContainer::Pointer  newParameters = kvl::AtlasMesh::PointDataContainer::New();
  collection->SetPointParameters( newParameters );
  int  newNumberOfLabels = 0;
  if ( mode == 0 )
  {
#if 1
    newNumberOfLabels = 4;
#else
    newNumberOfLabels = 5;
#endif
  }
  else
  {
    newNumberOfLabels = 2;
  }
  kvl::AtlasAlphasType  emptyAlphas( newNumberOfLabels );
  emptyAlphas.Fill( 0.0f );
  for ( kvl::AtlasMesh::PointDataContainer::ConstIterator it = oldParameters->Begin();
        it != oldParameters->End(); ++it )
  {
    kvl::PointParameters  params = it.Value();
    params.m_Alphas = emptyAlphas;
    newParameters->InsertElement( it.Index(), params );
  }


  // Loop over all entries in the lookup table.
  for ( kvl::CompressionLookupTable::CompressionLookupTableType::const_iterator  compressionIt
        = compressor->GetCompressionLookupTable().begin();
        compressionIt != compressor->GetCompressionLookupTable().end(); ++compressionIt )
  {
    // Determine the compressed label (i.e. the index in the old alphas)
    kvl::CompressionLookupTable::CompressedImageType::PixelType  oldIndex = ( *compressionIt ).second;

    // Determine which index in the new alphas the contribution should go to
    int  newIndex = 0;
    if ( mode == 0 )
    {
#if 0
      switch ( ( *compressionIt ).first )
      {
      case 24:  // CSF
        //case 122: // CSF-SA
      case 197: // Right-hippocampal_fissure
      case 215: // hippocampal_fissure
      case 505: // right_hippocampal_fissure
        //case 43:  // Right-Lateral-Ventricle
        //case 60:  // Right-VentralDC
        //case 76:  // Right-Lateral-Ventricles
        //case 221: // Inf_Lat_Vent
        //case 63:  // Right-choroid-plexus
        //case 44:  // Right-Inf-Lat-Vent
        newIndex = 0;
        break;
        //case 221: // Inf_Lat_Vent
      case 63:  // Right-choroid-plexus
      case 44:  // Right-Inf-Lat-Vent
        newIndex = 1;
        break;
      case 199: // Right-subiculum
      case 203: // parasubiculum
      case 204: // presubiculum
      case 205: // subiculum
      case 504: // right_presubiculum
      case 507: // right_subiculum
      case 198: // Right-CADG-head
      case 206: // CA1
      case 207: // CA2
      case 208: // CA3
      case 209: // CA4
      case 500: // right_CA2/3
      case 502: // right_CA1
      case 506: // right_CA4/DG
        //case 63:  // Right-choroid-plexus
      case 53:  // Right-Hippocampus
      case 42:  // Right-Cerebral-Cortex
        newIndex = 2;
        break;
      case 200: // Right-fimbria
      case 212: // fimbria
      case 503: // right_fimbria
      case 41:  // Right-Cerebral-White-Matter
      case 46:  // Right-Cerebellum-White-Matter
      case 219: // Cerebral_White_Matter
      case 223: // Cerebral_White_Matter_Edge
      case 5001: // Left-UnsegmentedWhiteMatter
      case 5002: // Right-UnsegmentedWhiteMatter
        //case 500: // right_CA2/3
        newIndex = 3;
        break;
      default :
        // We don't care about this label, put it in the last bag
        std::cerr << "Help: encountered a label that's not supposed to be here: " << ( *compressionIt ).first << std::endl;
        exit( -1 );
      }
#else

      //
      // About VentralDC: from http://www.cma.mgh.harvard.edu/manuals/segmentation
      //
      // The ventral diencephalon (VDC) is not an anatomical name for a single structure but a name given by the
      // CMA to a group of structures that generally cannot be distinguished from each other with standard MRI
      // images.  This "miscellaneous" area includes the hypothalamus, mammillary body, subthalamic nuclei,
      // substantia nigra, red nucleus, lateral geniculate nucleus (LGN), and medial geniculate nucleus (MGN).
      // White matter areas such as the zona incerta, cerebral peduncle (crus cerebri), and the lenticular fasciculus
      // are also included in this area.  The optic tract is included in this area in the most anterior extent.
      // Each structure fades in and out of the VDC at different times. Therefore, the VDC greatly varies from slice to slice.


      // Determine which index in the new alphas the contribution should go to
      switch ( ( *compressionIt ).first )
      {
      case 0: // Unknown
      case 30: // Left-vessel
      case 62: // Right-vessel
      case 29: // Left-undetermined
      case 61: // Right-undetermined
      case 85: // Optic-Chiasm
        newIndex = 0;
        break;
      case 4: // Left-Lateral-Ventricle
      case 43: // Right-Lateral-Ventricle
      case 5: // Left-Inf-Lat-Vent
      case 44: // Right-Inf-Lat-Vent
      case 14: // 3rd-Ventricle
      case 15: // 4th-Ventricle
      case 24: // CSF
      case 72: // 5th-Ventricle
        newIndex = 1;
        break;
      case 3: // Left-Cerebral-Cortex
      case 42: // Right-Cerebral-Cortex
      case 8: // Left-Cerebellum-Cortex
      case 47: // Right-Cerebellum-Cortex
      case 11: // Left-Caudate
      case 50: // Right-Caudate
      case 12: // Left-Putamen
      case 51: // Right-Putamen
      case 17: // Left-Hippocampus
      case 53: // Right-Hippocampus
      case 18: // Left-Amygdala
      case 54: // Right-Amygdala
      case 26: // Left-Accumbens-area
      case 58: // Right-Accumbens-area
      case 80: // non-WM-hypointensities
      case 10: // Left-Thalamus-Proper
      case 49: // Right-Thalamus-Proper
      case 31: // Left-choroid-plexus
      case 63: // Right-choroid-plexus
        newIndex = 2;
        break;
      case 2: // Left-Cerebral-White-Matter
      case 41: // Right-Cerebral-White-Matter
      case 7: // Left-Cerebellum-White-Matter
      case 46: // Right-Cerebellum-White-Matter
      case 16: // Brain-Stem
      case 28: // Left-VentralDC
      case 60: // Right-VentralDC
      case 13: // Left-Pallidum
      case 52: // Right-Pallidum
      case 77: // WM-hypointensities
        newIndex = 3;
        break;
      default :
        // We don't care about this label, put it in the last bag
        std::cerr << "Help: encountered a label that's not supposed to be here: " << ( *compressionIt ).first << std::endl;
        exit( -1 );
      }
#endif
    }
    else if ( mode == 1 )
    {
      // Mapping right hippcampal labels
      switch ( ( *compressionIt ).first )
      {
      case 53: // Right-Hippocampus
      case 500: // right_CA2-3
      case 502: // right_CA1
      case 503: // right_fimbria
      case 504: // right_presubiculum
      case 505: // right_hippocampal_fissure
      case 506: // right_CA4-DG
      case 507: // right_subiculum
        newIndex = 1;
        break;
      default :
        newIndex = 0;
      }

    }
    else
    {
      // Mapping left hippcampal labels
      switch ( ( *compressionIt ).first )
      {
      case 17: // Left-Hippocampus
      case 550: // left_CA2-3
      case 552: // left_CA1
      case 553: // left_fimbria
      case 554: // left_presubiculum
      case 555: // left_hippocampal_fissure
      case 556: // left_CA4-DG
      case 557: // left_subiculum
        newIndex = 1;
        break;
      default :
        newIndex = 0;
      }

    }

    // Now loop over all points, and add contributions of old alphas to new alphas
    kvl::AtlasMesh::PointDataContainer::ConstIterator  oldIt = oldParameters->Begin();
    kvl::AtlasMesh::PointDataContainer::Iterator  newIt = newParameters->Begin();
    for ( ; oldIt != oldParameters->End(); ++oldIt, ++newIt )
    {
      newIt.Value().m_Alphas[ newIndex ] += oldIt.Value().m_Alphas[ oldIndex ];
    }


  } // End loop over entries of the lookup table



  // Write out
  std::ostringstream  outputFileNameStream;
  outputFileNameStream << meshCollectionFileName << "_mapped.txt";
  if ( !collection->Write( outputFileNameStream.str().c_str() ) )
  {
    std::cerr << "Couldn't write to " <<  outputFileNameStream.str() << std::endl;
    exit( -1 );
  }
  std::cerr << "Just wrote to " <<  outputFileNameStream.str() << std::endl;



  return 0;
};

