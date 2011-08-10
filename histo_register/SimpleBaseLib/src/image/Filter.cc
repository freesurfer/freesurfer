// Licensed under MIT license; see license.txt.

#include <sbl/image/Filter.h>
#include <sbl/core/PathConfig.h>
#include <sbl/core/Command.h>
#include <sbl/math/MathUtil.h>
#include <sbl/image/Video.h>
namespace sbl {


// a generalized image filter
struct ColorImageFilter {
    String name;
    ColorImageFilterCallback callback;
};


// a generalized image filter
struct GrayImageFilter {
    String name;
    GrayImageFilterCallback callback;
};


// set of generalized filters for color images
Array<ColorImageFilter> &colorFilters() {
    static aptr<Array<ColorImageFilter> > s_colorImageFilters;
    if (s_colorImageFilters.get() == NULL) {
        s_colorImageFilters.reset( new Array<ColorImageFilter>() ); 
    }
    return *s_colorImageFilters;
}


// set of generalized filters for masks / grayscale images
Array<GrayImageFilter> &grayFilters() {
    static aptr<Array<GrayImageFilter> > s_grayImageFilters;
    if (s_grayImageFilters.get() == NULL) {
        s_grayImageFilters.reset( new Array<GrayImageFilter>() ); 
    }
    return *s_grayImageFilters;
}


/// register a generalized filter for color images
void registerFilterInternal( const String &name, ColorImageFilterCallback callback ) {
    ColorImageFilter *filter = new ColorImageFilter();
    filter->name = name;
    filter->callback = callback;
    colorFilters().append( filter );
}


/// register a generalized filters for masks / grayscale images
void registerFilterInternal( const String &name, GrayImageFilterCallback callback ) {
    GrayImageFilter *filter = new GrayImageFilter();
    filter->name = name;
    filter->callback = callback;
    grayFilters().append( filter );
}


//-------------------------------------------
// COMMANDS
//-------------------------------------------


// apply a generalized filter to a video
void filterVideo( Config &conf, bool gray ) {

    // get command parameters
    const String &filterName = conf.readString( "filterName" );
    String inputFileName = addDataPath( conf.readString( "inputFileName" ) );
    String outputFileName = addDataPath( conf.readString( "outputFileName" ) );
    int startFrameIndex = conf.readInt( "startFrameIndex", 0 );
    int endFrameIndex = conf.readInt( "endFrameIndex", 1000000000 );
    int frameStep = conf.readInt( "frameStep", 1 );
    if (conf.initialPass())
        return;

    // display parameters
    disp( 1, "filterName: %s", filterName.c_str() );
    disp( 1, "inputFileName: %s", inputFileName.c_str() );
    disp( 1, "outputFileName: %s", outputFileName.c_str() );
    disp( 1, "startFrameIndex: %d", startFrameIndex );
    disp( 1, "endFrameIndex: %d", endFrameIndex );
    disp( 1, "frameStep: %d", frameStep );

    // open input video
    InputVideo inputVideo( inputFileName );
    if (inputVideo.openSuccess() == false) {
        warning( "unable to open input video: %s", inputFileName.c_str() );
        return;
    }
    startFrameIndex = bound( startFrameIndex, 0, inputVideo.length() - 1 );
    endFrameIndex = bound( endFrameIndex, 0, inputVideo.length() - 1 );
    assertAlways( frameStep > 0 ); // fix(later): display warning, handle negative frameStep

    // open output video
    // fix(clean): would be nice to defer output size until have output image so allow resizing
    OutputVideo outputVideo( outputFileName, inputVideo.width(), inputVideo.height() );
    if (outputVideo.openSuccess() == false) {
        warning( "unable to open output video: %s", outputFileName.c_str() );
        return;
    }

    // get filter by filter name
    ColorImageFilterCallback colorFilterCallback = NULL;
    GrayImageFilterCallback grayFilterCallback = NULL;
    if (gray) {
        for (int i = 0; i < grayFilters().count(); i++) {
            if (grayFilters()[ i ].name == filterName)
                grayFilterCallback = grayFilters()[ i ].callback;
        }
        if (grayFilterCallback == NULL) {
            warning( "unable to find grayscale image filter: %s", filterName.c_str() );
            return;
        }
    } else {
        for (int i = 0; i < colorFilters().count(); i++) {
            if (colorFilters()[ i ].name == filterName)
                colorFilterCallback = colorFilters()[ i ].callback;
        }
        if (colorFilterCallback == NULL) {
            warning( "unable to find color image filter: %s", filterName.c_str() );
            return;
        }
    }

    // loop over input video
    for (int frameIndex = startFrameIndex; frameIndex <= endFrameIndex; frameIndex += frameStep) {
        if (gray) {
            aptr<ImageGrayU> image = inputVideo.frameGray( frameIndex );
            image = grayFilterCallback( *image, conf );
            outputVideo.append( *image );
        } else {
            aptr<ImageColorU> image = inputVideo.frame( frameIndex );
            image = colorFilterCallback( *image, conf );
            outputVideo.append( *image );
        }
    }
}


// apply a generalized filter to a video
void filterColorVideo( Config &conf ) {
    filterVideo( conf, false );    
}


// apply a generalized filter to a mask sequence
void filterGrayVideo( Config &conf ) {
    filterVideo( conf, true );    
}


//-------------------------------------------
// INIT / CLEAN-UP
//-------------------------------------------


// register commands, etc. defined in this module
void initFilter() {
    registerCommand( "vfiltcolor", filterColorVideo );
    registerCommand( "vfiltgray", filterColorVideo );
}


} // end namespace sbl

