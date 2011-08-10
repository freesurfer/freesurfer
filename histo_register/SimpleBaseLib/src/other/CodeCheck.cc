// Licensed under MIT license; see license.txt.

#include <sbl/other/CodeCheck.h>
#include <sbl/core/Command.h>
#include <sbl/core/Table.h>
#include <sbl/core/PathConfig.h>
#include <sbl/system/FileSystem.h>
namespace sbl {


/// get list of source files in given path and sub-dirs
Array<String> codeFileList( const String &path, const Array<String> &suffixes, bool excludeUnderscore ) {
    Array<String> fileNames;
    Array<String> dirList = dirFileList( path, "", "" );
    for (int i = 0; i < dirList.count(); i++) {
        const String &dirItem = dirList[ i ];

        // if no dot, try expanding as directory
        // fix(later): get dirs and files separately
        if (dirItem.contains( '.' ) == false && (excludeUnderscore == false || dirItem.startsWith( "_" ) == false)) {
            Array<String> subFileNames = codeFileList( path + "/" + dirItem, suffixes, excludeUnderscore );
            for (int j = 0; j < subFileNames.count(); j++) // fix(clean): replace with append array function 
                fileNames.appendCopy( subFileNames[ j ] );

        // if dot, assume file
        } else if (excludeUnderscore == false || dirItem.startsWith( "_" ) == false) {

            // if has one of the specified extensions
            for (int j = 0; j < suffixes.count(); j++) {
                if (dirItem.endsWith( suffixes[ j ] ))
                    fileNames.appendCopy( path + "/" + dirItem );
            }
        }
    }
    return fileNames;
}


/// check a source code file for coding conventions
void checkFile( const String &fileName, Table &table ) {
    bool header = fileName.endsWith( ".h" );

    // accumulate stats
    int lineCount = 0;
    int commentCount = 0;
    int missingCommentCount = 0;
    int namespaceCount = 0;
    int fixCleanerCount = 0;
    int fixFasterCount = 0;
    int fixSoonCount = 0;
    int fixOtherCount = 0;
    int badIncludeFormatCount = 0;
    int externalIncludeCount = 0;
    int badFormatCount = 0;
    int badTabCount = 0;
    int maxFuncLen = 0;

    // current state
    bool lastBlank = false;
    bool inFunction = 0;
    int funcLen = 0;

    // loop over the file
    File file( fileName, FILE_READ, FILE_TEXT );
    if (file.openSuccess()) {
        while (file.endOfFile() == false) {
            String line = file.readLine();
            String lineStrip = line.strip();
            lineCount++;

            // check for blank line
            if (lineStrip.length() == 0) {
                lastBlank = true;
                continue;
            }

            // check for namespace
            if (line.startsWith( "namespace " ))
                namespaceCount++;

            // check for comments
            bool comment = lineStrip.contains( "//" );
            if (comment) {
                commentCount++;
            } else if (lastBlank == true && lineStrip != "private:" && lineStrip != "protected:" ) {
                missingCommentCount++;
            }

            // check for fix messages
            if (comment) {
                if (lineStrip.contains( "// fix" )) {
                    if (lineStrip.contains( "fix(cleaner)" ) || lineStrip.contains( "fix(clean)" ))
                        fixCleanerCount++;
                    else if (lineStrip.contains( "fix(faster)" ))
                        fixFasterCount++;
                    else if (lineStrip.contains( "fix(soon)" ))
                        fixSoonCount++;
                    else 
                        fixOtherCount++; 
                }
            }

            // make sure using tabs instead of spaces (could also check inside line)
            if (line.get( 0 ) == ' ') {
                badTabCount++;
            }

            // check for function start/end (will also catch class/struct)
            if (header == false && line.get( 0 ) > ' ' && line.startsWith( "namespace" ) == false) {
                if (lineStrip.endsWith( "{" ))
                    inFunction = true;
                if (lineStrip.endsWith( "}" ))
                    inFunction = false;
            }

            // look for longest function
            if (inFunction) {
                funcLen++;
                if (funcLen > maxFuncLen) 
                    maxFuncLen = funcLen;
            } else {
                funcLen = 0;
            }

            // check includes
            if (lineStrip.startsWith( "#include" )) {
                String includeFile = line.rightOfFirst( ' ' );
                if (includeFile.startsWith( "\"" ))
                    badIncludeFormatCount++;
                if (includeFile.startsWith( "<sbl" ) == false && includeFile.startsWith( "<pvl" ) == false) {
                    externalIncludeCount++;
                }
            }
            
            // check statement formatting
            if (lineStrip.startsWith( "if" ) || lineStrip.startsWith( "for" ) || lineStrip.startsWith( "while" )) {
                if (lineStrip.startsWith( "if(" ) || lineStrip.startsWith( "if ( " )) {
                    disp( 1, "%s, line: %d; bad if", fileName.c_str(), lineCount );
                    badFormatCount++;
                }
                if (lineStrip.startsWith( "for(" ) || lineStrip.startsWith( "for ( " )) {
                    disp( 1, "%s, line: %d; bad for", fileName.c_str(), lineCount );
                    badFormatCount++;
                }
                if (lineStrip.startsWith( "while(" ) || lineStrip.startsWith( "while ( " )) {
                    disp( 1, "%s, line: %d; bad while", fileName.c_str(), lineCount );
                    badFormatCount++;
                }
                if (lineStrip.endsWith( " ) {" ) || lineStrip.endsWith( " )" )) {
                    disp( 1, "%s, line: %d; bad closing paren", fileName.c_str(), lineCount );
                    badFormatCount++;
                }
            }

            // MORE LATER:
            // - function name capitalization
            // - function arg spacing
            // - class method order: public / protected / private
            // - that destructors of base classes are virtual
            // - order of functions matches .cc and .h
            // - enum naming
            // - copy constructor and assignment operator defined for each class
            // - number of consecutive blank lines
            // - number of lines enclosed in /* */ comments

            // this line was not blank
            lastBlank = false;
        }
    }

    // add summary stats
    table.add( "file_name", fileName );
    table.add( "lines", lineCount );
    table.add( "comments", commentCount );
    table.add( "missing_comments", missingCommentCount );
    table.add( "fix_soon", fixSoonCount );
    table.add( "fix_other", fixOtherCount );
    table.add( "fix_cleaner", fixCleanerCount );
    table.add( "fix_faster", fixFasterCount );
    table.add( "namespaces", namespaceCount );
    table.add( "external_include", externalIncludeCount );
    table.add( "bad_include_format", badIncludeFormatCount );
    table.add( "bad_format", badFormatCount );
    table.add( "bad_tab", badTabCount );
    table.add( "max_func_len", maxFuncLen );
}


// check coding conventions for a set of source files
void checkCodeFormatting( Config &conf ) {

    // get command parameters
    const String &path = addDataPath( conf.readString( "path" ) );
    Array<String> suffixes = conf.readStrings( "suffixes" );
    bool excludeUnderscore = conf.readBool( "excludeUnderscore", true );
    if (conf.initialPass())
        return;
    
    // table of summary stats 
    Table table( "Code Formatting" );

    // loop over all source files
    Array<String> fileNames = codeFileList( path, suffixes, excludeUnderscore );
    disp( 1, "processing %d files...", fileNames.count() );
    for (int fileIndex = 0; fileIndex < fileNames.count(); fileIndex++) 
        checkFile( fileNames[ fileIndex ], table );

    // save the summary
    table.saveCSV( "codeFormatting.csv" );
}


// prepare code for distribution
void distributeCode( Config &conf ) {

    // get command parameters
    const String &sourcePath = addDataPath( conf.readString( "sourcePath" ) );
    const String &destPath = addDataPath( conf.readString( "destPath" ) );
    Array<String> suffixes = conf.readStrings( "suffixes" );
    bool excludeUnderscore = conf.readBool( "excludeUnderscore", true );
    if (conf.initialPass())
        return;

    // loop over all source files
    Array<String> fileNames = codeFileList( sourcePath, suffixes, excludeUnderscore );
    for (int fileIndex = 0; fileIndex < fileNames.count(); fileIndex++) {
        String inFileName = fileNames[ fileIndex ];
        String outFileName = destPath + inFileName.rightOf( sourcePath.length() - 1 );
        disp( 1, "%s->%s", inFileName.c_str(), outFileName.c_str() );
        File inFile( inFileName, FILE_READ, FILE_TEXT );
        File outFile( outFileName, FILE_WRITE, FILE_TEXT );

        // add license
        if (inFileName.contains( "external" ) == false) {
            File licenseFile( sourcePath + "../../license-short.txt", FILE_READ, FILE_TEXT ); // fix(later): path
            while (licenseFile.endOfFile() == false) 
                outFile.writeF( "%s\x0A", licenseFile.readLine().c_str() ); // write linefeed only at end of line
        }

        // process source file
        int lineCount = 0, tabCount = 0;
        while (inFile.endOfFile() == false) {
            String line = inFile.readLine();
            lineCount++;

            // if contains tab, expand tabs to spaces
            if (line.contains( '\t' )) {
                int len = line.length();
                String expandedLine( 0, len + 100 ); // fix(later): count tabs
                int pos = 0;
                for (int i = 0; i < len; i++) {
                    unsigned short c = line.get( i );
                    if (c == '\t') {
                        tabCount++;
                        expandedLine.set( pos++, ' ' ); // assume 4 spaces per tab
                        expandedLine.set( pos++, ' ' );
                        expandedLine.set( pos++, ' ' );
                        expandedLine.set( pos++, ' ' );
                    } else {
                        expandedLine.set( pos++, c );
                    }
                }
                expandedLine.set( pos, 0 ); // fix(later): should update length?
                line = expandedLine;
            }

            // write with line feed at end
            outFile.writeRawString( line.c_str() );
            outFile.writeRawString( "\x0A" );
        }

        // display diagnostic
        disp( 1, "inFileName: %s, lineCount: %d, tabCount: %d", inFileName.c_str(), lineCount, tabCount );
    }

    // remove the original files
    for (int fileIndex = 0; fileIndex < fileNames.count(); fileIndex++) {
        String fileName = fileNames[ fileIndex ];
        deleteFile( fileName + ".orig" );
    }
}


//-------------------------------------------
// INIT / CLEAN-UP
//-------------------------------------------


// register commands, etc. defined in this module
void initCodeCheck() {
    registerCommand( "codecheck", checkCodeFormatting );
    registerCommand( "codedist", distributeCode );
}


} // end namespace sbl

